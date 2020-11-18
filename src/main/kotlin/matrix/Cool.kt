package matrix

import hdf.hdf5lib.H5
import hdf.hdf5lib.HDF5Constants
import kotlin.reflect.full.memberProperties

data class Bins(
    val chrom: IntArray,
    val end: IntArray,
    val start: IntArray,
    val weight: DoubleArray
)

data class Chroms(
    val length: IntArray,
    val name: Array<String>
)

data class Indexes(
    val bin1_offset: LongArray,
    val chrom_offset: LongArray
)

data class Pixels(
    val bin1_id: LongArray,
    val bin2_id: LongArray,
    val count: IntArray
)

data class Cool(val groupId: Int) {
    val bins: Bins
    val chroms: Chroms
    val indexes: Indexes
    val pixels: Pixels

    init {
        val (chrom, end, start, weight) = getGroup(groupId, "bins")
        bins = Bins(chrom as IntArray, end as IntArray, start as IntArray, weight as DoubleArray)

        val (length, name) = getGroup(groupId, "chroms")
        chroms = Chroms(length as IntArray, name as Array<String>)

        val (bin1_offset, chrom_offset) = getGroup(groupId, "indexes")
        indexes = Indexes(bin1_offset as LongArray, chrom_offset as LongArray)

        val (bin1_id, bin2_id, count) = getGroup(groupId, "pixels")
        pixels = Pixels(bin1_id as LongArray, bin2_id as LongArray, count as IntArray)
    }
}

fun getDataset(datasetId: Int, datasetName: String): Any {
    val dspace = H5.H5Dget_space(datasetId)
    val rank = H5.H5Sget_simple_extent_ndims(dspace)

    val dims = longArrayOf(rank.toLong())
    H5.H5Sget_simple_extent_dims(dspace, dims, null)
    val size = dims[0].toInt()

    val tid = H5.H5Dget_type(datasetId)

    val dataRead: Any
    when (datasetName) {
        "weight" -> dataRead = DoubleArray(size)
        "name" -> {
            // Figure out the string size and number of strings
            val stringLength = H5.H5Tget_size(tid)
            val bufferSize = size * stringLength
            val byteBuff = ByteArray(bufferSize)

            // Read the string data into byte buff
            H5.H5Dread(
                datasetId, tid,
                HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT,
                byteBuff
            )

            // Convert byte array into string array
            dataRead = arrayOfNulls<String>(size)
            for (i in 0 until size) {
                dataRead[i] = String(byteBuff, i * stringLength, stringLength)
            }
            return dataRead
        }
        else -> {
            if (datasetName.startsWith("bin") || datasetName == "chrom_offset") {
                dataRead = LongArray(size)
            } else {
                dataRead = IntArray(size)
            }
        }

    }
    if (datasetId >= 0) H5.H5Dread(
        datasetId, tid,
        HDF5Constants.H5S_ALL, HDF5Constants.H5S_ALL, HDF5Constants.H5P_DEFAULT,
        dataRead
    )

    return dataRead
}

private fun getGroup(groupId: Int, groupName: String): ArrayList<Any> {
    var groupId = groupId
    if (groupId >= 0) groupId = H5.H5Gopen(groupId, groupName, HDF5Constants.H5P_DEFAULT)

    val members = H5.H5Gget_info(groupId)
    val nlinks = members.nlinks.toInt()
    val objNames = arrayOfNulls<String>(nlinks)
    val objTypes = IntArray(nlinks)
    val objRefs = LongArray(nlinks)

    var names_found = 0
    try {
        names_found = H5.H5Gget_obj_info_all(
            groupId, null, objNames,
            objTypes, null, objRefs, HDF5Constants.H5_INDEX_NAME
        )
    } catch (err: Throwable) {
        err.printStackTrace()
    }

    val datasetList = ArrayList<Any>()
    for (i in 0 until names_found) {
        if (objTypes[i] == HDF5Constants.H5O_TYPE_DATASET) {
            val datasetId = H5.H5Dopen(groupId, objNames[i], HDF5Constants.H5P_DEFAULT)

            datasetList.add(getDataset(datasetId, objNames[i].toString()))

            H5.H5Dclose(datasetId)
        }
    }
    return datasetList
}

fun writeCoolToFile(cool: Cool, fname: String, resolution: Int) {
    var fileId = -1
    var groupId1 = -1
    var groupId2 = -1
    // Create a new file using default properties.
    try {
        fileId = H5.H5Fcreate(
            fname, HDF5Constants.H5F_ACC_TRUNC,
            HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
        )
    } catch (e: Exception) {
        e.printStackTrace()
        System.err.println("Failed to create file:$fname")
        return
    }

    // Create a group in the file.
    try {
        if (fileId >= 0) {
            groupId1 = H5.H5Gcreate(
                fileId, "resolutions",
                HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
            )

            groupId2 = H5.H5Gcreate(
                groupId1, "$resolution",
                HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
            )
        }
    } catch (e: Exception) {
        e.printStackTrace()
    }

    val groups = arrayOf(cool.bins, cool.chroms, cool.indexes, cool.pixels)
    val groupNames = arrayOf("bins", "chroms", "indexes", "pixels")

    groups.forEachIndexed { index, group ->
        val gid = H5.H5Gcreate(
            groupId2, groupNames[index],
            HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
        )
        group::class.memberProperties.forEach { dataset ->
            // Create the data space for the 1D dataset
            var data = dataset.getter.call(group)
            val rank = 1


            when (dataset.name) {
                "weight" -> {
                    val dims = longArrayOf((data as DoubleArray).size.toLong())
                    val dataspaceId = H5.H5Screate_simple(rank, dims, null)

                    val datatype = HDF5Constants.H5T_NATIVE_DOUBLE
                    val datasetId = H5.H5Dcreate(
                        gid, dataset.name, datatype, dataspaceId,
                        HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
                    )

                    H5.H5Dwrite_double(
                        datasetId,
                        datatype,
                        HDF5Constants.H5S_ALL,
                        dataspaceId,
                        HDF5Constants.H5P_DEFAULT,
                        data
                    )

                    H5.H5Sclose(dataspaceId)
                    H5.H5Dclose(datasetId)
                }
                "name" -> {
                    val dims = longArrayOf((data as Array<String>).size.toLong())
                    val dataspaceId = H5.H5Screate_simple(rank, dims, null)

                    val strLength = 10
                    val datatype = H5.H5Tcopy(HDF5Constants.H5T_C_S1)
                    H5.H5Tset_size(datatype, strLength)

                    val datasetId = H5.H5Dcreate(
                        gid, dataset.name, datatype, dataspaceId,
                        HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
                    )

                    // TODO: 14.11.2020 expand with whitespaces to strLength
                    val byteBuff = data.joinToString("").toByteArray()
                    H5.H5Dwrite(
                        datasetId,
                        datatype,
                        HDF5Constants.H5S_ALL,
                        HDF5Constants.H5S_ALL,
                        HDF5Constants.H5P_DEFAULT,
                        byteBuff
                    );

                    H5.H5Sclose(dataspaceId)
                    H5.H5Dclose(datasetId)
                }
                else -> {
                    val dims: LongArray
                    val datatype: Int
                    if (dataset.name.startsWith("bin") || dataset.name == "chrom_offset") {
                        datatype = HDF5Constants.H5T_NATIVE_INT64
                        dims = longArrayOf((data as LongArray).size.toLong())
                    } else {
                        datatype = HDF5Constants.H5T_NATIVE_INT32
                        dims = longArrayOf((data as IntArray).size.toLong())
                    }
                    val dataspaceId = H5.H5Screate_simple(rank, dims, null)

                    val datasetId = H5.H5Dcreate(
                        gid, dataset.name, datatype, dataspaceId,
                        HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
                    )

                    H5.H5Dwrite(
                        datasetId,
                        datatype,
                        HDF5Constants.H5S_ALL,
                        dataspaceId,
                        HDF5Constants.H5P_DEFAULT,
                        data
                    )

                    H5.H5Dclose(datasetId)
                    H5.H5Sclose(dataspaceId)
                }
            }

        }
        if (gid >= 0) H5.H5Gclose(gid)

    }
    try {
        if (groupId2 >= 0) H5.H5Gclose(groupId2)
        if (groupId1 >= 0) H5.H5Gclose(groupId1)
        if (fileId >= 0) H5.H5Fclose(fileId)
    } catch (e: Exception) {
        e.printStackTrace()
    }
}

fun inverseWeight(cool: Cool, startChrom: Int, endChrom: Int, groupId: Int) {
    cool.bins.weight.reverse(startChrom, endChrom)

    val bGroupId = H5.H5Gopen(groupId, "bins", HDF5Constants.H5P_DEFAULT)

    val datasetId = H5.H5Dopen(bGroupId, "weight", HDF5Constants.H5P_DEFAULT)
    val tid = H5.H5Dget_type(datasetId)
    val dspace = H5.H5Dget_space(datasetId)

    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.bins.weight
    )

    H5.H5Dclose(datasetId)
}

fun inverseHorizontal(cool: Cool, startChrom: Int, endChrom: Int) {
    // TODO: 18.11.2020 think about collection. You can use array in first cycle
    //  and List in second
    val m = mutableMapOf<Int, Pair<MutableList<Long>, MutableList<Int>>>()

    val breakMap = mutableMapOf<Int, Int>()

    for (i in startChrom until endChrom) {
        val startIndex = cool.indexes.bin1_offset[i].toInt()

        var breakIndex = startIndex
        while (breakIndex < cool.pixels.bin1_id.size &&
            cool.pixels.bin1_id[startIndex] == cool.pixels.bin1_id[breakIndex] && cool.pixels.bin2_id[breakIndex] < endChrom
        ) breakIndex += 1
        breakMap[i] = breakIndex
//        val n
        val end = cool.indexes.bin1_offset[i + 1].toInt()

        val indices = breakIndex until end

        val bin2Array =
            cool.pixels.bin2_id.sliceArray(indices).toMutableList()
        val countArray =
            cool.pixels.count.sliceArray(indices).toMutableList()

        m[i] = Pair(bin2Array, countArray)
    }

    for (i in startChrom until endChrom) {
        val startIndex = cool.indexes.bin1_offset[i].toInt()
        val breakIndex = breakMap[i]

        println("row $i")
        for (j in startIndex until breakIndex!!) {
            val bin2Id = cool.pixels.bin2_id[j].toInt()
            println(bin2Id)

            m[bin2Id]!!.first.add(0, (endChrom - 1 - i + startChrom).toLong())
            m[bin2Id]!!.second.add(0, cool.pixels.count[j])
        }
    }

//    println("$startChrom $endChrom")

    var prevIndex = cool.indexes.bin1_offset[startChrom].toInt()
    for (i in endChrom - 1 downTo startChrom) {

        val length = m[i]!!.first.size

        for (j in 0 until length) {
            cool.pixels.bin1_id[prevIndex + j] = m.keys.sorted()[endChrom - 1 - i].toLong()
            cool.pixels.bin2_id[prevIndex + j] = m[i]!!.first[j]
            cool.pixels.count[prevIndex + j] = m[i]!!.second[j]
        }
        prevIndex += length
    }
}

fun inverseVertical(cool: Cool, startChrom: Int, endChrom: Int) {
    // TODO: 17.11.2020 think about multiple `indexOfFirst`
    var inverseIndex = -1
    var swapArray = mutableListOf<Pair<Long, Int>>()
    for (i in 0 until cool.indexes.bin1_offset[startChrom].toInt()) {
        if (inverseIndex == -1 && cool.pixels.bin2_id[i] >= startChrom && cool.pixels.bin2_id[i] < endChrom) {
            inverseIndex = i

            swapArray.add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
            continue
        }

        // TODO: 18.11.2020 check bin1_id changing
        // TODO: 18.11.2020 write logic in comment
        if (inverseIndex != -1) {
            if (cool.pixels.bin1_id[i - 1] == cool.pixels.bin1_id[i] && cool.pixels.bin2_id[i] < endChrom) {
                swapArray.add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
            } else {
                swapArray.asReversed().forEach {
                    cool.pixels.bin2_id[inverseIndex] = endChrom - 1 - (it.first - startChrom)
                    cool.pixels.count[inverseIndex] = it.second
                    inverseIndex += 1

                }
                inverseIndex = -1
                swapArray.clear()
            }
        }
    }
}

fun inverse(cool: Cool, groupId: Int) {
//    val binLength =
//        cool.indexes.chrom_offset.mapIndexed { index, l -> if (index > 0) l - cool.indexes.chrom_offset[index - 1] else l }
//    val w = binLength.indexOfFirst { it > 10 } - 1
    val w = 94

    val startChrom = cool.indexes.chrom_offset[w].toInt()
    val endChrom = cool.indexes.chrom_offset[w + 1].toInt()

    inverseWeight(cool, startChrom, endChrom, groupId)

    inverseHorizontal(cool, startChrom, endChrom)
    inverseVertical(cool, startChrom, endChrom)

    val pGroupId = H5.H5Gopen(groupId, "pixels", HDF5Constants.H5P_DEFAULT)

    var datasetId = H5.H5Dopen(pGroupId, "bin1_id", HDF5Constants.H5P_DEFAULT)
    var tid = H5.H5Dget_type(datasetId)
    var dspace = H5.H5Dget_space(datasetId)

    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.pixels.bin1_id
    )

    H5.H5Dclose(datasetId)

    datasetId = H5.H5Dopen(pGroupId, "bin2_id", HDF5Constants.H5P_DEFAULT)
    tid = H5.H5Dget_type(datasetId)
    dspace = H5.H5Dget_space(datasetId)

    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.pixels.bin2_id
    )

    H5.H5Dclose(datasetId)

    datasetId = H5.H5Dopen(pGroupId, "count", HDF5Constants.H5P_DEFAULT)
    tid = H5.H5Dget_type(datasetId)
    dspace = H5.H5Dget_space(datasetId)

    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.pixels.count
    )

    H5.H5Dclose(datasetId)

    if (pGroupId >= 0) H5.H5Gclose(pGroupId)
}

fun recalculateIndex(cool: Cool, groupId: Int) {
    var prevIndex = 0
    cool.pixels.bin1_id.forEachIndexed { index, l ->
        for (i in prevIndex..l.toInt()) {
            cool.indexes.bin1_offset[i] = index.toLong()
        }
        prevIndex = l.toInt() + 1
    }

    val groupId2 = H5.H5Gopen(groupId, "indexes", HDF5Constants.H5P_DEFAULT)

    val datasetId = H5.H5Dopen(groupId2, "bin1_offset", HDF5Constants.H5P_DEFAULT)
    val tid = H5.H5Dget_type(datasetId)
    val dspace = H5.H5Dget_space(datasetId)

    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.indexes.bin1_offset
    )

    H5.H5Dclose(datasetId)

    if (groupId2 >= 0) H5.H5Gclose(groupId2)
//    println(cool.indexes.bin1_offset.joinToString(separator = "\n"))
}


fun main() {
    val fname = "src/main/resources/2.cool"
    var fileId = -1
    var groupId = -1

    // Open file using the default properties
    try {
        fileId = H5.H5Fopen(fname, HDF5Constants.H5F_ACC_RDWR, HDF5Constants.H5P_DEFAULT)
    } catch (e: Exception) {
        e.printStackTrace()
    }

    // Open group using the default properties
    try {
        if (fileId >= 0) groupId = H5.H5Gopen(fileId, "resolutions/100000", HDF5Constants.H5P_DEFAULT)
    } catch (e: Exception) {
        e.printStackTrace()
    }

    val cool = Cool(groupId)

    inverse(cool, groupId)

    recalculateIndex(cool, groupId)

    // Close objects
    try {
        if (groupId >= 0) H5.H5Gclose(groupId)
        if (fileId >= 0) H5.H5Fclose(fileId)
    } catch (e: Exception) {
        e.printStackTrace()
    }

//    writeCoolToFile(cool, "src/main/resources/mat18_101k.cool", 100000)
}