package matrix

import hdf.hdf5lib.H5
import hdf.hdf5lib.HDF5Constants
import java.io.File
import kotlin.reflect.full.memberProperties

data class Bins(
    val chrom: MutableList<Int>,
    val end: IntArray,
    val start: IntArray,
    val weight: MutableList<Double>
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
    val bin1_id: MutableList<Long>,
    val bin2_id: MutableList<Long>,
    val count: MutableList<Int>
)

data class Cool(val groupId: Int) {
    val bins: Bins
    val chroms: Chroms
    val indexes: Indexes
    val pixels: Pixels

    init {
        val (chrom, end, start, weight) = getGroup(groupId, "bins")
        bins = Bins(chrom as MutableList<Int>, end as IntArray, start as IntArray, weight as MutableList<Double>)

        val (length, name) = getGroup(groupId, "chroms")
        chroms = Chroms(length as IntArray, name as Array<String>)

        val (bin1_offset, chrom_offset) = getGroup(groupId, "indexes")
        indexes = Indexes(bin1_offset as LongArray, chrom_offset as LongArray)

        val (bin1_id, bin2_id, count) = getGroup(groupId, "pixels")
        pixels = Pixels(bin1_id as MutableList<Long>, bin2_id as MutableList<Long>, count as MutableList<Int>)
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

    return when (datasetName) {
        "weight" -> {
            (dataRead as DoubleArray).toMutableList()
        }
        "bin1_id", "bin2_id" -> {
            (dataRead as LongArray).toMutableList()
        }
        "count", "chrom" -> {
            (dataRead as IntArray).toMutableList()
        }
        else -> dataRead
    }

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
    cool.bins.weight.subList(startChrom, endChrom).reverse()

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
        cool.bins.weight.toDoubleArray()
    )

    H5.H5Dclose(datasetId)
    H5.H5Gclose(bGroupId)
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

        val end = cool.indexes.bin1_offset[i + 1].toInt()

        val indices = breakIndex until end

        val bin2Array =
            cool.pixels.bin2_id.slice(indices)
        val countArray =
            cool.pixels.count.slice(indices)

        m[i] = Pair(bin2Array.toMutableList(), countArray.toMutableList())
    }

    for (i in startChrom until endChrom) {
        val startIndex = cool.indexes.bin1_offset[i].toInt()
        val breakIndex = breakMap[i]

//        println("row $i")
        for (j in startIndex until breakIndex!!) {
            val bin2Id = cool.pixels.bin2_id[j].toInt()
//            println(bin2Id)

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
        if (inverseIndex == -1 && cool.pixels.bin2_id[i] in startChrom until endChrom) {
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
    val w = 83

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
        cool.pixels.bin1_id.toLongArray()
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
        cool.pixels.bin2_id.toLongArray()
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
        cool.pixels.count.toIntArray()
    )

    H5.H5Dclose(datasetId)

    if (pGroupId >= 0) H5.H5Gclose(pGroupId)

    recalculateIndex(cool.pixels.bin1_id, cool, groupId)
}

fun recalculateIndex(newBin1: MutableList<Long>, cool: Cool, groupId: Int) {
    var prevIndex = 0

    // TODO: 25.11.2020 find built-in functions
    newBin1.forEachIndexed { index, l ->
        for (i in prevIndex..l.toInt()) {
            cool.indexes.bin1_offset[i] = index.toLong()
        }
        prevIndex = l.toInt() + 1
    }

    val groupId2 = H5.H5Gopen(groupId, "indexes", HDF5Constants.H5P_DEFAULT)

    var datasetId = H5.H5Dopen(groupId2, "bin1_offset", HDF5Constants.H5P_DEFAULT)
    var tid = H5.H5Dget_type(datasetId)
    var dspace = H5.H5Dget_space(datasetId)
    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.indexes.bin1_offset
    )
    H5.H5Dclose(datasetId)

    cool.bins.chrom.forEachIndexed { index, l ->
        cool.indexes.chrom_offset[l + 1] = index.toLong() + 1
    }

    datasetId = H5.H5Dopen(groupId2, "chrom_offset", HDF5Constants.H5P_DEFAULT)
    tid = H5.H5Dget_type(datasetId)
    dspace = H5.H5Dget_space(datasetId)
    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.indexes.chrom_offset
    )
    H5.H5Dclose(datasetId)

    if (groupId2 >= 0) H5.H5Gclose(groupId2)
//    println(cool.indexes.bin1_offset.joinToString(separator = "\n"))
}

fun moveWeight(cool: Cool, startChrom: Int, endChrom: Int, destination: Int, groupId: Int) {
    var buf1 = cool.bins.weight.slice(startChrom until endChrom)
    cool.bins.weight.subList(startChrom, endChrom).clear()
    cool.bins.weight.addAll(destination, buf1)

    val bGroupId = H5.H5Gopen(groupId, "bins", HDF5Constants.H5P_DEFAULT)

    var datasetId = H5.H5Dopen(bGroupId, "weight", HDF5Constants.H5P_DEFAULT)
    var tid = H5.H5Dget_type(datasetId)
    var dspace = H5.H5Dget_space(datasetId)

    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.bins.weight.toDoubleArray()
    )

    H5.H5Dclose(datasetId)

    val buf2 = cool.bins.chrom.slice(destination until endChrom)
    cool.bins.chrom.subList(destination, endChrom).clear()

    val movedBuf = arrayListOf<Int>()

    var m = buf2.groupingBy { it }.eachCount().toList()

    for (j in 0 until m.last().second) {
        movedBuf.add(m.first().first)
    }

    for (i in 0 until m.size - 1) {
        movedBuf.addAll(MutableList(m[i].second) { m[i].first + 1 })
    }

    cool.bins.chrom.addAll(destination, movedBuf)

    datasetId = H5.H5Dopen(bGroupId, "chrom", HDF5Constants.H5P_DEFAULT)
    tid = H5.H5Dget_type(datasetId)
    dspace = H5.H5Dget_space(datasetId)

    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.bins.chrom.toIntArray()
    )

    H5.H5Dclose(datasetId)
    H5.H5Gclose(bGroupId)
}

fun moveChrom(cool: Cool, destination: Int, endChrom: Int, groupId: Int) {
    val cGroupId = H5.H5Gopen(groupId, "chroms", HDF5Constants.H5P_DEFAULT)

    var datasetId = H5.H5Dopen(cGroupId, "length", HDF5Constants.H5P_DEFAULT)
    var tid = H5.H5Dget_type(datasetId)
    var dspace = H5.H5Dget_space(datasetId)

    cool.chroms.length.reverse(destination, endChrom)
    cool.chroms.length.reverse(destination + 1, endChrom)

    H5.H5Dwrite(
        datasetId,
        tid,
        HDF5Constants.H5S_ALL,
        dspace,
        HDF5Constants.H5P_DEFAULT,
        cool.chroms.length
    )

    H5.H5Dclose(datasetId)

//    datasetId = H5.H5Dopen(cGroupId, "name", HDF5Constants.H5P_DEFAULT)
//    tid = H5.H5Dget_type(datasetId)
//    dspace = H5.H5Dget_space(datasetId)
//
//    cool.chroms.name.reverse(destination, endChrom)
//    cool.chroms.name.reverse(destination + 1, endChrom)
//
//    H5.H5Dwrite(
//        datasetId,
//        tid,
//        HDF5Constants.H5S_ALL,
//        dspace,
//        HDF5Constants.H5P_DEFAULT,
//        cool.chroms.name
//    )
//
//    H5.H5Dclose(datasetId)

    H5.H5Gclose(cGroupId)
}

fun move(cool: Cool, groupId: Int) {
    val w = 94

    val startChrom = cool.indexes.chrom_offset[w].toInt()
    val endChrom = cool.indexes.chrom_offset[w + 1].toInt()

    val insertChrom = 80
    val destination = cool.indexes.chrom_offset[insertChrom].toInt()

    moveWeight(cool, startChrom, endChrom, destination, groupId)
    moveChrom(cool, insertChrom, w + 1, groupId)

    val newBin1 = LongArray(cool.pixels.bin1_id.size)
    val newBin2 = LongArray(cool.pixels.bin1_id.size)
    val newCount = IntArray(cool.pixels.bin1_id.size)

    var newIdx = 0

    newIdx = moveBeforeDst(cool, startChrom, endChrom, destination, newBin1, newBin2, newCount, newIdx)
    println(newIdx)
    newIdx = moveAfterDst(cool, startChrom, endChrom, destination, newBin1, newBin2, newCount, newIdx)
    println(newIdx)

    writeStuff(groupId, newBin1, newBin2, newCount)
    recalculateIndex(newBin1.toMutableList(), cool, groupId)
}

fun <A, B, C> partial2(f: (A, B) -> C, a: A): (B) -> C {
    return { b: B -> f(a, b) }
}

fun moveBeforeDst(
    cool: Cool,
    startChrom: Int,
    endChrom: Int,
    destination: Int,
    newBin1: LongArray,
    newBin2: LongArray,
    newCount: IntArray,
    newIdx: Int
): Int {
    var newIdx = newIdx

    val beforeDst = Array(destination) { ArrayList<Pair<Long, Int>>() }
    val afterDst = Array(destination) { ArrayList<Pair<Long, Int>>() }
    val source = Array(destination) { ArrayList<Pair<Long, Int>>() }
    val afterSource = Array(destination) { ArrayList<Pair<Long, Int>>() }

    for (i in 0 until cool.indexes.bin1_offset[destination].toInt()) {
        val row = cool.pixels.bin1_id[i].toInt()

        when {
            cool.pixels.bin2_id[i] < destination -> {
                beforeDst[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
            }
            cool.pixels.bin2_id[i] in destination until startChrom -> {
                val newBin2Id = cool.pixels.bin2_id[i] + endChrom - startChrom

                afterDst[row].add(Pair(newBin2Id, cool.pixels.count[i]))
            }
            cool.pixels.bin2_id[i] in startChrom until endChrom -> {
                val newBin2Id = destination + (cool.pixels.bin2_id[i] - startChrom)

                source[row].add(Pair(newBin2Id, cool.pixels.count[i]))
            }
            else -> {
                afterSource[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
            }
        }
    }

    fun action(i: Int, it: Pair<Long, Int>) {
        newBin1[newIdx] = i.toLong()
        newBin2[newIdx] = it.first
        newCount[newIdx] = it.second

        newIdx += 1
    }

    for (i in 0 until destination) {
        beforeDst[i].forEach(partial2(::action, i))
        afterDst[i].forEach(partial2(::action, i))
        source[i].forEach(partial2(::action, i))
        afterSource[i].forEach(partial2(::action, i))
    }

    return newIdx
}

fun moveAfterDst(
    cool: Cool,
    startChrom: Int,
    endChrom: Int,
    destination: Int,
    newBin1: LongArray,
    newBin2: LongArray,
    newCount: IntArray,
    newIdx: Int
): Int {
    var newIdx = newIdx
    val length = endChrom - startChrom

    val beforeSource = Array(startChrom - destination) { ArrayList<Pair<Long, Int>>() }
    val source = Array(length) { ArrayList<Pair<Long, Int>>() }
    val afterSource = Array(startChrom - destination) { ArrayList<Pair<Long, Int>>() }

    for (i in cool.indexes.bin1_offset[destination].toInt() until cool.indexes.bin1_offset[startChrom].toInt()) {
        val row = cool.pixels.bin1_id[i].toInt() - destination

        when {
            cool.pixels.bin2_id[i] < startChrom -> {
                val newBin2Id = cool.pixels.bin2_id[i] + length

                beforeSource[row].add(Pair(newBin2Id, cool.pixels.count[i]))
            }
            cool.pixels.bin2_id[i] in startChrom until endChrom -> {
                val newBin1Id = cool.pixels.bin1_id[i] + length

                source[cool.pixels.bin2_id[i].toInt() - startChrom].add(Pair(newBin1Id, cool.pixels.count[i]))
            }
            else -> {
                afterSource[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
            }
        }
    }

    val diagonal = Array(length) { ArrayList<Pair<Long, Int>>() }

    for (i in cool.indexes.bin1_offset[startChrom].toInt() until cool.indexes.bin1_offset[endChrom].toInt()) {
        val row = (cool.pixels.bin1_id[i] - startChrom).toInt()

        if (cool.pixels.bin2_id[i] in startChrom until endChrom) {
            val newBin2id = destination + (cool.pixels.bin2_id[i] - startChrom)

            diagonal[row].add(Pair(newBin2id, cool.pixels.count[i]))
        } else {
            source[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
        }
    }

//    println(diagonal.joinToString(separator = "\n"))

    fun action1(i: Int, it: Pair<Long, Int>) {
        newBin1[newIdx] = (i + destination).toLong()
        newBin2[newIdx] = it.first
        newCount[newIdx] = it.second

        newIdx += 1
    }

    for (i in 0 until length) {
        diagonal[i].forEach(partial2(::action1, i))
        source[i].forEach(partial2(::action1, i))
    }

    fun action2(i: Int, it: Pair<Long, Int>) {
        newBin1[newIdx] = (i + destination + length).toLong()
        newBin2[newIdx] = it.first
        newCount[newIdx] = it.second

        newIdx += 1
    }

    for (i in 0 until startChrom - destination) {
        beforeSource[i].forEach(partial2(::action2, i))
        afterSource[i].forEach(partial2(::action2, i))
    }

    for (i in cool.indexes.bin1_offset[endChrom].toInt() until cool.indexes.bin1_offset.last().toInt()) {
        newBin1[i] = cool.pixels.bin1_id[i]
        newBin2[i] = cool.pixels.bin2_id[i]
        newCount[i] = cool.pixels.count[i]
    }

    return newIdx
}


fun writeStuff(groupId: Int, newBin1: LongArray, newBin2: LongArray, newCount: IntArray) {
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
        newBin1
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
        newBin2
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
        newCount
    )

    H5.H5Dclose(datasetId)

    if (pGroupId >= 0) H5.H5Gclose(pGroupId)
}

fun main() {
    val fname = "src/main/resources/2.cool"
    File("src/main/resources/mat18_100k.cool").copyTo(File(fname), overwrite = true)

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

//    inverse(cool, groupId)
    move(cool, groupId)

    // Close objects
    try {
        if (groupId >= 0) H5.H5Gclose(groupId)
        if (fileId >= 0) H5.H5Fclose(fileId)
    } catch (e: Exception) {
        e.printStackTrace()
    }

//    writeCoolToFile(cool, "src/main/resources/mat18_101k.cool", 100000)
}