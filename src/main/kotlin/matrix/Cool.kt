package matrix

import hdf.hdf5lib.H5
import hdf.hdf5lib.HDF5Constants
import java.io.File
import kotlin.contracts.contract
import kotlin.reflect.full.memberProperties

data class Bins(
    val chrom: IntArray,
    val end: IntArray,
    val start: IntArray,
    var weight: DoubleArray
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
    var bin1_id: LongArray,
    var bin2_id: LongArray,
    var count: IntArray
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
            dataRead = if (datasetName.startsWith("bin") || datasetName == "chrom_offset") {
                LongArray(size)
            } else {
                IntArray(size)
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

fun getGroup(groupId: Int, groupName: String): ArrayList<Any> {
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

fun updateCoolFile(cool: Cool, groupId: Int) {
    /**
     * Update existing .cool file with new values
     */
    val groups = arrayOf(cool.bins, cool.chroms, cool.indexes, cool.pixels)
    val groupNames = arrayOf("bins", "chroms", "indexes", "pixels")

    groups.forEachIndexed { index, group ->
        val tableId = H5.H5Gopen(groupId, groupNames[index], HDF5Constants.H5P_DEFAULT)
        group::class.memberProperties.forEach { dataset ->
            var data = dataset.getter.call(group)
            val datasetId = H5.H5Dopen(tableId, dataset.name, HDF5Constants.H5P_DEFAULT)
            val tid = H5.H5Dget_type(datasetId)
            val dspace = H5.H5Dget_space(datasetId)

            when (dataset.name) {
                "weight" -> data = data as DoubleArray
                "name" -> {
//                    val dims = longArrayOf((data as Array<String>).size.toLong())
//                    val dataspaceId = H5.H5Screate_simple(rank, dims, null)
//
//                    val strLength = 10
//                    val datatype = H5.H5Tcopy(HDF5Constants.H5T_C_S1)
//                    H5.H5Tset_size(datatype, strLength)
//
//                    val datasetId = H5.H5Dcreate(
//                        gid, dataset.name, datatype, dataspaceId,
//                        HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
//                    )
//
//                    // TODO: 14.11.2020 expand with whitespaces to strLength
//                    val byteBuff = data.joinToString("").toByteArray()
//                    H5.H5Dwrite(
//                        datasetId,
//                        datatype,
//                        HDF5Constants.H5S_ALL,
//                        HDF5Constants.H5S_ALL,
//                        HDF5Constants.H5P_DEFAULT,
//                        byteBuff
//                    );
//
//                    H5.H5Sclose(dataspaceId)
//                    H5.H5Dclose(datasetId)
                }
                else -> {
                    data = if (dataset.name.startsWith("bin") || dataset.name == "chrom_offset") {
                        data as LongArray
                    } else {
                        data as IntArray
                    }
                }
            }

            if (dataset.name != "name") {
                H5.H5Dwrite(
                    datasetId,
                    tid,
                    HDF5Constants.H5S_ALL,
                    dspace,
                    HDF5Constants.H5P_DEFAULT,
                    data
                )
            }
            H5.H5Dclose(datasetId)
        }
        H5.H5Gclose(tableId)
    }
}

fun recalculateIndex(cool: Cool) {
    var prevIndex = 0

    // TODO: 25.11.2020 find built-in functions
    cool.pixels.bin1_id.forEachIndexed { index, l ->
        for (i in prevIndex..l.toInt()) {
            cool.indexes.bin1_offset[i] = index.toLong()
        }
        prevIndex = l.toInt() + 1
    }

    cool.bins.chrom.forEachIndexed { index, l ->
        cool.indexes.chrom_offset[l + 1] = index.toLong() + 1
    }
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

//    inverse(cool)
//    moveLeft(cool)
    moveRight(cool)

    updateCoolFile(cool, groupId)

    // Close objects
    try {
        if (groupId >= 0) H5.H5Gclose(groupId)
        if (fileId >= 0) H5.H5Fclose(fileId)
    } catch (e: Exception) {
        e.printStackTrace()
    }

//    writeCoolToFile(cool, "src/main/resources/mat18_101k.cool", 100000)
}