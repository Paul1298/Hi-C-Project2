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
    val bin1_offset: IntArray,
    val chrom_offset: IntArray
)

data class Pixels(
    val bin1_id: IntArray,
    val bin2_id: IntArray,
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
        indexes = Indexes(bin1_offset as IntArray, chrom_offset as IntArray)

        val (bin1_id, bin2_id, count) = getGroup(groupId, "pixels")
        pixels = Pixels(bin1_id as IntArray, bin2_id as IntArray, count as IntArray)
    }
}

fun getDataset(datasetId: Int): Any {
    val dspace = H5.H5Dget_space(datasetId)
    val rank = H5.H5Sget_simple_extent_ndims(dspace)

    val dims = longArrayOf(rank.toLong())
    H5.H5Sget_simple_extent_dims(dspace, dims, null)
    val size = dims[0].toInt()

    val tid = H5.H5Dget_type(datasetId)
    val class_name = H5.H5Tget_class_name(H5.H5Tget_class(tid))

    val dataRead: Any
    when (class_name) {
        "H5T_INTEGER" -> dataRead = IntArray(size)
        "H5T_FLOAT" -> dataRead = DoubleArray(size)
        "H5T_STRING" -> {
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
                dataRead[i] = String(byteBuff, i * stringLength, stringLength).trim { it <= ' ' }
            }
            return dataRead
        }
        else -> throw java.lang.Exception()
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

            datasetList.add(getDataset(datasetId))

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

    cool.javaClass.declaredFields.forEach { group ->
        val gname = group.name
        if (gname != "groupId") {
            val gid = H5.H5Gcreate(
                groupId2, gname,
                HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
            )
            println(gname)
            println(cool.bins.chrom.size)
            group.type.declaredFields.forEach { dataset ->
                val data =
                    cool::class.memberProperties.find { i ->
                        (i::class.memberProperties.find { it.name == dataset.name })?.name == gname
                    }
                print(data.toString())
            }
//                val data = groupObj::class.memberProperties.find { it.name == dataset.name }?.getter?.call(groupObj)
//                cool.bins.
//                val rank = 1
//                val dims = longArrayOf(dataset.)
//                val dspace = H5.H5Screate_simple(rank, )
//                val did = H5.H5Dcreate(
//                    gid, dataset.name.toLowerCase(),
//                    HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT, HDF5Constants.H5P_DEFAULT
//                )
            if (gid >= 0) H5.H5Gclose(gid)
        }

    }
    try {
        if (groupId2 >= 0) H5.H5Gclose(groupId2)
        if (groupId1 >= 0) H5.H5Gclose(groupId1)
        if (fileId >= 0) H5.H5Fclose(fileId)
    } catch (e: Exception) {
        e.printStackTrace()
    }
}

fun main() {
    val fname = "src/main/resources/mat18_100k.cool"
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

    // Close objects
    try {
        if (groupId >= 0) H5.H5Gclose(groupId)
        if (fileId >= 0) H5.H5Fclose(fileId)
    } catch (e: Exception) {
        e.printStackTrace()
    }

    writeCoolToFile(cool, "src/main/resources/mat18_101k.cool", 100000)
}