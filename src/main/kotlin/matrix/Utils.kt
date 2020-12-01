package matrix

import hdf.hdf5lib.H5
import hdf.hdf5lib.HDF5Constants
import kotlin.reflect.full.memberProperties

fun <A, B, C> partial2(f: (A, B) -> C, a: A): (B) -> C {
    return { b: B -> f(a, b) }
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
//                println((dataRead[i] as String).length)
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
    if (groupId < 0) throw Exception()
    val tableId = H5.H5Gopen(groupId, groupName, HDF5Constants.H5P_DEFAULT)

    val members = H5.H5Gget_info(tableId)
    val nlinks = members.nlinks.toInt()
    val objNames = arrayOfNulls<String>(nlinks)
    val objTypes = IntArray(nlinks)
    val objRefs = LongArray(nlinks)

    var names_found = 0
    try {
        names_found = H5.H5Gget_obj_info_all(
            tableId, null, objNames,
            objTypes, null, objRefs, HDF5Constants.H5_INDEX_NAME
        )
    } catch (err: Throwable) {
        err.printStackTrace()
    }

    val datasetList = ArrayList<Any>()
    for (i in 0 until names_found) {
        if (objTypes[i] == HDF5Constants.H5O_TYPE_DATASET) {
            val datasetId = H5.H5Dopen(tableId, objNames[i], HDF5Constants.H5P_DEFAULT)

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

            data = when (dataset.name) {
                "weight" -> data as DoubleArray
                "name" -> {
                    @Suppress("UNCHECKED_CAST")
                    (data as Array<String>).joinToString("").toByteArray()
                }
                else -> {
                    if (dataset.name.startsWith("bin") || dataset.name == "chrom_offset") {
                        data as LongArray
                    } else {
                        data as IntArray
                    }
                }
            }

            H5.H5Dwrite(
                datasetId,
                tid,
                HDF5Constants.H5S_ALL,
                dspace,
                HDF5Constants.H5P_DEFAULT,
                data
            )
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
