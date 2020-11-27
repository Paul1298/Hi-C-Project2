package matrix

import hdf.hdf5lib.H5
import hdf.hdf5lib.HDF5Constants

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