package matrix

import hdf.hdf5lib.H5
import hdf.hdf5lib.HDF5Constants

fun moveWeightRight(cool: Cool, startChrom: Int, endChrom: Int, destination: Int, groupId: Int) {
    val buf1 = cool.bins.weight.slice(startChrom until endChrom)
    cool.bins.weight.addAll(destination, buf1)
    cool.bins.weight.subList(startChrom, endChrom).clear()

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

    val buf2 = cool.bins.chrom.slice(startChrom until destination)

    val movedBuf = arrayListOf<Int>()

    // TODO: 27.11.2020 check isSorted
    val m = buf2.groupingBy { it }.eachCount().toList()

    movedBuf.addAll(MutableList(m.first().second) { m.first().first })

    // TODO: 27.11.2020 maybe forEach?
    for (i in 1 until m.size) {
        movedBuf.addAll(MutableList(m[i].second) { m[i].first - 1 })
    }

    cool.bins.chrom.addAll(destination, movedBuf)
    cool.bins.chrom.subList(startChrom, destination).clear()

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

fun moveChromRight(cool: Cool, startChrom: Int, destination: Int, groupId: Int) {
    val cGroupId = H5.H5Gopen(groupId, "chroms", HDF5Constants.H5P_DEFAULT)

    var datasetId = H5.H5Dopen(cGroupId, "length", HDF5Constants.H5P_DEFAULT)
    var tid = H5.H5Dget_type(datasetId)
    var dspace = H5.H5Dget_space(datasetId)

    cool.chroms.length.reverse(startChrom, destination)
    cool.chroms.length.reverse(startChrom, destination - 1)

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

fun moveRight(cool: Cool, groupId: Int) {
    val w = 93

    val startChrom = cool.indexes.chrom_offset[w].toInt()
    val endChrom = cool.indexes.chrom_offset[w + 1].toInt()

    val insertChrom = 96
    val destination = cool.indexes.chrom_offset[insertChrom].toInt()

    moveWeightRight(cool, startChrom, endChrom, destination, groupId)
    moveChromRight(cool, w, insertChrom, groupId)

    val newBin1 = LongArray(cool.pixels.bin1_id.size)
    val newBin2 = LongArray(cool.pixels.bin1_id.size)
    val newCount = IntArray(cool.pixels.bin1_id.size)

    var newIdx = 0

    newIdx = moveBeforeSrcRight(cool, startChrom, endChrom, destination, newBin1, newBin2, newCount, newIdx)
    newIdx = moveAfterSrcRight(cool, startChrom, endChrom, destination, newBin1, newBin2, newCount, newIdx)

    writeStuff(groupId, newBin1, newBin2, newCount)
    recalculateIndex(newBin1.toMutableList(), cool, groupId)
}

fun moveBeforeSrcRight(
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

    val beforeSrc = Array(startChrom) { ArrayList<Pair<Long, Int>>() }
    val source = Array(startChrom) { ArrayList<Pair<Long, Int>>() }
    val beforeDst = Array(startChrom) { ArrayList<Pair<Long, Int>>() }
    val afterDst = Array(startChrom) { ArrayList<Pair<Long, Int>>() }

    for (i in 0 until cool.indexes.bin1_offset[startChrom].toInt()) {
        val row = cool.pixels.bin1_id[i].toInt()

        when {
            cool.pixels.bin2_id[i] < startChrom -> {
                beforeSrc[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
            }
            cool.pixels.bin2_id[i] in startChrom until endChrom -> {
                val newBin2Id = destination - (endChrom - cool.pixels.bin2_id[i])

                source[row].add(Pair(newBin2Id, cool.pixels.count[i]))
            }
            cool.pixels.bin2_id[i] in endChrom until destination -> {
                val length = endChrom - startChrom
                val newBin2Id = cool.pixels.bin2_id[i] - length

                beforeDst[row].add(Pair(newBin2Id, cool.pixels.count[i]))
            }
            else -> {
                afterDst[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
            }
        }
    }

    fun action(i: Int, it: Pair<Long, Int>) {
        newBin1[newIdx] = i.toLong()
        newBin2[newIdx] = it.first
        newCount[newIdx] = it.second

        newIdx += 1
    }


    for (i in 0 until startChrom) {
        beforeSrc[i].forEach(partial2(::action, i))
        beforeDst[i].forEach(partial2(::action, i))
        source[i].forEach(partial2(::action, i))
        afterDst[i].forEach(partial2(::action, i))
    }

    return newIdx
}

fun moveAfterSrcRight(
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

    val diagonal = Array(length) { ArrayList<Pair<Long, Int>>() }
    val beforeDst = Array(destination - endChrom) { ArrayList<Pair<Long, Int>>() }
    val afterDst = Array(length) { ArrayList<Pair<Long, Int>>() }

    for (i in cool.indexes.bin1_offset[startChrom].toInt() until cool.indexes.bin1_offset[endChrom].toInt()) {
        val row = (cool.pixels.bin1_id[i] - startChrom).toInt()

        when {
            cool.pixels.bin2_id[i] in startChrom until endChrom -> {
                val newBin2Id = destination - (endChrom - cool.pixels.bin2_id[i])

                diagonal[row].add(Pair(newBin2Id, cool.pixels.count[i]))
            }
            cool.pixels.bin2_id[i] in endChrom until destination -> {
                val newBin1Id = destination - (endChrom - cool.pixels.bin1_id[i])
                val newBin2Id = (cool.pixels.bin2_id[i] - endChrom).toInt()

                beforeDst[newBin2Id].add(Pair(newBin1Id, cool.pixels.count[i]))
            }
            else -> {
                afterDst[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
            }
        }
    }

//    println(beforeDst.joinToString(separator = "\n"))

    val freeBefore = Array(destination - endChrom) { ArrayList<Pair<Long, Int>>() }
    val freeAfter = Array(destination - endChrom) { ArrayList<Pair<Long, Int>>() }

    for (i in cool.indexes.bin1_offset[endChrom].toInt() until cool.indexes.bin1_offset[destination].toInt()) {
        val row = cool.pixels.bin1_id[i].toInt() - endChrom

        if (cool.pixels.bin2_id[i] < destination) {
            val newBin2Id = cool.pixels.bin2_id[i] - length

            freeBefore[row].add(Pair(newBin2Id, cool.pixels.count[i]))
        }
        else {
            freeAfter[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
        }
    }

    fun action1(i: Int, it: Pair<Long, Int>) {
        newBin1[newIdx] = (startChrom + i).toLong()
        newBin2[newIdx] = it.first
        newCount[newIdx] = it.second
        if (newBin1[newIdx] == 113.toLong()) println("${newBin1[newIdx]} ${newBin2[newIdx]}")

        newIdx += 1
    }

    for (i in 0 until destination - endChrom) {
        freeBefore[i].forEach(partial2(::action1, i))
        beforeDst[i].forEach(partial2(::action1, i))
        freeAfter[i].forEach(partial2(::action1, i))
    }

    fun action2(i: Int, it: Pair<Long, Int>) {
        newBin1[newIdx] = (destination - length + i).toLong()
        newBin2[newIdx] = it.first
        newCount[newIdx] = it.second

        newIdx += 1
    }

    for (i in 0 until length) {
        diagonal[i].forEach(partial2(::action2, i))
        afterDst[i].forEach(partial2(::action2, i))
    }

    for (i in cool.indexes.bin1_offset[destination].toInt() until cool.indexes.bin1_offset.last().toInt()) {
        newBin1[i] = cool.pixels.bin1_id[i]
        newBin2[i] = cool.pixels.bin2_id[i]
        newCount[i] = cool.pixels.count[i]
    }

    return newIdx
}
