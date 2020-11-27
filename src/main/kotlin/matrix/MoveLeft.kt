package matrix

import hdf.hdf5lib.H5
import hdf.hdf5lib.HDF5Constants

fun moveWeight(cool: Cool, startChrom: Int, endChrom: Int, destination: Int, groupId: Int) {
    val buf1 = cool.bins.weight.slice(startChrom until endChrom)
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

    val m = buf2.groupingBy { it }.eachCount().toList()

    movedBuf.addAll(MutableList(m.last().second) { m.first().first })

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

fun moveLeft(cool: Cool, groupId: Int) {
    val w = 83

    val startChrom = cool.indexes.chrom_offset[w].toInt()
    val endChrom = cool.indexes.chrom_offset[w + 1].toInt()

    val insertChrom = 82
    val destination = cool.indexes.chrom_offset[insertChrom].toInt()

    moveWeight(cool, startChrom, endChrom, destination, groupId)
    moveChrom(cool, insertChrom, w + 1, groupId)

    val newBin1 = LongArray(cool.pixels.bin1_id.size)
    val newBin2 = LongArray(cool.pixels.bin1_id.size)
    val newCount = IntArray(cool.pixels.bin1_id.size)

    var newIdx = 0

    newIdx = moveBeforeDst(cool, startChrom, endChrom, destination, newBin1, newBin2, newCount, newIdx)
//    newIdx = moveAfterDst(cool, startChrom, endChrom, destination, newBin1, newBin2, newCount, newIdx)

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
                val length = endChrom - startChrom
                val newBin2Id = cool.pixels.bin2_id[i] + length

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
        source[i].forEach(partial2(::action, i))
        afterDst[i].forEach(partial2(::action, i))
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
            val newBin2Id = destination + (cool.pixels.bin2_id[i] - startChrom)

            diagonal[row].add(Pair(newBin2Id, cool.pixels.count[i]))
        } else {
            source[row].add(Pair(cool.pixels.bin2_id[i], cool.pixels.count[i]))
        }
    }

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