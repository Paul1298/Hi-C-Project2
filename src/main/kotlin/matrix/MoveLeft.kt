package matrix

fun moveWeight(cool: Cool, startChrom: Int, endChrom: Int, destination: Int) {
//    val buf1 = cool.bins.weight.slice(startChrom until endChrom)
//    cool.bins.weight.subList(startChrom, endChrom).clear()
//    cool.bins.weight.addAll(destination, buf1)
    cool.bins.weight.reverse(destination, startChrom)
    cool.bins.weight.reverse(startChrom, endChrom)

    cool.bins.weight.reverse(destination, endChrom)

    val buf = cool.bins.chrom.slice(destination until endChrom)
    val movedBuf = arrayListOf<Int>()
    val m = buf.groupingBy { it }.eachCount().toList()

    movedBuf.addAll(MutableList(m.last().second) { m.first().first })

    for (i in 0 until m.size - 1) {
        movedBuf.addAll(MutableList(m[i].second) { m[i].first + 1 })
    }

    for (i in destination until endChrom) {
        cool.bins.chrom[i] = buf[i - destination]
    }
}

fun moveChrom(cool: Cool, destination: Int, endChrom: Int) {
    cool.chroms.length.reverse(destination, endChrom)
    cool.chroms.length.reverse(destination + 1, endChrom)

    cool.chroms.name.reverse(destination, endChrom)
    cool.chroms.name.reverse(destination + 1, endChrom)
}

fun moveLeft(cool: Cool) {
    val w = 95

    val startChrom = cool.indexes.chrom_offset[w].toInt()
    val endChrom = cool.indexes.chrom_offset[w + 1].toInt()

    val insertChrom = 82
    val destination = cool.indexes.chrom_offset[insertChrom].toInt()

    moveWeight(cool, startChrom, endChrom, destination)
    moveChrom(cool, insertChrom, w + 1)

    val newBin1 = LongArray(cool.nnz)
    val newBin2 = LongArray(cool.nnz)
    val newCount = IntArray(cool.nnz)

    var newIdx = 0

    newIdx = moveBeforeDst(cool, startChrom, endChrom, destination, newBin1, newBin2, newCount, newIdx)
    moveAfterDst(cool, startChrom, endChrom, destination, newBin1, newBin2, newCount, newIdx)
    cool.pixels.bin1_id = newBin1
    cool.pixels.bin2_id = newBin2
    cool.pixels.count = newCount

    recalculateIndex(cool)
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
