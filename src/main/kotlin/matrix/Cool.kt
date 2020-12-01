package matrix

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

    val nnz: Int

    init {
        val (chrom, end, start, weight) = getGroup(groupId, "bins")
        bins = Bins(chrom as IntArray, end as IntArray, start as IntArray, weight as DoubleArray)

        val (length, name) = getGroup(groupId, "chroms")
        @Suppress("UNCHECKED_CAST")
        chroms = Chroms(length as IntArray, name as Array<String>)

        val (bin1_offset, chrom_offset) = getGroup(groupId, "indexes")
        indexes = Indexes(bin1_offset as LongArray, chrom_offset as LongArray)

        val (bin1_id, bin2_id, count) = getGroup(groupId, "pixels")
        pixels = Pixels(bin1_id as LongArray, bin2_id as LongArray, count as IntArray)

        nnz = bin1_id.size
    }
}