package matrix

import hdf.hdf5lib.H5
import hdf.hdf5lib.HDF5Constants
import java.io.File



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
}