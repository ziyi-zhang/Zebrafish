#pragma once

#include <zebrafish/Common.h>

#include <Eigen/Dense>
#include <vector>


namespace zebrafish {

bool GetDescription(const std::string &path, int &layerPerImg, int &numChannel);
/// Try to get the "slices=" & "channels=" info
/// Return true if successfully found

bool ReadTifFirstImg(const std::string &path, const int layerPerImg, const int numChannel, image_t &img, 
                     int r0 = -1, int c0 = -1, int r1 = -1, int c1 = -1);
/// This function only read channel 1 from the first 3D image

bool ReadTif(const std::string &path, const int layerPerImg, const std::vector<bool> &channelVec, const int targetNumImg, imageData_t &imgData, 
             int r0 = -1, int c0 = -1, int r1 = -1, int c1 = -1);

}
