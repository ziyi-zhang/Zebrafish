#pragma once

#include <zebrafish/Common.h>

#include <Eigen/Dense>
#include <vector>


namespace zebrafish {

bool ReadTifFirstImg(const std::string &path, const int layerPerImg, const int numChannel, image_t &img);
/// This function only read channel 1 from the first 3D image

bool ReadTif(const std::string &path, const int layerPerImg, const std::vector<bool> &channelVec, const int targetNumImg, imageData_t &imgData);

}
