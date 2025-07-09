#pragma once

#include <zebrafish/Common.h>

#include <Eigen/Dense>
#include <vector>

namespace zebrafish {

bool GetImageDescription(const std::string &path, int &layerPerImg,
                         int &numChannel, int &ttlFrames);

bool ReadImageFirstFrame(const std::string &path, const int layerPerImg,
                         const int numChannel, image_t &img, int r0 = -1,
                         int c0 = -1, int r1 = -1, int c1 = -1,
                         int channelToLoad = 0);

bool ReadImage(const std::string &path, const int layerPerImg,
               const std::vector<bool> &channelVec, const int targetNumImg,
               imageData_t &imgData, int r0 = -1, int c0 = -1, int r1 = -1,
               int c1 = -1);

bool WriteTif(const std::string &path, const image_t &imgData, int sliceBegin,
              int sliceEnd);

} // namespace zebrafish
