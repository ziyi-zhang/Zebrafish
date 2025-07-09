#include <zebrafish/TiffReader.h>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace zebrafish;

int main(int argc, char const *argv[]) {
  const std::string path = argv[1];

  int layerPerImg, channelPerSlice, ttlFrames;
  bool ok = GetImageDescription(path, layerPerImg, channelPerSlice, ttlFrames);
  imageData_t imgData;
  imgData.resize(ttlFrames);
  std::vector<bool> channelVec(channelPerSlice, false);
  channelVec[0] = true;

  ReadImage(path, layerPerImg, channelVec, ttlFrames, imgData);

  const int num_zero = ceil(log(ttlFrames) / log(10));
  std::cout << "done reading!" << std::endl;

  for (int i = 0; i < ttlFrames; ++i) {
    std::ostringstream ss;
    ss << std::setw(num_zero) << std::setfill('0') << i;

    std::string out_path =
        path.substr(0, path.size() - 4) + "_" + ss.str() + ".tif";
    WriteTif(out_path, imgData[i], 0, layerPerImg - 1);
  }

  return 0;
}
