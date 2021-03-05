#include <zebrafish/TiffReader.h>

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace zebrafish;

int main(int argc, char const *argv[])
{
    const std::string path = "/Users/teseo/Downloads/C2-em3_13nov2020-3.tif";

    int layerPerImg, channelPerSlice, ttlFrames;
    bool ok = GetDescription(path, layerPerImg, channelPerSlice, ttlFrames);
    imageData_t imgData;
    imgData.resize(ttlFrames);
    std::vector<bool> channelVec(channelPerSlice, false);
    channelVec[0] = true;

    ReadTif(path, layerPerImg, channelVec, ttlFrames, imgData);

    const int num_zero = ceil(log(ttlFrames) / log(10));
    std::cout << "done reading!" << std::endl;

    for (int i = 0; i < ttlFrames; ++i)
    {
        std::ostringstream ss;
        ss << std::setw(num_zero) << std::setfill('0') << i;

        std::string out_path = path.substr(0, path.size() - 4) + "_" + ss.str() + ".tif";
        WriteTif(out_path, imgData[i], 0, layerPerImg - 1);
    }

    return 0;
}
