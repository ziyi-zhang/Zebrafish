#include <zebrafish/Common.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/Logger.hpp>

#include <tinytiffreader.h>
#include <iostream>
#include <string>
#include <vector>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

namespace {
    template <typename T>
    bool read_mem_to_eigen(TinyTIFFReaderFile *tiffr, T *image, Eigen::MatrixXd &img) {
        uint32_t width = TinyTIFFReader_getWidth(tiffr);
        uint32_t height = TinyTIFFReader_getHeight(tiffr);

        if(!image)
            image = (T *)calloc(width * height, sizeof(T));
        bool ok = TinyTIFFReader_getSampleData(tiffr, image, 0);

        //image in (ROW-MAJOR!)
        if (ok) {
            static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp;
            tmp = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(image, width, height);

            // normalize img
            // double min = tmp.minCoeff();
            // double max = tmp.maxCoeff();

            // img = ((tmp.template cast<double>().array() - min) / (max - min)).rowwise().reverse();
            tmp.transposeInPlace();
            img = tmp.template cast<double>();
        }

        return ok;
    }
}  // anonymous namespace

////////////////////////////////////////////////////////////

bool ReadTifFirstImg(const std::string &path, const int layerPerImg, const int numChannel, image_t &img) {
/// This function only read channel 1 from the first 3D image

    std::vector<bool> channelVec(numChannel, false);
    channelVec[0] = true;

    imageData_t imgData;
    bool ok;

    ok = ReadTif(path, layerPerImg, channelVec, 1, imgData);
    
    if (ok) img = imgData[0];
    return ok;
}


bool ReadTif(const std::string &path, const int layerPerImg, const std::vector<bool> &channelVec, const int targetNumImg, imageData_t &imgData) {

    bool ok = false;
    TinyTIFFReaderFile *tiffr = NULL;
    tiffr = TinyTIFFReader_open(path.c_str());

    if (!tiffr) {

        logger().error("ERROR reading (not existent, not accessible or no TIFF file)");
        std::cerr << "ERROR reading (not existent, not accessible or no TIFF file)" << std::endl;
    } else {

        logger().info(path);
        logger().info(TinyTIFFReader_getImageDescription(tiffr));

        const uint32_t ttlSlices = TinyTIFFReader_countFrames(tiffr);
        const uint32_t channelNum = channelVec.size();
        const uint32_t slicePerImg = channelNum * layerPerImg;
        const uint32_t targetSlices = slicePerImg * targetNumImg;
            // We want to read the first "targetNumImg" 3D images
            // Each 3D image has "layerPerImg" layers
            // Each layer (a 2D matrix) has "channelNum" copies with different channels

        uint32_t countSlice = 0;

        image_t img;  // vector of sliceMat
        Eigen::MatrixXd sliceMat;

        uint8_t *buffer8 = nullptr;
        uint16_t *buffer16 = nullptr;
        uint32_t *buffer32 = nullptr;

        do {
            // prepare to read this slice
            ok = true;
            ++countSlice;

            // Is this a desired channel?
            if (channelVec[(countSlice-1) % channelNum]) {

                // read data to "sliceMat"
                uint16_t samples = TinyTIFFReader_getBitsPerSample(tiffr);
                if(samples == 8) {
                    ok = read_mem_to_eigen<uint8_t>(tiffr, buffer8, sliceMat);
                }
                else if (samples == 16) {
                    ok = read_mem_to_eigen<uint16_t>(tiffr, buffer16, sliceMat);
                    // scale to unsigned char for visualization
                    sliceMat /= double((1 << 8) + 1);
                }
                else if (samples == 32) {
                    ok = read_mem_to_eigen<uint32_t>(tiffr, buffer32, sliceMat);
                    // scale to unsigned char for visualization
                    sliceMat /= double((1 << 24) + 1);
                } 
                else {
                    logger().error("ERROR: TinyTIFFReader_getBitsPerSample wrong format");
                    std::cerr << "ERROR: TinyTIFFReader_getBitsPerSample wrong format" << std::endl;
                    ok = false;
                    break;
                }
                std::string descStr = TinyTIFFReader_getImageDescription(tiffr);

                // push a slice to "img"
                img.push_back(sliceMat);
            }

            // This needs to be executed even if this is not a desired channel
            // if "img" is a full 3D image, push to "imgData"
            if (countSlice % slicePerImg == 0) {
                imgData.push_back(img);
                img.clear();
            }

            // progress
            if (countSlice % 20 == 0) {
                logger().trace("Processed {} slices / {} slices (total {})", countSlice, targetSlices, ttlSlices);
            }

        } while (TinyTIFFReader_readNext(tiffr) && (countSlice < targetSlices));

        free(buffer8);
        free(buffer16);
        free(buffer32);

        if (countSlice < targetSlices) {
        // The tiff image is very small, even small than the default guess

            if (img.empty()) 
                ok = false;
            else {
                imgData.push_back(img);
                img.clear();
            }
        }
    }
    
    TinyTIFFReader_close(tiffr);

    return ok;
}

}  // namespace zebrafish
