#include <zebrafish/Common.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/Logger.hpp>

#include <tinytiffreader.h>
#include <tinytiffwriter.h>
#include <iostream>
#include <string>
#include <vector>
#include <regex>

////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

namespace {
    template <typename T>
    bool read_mem_to_eigen(TinyTIFFReaderFile *tiffr, T *image, Eigen::MatrixXd &img) {
        uint32_t width = TinyTIFFReader_getWidth(tiffr);
        uint32_t height = TinyTIFFReader_getHeight(tiffr);

        logger().debug("L22 width={} height={}", width, height);
        if(!image)
            image = (T *)calloc(width * height, sizeof(T));
        bool ok = TinyTIFFReader_getSampleData(tiffr, image, 0);

        //image in (ROW-MAJOR!)
        if (ok) {
            logger().debug("L29 read_mem_to_eigen OK");
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> tmp;
            tmp = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >(image, width, height);

            // normalize img
            // double min = tmp.minCoeff();
            // double max = tmp.maxCoeff();

            // img = ((tmp.template cast<double>().array() - min) / (max - min)).rowwise().reverse();
            logger().debug("L38 tmp mapped");
            tmp.transposeInPlace();
            img = tmp.template cast<double>();
        }

        return ok;
    }
}  // anonymous namespace

////////////////////////////////////////////////////////////

bool GetDescription(const std::string &path, int &layerPerImg, int &numChannel, int &ttlFrames) {

    TinyTIFFReaderFile *tiffr = NULL;
    tiffr = TinyTIFFReader_open(path.c_str());
    logger().info("Trying to phase image description with path {}", path);

    if (tiffr) {
        // get a random description
        std::string str = TinyTIFFReader_getImageDescription(tiffr);
        // search for "slices=" & "channels=" & "frames="
        std::smatch m;
        std::regex e1("(slices=)[0-9]+");
        std::regex_search(str, m, e1);
        if (m.size() == 0) {
            // did not find "slices="
            return false;  
        } else {
            std::string temp = m[0].str();
            int num = std::stoi(temp.substr(7));
            if (num > 0 && num < 1000) {
                layerPerImg = num;
            } else
                return false;
        }

        std::regex e2("(channels=)[0-9]+");
        std::regex_search(str, m, e2);
        if (m.size() == 0) {
            // did not find "channels="
            return false;  
        } else {
            std::string temp = m[0].str();
            int num = std::stoi(temp.substr(9));
            if (num > 0 && num < 10)
                numChannel = num;
            else
                return false;
        }

        std::regex e3("(frames=)[0-9]+");
        std::regex_search(str, m, e3);
        if (m.size() == 0) {
            // did not find "frames="
            return false;  
        } else {
            std::string temp = m[0].str();
            int num = std::stoi(temp.substr(7));
            if (num > 0 && num < 1000)
                ttlFrames = num;
            else
                return false;
        }

        logger().info("Image description successfully phased: layerPerImg={} numChannel={} frames={}", layerPerImg, numChannel, ttlFrames);
        return true;
    }

    return false;
}


bool ReadTifFirstFrame(const std::string &path, const int layerPerImg, const int numChannel, image_t &img, int r0, int c0, int r1, int c1, int channelToLoad) {
/// This function only reads one channel from the first frame 3D image

    std::vector<bool> channelVec(numChannel, false);
    assert(channelToLoad < numChannel);
    channelVec[channelToLoad] = true;

    imageData_t imgData_;
    bool ok;

    ok = ReadTif(path, layerPerImg, channelVec, 1, imgData_, r0, c0, r1, c1);

    logger().info(" exiting ReadTifFirstFrame() with status {}", ok);
    if (ok) img = imgData_[0];
    return ok;
}


bool ReadTif(const std::string &path, const int layerPerImg, const std::vector<bool> &channelVec, const int targetNumImg, imageData_t &imgData, int r0, int c0, int r1, int c1) {

    bool ok = false;
    TinyTIFFReaderFile *tiffr = NULL;
    tiffr = TinyTIFFReader_open(path.c_str());

    if (!tiffr) {

        logger().error("ERROR reading (not existent, not accessible or no TIFF file)");
        std::cerr << "ERROR reading (not existent, not accessible or no TIFF file)" << std::endl;
    } else {

        logger().info(path);
        logger().info(TinyTIFFReader_getImageDescription(tiffr));

        imgData.resize(targetNumImg);
        int currentImgCount = 0;

        logger().debug("L146");

        const uint32_t ttlSlices = TinyTIFFReader_countFrames(tiffr);
        const uint32_t channelNum = channelVec.size();
        const uint32_t slicePerImg = channelNum * layerPerImg;
        const uint32_t targetSlices = slicePerImg * targetNumImg;
            // We want to read the first "targetNumImg" 3D images
            // Each 3D image has "layerPerImg" layers
            // Each layer (a 2D matrix) has "channelNum" copies with different channels

        logger().debug("L156");
        uint32_t countSlice = 0;

        image_t img;  // vector of sliceMat
        Eigen::MatrixXd sliceMat;

        uint8_t *buffer8 = nullptr;
        uint16_t *buffer16 = nullptr;
        uint32_t *buffer32 = nullptr;

        do {
            logger().debug("L167 do while start");
            // prepare to read this slice
            ok = true;
            ++countSlice;

            // Is this a desired channel?
            if (channelVec[(countSlice-1) % channelNum]) {

                logger().debug("L175 Desired channel");
                // read data to "sliceMat"
                uint16_t samples = TinyTIFFReader_getBitsPerSample(tiffr);
                if(samples == 8) {
                    ok = read_mem_to_eigen<uint8_t>(tiffr, buffer8, sliceMat);
                    // scale to 0~1 for visualization
                    sliceMat /= double((1 << 8) - 1);
                }
                else if (samples == 16) {
                    ok = read_mem_to_eigen<uint16_t>(tiffr, buffer16, sliceMat);
                    // scale to 0~1 for visualization
                    sliceMat /= double((1 << 16) - 1);
                }
                else if (samples == 32) {
                    ok = read_mem_to_eigen<uint32_t>(tiffr, buffer32, sliceMat);
                    // scale to 0~1 for visualization
                    sliceMat /= double((1ll << 32) - 1);
                } 
                else {
                    logger().error("ERROR: TinyTIFFReader_getBitsPerSample wrong format");
                    std::cerr << "ERROR: TinyTIFFReader_getBitsPerSample wrong format" << std::endl;
                    ok = false;
                    break;
                }

                logger().debug("L200");
                // push a slice to "img"
                if (r0 >=0 && c0 >= 0 && r1 >=0 && c1 >= 0)
                    img.push_back(sliceMat.block(r0, c0, r1-r0+1, c1-c0+1));
                else
                    img.push_back(sliceMat);
            }

            // This needs to be executed even if this is not a desired channel
            // if "img" is a full 3D image, push to "imgData"
            logger().debug("L210");
            if (countSlice % slicePerImg == 0) {
                imgData[currentImgCount++] = img;
                img.clear();
            }

            logger().debug("L216");
            // progress
            if (countSlice % 50 == 0) {
                logger().trace("Processed {} slices / {} slices (total {})", countSlice, targetSlices, ttlSlices);
            }

        } while (TinyTIFFReader_readNext(tiffr) && (countSlice < targetSlices));

        logger().debug("L224");
        free(buffer8);
        free(buffer16);
        free(buffer32);

        if (countSlice < targetSlices) {
        // The tiff image is very small, even smaller than the default guess

            if (img.empty()) 
                ok = false;
            else {
                logger().debug("L235 TOO SMALL");
                imgData[currentImgCount++] = img;
                img.clear();
            }
        }
    }
    
    logger().debug("L242 close");
    TinyTIFFReader_close(tiffr);

    return ok;
}


bool WriteTif(const std::string &path, const image_t &imgData, int sliceBegin, int sliceEnd) {

    bool ok = true;
    TinyTIFFFile *tiffr = NULL;
    const uint16_t samples = 8;
    assert(!imgData.empty());
    const int rows = imgData[0].rows();
    const int cols = imgData[0].cols();
    const int slices = imgData.size();

    assert(sliceBegin >= 0 && sliceBegin < slices);
    assert(sliceEnd >= 0 && sliceEnd < slices);

    tiffr = TinyTIFFWriter_open(path.c_str(), samples, cols, rows);

    if (!tiffr) {

        logger().error("ERROR writing TIFF %s", path);
        std::cerr << "ERROR writing TIFF " << path << std::endl;
        ok = false;
    } else {

        for (uint16_t slice=sliceBegin; slice<sliceEnd; slice++) {

            Eigen::MatrixXd tempDouble = (imgData[slice].array() * 255.0).transpose();
            Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> temp = tempDouble.rowwise().reverse().cast<uint8_t>();
            uint8_t* data = temp.data();
            TinyTIFFWriter_writeImage(tiffr, data);
        }
    }

    TinyTIFFWriter_close(tiffr);
    return ok;
}


}  // namespace zebrafish
