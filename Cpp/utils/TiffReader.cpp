#include <zebrafish/TiffReader.h>

#include <tinytiffreader.h>
#include <iostream>
#include <string>
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


    bool read_tif_image(const std::string &path, std::vector<Eigen::MatrixXd> &img) {

        bool ok = false;
        TinyTIFFReaderFile *tiffr = NULL;
        tiffr = TinyTIFFReader_open(path.c_str());
        if (!tiffr) {
            std::cerr << "    ERROR reading (not existent, not accessible or no TIFF file)" << std::endl;
        } else {
            std::cout << TinyTIFFReader_getImageDescription(tiffr) << std::endl;
            const uint32_t frames = TinyTIFFReader_countFrames(tiffr);

            uint32_t targetFrames = 82;  // only read first frame
            uint32_t countFrame = 0;
            bool channel1 = true;  // only read from channel 1

            Eigen::MatrixXd frameMat;

            uint32_t frame = 0;

            uint8_t *buffer8 = nullptr;
            uint16_t *buffer16 = nullptr;
            uint32_t *buffer32 = nullptr;

            do {
                if (channel1 && (countFrame&1)) continue;  // skip channel 2
                uint16_t samples = TinyTIFFReader_getBitsPerSample(tiffr);

                ok = true;
                ++frame;

                if(samples == 8)
                {
                    ok = read_mem_to_eigen<uint8_t>(tiffr, buffer8, frameMat);
                }
                else if (samples == 16)
                {
                    ok = read_mem_to_eigen<uint16_t>(tiffr, buffer16, frameMat);
                }
                else if (samples == 32)
                {
                    ok = read_mem_to_eigen<uint32_t>(tiffr, buffer32, frameMat);
                }
                else
                {
                    std::cerr << "    ERROR wrong format" << std::endl;
                    ok = false;
                    break;
                }
                std::string str = TinyTIFFReader_getImageDescription(tiffr);

                if (frame % 20 == 0) {
                    std::cout << "Processed " << frame << "/" << frames << std::endl;
                }

                img.push_back(frameMat);
            } while (TinyTIFFReader_readNext(tiffr) && (++countFrame < targetFrames));

            free(buffer8);
            free(buffer16);
            free(buffer32);
        }
        TinyTIFFReader_close(tiffr);
        return ok;
    }
}
