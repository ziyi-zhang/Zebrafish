#include <zebrafish/TiffReader.h>

#include <tinytiffreader.h>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish
{
    namespace
    {
        template <typename T>
        bool read_mem_to_eigen(TinyTIFFReaderFile *tiffr, T *image, Eigen::MatrixXd &img)
        {
            uint32_t width = TinyTIFFReader_getWidth(tiffr);
            uint32_t height = TinyTIFFReader_getHeight(tiffr);

            if(!image)
                image = (T *)calloc(width * height, sizeof(T));
            bool ok = TinyTIFFReader_getSampleData(tiffr, image, 0);

            //image in (ROW-MAJOR!)
            if (ok)
            {
                static Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> tmp;
                tmp = Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>(image, width, height);

                //normalize img
                double min = tmp.minCoeff();
                double max = tmp.maxCoeff();

                img = ((tmp.template cast<double>().array() - min) / (max - min)).rowwise().reverse();
            }

            return ok;
        }
    } // namespace

    bool read_tif_image(const std::string &path, std::vector<Eigen::MatrixXd> &img)
    {
        bool ok = false;
        TinyTIFFReaderFile *tiffr = NULL;
        tiffr = TinyTIFFReader_open(path.c_str());
        if (!tiffr)
        {
            std::cerr << "    ERROR reading (not existent, not accessible or no TIFF file)" << std::endl;
        }
        else
        {
            std::cout << TinyTIFFReader_getImageDescription(tiffr) << std::endl;
            const uint32_t frames = TinyTIFFReader_countFrames(tiffr);

            Eigen::MatrixXd tmp;

            uint32_t frame = 0;

            uint8_t *buffer8 = nullptr;
            uint16_t *buffer16 = nullptr;
            uint32_t *buffer32 = nullptr;
            do{

                uint16_t samples = TinyTIFFReader_getBitsPerSample(tiffr);

                ok = true;
                ++frame;

                if(samples == 8)
                {
                    ok = read_mem_to_eigen<uint8_t>(tiffr, buffer8, tmp);
                }
                else if (samples == 16)
                {
                    ok = read_mem_to_eigen<uint16_t>(tiffr, buffer16, tmp);
                }
                else if (samples == 32)
                {
                    ok = read_mem_to_eigen<uint32_t>(tiffr, buffer32, tmp);
                }
                else
                {
                    std::cerr << "    ERROR wrong format" << std::endl;
                    ok = false;
                    break;
                }

                if(frame%50 == 0)
                    std::cout<<"Processed "<< frame <<"/"<< frames<<std::endl;

                img.push_back(tmp);
            }
            while (TinyTIFFReader_readNext(tiffr));

            free(buffer8);
            free(buffer16);
            free(buffer32);
        }
        TinyTIFFReader_close(tiffr);
        return ok;
    }

}

