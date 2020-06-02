#include <zebrafish/TiffReader.h>
#include <zebrafish/common.h>
#include <zebrafish/Bspline.h>

#include <CLI/CLI.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
using namespace zebrafish;

int main(int argc, char **argv)
{
    std::string image_path = "";

    CLI::App command_line{"ZebraFish"};
    command_line.add_option("-i,--img", image_path, "Input TIFF image to process")->check(CLI::ExistingFile);

    try
    {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError &e)
    {
        return command_line.exit(e);
    }

    image_t image;
    read_tif_image(image_path, image);
    
    // clip image
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &img = *it;
        // img = img.block(305, 333, 638-306, 717-334);
        img = img.block(305, 333, 20, 20);
    }


    cout << "Total number of frames picked = " << image.size() << endl;
    // ofstream file("output.log");
    // Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
    // file << image[0].format(OctaveFmt);
    /*
    Eigen::MatrixXd frame1 = image[0];
    Eigen::MatrixXd frame2 = image[1];
    file << frame1 << std::endl;
    file << frame2;
    */

    zebrafish::bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 1, 1, 1);

    Eigen::VectorXd t(10);
    Eigen::Matrix<double, 10, 3> V;
    V << 5.0665,    8.3141,   12.5303,
   13.6190,   13.1065,   10.2206,
    4.0463,    5.8185,    7.5095,
   11.7491,    6.6380,    9.1325,
   12.1730,    5.4554,    8.0181,
   12.6869,    5.3607,    4.7597,
    4.8444,   12.6929,    6.3992,
    7.9978,    9.7970,    5.2332,
    6.5987,    9.4986,    5.8391,
   12.0007,    5.4495,    6.3995;
    bsplineSolver.Interp3D(V, t);
    cout << "Interpolation result:" << endl << V << endl;
    cout << t << endl;

    return EXIT_SUCCESS;
}
