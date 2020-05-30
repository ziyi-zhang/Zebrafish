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

    Eigen::VectorXd t(1);
    Eigen::Matrix<double, 1, 3> V;
    V << 23.2, 45.3, 10.9;
    // zebrafish::Interp3D(image, V, t);

    zebrafish::bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 1, 1, 1);
    // cout << t(0) << endl;

    return EXIT_SUCCESS;
}
