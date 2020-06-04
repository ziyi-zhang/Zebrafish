#include <zebrafish/TiffReader.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Cylinder.h>

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
        img = img.block(305, 333, 100, 100);
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
    bsplineSolver.CalcControlPts(image, 0.3, 0.3, 1);

    zebrafish::cylinder cylinder;
    if (!cylinder.SampleCylinder(image, bsplineSolver, 25, 35, 16, 5, 4)) {
        cerr << "Invalid cylinder to interpolate!";
    }
    cout << "Evaluated result: " << cylinder.EvaluateCylinder(image, bsplineSolver) << endl;

    return EXIT_SUCCESS;
}
