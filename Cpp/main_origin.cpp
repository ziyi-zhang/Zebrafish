#include <zebrafish/TiffReader.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/Quantile.h>

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
    cout << "====================================================" << endl;
    read_tif_image(image_path, image);
    cout << "Total number of frames picked = " << image.size() << endl;
    
    // clip image
    double maxPixel = 0, tempMaxPixel;
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &img = *it;
        // img = img.block(305, 333, 638-306, 717-334);
        img = img.block(305, 333, 100, 100);
    }
    cout << "Each layer clipped to be " << image[0].rows() << " x " << image[0].cols() << endl;
    // normalize all layers
    double quantile = zebrafish::QuantileImage(image, pixelQuantile);
    cout << "Quantile of image with q=" << pixelQuantile << " is " << quantile << endl;
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &img = *it;
        // img.array() /= quantile;
    }
    cout << "Image normalized: most pixels will have value between 0 and 1" << endl;

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
    bsplineSolver.CalcControlPts(image, 0.5, 0.5, 1);

    zebrafish::cylinder cylinder;
    if (!cylinder.SampleCylinder(image, bsplineSolver, 66, 54, 16, 5, 4)) {
        cerr << "Invalid cylinder" << endl;
    }
    cout << "Evaluated result: " << cylinder.EvaluateCylinder(image, bsplineSolver) << endl;

    return EXIT_SUCCESS;
}
