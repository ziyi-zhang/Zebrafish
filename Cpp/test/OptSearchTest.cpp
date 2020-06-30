// Search Optimization space test (image)
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/LBFGS.h>
#include <zebrafish/autodiff.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/TiffReader.h>
#include <zebrafish/Quantile.h>
#include <math.h>
#include <CLI/CLI.hpp>
#include <string>
#include <igl/png/writePNG.h>

using namespace std;
using namespace Eigen;
using namespace zebrafish;
using namespace LBFGSpp;


DECLARE_DIFFSCALAR_BASE();

int main(int argc, char **argv) {

    // logger
    bool is_quiet = false;
    std::string log_file = "";
    int log_level = 0;
    Logger::init(!is_quiet, log_file);
    log_level = std::max(0, std::min(6, log_level));
    spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
    spdlog::flush_every(std::chrono::seconds(3));

    // read in
    std::string image_path = "";
    CLI::App command_line{"ZebraFish"};
    command_line.add_option("-i,--img", image_path, "Input TIFF image to process")->check(CLI::ExistingFile);

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
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
        img = img.block(305-5, 333-5, 37, 80);  // pt 1234
    }
    cout << "Each layer clipped to be " << image[0].rows() << " x " << image[0].cols() << endl;
    // normalize all layers
    double quantile = zebrafish::QuantileImage(image, 0.995);
    cout << "Quantile of image with q=0.995 is " << quantile << endl;
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &img = *it;
        img.array() /= quantile;
    }
    cout << "Image normalized: most pixels will have value between 0 and 1" << endl;

    ///////////////////////////////////////////////////////////////////////////////////////////
    // main
    const int xmax = image[0].rows();
    const int ymax = image[0].cols();
    const int zmax = image.size();

    // prepare B-spline
    bspline bsplineSolver;
    const int bsplineDegree = 2;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

    // prepare cylinder
    cylinder cylinder;
    cylinder.UpdateBoundary(image);

    int ir;
    double x, y, z, r, ans;
    const int rSize = 4;
    double rArray[4] = {2.0, 2.8, 3.6, 4.4};
    for (x=0; x<xmax; x++)
        for (y=0; y<ymax; y++)
            for (z=0; z<zmax; z++) 
                for (ir=0; ir<rSize; ir++) {

                    r = rArray[ir];
                    if (!cylinder.EvaluateCylinder(bsplineSolver, x, y, z, r, 3, ans))
                        continue;
                    cout << x << " " << y << " " << z << " " << r << " " << ans << endl;
                }

    logger().info("#cylidners = {}", xmax * ymax * zmax * rSize);
}
