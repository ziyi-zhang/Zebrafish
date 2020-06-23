// Cylinder TEST (image)
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
        img = img.block(305, 333, 30-5, 30-5);
    }
    cout << "Each layer clipped to be " << image[0].rows() << " x " << image[0].cols() << endl;
    // normalize all layers
    double quantile = zebrafish::QuantileImage(image, pixelQuantile);
    cout << "Quantile of image with q=" << pixelQuantile << " is " << quantile << endl;
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &img = *it;
        img.array() /= quantile;
    }
    cout << "Image normalized: most pixels will have value between 0 and 1" << endl;

    ///////////////////////////////////////////////////////////////////////////////////////////
    // main
    int quadArray[] = {99, 0, 1, 2, 3, 4, 5, 6};

    // prepare B-spline
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7);

    // prepare cylinder
    cylinder cylinder;
    double xx, yy, rr;
    MatrixXd query;
    int numQuery = 6;
    query.resize(numQuery, 3);
    query << 12.3524, 11.3772, 2.80271, 
             13, 14, 3, 
             10, 15, 3, 
             8, 12, 4,
             13, 10, 4,
             14, 13, 4;
    for (int j=0; j<numQuery; j++) {

        xx = query(j, 0);
        yy = query(j, 1);
        rr = query(j, 2);
        logger().info("===============");

        double refEnergy;
        for (int i=0; i<8; i++) {

            diskQuadMethod = quadArray[i];
            cylinder.LoadQuadParas();
            if (!cylinder.SampleCylinder(image, bsplineSolver, xx, yy, 32, rr, 3))
                cerr << "Invalid cylinder 1" << endl;
            DScalar ans1 = cylinder.EvaluateCylinder(image, bsplineSolver);

            if (i == 0) refEnergy = ans1.getValue();
            logger().info("QuadMethod = {}, eval = {}, diff = {}", diskQuadMethod, ans1.getValue(), refEnergy-ans1.getValue());
        }
    }
}
