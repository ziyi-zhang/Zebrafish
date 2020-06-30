// interp visualization (image)
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
        img = img.block(305, 333, 24, 90);  // pt1, 2, 3, 4
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
    // prepare B-spline
    bspline bsplineSolver;
    const int bsplineDegree = 2;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

    Eigen::Matrix<double, Eigen::Dynamic, 2> sample;
    Eigen::Matrix<double, Eigen::Dynamic, 1> res;
    sample.resize((image[0].rows()-1) * (image[0].cols()-1), 2);
    double x, y;
    int count = 0;
    for (x=0; x<image[0].rows()-1; x++)
        for (y=0; y<image[0].cols()-1; y++) {
            sample(count, 0) = x;
            sample(count, 1) = y;
            count++;
        }
    bsplineSolver.Interp3D(sample, 33, res);

    // cout
    int i;
    printf(">>>>>>>>>>>\n");
    for (i=0; i<res.rows(); i++) {
        cout << sample(i, 0) << " " << sample(i, 1) << " " << res(i, 0) << endl;
    }
}
