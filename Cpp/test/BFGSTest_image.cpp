// BFGS Test (real image space)
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

class CylinderEnergy {
private:
    bspline bsplineSolver;
    cylinder cylinder;
    int evalCount;

public:
    // constructor
    CylinderEnergy(const image_t &image) {

        // prepare B-spline
        const int bsplineDegree = 2;
        bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

        // prepare cylinder
        cylinder.UpdateBoundary(image);
    }

    // evaluate
    double operator()(const VectorXd& x, VectorXd& grad) {

        DScalar ans;

        if (!cylinder.EvaluateCylinder(bsplineSolver, DScalar(0, x(0)), DScalar(1, x(1)), 32, DScalar(2, x(2)), 3, ans)) {
            // cout << "Invalid cylinder - ";
            grad.setZero();
            return 1.0;
        }

        grad.resize(3, 1);
        grad = ans.getGradient();
        cout << evalCount++ << " " << x(0) << " " << x(1) << " " << x(2) << " " << ans.getValue() << endl;
        cout << "Grad = " << ans.getGradient().transpose() << endl;
        return ans.getValue();
    }

    void ResetCount() {
        evalCount = 0;
    }
};


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
        img = img.block(305, 333, 25, 25);  // pt1
        // img = img.block(305, 350, 25, 25);  // pt2
        // img = img.block(305, 370, 25, 25);  // pt3
        // img = img.block(305, 389, 25, 25);  // pt4
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
    CylinderEnergy energy(image);

    LBFGSParam<double> param;
    param.epsilon = 1e-4;
    param.max_iterations = 15;

    LBFGSSolver<double> solver(param);
    VectorXd xx(3, 1);

    double delta = 5;
    int i, j, it, x, y;
    double res;
    int x0 = 12, y0 = 11;
    // pt1: 12 11
    // pt2: 11 9
    // pt3: 13 8
    // pt4: 13 8

    for (i=-delta; i<=delta; i++)
        for (j=-delta; j<=delta; j++) {

            x = x0 + i;
            y = y0 + j;
            xx(0) = x;  // starting point
            xx(1) = y;
            xx(2) = 4;

            // call optimizer
            printf(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
            energy.ResetCount();
            it = solver.minimize(energy, xx, res);

            printf("<<<<< Summary <<<<<\n");
            cout << "s_start = " << x << " " << y << endl;
            cout << "xres = " << xx.transpose() << endl;
            cout << "iterations = " << it << endl;
            logger().info("Timer");
        }

    logger().info("Timer");
}
