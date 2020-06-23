// interp visualization (test space)
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


double func(double x, double y, double z) {

    // return x + y;

    // x^2 + y^2
    // return (x-14.5)*(x-14.5) + (y-14.5)*(y-14.5);

    // (x^2 + y^2)^(3/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 1.5);

    // x^4 + y^4 + 2 * x^2 * y^2
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5) +
      //   2 * (y-14.5)*(y-14.5)* (x-14.5)*(x-14.5);

    // x^5 + y^5
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5)*(y-14.5);

    // (x^2 + y^2)^(5/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 2.5);

    // LBFGS test
    if ((x-14.25)*(x-14.25) + (y-14.85)*(y-14.85) > 16)
        return 1;
    else
        return 0;
}


DECLARE_DIFFSCALAR_BASE();

class CylinderEnergy {
private:
    const image_t &image;
    bspline bsplineSolver;
    cylinder cylinder;
    int evalCount;

public:
    // constructor
    CylinderEnergy(const image_t &image_) : image(image_) {

        // prepare B-spline
        bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7);
    }

    // evaluate
    double operator()(const VectorXd& x, VectorXd& grad) {

        if (!cylinder.SampleCylinder(image, bsplineSolver, x(0), x(1), 3, x(2), 3)) {
            // cout << "Invalid cylinder - ";
            grad.setZero();
            return 1000;
        }

        DScalar ans = cylinder.EvaluateCylinder(image, bsplineSolver);
        grad.resize(3, 1);
        grad = ans.getGradient();
        cout << evalCount++ << " " << x(0) << " " << x(1) << " " << x(2) << " " << ans.getValue() << endl;
        cout << "grad = " << grad.transpose() << endl;
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

    image_t image;  // 30 * 30 * 10
    int sizeX, sizeY, sizeZ;
    double x, y, z;

    sizeX = 30;  // 0, 1, ..., 29
    sizeY = 30;
    sizeZ = 10;

    // generate sample grid (3D)
    double maxPixel = 0;
    for (z=0; z<sizeZ; z++) {
        
        MatrixXd layer(sizeX, sizeY);
        for (x=0; x<sizeX; x++)
            for (y=0; y<sizeY; y++) {
                layer(x, y) = func(x, y, z);
            }

        image.push_back(layer);
        if (layer.maxCoeff() > maxPixel) maxPixel = layer.maxCoeff();
    }   

    ///////////////////////////////////////////////////////////////////////////////////////////
    // main
    // prepare B-spline
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7);

    Eigen::Matrix<DScalar, Eigen::Dynamic, 2> sample;
    Eigen::Matrix<DScalar, Eigen::Dynamic, 1> res;
    sample.resize(image[0].rows() * image[0].cols(), 2);
    int count = 0;
    for (x=0; x<image[0].rows(); x++)
        for (y=0; y<image[0].cols(); y++) {
            sample(count, 0) = DScalar(x);
            sample(count, 1) = DScalar(y);
            count++;
        }
    bsplineSolver.Interp3D(sample, DScalar(3), res);

    // cout
    int i;
    cout << ">>>>>>>>>>>" << endl;
    for (i=0; i<res.rows(); i++) {
        cout << sample(i, 0).getValue() << " " << sample(i, 1).getValue() << " " << res(i, 0).getValue() << endl;
    }
}
