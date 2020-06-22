// BFGS Test (test space)
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/LBFGS.h>
#include <zebrafish/autodiff.h>
#include <zebrafish/Logger.hpp>
#include <math.h>

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
    if ((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5) > 9)
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

public:
    // constructor
    CylinderEnergy(const image_t &image_) : image(image_) {

        // prepare B-spline
        bsplineSolver.CalcControlPts(image, 1, 1, 1);
    }

    // evaluate
    double operator()(const VectorXd& x, VectorXd& grad) {

        cout << "=============================" << endl;
        cout << "F(" << x.transpose() << ")" << endl;

        if (!cylinder.SampleCylinder(image, bsplineSolver, x(0), x(1), 32, x(2), 3)) {
            cout << "Invalid cylinder - ";
            grad.setZero();
            return 1000;
        }

        DScalar ans = cylinder.EvaluateCylinder(image, bsplineSolver);
        grad.resize(3, 1);
        grad = ans.getGradient();
        cout << "F(" << x.transpose() << ") = " << ans.getValue() << endl;
        cout << "    Grad = " << grad.transpose() << endl;
        return ans.getValue();
    }
};


int main() {

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
    int x, y, z;

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

    // normalize it
    /*
    for (z=0; z<sizeZ; z++) {
        
        MatrixXd &layer = image[z];
        layer.array() /= maxPixel;
    }
    */

    CylinderEnergy energy(image);

    LBFGSParam<double> param;
    param.epsilon = 1e-3;
    param.max_iterations = 10;

    LBFGSSolver<double> solver(param);
    VectorXd xx(3, 1);
    xx << 13, 15, 4;
    int it;
    double res;
    it = solver.minimize(energy, xx, res);

    cout << "iterations = " << it << endl;
    cout << "Fmin = " << res << endl;
}
