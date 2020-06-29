// Cylinder & autodiff TEST
// Standard gradient correctness test
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>
#include <math.h>

using namespace std;
using namespace Eigen;
using namespace zebrafish;

double func(double x, double y, double z) {

    // return x + y;

    // x^2 + y^2
    // return (x-14.5)*(x-14.5) + (y-14.5)*(y-14.5);

    // (x^2 + y^2)^(3/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 1.5);
    return (x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5);

    // x^4 + y^4 + 2 * x^2 * y^2
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5) +
      //   2 * (y-14.5)*(y-14.5)* (x-14.5)*(x-14.5);

    // x^5 + y^5
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5)*(y-14.5);

    // (x^2 + y^2)^(5/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 2.5);
}

DECLARE_DIFFSCALAR_BASE();

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
    // prepare B-spline
    const int bsplineDegree = 2;
    bspline bsplineSolver;
    bsplineSolver.SetResolution(0.325, 0.325, 0.5);
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

    // prepare cylinder
    cylinder cylinder;
    cylinder.UpdateBoundary(image);
    double eps = 1e-4;
    DScalar xx = DScalar(0, 13.28), yy = DScalar(1, 16.52), rr = DScalar(2, 4.33);
    DScalar ans1, ans2, ans3, ans4;

    // evaluate
    if (!cylinder.EvaluateCylinder(bsplineSolver, xx, yy, 4, rr, 3, ans1))
            cerr << "Invalid cylinder 1" << endl;
    if (!cylinder.EvaluateCylinder(bsplineSolver, xx+eps, yy, 4, rr, 3, ans2))
            cerr << "Invalid cylinder 2" << endl;
    if (!cylinder.EvaluateCylinder(bsplineSolver, xx, yy+eps, 4, rr, 3, ans3))
            cerr << "Invalid cylinder 3" << endl;
    if (!cylinder.EvaluateCylinder(bsplineSolver, xx, yy, 4, rr+eps, 3, ans4))
            cerr << "Invalid cylinder 4" << endl;

    // report result
    cout.precision(12);
    cout << "E(t) = " << ans1 << endl;
    cout << "E(t)+dx*G = " << ans1.getValue()+eps*ans1.getGradient()(0) << "    E(t+dx) = " << ans2.getValue() << endl;
    cout << "E(t)+dy*G = " << ans1.getValue()+eps*ans1.getGradient()(1) << "    E(t+dy) = " << ans3.getValue() << endl;
    cout << "E(t)+dr*G = " << ans1.getValue()+eps*ans1.getGradient()(2) << "    E(t+dr) = " << ans4.getValue() << endl;
}
