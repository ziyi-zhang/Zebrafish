// Cylinder & autodiff TEST
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>

#include <igl/png/writePNG.h>

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

    // user input
    resolutionX = 0.325;
    resolutionY = 0.325;
    resolutionZ = 0.5;

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

    // prepare B-spline
    bspline bsplineSolver;
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 1);

    // prepare cylinder
    cylinder cylinder;
    double eps = 1e-4;
    double xx = 13.28, yy = 16.52, rr = 4.33;
    if (!cylinder.SampleCylinder(image, bsplineSolver, xx, yy, 4, rr, 3))
        cerr << "Invalid cylinder 1" << endl;
    DScalar ans1 = cylinder.EvaluateCylinder(image, bsplineSolver);
    if (!cylinder.SampleCylinder(image, bsplineSolver, xx+eps, yy, 4, rr, 3))
        cerr << "Invalid cylinder 2" << endl;
    DScalar ans2 = cylinder.EvaluateCylinder(image, bsplineSolver);
    if (!cylinder.SampleCylinder(image, bsplineSolver, xx, yy+eps, 4, rr, 3))
        cerr << "Invalid cylinder 3" << endl;
    DScalar ans3 = cylinder.EvaluateCylinder(image, bsplineSolver);
    if (!cylinder.SampleCylinder(image, bsplineSolver, xx, yy, 4, rr+eps, 3))
        cerr << "Invalid cylinder 4" << endl;
    DScalar ans4 = cylinder.EvaluateCylinder(image, bsplineSolver);

    cout.precision(10);
    cout << "E(t) = " << ans1 << endl;
    cout << "E(t)+dx*G = " << ans1.getValue()+eps*ans1.getGradient()(0) << "    E(t+dx) = " << ans2.getValue() << endl;
    cout << "E(t)+dy*G = " << ans1.getValue()+eps*ans1.getGradient()(1) << "    E(t+dy) = " << ans3.getValue() << endl;
    cout << "E(t)+dr*G = " << ans1.getValue()+eps*ans1.getGradient()(2) << "    E(t+dr) = " << ans4.getValue() << endl;

    /*
    cout << "Evaluated result: " << ans.getValue() << endl;
    cout << "Gradient: " << ans << endl;
    cout << "maxPixel = " << maxPixel << "  normalized res = " << ans.getValue() / maxPixel << endl;
    */
}
