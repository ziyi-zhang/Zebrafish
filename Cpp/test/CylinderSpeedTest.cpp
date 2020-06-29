// Cylinder speed TEST
// Eval a lot of cylinders to test speed
// do not care about result
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>
#include <math.h>
#include <random>

using namespace std;
using namespace Eigen;
using namespace zebrafish;

double func(double x, double y, double z) {

    // return x + y;

    // x^2 + y^2
    return (x-14.5)*(x-14.5) + (y-14.5)*(y-14.5);

    // (x^2 + y^2)^(3/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 1.5);

    // x^4 + y^4 + 2 * x^2 * y^2
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5) +
    //      2 * (y-14.5)*(y-14.5)* (x-14.5)*(x-14.5);

    // x^5 + y^5
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5)*(y-14.5);

    // (x^2 + y^2)^(5/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 2.5);

    // (x^2 + y^2)^3
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 3.0);
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
    sizeZ = 30;

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

    // random
    srand(time(NULL));
    uniform_real_distribution<double> unif(8, sizeX-1-8);
    default_random_engine re(0);

    // eval
    const int trialNum = 2e4;
    double xx, yy, rr = 4;
    double ans;
    DScalar xxDS, yyDS, rrDS = DScalar(2, 4.0);
    DScalar ansDS;
    logger().info("Before Eval");
    for (int i=0; i<trialNum; i++) {

        xx = unif(re);
        yy = unif(re);
        xxDS = DScalar(0, xx);
        yyDS = DScalar(1, yy);

        if (!cylinder.EvaluateCylinder(bsplineSolver, xx, yy, 4, rr, 3, ans))
        //if (!cylinder.EvaluateCylinder(bsplineSolver, xxDS, yyDS, 4, rrDS, 3, ansDS))
            cerr << "Invalid cylinder" << endl;
    }
    logger().info("After Eval");

    // report
    cout << flush;
    cout << "#cylinders = " << trialNum << endl;
    cout << "#Interps = " << trialNum * 2 * 4 * 57 << endl;
}
