// Cylinder TEST Parallel
// Care about correctness
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <math.h>

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

int main(int argc, char* argv[]) {

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

    // TBB
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    unsigned int num_threads = atoi(argv[1]); // std::max(1u, std::thread::hardware_concurrency());
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
    logger().info("Desired #threads = {}", num_threads);

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

    const int bsplineDegree = 2;
    for (int i=0; i<7; i++) {

        // prepare B-spline
        struct quadrature quad;
        double size = sizeArray[i];
        bspline bsplineSolver(quad);
        bsplineSolver.SetResolution(0.325, 0.325, 0.5);
        bsplineSolver.CalcControlPts(image, size, size, size, bsplineDegree);

        // double xx = 14.5, yy = 14.5, rr = 5;
        // double ans;
        DScalar xx = DScalar(14.5), yy = DScalar(14.5), rr = DScalar(5);
        DScalar ans;

        if (!cylinder::IsValid(bsplineSolver, xx, yy, 4, rr, 3))
            cerr << "Invalid cylinder" << endl;
        cylinder::EvaluateCylinder(bsplineSolver, xx, yy, 4, rr, 3, ans);

        cout.precision(10);

        double theory = -25.0;
        cout << "Evaluated result: " << ans << " Error = " << ans - theory << endl;
        cout << "Degree = " << bsplineSolver.Get_degree() << endl;
        cout << "control size = " << size << endl;
        cout << "maxPixel = " << maxPixel << "  normalized res = " << ans / maxPixel << endl;
        cout << endl << flush;
    }
}