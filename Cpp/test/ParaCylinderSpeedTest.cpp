// Cylinder speed TEST parallel
// Eval a lot of cylinders to test speed
// do not care about result
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

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

    sizeX = 100;  // 0, 1, ..., 29
    sizeY = 100;
    sizeZ = 30;

    if (argc < 2) {
        cerr << "Not enough input arguments." << endl;
        return 0;
    }

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

    // prepare B-spline
    struct quadrature quad;
    const int bsplineDegree = 2;
    bspline bsplineSolver(quad);
    bsplineSolver.SetResolution(0.325, 0.325, 0.5);
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

    // random
    srand(time(NULL));
    uniform_real_distribution<double> unif(8, sizeX-1-8);
    default_random_engine re(0);

    // eval
    const int trialNum = 1e6;
    double xx, yy, rr = 4;
    double ans;
    DScalar xxDS, yyDS, rrDS = DScalar(2, 4.0);
    DScalar ansDS;
    // create input/output array
    vector<double> sampleOutput;
    sampleOutput.resize(trialNum);
    MatrixXd sampleInput;
    sampleInput.resize(trialNum, 2);
    for (int i=0; i<trialNum; i++) {
        xx = unif(re);
        yy = unif(re);
        sampleInput(i, 0) = xx;
        sampleInput(i, 1) = yy;
    }

    logger().info("Before Eval");
    tbb::parallel_for( tbb::blocked_range<int>(0, trialNum),
        [&](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {

                cylinder::EvaluateCylinder(bsplineSolver, sampleInput(ii, 0), sampleInput(ii, 1), 4, rr, 3, sampleOutput[ii]);
            }
        });
    logger().info("After Eval");

    // report
    cout << flush;
    cout << "#cylinders = " << trialNum << endl;
    cout << "#Interps = trialNum * 2 * 4 * 57 = " << trialNum * 2 * 4 * 57 << endl;
}
