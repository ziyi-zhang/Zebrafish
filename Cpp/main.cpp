// Interp TEST
// Test interpolation of millions of random points
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Logger.hpp>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/enumerable_thread_specific.h>

#include <math.h>
#include <stdlib.h>
#include <random>
#include <queue>
#include <vector>
#include <chrono>

using namespace std;
using namespace Eigen;
using namespace zebrafish;

double func(double x, double y, double z) {

    // return x + y;

    // x^2 + y^2
    return (x-14.5)*(x-14.5) + (y-14.5)*(y-14.5);

    // (x^2 + y^2)^(3/2)
    // return pow((x-14.5)*(x-14.5) + (y-14.5)*(y-14.5), 1.5);
    // return (x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5);

    // x^4 + y^4 + 2 * x^2 * y^2
    // return (x-14.5)*(x-14.5)*(x-14.5)*(x-14.5) + (y-14.5)*(y-14.5)*(y-14.5)*(y-14.5) +
    //     2 * (y-14.5)*(y-14.5)* (x-14.5)*(x-14.5);

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
    int sizeX, sizeY, sizeZ, i;
    double x, y, z;

    sizeX = 100;  // 0, 1, ..., 29
    sizeY = 100;
    sizeZ = 30;

    // TBB
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
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
    bsplineSolver.CalcControlPts_um(image, 0.6, 0.6, 0.6, bsplineDegree);

    // random
    srand(time(NULL));
    uniform_real_distribution<double> unif(0, sizeX-1);
    default_random_engine re(0);

    // moving median
    double tt, l_;
    priority_queue<double, vector<double>, greater<double> > u;
    priority_queue<double, vector<double>, less<double> > l;

    // interp test
    double err, sumerr = 0, minerr = 1.0, maxerr = 0.0;
    const int trialNum = 1e8;
    Matrix<double, Dynamic, 2> sampleInput;
    Matrix<double, Dynamic, 1> sampleOutput;
    sampleInput.resize(trialNum, 2);
    sampleOutput.resize(trialNum, 1);

    for (i = 0; i<trialNum; i++) {

        x = unif(re);
        y = unif(re);
        sampleInput(i, 0) = x;
        sampleInput(i, 1) = y;
    }

    for (z = 5; z <= 5; z++) {

        logger().info("Before Interp");

        tbb::parallel_for( tbb::blocked_range<int>(0, trialNum),
        [&](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {

                sampleOutput(ii) = bsplineSolver.Interp3D(sampleInput(ii, 0), sampleInput(ii, 1), z);
                // printf("%d ", ii);
            }
        });

        logger().info("After Interp");

        // calculate theoretical output
        for (i = 0; i<trialNum; i++) {

            err = func(sampleInput(i, 0), sampleInput(i, 1), z) - sampleOutput(i);
            err = fabs(err);

            // stats
            sumerr += err;
            if (err > maxerr) maxerr = err;
            if (err < minerr) minerr = err;

            // moving median
            if (l.empty()) {l.push(err); continue;}
            l_ = l.top();
            if (err >= l_) {u.push(err);} else {l.push(err);}
            if (u.size() > l.size()) {
                tt = u.top();
                u.pop();
                l.push(tt);
            }
            if (l.size() > u.size() + 1) {
                tt = l.top();
                l.pop();
                u.push(tt);
            }
        }
    }

    // output the first several samples
    cout << "==============================" << endl;
    cout << "First 5 sample points x, y = " << endl;
    for (i=0; i<5; i++) {
        cout << sampleInput(i, 0) << " " << sampleInput(i, 1) << endl;
    }

    // interp test report
    cout << "==============================" << endl;
    cout << "Degree = " << bsplineSolver.Get_degree() << endl;
    cout << "Interp trial number = " << trialNum << endl;
    cout << "Mean error = " << sumerr / double(trialNum) << endl;
    cout << "Median error = " << l.top() << endl;
    cout << "Min  error = " << minerr << endl;
    cout << "Max  error = " << maxerr << endl;
}
