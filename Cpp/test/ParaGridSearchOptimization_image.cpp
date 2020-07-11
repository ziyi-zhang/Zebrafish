// Grid Search + BFGS (image) Parallel
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/autodiff.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/TiffReader.h>
#include <zebrafish/Quantile.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <math.h>
#include <CLI/CLI.hpp>
#include <LBFGS.h>
#include <string>
#include <igl/png/writePNG.h>

using namespace std;
using namespace Eigen;
using namespace zebrafish;
using namespace LBFGSpp;


DECLARE_DIFFSCALAR_BASE();



bool ValidStartingPoint(const image_t &image, const bspline &bsp, double x, double y, double z, double r) {

    static const double thres = QuantileImage(image, 0.85);

    // valid cylinder
    if (!cylinder::IsValid(bsp, x, y, z, r, 3.0)) return false;
    // membrane
    const double d1 = round(r * 1.2), d2 = round(r * 0.85);
    const MatrixXd &layer = image[round(z)];
    int count = 0;
    if (layer(x-d2, y-d2) > thres) count++;
    if (layer(x-d1, y   ) > thres) count++;
    if (layer(x-d2, y+d2) > thres) count++;
    if (layer(x   , y+d1) > thres) count++;
    if (layer(x+d2, y+d2) > thres) count++;
    if (layer(x+d1, y   ) > thres) count++;
    if (layer(x+d2, y-d2) > thres) count++;
    if (layer(x   , y-d1) > thres) count++;
    if (count < 7) return false;

    return true;
}


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
    unsigned int num_threads = 32; // std::max(1u, std::thread::hardware_concurrency());
    int lsMethod = 2;
    CLI::App command_line{"ZebraFish"};
    command_line.add_option("-i,--img", image_path, "Input TIFF image to process")->check(CLI::ExistingFile);
    command_line.add_option("-n", num_threads, "Input number of threads");
    command_line.add_option("-l", lsMethod, "Input least square solver method");

    try {
        command_line.parse(argc, argv);
    }
    catch (const CLI::ParseError &e) {
        return command_line.exit(e);
    }

    // TBB
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
    logger().info("Desired #threads = {}", num_threads);

    // read image
    image_t image;
    cout << "====================================================" << endl;
    read_tif_image(image_path, image);
    cout << "Total number of frames picked = " << image.size() << endl;

    // clip image
    for (auto it=image.begin(); it!=image.end(); it++) {
        Eigen::MatrixXd &img = *it;
        static Eigen::MatrixXd tmp;
        // tmp = img.block(377, 304, 200, 200);  // DEBUG only
        // tmp = img.block(305, 333, 638-306, 717-334);  // for 6Jan2020
        tmp = img.block(377, 304, 696-377, 684-304);  // for 6June_em1
        img = tmp;
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
    quadrature quad;
    bspline bsplineSolver(quad);
    const int bsplineDegree = 2;
    bsplineSolver.Set_leastSquareMethod(lsMethod);
    bsplineSolver.SetResolution(0.325, 0.325, 0.5);
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);
    return 0;
    /////////////////////////////////////////////////
    // prepare sampleInput & sampleOutput for grid search
    MatrixXd sampleInput, sampleOutput;
    const int Nx = bsplineSolver.Get_Nx();
    const int Ny = bsplineSolver.Get_Ny();
    const int Nz = bsplineSolver.Get_Nz();
    vector<double> rArray = {3, 4, 5, 6};
    const int Nr = rArray.size();
    const double gapX = 1.0;
    const double gapY = 1.0;
    const double gapZ = 1.0;
    double xx, yy, zz, rr;
    int sampleCount = 0;
    sampleInput.resize(Nx*Ny*Nz*Nr, 4);
    for (int ix=0; ix<floor(Nx/gapX); ix++)
        for (int iy=0; iy<floor(Ny/gapY); iy++)
            for (int iz=0; iz<floor(Nz/gapZ); iz++)
                for (int ir=0; ir<Nr; ir++) {

                    xx = ix * gapX;
                    yy = iy * gapY;
                    zz = iz * gapZ;
                    rr = rArray[ir];
                    if (!ValidStartingPoint(image, bsplineSolver, xx, yy, zz, rr)) continue;

                    sampleInput(sampleCount, 0) = xx;
                    sampleInput(sampleCount, 1) = yy;
                    sampleInput(sampleCount, 2) = zz;
                    sampleInput(sampleCount, 3) = rr;
                    sampleCount++;
                }
    sampleOutput.resize(sampleCount, 1);
    logger().info("Grid search #starting points = {}", sampleCount);
    logger().info("Grid search samples prepared...");

    // Search
    logger().info("Before search");
    tbb::parallel_for( tbb::blocked_range<int>(0, sampleCount),
        [&sampleInput, &sampleOutput, &bsplineSolver](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {
                cylinder::EvaluateCylinder(bsplineSolver, sampleInput(ii, 0), sampleInput(ii, 1), sampleInput(ii, 2), sampleInput(ii, 3), 3, sampleOutput(ii));
            }
        });
    logger().info("After search");

    /////////////////////////////////////////////////
    // prepare LBFGS
    LBFGSParam<double> param;
    param.epsilon = 1e-4;
    param.max_iterations = 15;

    // prepare sampleInput & sampleOutput for Newtons
    MatrixXd sampleInput_Newton, sampleOutput_Newton;
    sampleInput_Newton.resize(sampleCount, 4);
    int sampleCountNewton = 0;
    const double energyThres = -0.05;
    for (int i=0; i<sampleCount; i++) {
        if (sampleOutput(i) > energyThres) continue;
        
        sampleInput_Newton(sampleCountNewton, 0) = sampleInput(i, 0);  // x
        sampleInput_Newton(sampleCountNewton, 1) = sampleInput(i, 1);  // y
        sampleInput_Newton(sampleCountNewton, 2) = sampleInput(i, 2);  // z
        sampleInput_Newton(sampleCountNewton, 3) = sampleInput(i, 3);  // r
        sampleCountNewton++;
    }
    sampleOutput_Newton.resize(sampleCountNewton, 6);  // x, y, z, r, energy, iteration
    logger().info("Optimization #starting points = {}", sampleCountNewton);

    // Optimization
    logger().info("Before optimization");
    tbb::parallel_for( tbb::blocked_range<int>(0, sampleCountNewton),
        //////////////////////////////////////
        // lambda function for parallel_for //
        //////////////////////////////////////
        [&sampleInput_Newton, &sampleOutput_Newton, &bsplineSolver, &param]
        (const tbb::blocked_range<int> &r) {
        
        // NOTE: LBFGSSolver is NOT thread safe. This must be 
        LBFGSSolver<double> solver(param);

        // NOTE: the "variable count" used by "Autodiff" will be stored in 
        //       thread-local memory, so this must be set for every thread
        DiffScalarBase::setVariableCount(3);

        for (int ii = r.begin(); ii != r.end(); ++ii) {    

                Eigen::VectorXd vec(3, 1);
                vec(0) = sampleInput_Newton(ii, 0);  // x
                vec(1) = sampleInput_Newton(ii, 1);  // y
                vec(2) = sampleInput_Newton(ii, 3);  // r
                double res;

                ///////////////////////////////////
                // lambda function for optimizer //
                ///////////////////////////////////
                auto func = [&sampleInput_Newton, &bsplineSolver, ii]
                (const VectorXd& x, VectorXd& grad) {

                        DScalar ans;

                        if (!cylinder::IsValid(bsplineSolver, x(0), x(1), sampleInput_Newton(ii, 2), x(2), 3)) {
                            grad.setZero();
                            return 1.0;
                        }
                        cylinder::EvaluateCylinder(bsplineSolver, DScalar(0, x(0)), DScalar(1, x(1)), sampleInput_Newton(ii, 2), DScalar(2, x(2)), 3, ans);
                        grad.resize(3, 1);
                        grad = ans.getGradient();
                        return ans.getValue();
                    };
                // NOTE: the template of "solver.minimize" does not accept a temprary variable (due to non-const argument)
                //       so we define a "func" and pass it in every time
                int it = solver.minimize(func, vec, res);
                ///////////////////////////////////

                sampleOutput_Newton(ii, 0) = vec(0);  // x
                sampleOutput_Newton(ii, 1) = vec(1);  // y
                sampleOutput_Newton(ii, 2) = sampleInput_Newton(ii, 2);  // z
                sampleOutput_Newton(ii, 3) = vec(2);  // r
                sampleOutput_Newton(ii, 4) = res;     // energy
                sampleOutput_Newton(ii, 5) = it;      // iteration
        }
    });
    logger().info("After optimization");

    // print result
    for (int i=0; i<sampleCountNewton; i++) {
        // x y z r energy iter   x_start y_start z_start r_start
        printf("%.3f %.3f %.1f %.3f %f %.0f   %.3f %.3f %.1f %.3f\n", sampleOutput_Newton(i, 0), sampleOutput_Newton(i, 1), sampleOutput_Newton(i, 2), 
               sampleOutput_Newton(i, 3), sampleOutput_Newton(i, 4), sampleOutput_Newton(i, 5), sampleInput_Newton(i, 0), sampleInput_Newton(i, 1), 
               sampleInput_Newton(i, 2), sampleInput_Newton(i, 3));
    }
}
