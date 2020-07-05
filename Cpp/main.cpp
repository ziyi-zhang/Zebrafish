// Cylinder Grid Search (image) Parallel
#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
// #include <zebrafish/LBFGS.h>
#include <zebrafish/autodiff.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/TiffReader.h>
#include <zebrafish/Quantile.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <math.h>
#include <CLI/CLI.hpp>
#include <string>
#include <igl/png/writePNG.h>

using namespace std;
using namespace Eigen;
using namespace zebrafish;
// using namespace LBFGSpp;


DECLARE_DIFFSCALAR_BASE();


bool ValidStartingPoint(const image_t &image, const bspline &bsp, double x, double y, double z, double r) {

    static const double thres = QuantileImage(image, 0.7);
    
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
    if (count < 8) return false;

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
    CLI::App command_line{"ZebraFish"};
    command_line.add_option("-i,--img", image_path, "Input TIFF image to process")->check(CLI::ExistingFile);
    command_line.add_option("-n", num_threads, "Input number of threads");

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
        // img = img.block(305, 333, 638-306, 717-334);
        tmp = img.block(305, 333, 638-306, 717-334);
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
    bsplineSolver.CalcControlPts(image, 0.7, 0.7, 0.7, bsplineDegree);

    // prepare sampleInput & sampleOutput
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
        [&](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {
                cylinder::EvaluateCylinder(bsplineSolver, sampleInput(ii, 0), sampleInput(ii, 1), sampleInput(ii, 2), sampleInput(ii, 3), 3, sampleOutput(ii));
            }
        });
    logger().info("After search");

    // store results
    const double energyThres = -0.05;
    for (int i=0; i<sampleCount; i++) {
        if (sampleOutput(i) > energyThres) continue;
        printf("%.2f %.2f %.2f %.2f %f\n", sampleInput(i, 0), sampleInput(i, 1), sampleInput(i, 2), sampleInput(i, 3), sampleOutput(i));
    }
}
