#include <zebrafish/GUI.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

namespace zebrafish {

namespace {

bool ValidStartingPoint(const image_t &image, const bspline &bsp, double x, double y, double z, double r) {

    static const double thres = QuantileImage(image, 0.85);

    // valid cylinder
    if (!cylinder::IsValid(bsp, x, y, z, r, 3.0)) return false;
    // membrane
    const double d1 = round(r * 1.2), d2 = round(r * 0.85);
    const Eigen::MatrixXd &layer = image[round(z)];
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

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 1: Grid Search & Optimization

void GUI::DrawStage1() {

    if (ImGui::Button("Compute Bspline")) {

        logger().info("Computing Bspine for the first frame");
        const int bsplineDegree = 2;
        bsplineSolver.SetResolution(resolutionX, resolutionY, resolutionZ);
        bsplineSolver.CalcControlPts(img, 0.7, 0.7, 0.7, bsplineDegree);
    }

    ImGui::Separator();

    ImGui::Text("Grid Search");
    if (ImGui::Button("Grid Search")) {
        
        GridSearch();
    }

    ImGui::Separator();

    ImGui::Text("Optimization");
    if (ImGui::Button("Optimization")) {
        
        Optimization();
    }

    ImGui::Text("Stage 1: Grid Search & Optimization");
}


////////////////////////////////////////////////////////////////////////////////////////
// Grid Search


void GUI::GridSearch() {

    // prepare gridSampleInput & gridSampleOutput for grid search
    const int Nx = bsplineSolver.Get_Nx();
    const int Ny = bsplineSolver.Get_Ny();
    const int Nz = bsplineSolver.Get_Nz();
    std::vector<double> rArray = {3, 4, 5, 6, 7};
    const int Nr = rArray.size();
    const double gapX = 1.0;
    const double gapY = 1.0;
    const double gapZ = 1.0;
    double xx, yy, zz, rr;
    int sampleCount = 0;
    gridSampleInput.resize(Nx*Ny*Nz*Nr, 4);
    for (int ix=0; ix<floor(Nx/gapX); ix++)
        for (int iy=0; iy<floor(Ny/gapY); iy++)
            for (int iz=0; iz<floor(Nz/gapZ); iz++)
                for (int ir=0; ir<Nr; ir++) {

                    xx = ix * gapX;
                    yy = iy * gapY;
                    zz = iz * gapZ;
                    rr = rArray[ir];
                    if (!ValidStartingPoint(img, bsplineSolver, xx, yy, zz, rr)) continue;

                    gridSampleInput(sampleCount, 0) = xx;
                    gridSampleInput(sampleCount, 1) = yy;
                    gridSampleInput(sampleCount, 2) = zz;
                    gridSampleInput(sampleCount, 3) = rr;
                    sampleCount++;
                }
    gridSampleOutput.resize(sampleCount, 1);
    logger().info("Grid search #starting points = {}", sampleCount);
    logger().info("Grid search samples prepared...");

    // Search
    logger().info(">>>>>>>>>> Before grid search >>>>>>>>>>");
    tbb::parallel_for( tbb::blocked_range<int>(0, sampleCount),
        [this/*.bsplineSolver, .gridSampleInput, .gridSampleOutput*/](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {
                cylinder::EvaluateCylinder(bsplineSolver, gridSampleInput(ii, 0), gridSampleInput(ii, 1), gridSampleInput(ii, 2), gridSampleInput(ii, 3), 3, gridSampleOutput(ii));
            }
        });
    logger().info("<<<<<<<<<< After grid search <<<<<<<<<<");
}


////////////////////////////////////////////////////////////////////////////////////////
// Grid Search


void GUI::Optimization() {

    logger().info("Not implemented");
}

}  // namespace zebrafish
