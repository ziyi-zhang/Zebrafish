#include <zebrafish/GUI.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/FileDialog.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Quantile.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <sstream>

namespace zebrafish {

namespace {

bool ValidGridSearchPoint(const image_t &image, const bspline &bsp, double x, double y, double z, double r) {

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


std::string Vec2Str(const Eigen::VectorXd &vec) {

    std::stringstream sstr;
    sstr << vec;
    return sstr.str();
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 3: Grid Search

void GUI::DrawStage3() {

    ImGui::Separator(); ////////////////////////

    if (ImGui::CollapsingHeader("Grid Search", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::Button("Run Grid Search")) {
            
            GridSearch();
        }
    }

    ImGui::Separator(); ////////////////////////

    if (ImGui::CollapsingHeader("Grid Search Result", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::Text("Histogram of energy");


        ImGui::PushItemWidth(zebrafishWidth / 3.0);
        ImGui::InputDouble("Grid Search Energy Threshold", &gridEnergyThres);
        ImGui::PopItemWidth();

        if (ImGui::Button("Register Promising Results")) {

            // put into effect

        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Text("Stage 3: Grid Search");
}


////////////////////////////////////////////////////////////////////////////////////////
// Grid Search


void GUI::GridSearch() {

    // prepare gridSampleInput & gridSampleOutput for grid search
    Eigen::MatrixXd gridSampleInput, gridSampleOutput;
    const int Nx = bsplineSolver.Get_Nx();
    const int Ny = bsplineSolver.Get_Ny();
    const int Nz = bsplineSolver.Get_Nz();
    const int Nr = std::round((rArrayMax_grid - rArrayMin_grid) / rArrayGap_grid);

    Eigen::VectorXd rArray;
    rArray.resize(Nr);
    for (int i=0; i<Nr; i++)
        rArray(i) = rArrayMin_grid + i * rArrayGap_grid;
    if (rArray(Nr-1) > rArrayMax_grid) rArray(Nr-1) = rArrayMax_grid;
    
    double xx, yy, zz, rr;
    int sampleCount = 0;
    gridSampleInput.resize(Nx*Ny*Nz*Nr, 4);
    for (int ix=0; ix<floor(Nx/gapX_grid); ix++)
        for (int iy=0; iy<floor(Ny/gapY_grid); iy++)
            for (int iz=0; iz<floor(Nz/gapZ_grid); iz++)
                for (int ir=0; ir<Nr; ir++) {

                    xx = ix * gapX_grid;
                    yy = iy * gapY_grid;
                    zz = iz * gapZ_grid;
                    rr = rArray(ir);
                    if (!ValidGridSearchPoint(img, bsplineSolver, xx, yy, zz, rr)) continue;

                    gridSampleInput(sampleCount, 0) = xx;
                    gridSampleInput(sampleCount, 1) = yy;
                    gridSampleInput(sampleCount, 2) = zz;
                    gridSampleInput(sampleCount, 3) = rr;
                    sampleCount++;
                }
    gridSampleOutput.resize(sampleCount, 1);
    logger().info("Grid search samples prepared...");
    logger().debug("gapX = {:.2f}  gapY = {:.2f}  gapZ = {:.2f}", gapX_grid, gapY_grid, gapZ_grid);
    logger().debug("Calculated rArray = {}", Vec2Str(rArray));

    // Parallel Grid Search
    logger().info(">>>>>>>>>> Before grid search >>>>>>>>>>");
    logger().info("Grid search #points = {}", sampleCount);
    tbb::parallel_for( tbb::blocked_range<int>(0, sampleCount),
        [&gridSampleInput, &gridSampleOutput, this/*.bsplineSolver*/](const tbb::blocked_range<int> &r) {

            const double cylinderHeight = 3.0;
            for (int ii = r.begin(); ii != r.end(); ++ii) {
                cylinder::EvaluateCylinder(bsplineSolver, gridSampleInput(ii, 0), gridSampleInput(ii, 1), gridSampleInput(ii, 2), gridSampleInput(ii, 3), cylinderHeight, gridSampleOutput(ii));
            }
        });
    logger().info("<<<<<<<<<< After grid search <<<<<<<<<<");

    // register results to "pointRecord"
    UpdateSampleNewton(gridSampleInput, gridSampleOutput);
}


}  // namespace zebrafish
