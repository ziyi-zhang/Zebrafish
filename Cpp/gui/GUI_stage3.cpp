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

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 3: Grid Search & Optimization

void GUI::DrawStage3() {

    ImGui::Separator(); ////////////////////////

    ImGui::Text("Grid Search");
    if (ImGui::Button("Grid Search")) {
        
        GridSearch();
    }

    ImGui::Text("Histogram of energy distribution");
    ImGui::Text("Not implemented yet");

    ImGui::PushItemWidth(zebrafishWidth / 3.0);
    ImGui::InputDouble("Grid Search Energy Threshold", &gridEnergyThres);
    ImGui::PopItemWidth();

    ImGui::Separator(); ////////////////////////

    ImGui::Text("Optimization");
    if (ImGui::Button("Optimization")) {
        
        Optimization();
    }

    ImGui::Text("Stage 3: Grid Search & Optimization");
}


////////////////////////////////////////////////////////////////////////////////////////
// Grid Search


void GUI::GridSearch() {

    // prepare gridSampleInput & gridSampleOutput for grid search
    Eigen::MatrixXd gridSampleInput, gridSampleOutput;
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
                    if (!ValidGridSearchPoint(img, bsplineSolver, xx, yy, zz, rr)) continue;

                    gridSampleInput(sampleCount, 0) = xx;
                    gridSampleInput(sampleCount, 1) = yy;
                    gridSampleInput(sampleCount, 2) = zz;
                    gridSampleInput(sampleCount, 3) = rr;
                    sampleCount++;
                }
    gridSampleOutput.resize(sampleCount, 1);
    logger().info("Grid search samples prepared...");

    // Parallel Grid Search
    logger().info(">>>>>>>>>> Before grid search >>>>>>>>>>");
    logger().info("Grid search #points = {}", sampleCount);
    tbb::parallel_for( tbb::blocked_range<int>(0, sampleCount),
        [&gridSampleInput, &gridSampleOutput, this/*.bsplineSolver*/](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {
                cylinder::EvaluateCylinder(bsplineSolver, gridSampleInput(ii, 0), gridSampleInput(ii, 1), gridSampleInput(ii, 2), gridSampleInput(ii, 3), 3, gridSampleOutput(ii));
            }
        });
    logger().info("<<<<<<<<<< After grid search <<<<<<<<<<");

    // register results to "pointRecord"
    UpdateSampleNewton(gridSampleInput, gridSampleOutput);
}


////////////////////////////////////////////////////////////////////////////////////////
// Optimization


void GUI::Optimization() {

    logger().info("Not implemented");
}


void GUI::UpdateSampleNewton(const Eigen::MatrixXd &gridSampleInput, const Eigen::MatrixXd &gridSampleOutput) {

    assert(gridSampleInput.rows() == gridSampleOutput.rows());
    const int N = gridSampleOutput.rows();
    int M = 0;
    int i, count = 0;

    for (i=0; i<N; i++)
        if (gridSampleOutput(i) < gridEnergyThres) M++;  // count to reserve space

    pointRecord.num = M;
    pointRecord.alive.resize(M, Eigen::NoChange);
    pointRecord.grid_search.resize(M, Eigen::NoChange);
    pointRecord.optimization.resize(M, Eigen::NoChange);

    for (i=0; i<N; i++) {
        if (gridSampleOutput(i) > gridEnergyThres) continue;

        // alive
        pointRecord.alive(count) = true;

        // grid search
        pointRecord.grid_search(count, 0) = gridSampleInput(i, 0);
        pointRecord.grid_search(count, 1) = gridSampleInput(i, 1);
        pointRecord.grid_search(count, 2) = gridSampleInput(i, 2);
        pointRecord.grid_search(count, 3) = gridSampleInput(i, 3);
        pointRecord.grid_search(count, 4) = gridSampleOutput(i);

        // optimization
        pointRecord.optimization(count, 0) = 0;
        pointRecord.optimization(count, 1) = 0;
        pointRecord.optimization(count, 2) = 0;
        pointRecord.optimization(count, 3) = 0;
        pointRecord.optimization(count, 4) = 0;
        pointRecord.optimization(count, 5) = 0;

        count++;
    }
}

}  // namespace zebrafish
