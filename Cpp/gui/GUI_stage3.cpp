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
    sstr << vec.transpose();
    return sstr.str();
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 3: Grid Search

void GUI::DrawStage3() {

    // Visualize promising starting points
    if (showPromisingPoints) {

        viewer.data().point_size = 8.0f;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.87, 0.33, 0.33;

        static Eigen::MatrixXd locations;
        locations.resize(pointRecord.num, 3);
        locations.col(0) = pointRecord.grid_search.col(1);
        locations.col(1) = imgRows - pointRecord.grid_search.col(0).array();
        locations.col(2) = pointRecord.grid_search.col(2);

        viewer.data().add_points(
            locations,
            pointColor
        );
    }

    ImGui::Separator(); ////////////////////////

    if (ImGui::CollapsingHeader("Grid Search", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::TreeNode("Advanced search range")) {
            
            const float inputWidth = ImGui::GetWindowWidth() / 3.0;
            ImGui::PushItemWidth(inputWidth);
            ImGui::InputDouble("Gap X", &gapX_grid);
            ImGui::InputDouble("Gap Y", &gapY_grid);
            ImGui::InputDouble("Gap Z", &gapZ_grid);

            ImGui::InputDouble("rArray min", &rArrayMin_grid);
            ImGui::InputDouble("rArray max", &rArrayMax_grid);
            ImGui::InputDouble("rArray gap", &rArrayGap_grid);
            ImGui::PopItemWidth();
            
            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ImGui::Button("Run Grid Search")) {
            
            GridSearch();
        }
    }

    ImGui::Separator(); ////////////////////////

    if (ImGui::CollapsingHeader("Grid Search Result", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImGui::Checkbox("Show promising locations", &showPromisingPoints);

        ImGui::Text("Histogram of energy");


        const float inputWidth = ImGui::GetWindowWidth() / 3.0;
        ImGui::PushItemWidth(inputWidth);
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
    const int Nx = bsplineSolver.Get_Nx();
    const int Ny = bsplineSolver.Get_Ny();
    const int Nz = bsplineSolver.Get_Nz();
    const int Nr = std::round((rArrayMax_grid - rArrayMin_grid) / rArrayGap_grid + 1);

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
    logger().debug("Calculated rArray = [ {} ]", Vec2Str(rArray));

    // Parallel Grid Search
    logger().info(">>>>>>>>>> Before grid search >>>>>>>>>>");
    logger().info("Grid search #points = {}", sampleCount);
    tbb::parallel_for( tbb::blocked_range<int>(0, sampleCount),
        [this/*.bsplineSolver, .gridSampleInput, .gridSampleOutput*/](const tbb::blocked_range<int> &r) {

            const double cylinderHeight = 3.0;
            for (int ii = r.begin(); ii != r.end(); ++ii) {
                cylinder::EvaluateCylinder(bsplineSolver, gridSampleInput(ii, 0), gridSampleInput(ii, 1), gridSampleInput(ii, 2), gridSampleInput(ii, 3), cylinderHeight, gridSampleOutput(ii));
            }
        });
    logger().info("<<<<<<<<<< After grid search <<<<<<<<<<");

    // register results to "pointRecord"
    UpdateSampleNewton(gridSampleInput, gridSampleOutput);
}


void GUI::UpdateSampleNewton(const Eigen::MatrixXd &gridSampleInput, const Eigen::MatrixXd &gridSampleOutput) {
/// This function will clear "pointRecord" and intialize it with grid search results

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
