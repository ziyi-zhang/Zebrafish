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

std::string Vec2Str(const Eigen::VectorXd &vec) {

    std::stringstream sstr;
    sstr << vec.transpose();
    return sstr.str();
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 3: Grid Search

void GUI::DrawStage3() {

    if (stage2to3Flag) {

        membraneThres = QuantileImage(imgData[0], 0.85, layerBegin, layerEnd);
        stage2to3Flag = false;
    }

    // Visualize promising starting points
    static int gridSearchPointSize = 7;
    if (showPromisingPoints) {

        static double gridEnergyThres_cache = -1;
        if (gridEnergyThres != gridEnergyThres_cache) {
            // Update the points to visualize when "gridEnergyThres" has changed or 
            // "Grid Search" has been re-computed
            UpdatePromisingPointLoc();
            gridEnergyThres_cache = gridEnergyThres;
        }

        viewer.data().point_size = gridSearchPointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.87, 0.33, 0.33;

        if (promisingPointLoc.rows() > 0) {
            // show promising points
            viewer.data().add_points(
                promisingPointLoc,
                pointColor
            );
        }

        ////// DEBUG ONLY //////
        Eigen::MatrixXd tempLoc;
        tempLoc.resize(3, 3);
        tempLoc << 0, 0, 1, 
                   imgCols, imgRows, layerPerImg, 
                   imgCols-1, imgRows-1, layerPerImg-1;
        Eigen::MatrixXd debugPointColor(1, 3);
        debugPointColor << 0.33, 0.83, 0.33;
        viewer.data().add_points(tempLoc, debugPointColor);
        ////// DEBUG ONLY //////
    }

    ImGui::Separator(); ////////////////////////

    if (ImGui::CollapsingHeader("Grid Search", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::TreeNode("Advanced search range")) {
            
            const float inputWidth = ImGui::GetWindowWidth() / 3.0;
            ImGui::PushItemWidth(inputWidth);
            ImGui::InputDouble("Gap X", &gapX_grid, 0.0, 0.0, "%.2f");
            ImGui::InputDouble("Gap Y", &gapY_grid, 0.0, 0.0, "%.2f");
            ImGui::InputDouble("Gap Z", &gapZ_grid, 0.0, 0.0, "%.2f");

            ImGui::InputDouble("Radius min", &rArrayMin_grid, 0.0, 0.0, "%.2f");
            ImGui::InputDouble("Radius max", &rArrayMax_grid, 0.0, 0.0, "%.2f");
            ImGui::InputDouble("Radius gap", &rArrayGap_grid, 0.0, 0.0, "%.2f");
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("The following array of radii will be searched\n[R_min R_min+R_gap R_min+2*R_gap ... R_max]");
            }

            ImGui::Checkbox("Reverse color", &reverseColor);
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("By default the markers should have darker color. If checked, the program will search for light color markers.");
            }
            ImGui::Checkbox("Skip membrane area check", &skipMembrane);
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("By default only cylinders in the membrane area will be used.");
            }
            ImGui::PopItemWidth();
            
            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ImGui::Button("Start Grid Search")) {

            GridSearch();
            // register results to "pointRecord"
            UpdateSampleNewton(gridSampleInput, gridSampleOutput);
            UpdatePromisingPointLoc();
            UpdateGridEnergyHist();
            stage3Lock = true;

            logger().debug("   <button> Start Grid Search");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Start grid search with the specified resolution. This may take some time.");
        }
    }

    ImGui::Separator(); ////////////////////////

    if (ImGui::CollapsingHeader("Grid Search Result", ImGuiTreeNodeFlags_DefaultOpen)) {

        // Histogram
        ImGui::Text("Histogram of grid search energy");

        const float width = ImGui::GetWindowWidth() * 0.75f - 2;
        ImGui::PushItemWidth(width + 2);

        ImVec2 before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", gridEnergyHist.hist.data(), gridEnergyHist.hist.size(), 0, NULL, 0, gridEnergyHist.hist.maxCoeff(), ImVec2(0, 80));
        ImVec2 after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        float ratio = ((gridEnergyThres-gridEnergyHist.minValue)/(gridEnergyHist.maxValue-gridEnergyHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        ImDrawList *drawList = ImGui::GetWindowDrawList();
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();

        ImGui::PushItemWidth(zebrafishWidth * 0.75);
        ImGui::Text("Grid search energy threshold");
        if (ImGui::SliderFloat("", &gridEnergyThres, gridEnergyHist.minValue, gridEnergyHist.maxValue, "%.3f")) {

            UpdateSampleNewton(gridSampleInput, gridSampleOutput);
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Set the grid search energy threshold. The next stage only optimizes points with smaller energy.\nA large threshold makes the detection more stable.\nA small threshold makes optimization faster.");
        }
        ImGui::PopItemWidth();

        if (ImGui::TreeNode("Advanced visualization")) {

            const float inputWidth = ImGui::GetWindowWidth() / 3.0;
            ImGui::PushItemWidth(inputWidth);
            ImGui::Checkbox("Show grid search promising location", &showPromisingPoints);
            ImGui::SliderInt("Point Size", &gridSearchPointSize, 1, 30);
            ImGui::PopItemWidth();
            
            ImGui::TreePop();
            ImGui::Separator();
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// Grid Search


void GUI::GridSearch() {

    // prepare gridSampleInput & gridSampleOutput for grid search
    const int Nx = bsplineArray[0].Get_Nx();
    const int Ny = bsplineArray[0].Get_Ny();
    const int Nz = bsplineArray[0].Get_Nz();
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
                    if (!ValidGridSearchPoint(imgData[0], bsplineArray[0], skipMembrane, membraneThres, xx, yy, zz, rr)) continue;

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
        [this/*.bsplineArray[0], .gridSampleInput, .gridSampleOutput, .cylinderHeight*/](const tbb::blocked_range<int> &r) {

            for (int ii = r.begin(); ii != r.end(); ++ii) {
                cylinder::EvaluateCylinder(bsplineArray[0], gridSampleInput(ii, 0), gridSampleInput(ii, 1), gridSampleInput(ii, 2), gridSampleInput(ii, 3), cylinderHeight, gridSampleOutput(ii), reverseColor);
            }
        });
    logger().info("<<<<<<<<<< After grid search <<<<<<<<<<");
}


bool GUI::InMembraneArea(const image_t &image, const double thres, double x, double y, double z, double r) {
// Is the cylinder in the membrane area?
// Check this by sampling 8 nearby locations. They must have light color.

    // assume cylinder height 3, peripheral sqrt(2)*r
    const double d1 = round(r * 1.2), d2 = round(r * 0.85);
    const Eigen::MatrixXd &layer = image[round(z + 1.5)];
    // if (layer(x, y) > thres) return false;  // center must be dark
    int count = 0;
    if (layer(x-d2, y-d2) > thres) count++;
    if (layer(x-d1, y   ) > thres) count++;
    if (layer(x-d2, y+d2) > thres) count++;
    if (layer(x   , y+d1) > thres) count++;
    if (layer(x+d2, y+d2) > thres) count++;
    if (layer(x+d1, y   ) > thres) count++;
    if (layer(x+d2, y-d2) > thres) count++;
    if (layer(x   , y-d1) > thres) count++;
    
    if (count < 7) 
        return false;
    else
        return true;
}


bool GUI::ValidGridSearchPoint(const image_t &image, const bspline &bsp, bool skipMembrane, double membraneThres, double x, double y, double z, double r) {

    // valid cylinder
    if (!cylinder::IsValid(bsp, x, y, z, r, cylinderHeight)) return false;
    // whether in membrane area
    if (!skipMembrane)
        if (!InMembraneArea(image, membraneThres, x, y, z, r)) return false;

    return true;
}


void GUI::UpdateSampleNewton(const Eigen::MatrixXd &gridSampleInput, const Eigen::MatrixXd &gridSampleOutput) {
/// This function will clear "pointRecord" and intialize it with grid search results

    // Note: gridSampleInput has a much larger reserved space, but the remaining should be empty
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


void GUI::UpdatePromisingPointLoc() {

    const int N = pointRecord.num;
    if (N == 0) return;

    promisingPointLoc.resize(N, 3);

    promisingPointLoc.col(0) = pointRecord.grid_search.col(1).array() + 0.5;
    promisingPointLoc.col(1) = (imgRows-0.5) - pointRecord.grid_search.col(0).array();
    promisingPointLoc.col(2) = pointRecord.grid_search.col(2);

    logger().info("   [Visualization] Grid search promising points updated: total number = {}", N);
}


void GUI::UpdateGridEnergyHist() {

    double minValue = gridSampleOutput.minCoeff();
    double maxValue = gridSampleOutput.maxCoeff();
    const int N = gridSampleOutput.rows();
    assert(N > 0);

    const double epsilon = 0.0001;  // to make sure every number lies inside
    maxValue += epsilon;
    minValue -= epsilon;
    const double gap = (maxValue - minValue) / double(histBars);

    gridEnergyHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    gridEnergyHist.minValue = minValue;
    gridEnergyHist.maxValue = maxValue;

    for (int i=0; i<N; i++) {
        gridEnergyHist.hist( std::floor((gridSampleOutput(i) - minValue)/gap) )++;
    }
}

}  // namespace zebrafish
