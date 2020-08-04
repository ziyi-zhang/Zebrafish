#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>

#include <igl/unproject_onto_mesh.h>
#include <vector>
#include <set>


namespace zebrafish {

// simple implementation of union-find set
// FIXME: standardize this
int query(std::vector<int> &a, int i) {

    if (a[i] != i)
        a[i] = query(a, a[i]);
    return a[i];
}


void join(std::vector<int> &a, std::vector<int> &rank, int t1, int t2) {

    int t1_ = query(a, t1);
    int t2_ = query(a, t2);
    
    if (rank[t1_] < rank[t2_])
        a[t1_] = t2_;
    else if (rank[t1_] > rank[t2_])
        a[t2_] = t1_;
    else {
        a[t1_] = t2_;
        rank[t2_]++;
    }
}

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 5: Filter & Cluster

void GUI::DrawStage5() {

    if (stage4to5Flag) {

        UpdateCylEnergyHist();
        UpdateCylRadiusHist();
        UpdateCylIterHist();
        stage4to5Flag = false;
    }

    // Visualize filtered cylinder points
    static int filterPointSize = 7;
    if (showCylFilterPoints) {

        static float cylinderEnergyThres_cache = -1;
        static float cylinderRadiusThres_cache = -1;
        static int   cylinderIterThres_cache = -1;
        if (cylinderEnergyThres != cylinderEnergyThres_cache ||
            cylinderRadiusThres != cylinderRadiusThres_cache ||
            cylinderIterThres   != cylinderIterThres_cache) {
            // Update the points to visualize when the three thresholds have changed

            CylinderFilter();
            UpdateCylPointLoc();

            cylinderEnergyThres_cache = cylinderEnergyThres;
            cylinderRadiusThres_cache = cylinderRadiusThres;
            cylinderIterThres_cache   = cylinderIterThres;
        }

        viewer.data().point_size = filterPointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.87, 0.33, 0.33;

        if (cylPointLoc.rows() > 0) {
            // show optimized points
            viewer.data().add_points(
                cylPointLoc,
                pointColor
            );
        }

        ////// DEBUG ONLY //////
        Eigen::MatrixXd tempLoc;
        tempLoc.resize(3, 3);
        tempLoc << 0, 0, 1, 
                   imgCols, imgRows, 1, 
                   imgCols-1, imgRows-1, 1;
        Eigen::MatrixXd debugPointColor(1, 3);
        debugPointColor << 0.33, 0.83, 0.33;
        viewer.data().add_points(tempLoc, debugPointColor);
        ////// DEBUG ONLY //////
    }

    // Visualize filtered cluster points
    if (showClusterFilterPoints) {

        static int clusterSizeThres_cache = -1;
        if (clusterSizeThres != clusterSizeThres_cache) {
            // Update the cluster points to visualize when the threshold have changed

            ClusterFilter();
            UpdateClusterPointLoc();

            clusterSizeThres_cache = clusterSizeThres;
        }

        viewer.data().point_size = filterPointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.99, 0.41, 0.01;

        if (clusterPointLoc.rows() > 0) {
            // show optimized cluster points
            viewer.data().add_points(
                clusterPointLoc,
                pointColor
            );
        }

        ////// DEBUG ONLY //////
        Eigen::MatrixXd tempLoc;
        tempLoc.resize(3, 3);
        tempLoc << 0, 0, 1, 
                   imgCols, imgRows, 1, 
                   imgCols-1, imgRows-1, 1;
        Eigen::MatrixXd debugPointColor(1, 3);
        debugPointColor << 0.33, 0.83, 0.33;
        viewer.data().add_points(tempLoc, debugPointColor);
        ////// DEBUG ONLY //////
    }

    // Visualize mouse picked cluster
    if (rejectActive && rejectHit) {

        Eigen::MatrixXd rejectLoc;
        static const double deltaZ = 0.3;
        int rejectHitNum = rejectHitIndex.rows();
        rejectLoc.resize(rejectHitNum, 3);
        for (int i=0; i<rejectHitNum; i++) {
            rejectLoc(i, 0) = clusterRecord.loc(rejectHitIndex(i), 1) + 0.5;
            rejectLoc(i, 1) = (imgRows-0.5) - clusterRecord.loc(rejectHitIndex(i), 0);
            rejectLoc(i, 2) = clusterRecord.loc(rejectHitIndex(i), 2) + deltaZ;
        }

        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.0, 1.0, 1.0;
        viewer.data().add_points(rejectLoc, pointColor);
    }

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::Checkbox("Show cylinder filtered locations", &showCylFilterPoints);
    ImGui::Checkbox("Show cluster filtered locations", &showClusterFilterPoints);

    ImGui::Separator(); /////////////////////////////////////////

    ImGui::PushItemWidth(zebrafishWidth / 3.0);
    ImGui::SliderInt("Point Size", &filterPointSize, 1, 30);
    ImGui::PopItemWidth();

    ImGui::Separator(); /////////////////////////////////////////
    
    // ----------------------------------------------------------

    if (ImGui::CollapsingHeader("Cylinder Filter", ImGuiTreeNodeFlags_DefaultOpen)) {

        ImVec2 before, after;
        float ratio;
        ImDrawList *drawList = ImGui::GetWindowDrawList();
        // Histogram of cylinder energy
        ImGui::Text("Histogram of energy");

        const float width = ImGui::GetWindowWidth() * 0.75f - 2;
        ImGui::PushItemWidth(width + 2);

        before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", cylEnergyHist.hist.data(), cylEnergyHist.hist.size(), 0, NULL, 0, cylEnergyHist.hist.maxCoeff(), ImVec2(0, 80));
        after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        ratio = ((cylinderEnergyThres-cylEnergyHist.minValue)/(cylEnergyHist.maxValue-cylEnergyHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();
        ImGui::SliderFloat("energy threshold", &cylinderEnergyThres, cylEnergyHist.minValue, cylEnergyHist.maxValue);

        ImGui::Separator(); /////////////////////////////////////////

        // Histogram of cylinder radius
        ImGui::Text("Histogram of cylinder radius");

        ImGui::PushItemWidth(width + 2);

        before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", cylRadiusHist.hist.data(), cylRadiusHist.hist.size(), 0, NULL, 0, cylRadiusHist.hist.maxCoeff(), ImVec2(0, 80));
        after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        ratio = ((cylinderRadiusThres-cylRadiusHist.minValue)/(cylRadiusHist.maxValue-cylRadiusHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();
        ImGui::SliderFloat("maximum radius", &cylinderRadiusThres, cylRadiusHist.minValue, cylRadiusHist.maxValue, "%.2f pixels");
        
        ImGui::Separator(); /////////////////////////////////////////

        // Histogram of cylinder radius
        ImGui::Text("Histogram of optimization iterations");

        ImGui::PushItemWidth(width + 2);

        before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", cylIterHist.hist.data(), cylIterHist.hist.size(), 0, NULL, 0, cylIterHist.hist.maxCoeff(), ImVec2(0, 80));
        after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        ratio = ((cylinderIterThres-cylIterHist.minValue)/(cylIterHist.maxValue-cylIterHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();
        ImGui::SliderInt("maximum iteration", &cylinderIterThres, 1, optimMaxIt, "%d iterations");
    }

    ImGui::Separator(); /////////////////////////////////////////

    // ----------------------------------------------------------

    if (ImGui::CollapsingHeader("Cluster Filter", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::TreeNode("Advanced cluster config")) {

            ImGui::InputFloat("Cluster dist thres", &clusterDistThres);

            ImGui::TreePop();
            ImGui::Separator();
        }

        if (ImGui::Button("Cluster")) {
            
            Cluster();
            UpdateClusterPointLoc();
            UpdateClusterSizeHist();

            propertyListType = 1;
            // update visualized points
            showCylFilterPoints = false;
            showClusterFilterPoints = true;
            logger().debug("   <button> Cluster");
        }
        
        ImGui::Separator(); /////////////////////////////////////////

        // Histogram of cluster size
        ImVec2 before, after;
        ImDrawList *drawList = ImGui::GetWindowDrawList();
        ImGui::Text("Histogram of cluster size");

        const float width = ImGui::GetWindowWidth() * 0.75f - 2;
        ImGui::PushItemWidth(width + 2);

        before = ImGui::GetCursorScreenPos();
        ImGui::PlotHistogram("", clusterSizeHist.hist.data(), clusterSizeHist.hist.size(), 0, NULL, 0, clusterSizeHist.hist.maxCoeff(), ImVec2(0, 80));
        after = ImGui::GetCursorScreenPos();
        after.y -= ImGui::GetStyle().ItemSpacing.y;

        float ratio = ((float)(clusterSizeThres-clusterSizeHist.minValue)/(float)(clusterSizeHist.maxValue-clusterSizeHist.minValue));
        ratio = std::min(1.0f, std::max(0.0f, ratio));
        drawList->PushClipRectFullScreen();
        drawList->AddLine(
            ImVec2(before.x + width * ratio, before.y), 
            ImVec2(before.x + width * ratio, after.y), 
            IM_COL32(50, 205, 50, 255), 
            2.0f
        );
        drawList->PopClipRect();
        ImGui::PopItemWidth();
        ImGui::SliderInt("Minimal cluster size", &clusterSizeThres, clusterSizeHist.minValue, clusterSizeHist.maxValue);
        
        ImGui::Separator();
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Mouse Reject", ImGuiTreeNodeFlags_DefaultOpen)) {
    
        std::vector<std::string> typeName{"Single cluster", "Area clusters"};
        ImGui::Combo("Reject mode", &rejectMode, typeName);
        ImGui::Checkbox("Mouse reject", &rejectActive);

        if (ImGui::TreeNode("Advanced mouse pick")) {

            ImGui::InputDouble("Mouse pick radius square", &mousePickDistSquareThres);
            ImGui::TreePop();
            ImGui::Separator();
        }

        ImGui::Separator();
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::TreeNode("Advanced finalize")) {
    
        ImGui::InputFloat("Finalize cluster dist thres", &finalizeClusterDistThres);
        if (ImGui::Button("Finalize cluster locations")) {
            FinalizeClusterLoc();
            propertyListType = 2;
            logger().debug("   <button> Finalize cluster locations");
        }

        ImGui::TreePop();
        ImGui::Separator();
    }

    ImGui::Text("Stage 5: Filter & Cluster");
}


////////////////////////////////////////////////////////////////////////////////////////
// Cylidner Filter

void GUI::CylinderFilter() {

    const int N = pointRecord.num;

    for (int i=0; i<N; i++) {

        pointRecord.alive(i) = 
            (pointRecord.optimization(i, 4) < cylinderEnergyThres) &&
            (pointRecord.optimization(i, 3) < cylinderRadiusThres) &&
            (pointRecord.optimization(i, 5) < cylinderIterThres);
    }
}


void GUI::UpdateCylPointLoc() {

    const int N = pointRecord.num;
    int M = 0, i, count;
    Eigen::MatrixXd tempLoc;

    for (i=0; i<N; i++)
        if (pointRecord.alive(i)) M++;
    if (M == 0) return;

    cylPointLoc.resize(M, 3);
    tempLoc.resize(M, 3);

    count = 0;
    for (i=0; i<N; i++)
        if (pointRecord.alive(i)) {
            tempLoc(count, 0) = pointRecord.optimization(i, 0);
            tempLoc(count, 1) = pointRecord.optimization(i, 1);
            tempLoc(count, 2) = pointRecord.optimization(i, 2);
            count++;
        }

    assert(count == M);

    cylPointLoc.col(0) = tempLoc.col(1).array() + 0.5;
    cylPointLoc.col(1) = (imgRows-0.5) - tempLoc.col(0).array();
    cylPointLoc.col(2) = tempLoc.col(2);

    logger().info("   [Visualization] Filtered points updated: total number = {}", M);
}


void GUI::UpdateCylEnergyHist() {
// Note: only count cylinders with energy smaller than 0
//       some cylinders failed to converge and thus have energy 1.0

    Eigen::MatrixXd energyCol = pointRecord.optimization.col(4);
    double minValue = energyCol.minCoeff();
    double maxValue = 0.0;
    const int N = pointRecord.num;
    assert(N > 0);

    const double epsilon = 0.0001;  // to make sure every number lies inside
    minValue -= epsilon;
    const double gap = (maxValue - minValue) / double(histBars);

    cylEnergyHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylEnergyHist.minValue = minValue;
    cylEnergyHist.maxValue = maxValue;

    int idx;
    for (int i=0; i<N; i++) {
        idx = std::floor((energyCol(i) - minValue)/gap);
        if (idx >= histBars) continue;  // energy > 0.0
        cylEnergyHist.hist(idx)++;
    }
}


void GUI::UpdateCylRadiusHist() {

    Eigen::MatrixXd radiusCol = pointRecord.optimization.col(3);
    double minValue = 0.0;
    double maxValue = 12.0;
    const int N = pointRecord.num;
    assert(N > 0);

    const double gap = (maxValue - minValue) / double(histBars);

    cylRadiusHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylRadiusHist.minValue = minValue;
    cylRadiusHist.maxValue = maxValue;

    int idx;
    for (int i=0; i<N; i++) {
        idx = std::floor((radiusCol(i) - minValue)/gap);
        if (idx >= histBars) idx = histBars - 1;
        if (idx < 0) idx = 0;
        cylRadiusHist.hist(idx)++;
    }
}


void GUI::UpdateCylIterHist() {

    Eigen::MatrixXd iterCol = pointRecord.optimization.col(5);
    double minValue = iterCol.minCoeff();
    double maxValue = iterCol.maxCoeff();
    const int N = pointRecord.num;
    assert(N > 0);

    const double epsilon = 0.0001;  // to make sure every number lies inside
    maxValue += epsilon;
    minValue -= epsilon;
    const double gap = (maxValue - minValue) / double(histBars);

    cylIterHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    cylIterHist.minValue = minValue;
    cylIterHist.maxValue = maxValue;

    for (int i=0; i<N; i++) {
        cylIterHist.hist( std::floor((iterCol(i) - minValue)/gap) )++;
    }
}

// ---------------------------------------------------------------
// cluster

void GUI::Cluster() {
// Note: This step may take a few seconds

    const int N = pointRecord.num;  // total number of cylinders
    int M = 0;  // number of alive cylinders
    int i, j, count;
    Eigen::MatrixXd opt_temp, x_temp, x_sorted, sortIdx, belongIdx;

    // determine the number of alive cylinders
    for (i=0; i<N; i++)
        if (pointRecord.alive(i))
            M++;
    opt_temp.resize(M, 5);
    x_temp.resize(M, 1);

    // copy xyzr-energy sub-matrix
    count = 0;
    for (i=0; i<N; i++)
        if (pointRecord.alive(i)) {
            x_temp(count) = pointRecord.optimization(i, 0);  // x
            opt_temp(count, 0) = pointRecord.optimization(i, 0);  // x
            opt_temp(count, 1) = pointRecord.optimization(i, 1);  // y
            opt_temp(count, 2) = pointRecord.optimization(i, 2);  // z
            opt_temp(count, 3) = pointRecord.optimization(i, 3);  // r
            opt_temp(count, 4) = pointRecord.optimization(i, 4);  // energy
            count++;
        }

    // sort x
    igl::sort(x_temp, 1, true, x_sorted, sortIdx);

    // prepare union-find set
    /// FIXME: standardize this
    std::vector<int> clusterSet, setRank;
    clusterSet.reserve(M);
    setRank.reserve(M);
    for (i=0; i<M; i++) {
        clusterSet[i] = i;
        setRank[i] = 0;
    }

    // cluster
    double dist_square;
    for (i=0; i<M; i++)
        for (j=i+1; j<M; j++) {

            // early stop
            if (x_sorted(i) - x_sorted(j) > clusterDistThres) continue;
            // if not same layer
            if (opt_temp(sortIdx(i), 2) != opt_temp(sortIdx(j), 2)) continue;

            // dist_square = x^2 + y^2
            dist_square = (opt_temp(sortIdx(i), 0) - opt_temp(sortIdx(j), 0))*(opt_temp(sortIdx(i), 0) - opt_temp(sortIdx(j), 0)) +
                          (opt_temp(sortIdx(i), 1) - opt_temp(sortIdx(j), 1))*(opt_temp(sortIdx(i), 1) - opt_temp(sortIdx(j), 1));

            if (dist_square < clusterDistThres*clusterDistThres)
                join(clusterSet, setRank, sortIdx(i), sortIdx(j));
        }

    // calculate belongIdx
    int numClusters;
    std::set<int> clusters;
    belongIdx.resize(M, 1);
    for (i=0; i<M; i++) {
        belongIdx(i) = query(clusterSet, i);
        clusters.insert(belongIdx(i));
    }
    numClusters = clusters.size();

    // prepare clusterRecord
    clusterRecord.num = numClusters;
    clusterRecord.alive.resize(numClusters, 1);
    clusterRecord.loc = Eigen::MatrixXd::Zero(numClusters, 4);
    clusterRecord.energy = Eigen::MatrixXd::Zero(numClusters, 1);
    clusterRecord.size = Eigen::MatrixXi::Zero(numClusters, 1);
    for (i=0; i<numClusters; i++)
        clusterRecord.alive(i) = true;

    // fill in data
    count = 0;
    for (int clusterAlias : clusters) {

        for (i=0; i<M; i++) {

            if (belongIdx(i) != clusterAlias) continue;
            clusterRecord.loc(count, 0) += opt_temp(i, 0);  // x
            clusterRecord.loc(count, 1) += opt_temp(i, 1);  // y
            clusterRecord.loc(count, 2) += opt_temp(i, 2);  // z
            clusterRecord.loc(count, 3) += opt_temp(i, 3);  // r
            clusterRecord.energy(count)    += opt_temp(i, 4);  // energy
            clusterRecord.size(count)++;  // size
        }
        // next cluster
        count++;
    }
    for (i=0; i<numClusters; i++) {

        clusterRecord.loc(i, 0) /= clusterRecord.size(i);  // x
        clusterRecord.loc(i, 1) /= clusterRecord.size(i);  // y
        clusterRecord.loc(i, 2) /= clusterRecord.size(i);  // z
        clusterRecord.loc(i, 3) /= clusterRecord.size(i);  // r
        clusterRecord.energy(i) /= clusterRecord.size(i);  // energy
    }

    logger().info("[Cluster] #(Alive Cylinders) = {} | #Cluster = {}", M, numClusters);
}


void GUI::ClusterFilter() {

    const int N = clusterRecord.num;

    for (int i=0; i<N; i++) {

        clusterRecord.alive(i) = 
            clusterRecord.size(i) >= clusterSizeThres;
    }
}


void GUI::UpdateClusterPointLoc() {

    const int N = clusterRecord.num;
    int M = 0, i, count;
    Eigen::MatrixXd tempLoc;

    for (i=0; i<N; i++)
        if (clusterRecord.alive(i)) M++;
    if (M == 0) return;

    clusterPointLoc.resize(M, 3);
    tempLoc.resize(M, 3);

    count = 0;
    for (i=0; i<N; i++)
        if (clusterRecord.alive(i)) {
            tempLoc(count, 0) = clusterRecord.loc(i, 0);
            tempLoc(count, 1) = clusterRecord.loc(i, 1);
            tempLoc(count, 2) = clusterRecord.loc(i, 2);
            count++;
        }

    assert(count == M);

    clusterPointLoc.col(0) = tempLoc.col(1).array() + 0.5;
    clusterPointLoc.col(1) = (imgRows-0.5) - tempLoc.col(0).array();
    clusterPointLoc.col(2) = tempLoc.col(2);

    logger().info("   [Visualization] Filtered clusters updated: total number = {}", M);
}


void GUI::UpdateClusterSizeHist() {

    Eigen::Matrix<int, Eigen::Dynamic, 1> sizeCol = clusterRecord.size;
    int minValue = 0;
    int maxValue = sizeCol.maxCoeff();
    const int N = clusterRecord.num;
    assert(N > 0);

    const double gap = (double)(maxValue - minValue) / double(histBars);

    clusterSizeHist.hist = Eigen::MatrixXf::Zero(histBars, 1);
    clusterSizeHist.minValue = minValue;
    clusterSizeHist.maxValue = maxValue;

    int idx;
    for (int i=0; i<N; i++) {
        idx = std::floor((double)(sizeCol(i) - minValue)/gap);
        if (idx >= histBars) idx = histBars - 1;
        if (idx < 0) idx = 0;
        clusterSizeHist.hist(idx)++;
    }
}


void GUI::FinalizeClusterLoc() {
// Do another round of cluster to get "markerRecord" from "clusterRecord"
/// FIXME: the code is very similar to "Cluster" but not identical
/// FIXME: The variables in this function is NOT properly named (copied from "Cluster")

    const int N = clusterRecord.num;  // total number of clusters
    int M = 0;  // number of alive clusters
    int i, j, count;
    Eigen::MatrixXd opt_temp, x_temp, x_sorted, sortIdx, belongIdx;

    // determine the number of alive clusters
    for (i=0; i<N; i++)
        if (clusterRecord.alive(i))
            M++;
    opt_temp.resize(M, 6);
    x_temp.resize(M, 1);

    // copy xyzr-energy-size sub-matrix
    count = 0;
    for (i=0; i<N; i++)
        if (clusterRecord.alive(i)) {
            x_temp(count) = clusterRecord.loc(i, 0);  // x
            opt_temp(count, 0) = clusterRecord.loc(i, 0);  // x
            opt_temp(count, 1) = clusterRecord.loc(i, 1);  // y
            opt_temp(count, 2) = clusterRecord.loc(i, 2);  // z
            opt_temp(count, 3) = clusterRecord.loc(i, 3);  // r
            opt_temp(count, 4) = clusterRecord.energy(i);  // energy
            opt_temp(count, 5) = clusterRecord.size(i);  // size
            count++;
        }

    // sort x
    igl::sort(x_temp, 1, true, x_sorted, sortIdx);

    // prepare union-find set
    /// FIXME: standardize this
    std::vector<int> clusterSet, setRank;
    clusterSet.reserve(M);
    setRank.reserve(M);
    for (i=0; i<M; i++) {
        clusterSet[i] = i;
        setRank[i] = 0;
    }

    // cluster
    double dist_square;
    for (i=0; i<M; i++)
        for (j=i+1; j<M; j++) {

            // early stop
            if (x_sorted(i) - x_sorted(j) > finalizeClusterDistThres) continue;
            // dont care about depth layer

            // dist_square = x^2 + y^2
            dist_square = (opt_temp(sortIdx(i), 0) - opt_temp(sortIdx(j), 0))*(opt_temp(sortIdx(i), 0) - opt_temp(sortIdx(j), 0)) +
                          (opt_temp(sortIdx(i), 1) - opt_temp(sortIdx(j), 1))*(opt_temp(sortIdx(i), 1) - opt_temp(sortIdx(j), 1));

            if (dist_square < finalizeClusterDistThres*finalizeClusterDistThres)
                join(clusterSet, setRank, sortIdx(i), sortIdx(j));
        }

    // calculate belongIdx
    int numClusters;
    std::set<int> clusters;
    belongIdx.resize(M, 1);
    for (i=0; i<M; i++) {
        belongIdx(i) = query(clusterSet, i);
        clusters.insert(belongIdx(i));
    }
    numClusters = clusters.size();

    // prepare markerRecord
    markerRecord_t markerFirstFrame;
    markerFirstFrame.num = numClusters;
    markerFirstFrame.loc = Eigen::MatrixXd::Zero(numClusters, 4);
    markerFirstFrame.energy = Eigen::MatrixXd::Zero(numClusters, 1);
    markerFirstFrame.size = Eigen::MatrixXi::Zero(numClusters, 1);

    // fill in data
    count = 0;
    for (int clusterAlias : clusters) {

        for (i=0; i<M; i++) {

            if (belongIdx(i) != clusterAlias) continue;
            int clusterSize = opt_temp(i, 5);
            markerFirstFrame.loc(count, 0) += opt_temp(i, 0) * clusterSize;  // x
            markerFirstFrame.loc(count, 1) += opt_temp(i, 1) * clusterSize;  // y
            markerFirstFrame.loc(count, 2) += opt_temp(i, 2) * clusterSize;  // z
            markerFirstFrame.loc(count, 3) += opt_temp(i, 3) * clusterSize;  // r
            markerFirstFrame.size(count)   += clusterSize;  // size
        }
        // next cluster
        count++;
    }
    for (i=0; i<numClusters; i++) {

        markerFirstFrame.loc(i, 0) /= markerFirstFrame.size(i);  // x
        markerFirstFrame.loc(i, 1) /= markerFirstFrame.size(i);  // y
        markerFirstFrame.loc(i, 2) /= markerFirstFrame.size(i);  // z
        markerFirstFrame.loc(i, 3) /= markerFirstFrame.size(i);  // r

        // calculate energy (directly)
        cylinder::EvaluateCylinder(bsplineArray[0], markerFirstFrame.loc(i, 0), markerFirstFrame.loc(i, 1), markerFirstFrame.loc(i, 2), markerFirstFrame.loc(i, 3), 3.0, markerFirstFrame.energy(i));
    }

    // push to "markerArray"
    markerArray.clear();
    markerArray.push_back(markerFirstFrame);

    logger().info("[Finalize Cluster] #(Alive clusters) = {} | #Markers = {}", M, numClusters);
}

////////////////////////////////////////////////////////////////////////////////////////
// mouse reject

void GUI::MouseSelectCluster(const Eigen::Vector2f &mouse) {
/// called by "MouseMoveCallback"

    static int fid;
    static Eigen::Vector3f bc;
    static double x, y;
    static bool hit;

    hit = igl::unproject_onto_mesh(
        Eigen::Vector2f(mouse(0), viewer.core().viewport(3)-mouse(1)), 
        viewer.core().view, 
        viewer.core().proj,
        viewer.core().viewport, 
        V, 
        F, 
        fid, 
        bc);

    if (hit) {
        // has hit

        // FIXME: why minus 0.5?
        y =                 (V(F(fid, 0), 0) * bc(0) + V(F(fid, 1), 0) * bc(1) + V(F(fid, 2), 0) * bc(2)) - 0.5;
        x = imgRows - 0.5 - (V(F(fid, 0), 1) * bc(0) + V(F(fid, 1), 1) * bc(1) + V(F(fid, 2), 1) * bc(2));

        if (rejectMode == REJECT_SINGLE) {
            // search for the nearest cluster
            double dist_square, minDistSquare;
            minDistSquare = mousePickDistSquareThres;  // reset
            for (int i=0; i<clusterRecord.num; i++) {

                if (!clusterRecord.alive(i)) continue;
                dist_square = (clusterRecord.loc(i, 0)-x)*(clusterRecord.loc(i, 0)-x) + (clusterRecord.loc(i, 1)-y)*(clusterRecord.loc(i, 1)-y);
                if (dist_square >= mousePickDistSquareThres) continue;
                if (dist_square < minDistSquare) {
                    // find a new nearest cluster

                    minDistSquare = dist_square;
                    rejectHitIndex.resize(1, 1);
                    rejectHitIndex(0) = i;
                }
            }
            // update rejectHit
            rejectHit = minDistSquare < mousePickDistSquareThres;
        } else if (rejectMode == REJECT_AREA) {
            // search for all near clusters
            Eigen::VectorXi tempIndex(clusterRecord.num, 1);
            double dist_square;
            int hitCount = 0;
            for (int i=0; i<clusterRecord.num; i++) {

                if (!clusterRecord.alive(i)) continue;
                dist_square = (clusterRecord.loc(i, 0)-x)*(clusterRecord.loc(i, 0)-x) + (clusterRecord.loc(i, 1)-y)*(clusterRecord.loc(i, 1)-y);
                if (dist_square >= mousePickDistSquareThres) continue;
                // find a cluster that is close enough

                tempIndex(hitCount) = i;
                hitCount++;
            }
            rejectHitIndex = tempIndex.block(0, 0, hitCount, 1);
            // update rejectHit
            rejectHit = hitCount > 0;
        }
        
            // DEBUG PURPOSE
            // if (rejectHit) {
            //     logger().debug("reject hit: {} {} {}", clusterRecord.loc(rejectHitIndex(0), 0), clusterRecord.loc(rejectHitIndex(0), 1), clusterRecord.loc(rejectHitIndex(0), 2));
            //     logger().debug("minDist = {}", minDistSquare);
            // }
    } else {

        // update rejectHit
        rejectHit = false;
    }
}


void GUI::MouseRejectCluster() {
/// called by "MouseDownCallback"

    if (!rejectHit) return;

    const int hitCount = rejectHitIndex.rows();
    for (int i=0; i<hitCount; i++) {
        clusterRecord.alive(rejectHitIndex(i)) = false;
    }

    rejectHit = false;
    UpdateClusterPointLoc();
    logger().info("Mouse pick rejected {} clusters", hitCount);
}

}  // namespace zebrafish
