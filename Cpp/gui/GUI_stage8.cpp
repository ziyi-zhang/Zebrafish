#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/VTUwriter.h>
#include <zebrafish/TiffReader.h>
#include <zebrafish/zebra-analysis.hpp>
#include <zebrafish/Padding.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <LBFGS.h>
#include <string>
#include <ctime>
#include <stdio.h>

namespace zebrafish
{

    namespace
    {

        using std::string;
        string GetFileName(const string &imagePath, int frameindex, const string &ext, const string &extra)
        {

            string fileName;
            // erase extension (.tif / .tiff)
            size_t lastdot = imagePath.find_last_of(".");
            if (lastdot == string::npos)
                fileName = imagePath;
            else
                fileName = imagePath.substr(0, lastdot);
            // erase path
            //size_t lastDelimiter = fileName.find_last_of("/");
            //if (lastDelimiter != string::npos) fileName.erase(0, lastDelimiter + 1);
            //lastDelimiter = fileName.find_last_of("\\");
            //if (lastDelimiter != string::npos) fileName.erase(0, lastDelimiter + 1);
            // add extra
            fileName += "-";
            fileName += extra;
            // add time
            char buffer[80];
            time_t rawTime;
            struct tm *timeInfo;
            time(&rawTime);
            timeInfo = localtime(&rawTime);
            strftime(buffer, 80, "%m_%dT%H_%M", timeInfo);
            fileName += "-";
            fileName += buffer;
            // add index
            if (frameindex >= 0)
            {
                fileName += "-frame";
                fileName += std::to_string(frameindex);
            }

            // add extension
            fileName += ext;

            return fileName;
        }

        void GetDisplacement(const std::vector<Eigen::MatrixXd> &markerPointLocArray, int currFrameIdx, bool incremental, Eigen::MatrixXd &disp)
        {

            disp.resizeLike(markerPointLocArray[currFrameIdx]);

            // loc[currFrame] - loc[lastFrame]
            if (currFrameIdx == 0)
            {
                disp.setZero();
            }
            else
            {
                if (incremental)
                    disp = markerPointLocArray[currFrameIdx] - markerPointLocArray[currFrameIdx - 1];
                else // accumulative
                    disp = markerPointLocArray[currFrameIdx] - markerPointLocArray[0];
            }
        }

        double GetTriArea(const Eigen::MatrixXd &meshPoint, int v1, int v2, int v3)
        {
            // should have used cross product...

            double a2 = (meshPoint(v1, 0) - meshPoint(v2, 0)) * (meshPoint(v1, 0) - meshPoint(v2, 0)) +
                        (meshPoint(v1, 1) - meshPoint(v2, 1)) * (meshPoint(v1, 1) - meshPoint(v2, 1)) +
                        (meshPoint(v1, 2) - meshPoint(v2, 2)) * (meshPoint(v1, 2) - meshPoint(v2, 2));
            double b2 = (meshPoint(v1, 0) - meshPoint(v3, 0)) * (meshPoint(v1, 0) - meshPoint(v3, 0)) +
                        (meshPoint(v1, 1) - meshPoint(v3, 1)) * (meshPoint(v1, 1) - meshPoint(v3, 1)) +
                        (meshPoint(v1, 2) - meshPoint(v3, 2)) * (meshPoint(v1, 2) - meshPoint(v3, 2));
            double c2 = (meshPoint(v3, 0) - meshPoint(v2, 0)) * (meshPoint(v3, 0) - meshPoint(v2, 0)) +
                        (meshPoint(v3, 1) - meshPoint(v2, 1)) * (meshPoint(v3, 1) - meshPoint(v2, 1)) +
                        (meshPoint(v3, 2) - meshPoint(v2, 2)) * (meshPoint(v3, 2) - meshPoint(v2, 2));

            return 0.25 * std::sqrt(4.0 * a2 * b2 - (a2 + b2 - c2) * (a2 + b2 - c2)); // Helen's formula
        }

        void GetPerpendicularVector(const Eigen::MatrixXd &meshPoint, int v1, int v2, int v3, Eigen::Vector3d &V31_perp)
        {

            using Eigen::Vector3d;

            Vector3d Vik, Vij, Vip, Vpj;
            Vij << meshPoint(v2, 0) - meshPoint(v1, 0), meshPoint(v2, 1) - meshPoint(v1, 1), meshPoint(v2, 2) - meshPoint(v1, 2);
            Vik << meshPoint(v3, 0) - meshPoint(v1, 0), meshPoint(v3, 1) - meshPoint(v1, 1), meshPoint(v3, 2) - meshPoint(v1, 2);
            Vip = Vik.normalized() * (Vij.dot(Vik.normalized()));
            Vpj = Vij - Vip;
            V31_perp = Vpj.normalized() * Vik.norm();
        }

        void GetJacobian(const Eigen::MatrixXd &meshPoint, const Eigen::MatrixXi &markerMeshArray, const Eigen::MatrixXd &displacement, Eigen::MatrixXd &jacobian)
        {
            // Calculate gradient on a mesh for x, y, z
            // Ref: Geometric Modeling - Daniele Panozzo - Single Patch Parametrization

            const int N = markerMeshArray.rows(); // number of triangles
            int v1, v2, v3;
            double area;
            Eigen::Vector3d Vki_perp, Vij_perp;
            jacobian.resize(N, 3);

            for (int i = 0; i < N; i++)
            { // i-th triangle

                v1 = markerMeshArray(i, 0);
                v2 = markerMeshArray(i, 1);
                v3 = markerMeshArray(i, 2);
                GetPerpendicularVector(meshPoint, v1, v2, v3, Vki_perp);
                GetPerpendicularVector(meshPoint, v2, v3, v1, Vij_perp);
                area = GetTriArea(meshPoint, v1, v2, v3);
                // DEBUG PURPOSE
                /*
        using namespace std;
        cout << "-->>" << endl;
        cout << meshPoint.row(v1) << endl;
        cout << meshPoint.row(v2) << endl;
        cout << meshPoint.row(v3) << endl;
        cout << Vki_perp.transpose() << endl;
        cout << Vij_perp.transpose() << endl;
        cout << area << endl;
        */
                // DEBUG PURPOSE

                for (int dim = 0; dim < 3; dim++)
                { // x, y, z

                    jacobian(i, dim) = (displacement(v2, dim) - displacement(v1, dim)) * Vki_perp(dim) +
                                       (displacement(v3, dim) - displacement(v1, dim)) * Vij_perp(dim);
                    jacobian(i, dim) /= (2.0 * area);
                }
            }
        }

        void GetPhysicalLocation(const std::vector<Eigen::MatrixXd> &locArray, double resolutionX, double resolutionY, double resolutionZ, std::vector<Eigen::MatrixXd> &locArray_out)
        {

            const int N = locArray.size();
            locArray_out.resize(N);
            bool physical;

            if (resolutionX > 0 && resolutionY > 0 && resolutionZ > 0)
            {
                physical = true;
                logger().info("Using physical unit. Row dist = {}, col dist = {}, depth dist = {}", resolutionX, resolutionY, resolutionZ);
            }
            else
            {
                physical = false; // the user does not input valid physical resolution
                logger().info("Physical resolution invalid. Using pixel as unit.");
            }

            for (int i = 0; i < locArray.size(); i++)
            {

                // Note: locArray has inverted XY
                locArray_out[i] = locArray[i];
                if (physical)
                {
                    locArray_out[i].col(0) *= resolutionY;
                    locArray_out[i].col(1) *= resolutionX;
                    locArray_out[i].col(2) *= resolutionZ;
                }
            }
        }

    } // anonymous namespace

    ////////////////////////////////////////////////////////////////////////////////////////
    // Stage 8: Dispalcement & Export

    void GUI::DrawStage8()
    {

        // cropActive
        if (meanCrop.cropActive && meanCrop.downClicked)
        {
            // crop activated & has been updated (not default value)
            Eigen::MatrixXd lineColor(1, 3);
            lineColor << 0.77, 0.28, 0.24;
            viewer.data().line_width = 2.0f;

            // upper-left corner (x0, y0)
            // lower-right corner (x1, y1)
            float x0 = meanCrop.baseLoc(0);
            float x1 = std::max(x0, meanCrop.currentLoc(0));
            float y0 = meanCrop.baseLoc(1);
            float y1 = std::min(y0, meanCrop.currentLoc(1));
            DrawRect(x0, y0, x1, y1, lineColor);
        }

        // showCropArea
        if (meanCrop.showCropArea && currentLoadedFrames > 0)
        {
            // show the area specified by current [r0, c0] x [r1, c1]
            Eigen::MatrixXd lineColor(1, 3);
            lineColor << 0.77, 0.28, 0.24;
            viewer.data().line_width = 2.0f;

            // upper-left corner (x0, y0)
            // lower-right corner (x1, y1)
            float x0 = (meanCrop.c0 == -1) ? 0 : meanCrop.c0;
            float x1 = (meanCrop.c1 == -1) ? imgCols : meanCrop.c1;
            float y0 = (meanCrop.r0 == -1) ? imgRows : imgRows - meanCrop.r0;
            float y1 = (meanCrop.r1 == -1) ? 0 : imgRows - meanCrop.r1;
            DrawRect(x0, y0, x1, y1, lineColor);
        }

        ImGui::Separator(); /////////////////////////////////////////

        // Visualize marker cluster points
        static int pointSize = 7;
        if (showMarkerPoints) {

            viewer.data().point_size = pointSize;

            if (!markerPointLocArray.empty()) {
                // show optimized markers
                if (!manualOverrideMarkerVis) {
                    // show markers in the frame that is currently focused
                    viewer.data().add_points(
                        markerPointLocArray[frameToShow],
                        markerPointColor);
                } else {
                    // show markers in manually selected frames
                    for (int i = 0; i < markerPointStatusArray.rows(); i++) {
                        if (!markerPointStatusArray(i))
                            continue;
                        viewer.data().add_points(
                            markerPointLocArray[i],
                            markerPointColor);
                    }
                }
            }
        }

        // Visualize unsuccessfully optimized points
        PlotBadDCPoints();

        // Visualize meshes
        DrawMarkerMesh();

        // visualize optical flow points
        // FIXME: DO NOT DO THIS
        if (showOpticalFlow)
        {

            viewer.data().point_size = pointSize;
            Eigen::MatrixXd pointColor(1, 3);
            pointColor << 0.91, 0.69, 0.01;

            if (!opticalFlowCorrection.empty() && frameToShow < opticalFlowCorrection.size())
            {
                // show optical flow corrected location in the frame that is currently focused
                Eigen::MatrixXd temp(markerPointLocArray[frameToShow].rows(), markerPointLocArray[frameToShow].cols());
                temp.col(0) = markerPointLocArray[frameToShow].col(0) + opticalFlowCorrection[frameToShow].col(1);
                temp.col(1) = markerPointLocArray[frameToShow].col(1) - opticalFlowCorrection[frameToShow].col(0);
                temp.col(2) = markerPointLocArray[frameToShow].col(2) + opticalFlowCorrection[frameToShow].col(2);
                viewer.data().add_points(
                    temp,
                    pointColor);
            }
        }

        // Visualize manual marker drag code
        MarkerDragVisualization();

        DrawReferenceDots();

        ImGui::Separator(); /////////////////////////////////////////

        if (analysisInputPath.empty() && ImGui::CollapsingHeader("Displacement", ImGuiTreeNodeFlags_DefaultOpen)) {  // do not draw if in re-analysis mode

            static bool logEnergy = false;
            static std::string calcDispStr = "";
            if (ImGui::TreeNode("Advanced depth correction")) {

                const float inputWidth = ImGui::GetWindowWidth() / 3.0;
                ImGui::PushItemWidth(inputWidth);

                if (ImGui::SliderFloat("DC gap", &depthCorrectionGap, 0, 0.3, "%.3f pixels")) {
                    calcDispStr = "";
                }
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Depth search gap in pixels");
                }
                if (ImGui::SliderInt("DC trial numbers", &depthCorrectionNum, 0, 50, "%d * gap")) {
                    calcDispStr = "";
                }
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Depth search trial numbers. A vertical interval of length [2*num+1]x[gap] pixels will be searched to determine whether the depth should be modified.");
                }
                if (ImGui::SliderFloat("Max XY displacement", &optimMaxXYDisp, 0.5, 9, "%.2f pixels")) {
                    calcDispStr = "";
                }
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Maximum displacement in XY plane during depth correction");
                }
                ImGui::Checkbox("Second round DC", &secondRoundDepthCorrection);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Increase depth search precision but time-consuming");
                }
                ImGui::Checkbox("Log energy matrix", &logEnergy);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("[Debug purpose] whether log the energy for all depth trials");
                }

                ImGui::PopItemWidth();

                ImGui::TreePop();
                ImGui::Separator();
            }

            if (ImGui::Button("Calculate Displacement")) {
                try {
                    if (OptimizeAllFrames(logEnergy))
                        calcDispStr = "Successful";
                    else
                        calcDispStr = "exception: see log";
                    logger().debug("   <button> Calculate Displacement");
                } catch (const std::exception &e) {
                    logger().error("   <button> [Calculate Displacement] Fatal error encountered.");
                    std::cerr << "   <button> [Calculate Displacement] Fatal error encountered." << std::endl;
                    calcDispStr = "Fatal error";
                }
            }
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Calculate displacement for all loaded frames");
            }
            ImGui::SameLine();
            ImGui::Text("%s", calcDispStr.c_str());

            ImGui::Separator(); /////////////////////////////////////////

            if (ImGui::TreeNode("Advanced visualization     ")) {

                const float inputWidth = ImGui::GetWindowWidth() / 3.0;
                ImGui::PushItemWidth(inputWidth);
                ImGui::Checkbox("Show background image", &showBackgroundImage);
                ImGui::Checkbox("Show marker centers", &showMarkerPoints);
                ImGui::Checkbox("Show mesh", &showMarkerMesh);
                ImGui::SliderInt("Point size", &pointSize, 1, 30);
                ImGui::SliderFloat("Line width", &lineWidth, 1, 16);
                ImGui::PopItemWidth();

                ImGui::TreePop();
            }

            // this loads the code to render the GUI about mouse draging
            ImGui::Separator(); /////////////////////////////////////////
            RenderMarkerDragGUI();
        }

        ImGui::Separator(); /////////////////////////////////////////

        if (ImGui::CollapsingHeader("Analysis", ImGuiTreeNodeFlags_DefaultOpen)) {

            if (ImGui::Button("Padding test")) {
                Eigen::MatrixXd appendV;
                Eigen::MatrixXi appendF;                
                padding::ComputeOneRing(analysisPara.V[0], analysisPara.F, analysisPara.markerRCMap, appendV, appendF);
                std::cerr << "markerMeshArray\n" << analysisPara.F << std::endl;
                padding::AddOneRing<Eigen::MatrixXd>(appendV, appendF, analysisPara.V[0], analysisPara.F);
                UpdateAnalysisPointLocArray();
                std::cerr << "markerMeshArray\n" << analysisPara.F << std::endl;
            }

            // average displacement area crop
            if (ImGui::Checkbox("[Mouse] global disp area", &meanCrop.cropActive)) {
                if (!meanCrop.cropActive)
                    logger().debug("[Mouse] global disp area: de-activated.");
                else {
                    logger().debug("[Mouse] global disp area: activated.");
                    meanCrop.showCropArea = true;
                }
            }
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Crop an area that will be used to estimate the global displacement.\nGlobal displacement will be subtracted in analysis. By default we use the entire image.");
            }

            ImGui::Separator(); /////////////////////////////////////////

            if (ImGui::TreeNode("Advanced analysis")) {

                const float inputWidth = ImGui::GetWindowWidth() / 3.0;
                ImGui::PushItemWidth(inputWidth);
                ImGui::InputDouble("offset", &analysisPara.offset);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Diagonal multiplier for box mesh");
                }
                ImGui::InputDouble("radius-edge ratio", &analysisPara.radius_edge_ratio);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Radius-edge ratio in TetGen. Cannot be smaller than 0.707.");
                }
                ImGui::InputDouble("max tetrahedral volume", &analysisPara.max_tet_vol);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Maximum tetrahedral volume in TetGen");
                }
                ImGui::InputDouble("E (Young's modulus)", &analysisPara.E);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Young's Modulus [Pascal]");
                }
                ImGui::InputDouble("nu (Poisson's ratio)", &analysisPara.nu);
                ImGui::Checkbox("linear material", &analysisPara.is_linear);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Use non-linear material");
                }
                ImGui::InputInt("discretization order", &analysisPara.discr_order);
                ImGui::InputInt("#refine", &analysisPara.n_refs);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Number of mesh uniform refinements");
                }
                ImGui::InputDouble("vismesh_rel_area", &analysisPara.vismesh_rel_area);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Desnsity of the output visualization");
                }
                ImGui::InputInt("upsample", &analysisPara.upsample);
                if (showTooltip && ImGui::IsItemHovered()) {
                    ImGui::SetTooltip("Recursively upsample to get a finer mesh");
                }
                ImGui::PopItemWidth();

                ImGui::TreePop();
            }

            static std::string runAnalysisStr = "";
            if (ImGui::Button("Run analysis")) {
                try {
                    if (!analysisInputPath.empty()) {
                        // re-run a previous experiment result
                        compute_analysis(
                            analysisPara.V,
                            analysisPara.F,
                            analysisInputPath,
                            analysisPara.E,
                            analysisPara.nu,
                            analysisPara.offset,
                            analysisPara.radius_edge_ratio,
                            analysisPara.max_tet_vol,
                            analysisPara.discr_order,
                            analysisPara.is_linear,
                            analysisPara.n_refs,
                            analysisPara.vismesh_rel_area,
                            analysisPara.upsample,
                            analysisPara.markerRCMap,
                            false);
                    }
                    else {

                        // prepare the displacement
                        std::vector<Eigen::MatrixXd> analysisDisplacementVec;
                        std::vector<Eigen::MatrixXd> markerPointLocArray_phy;
                        GetPhysicalLocation(markerPointLocArray, resolutionX, resolutionY, resolutionZ, markerPointLocArray_phy);
                        // determine markers that are used to compute global displacement
                        std::vector<bool> markerInAvgDispArea;
                        GetMarkersInAvgDispArea(markerInAvgDispArea);
                        // remove global displacement
                        for (int i = 0; i < currentLoadedFrames; i++) {

                            Eigen::MatrixXd V_analysis = markerPointLocArray_phy[i];
                            // remove the global movement
                            Eigen::RowVectorXd meanV(1, 3); // (V_analysis - markerPointLocArray_phy[0]).colwise().mean();
                            meanV << 0.0f, 0.0f, 0.0f;
                            int count = 0;
                            for (int j = 0; j < markerInAvgDispArea.size(); j++)
                            {
                                if (!markerInAvgDispArea[j])
                                    continue;
                                meanV += V_analysis.row(j) - markerPointLocArray_phy[0].row(j);
                                count += 1;
                            }
                            meanV /= double(count);
                            logger().debug("Avg global displacement in frame {}: {} {} {}", i, meanV(0), meanV(1), meanV(2));

                            V_analysis.rowwise() -= meanV;
                            // push to new analysis vector
                            analysisDisplacementVec.push_back(V_analysis); // want accumulativeDisplacement (relative)
                        }

                        // run analysis
                        std::string path = GetFileName(imagePath, -1, "", "analysis");
                        // [NOTE]: E is in unit of [Pascal], displacement is in unit of [um]. So multiply by 1e-6
                        compute_analysis(
                            analysisDisplacementVec, 
                            markerMeshArray, 
                            path, 
                            analysisPara.E, 
                            analysisPara.nu, 
                            analysisPara.offset, 
                            analysisPara.radius_edge_ratio,
                            analysisPara.max_tet_vol, 
                            analysisPara.discr_order, 
                            analysisPara.is_linear, 
                            analysisPara.n_refs, 
                            analysisPara.vismesh_rel_area, 
                            analysisPara.upsample, 
                            analysisPara.markerRCMap,
                            true);
                    }

                    runAnalysisStr = "Done";
                }
                catch (const std::exception &e) {
                    logger().error("   <button> [Run analysis] Fatal error encountered.");
                    std::cerr << "   <button> [Run analysis] Fatal error encountered." << std::endl;
                    runAnalysisStr = "Fatal error";
                }
            }
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Run simulation and save result to analysis VTU file");
            }
            ImGui::SameLine();
            ImGui::Text("%s", runAnalysisStr.c_str());
        }

        ImGui::Separator(); /////////////////////////////////////////

        static bool onlySaveFirstFrameMesh = true;
        static bool saveAccumulativeDisplacement = true;
        static bool saveAccumulativeDisplacement_relative = true;
        static bool saveIncrementalDisplacement = true;
        static bool saveIncrementalDisplacement_relative = true;
        static bool saveMeshPoint = false;  // vtu
        static bool saveMeshCell = false;   // vtu
        static bool saveMarkerImage = true; // tif
        static bool saveCellImage = true;   // tif
        static int cellChannel = -1;
        if (cellChannel < 0)
        {
            // make a guess of which channle is the cell channel
            if (channelPerSlice < 2)
            {
                cellChannel = 0;
                saveCellImage = false;
            }
            else if (channelToLoad == 0)
            {
                cellChannel = 1;
            }
            else
            {
                cellChannel = 0;
            }
        }

        // static bool saveEntireImage = false;
        if (analysisInputPath.empty() && ImGui::CollapsingHeader("Save & Export", ImGuiTreeNodeFlags_DefaultOpen)) {

            static std::string saveStr;
            static bool meshPointSaveFlag, meshCellSaveFlag, imageSaveFlag;
            if (ImGui::Button("Export TIFF (and VTU)")) {

                meshPointSaveFlag = true;
                meshCellSaveFlag = true;
                imageSaveFlag = true;

                // Save mesh to VTU (point data)
                if (saveMeshPoint)
                    meshPointSaveFlag = SaveMeshToVTU_point(onlySaveFirstFrameMesh, saveAccumulativeDisplacement, saveAccumulativeDisplacement_relative, saveIncrementalDisplacement, saveIncrementalDisplacement_relative);
                // Save mesh to VTU (cell data)
                if (saveMeshCell)
                    meshCellSaveFlag = SaveMeshToVTU_cell(onlySaveFirstFrameMesh);
                // Save image to TIFF
                imageSaveFlag = SaveImageToTIFF(saveMarkerImage, saveCellImage, cellChannel);

                saveStr = (meshPointSaveFlag && meshCellSaveFlag && imageSaveFlag) ? "Data saved" : "Save failed";

                logger().debug("   <button> Export VTU & TIFF");
            }
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Export the cropped area as new TIFF images. (Optional) Save the raw displacements and the moving mesh. ");
            }
            ImGui::SameLine();
            ImGui::Text("%s", saveStr.c_str());

            if (ImGui::TreeNode("Advanced export")) {

                ImGui::Checkbox("Save cropped image (marker channel)", &saveMarkerImage);
                ImGui::Checkbox("Save cropped image (cell channel)", &saveCellImage);
                ImGui::SliderInt("Cell channel", &cellChannel, 0, channelPerSlice - 1);
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("The channel with cell image.\nNote: the loaded channel is the one with marker information.");
                }

                ImGui::Separator();
                /////////////////////////////////////

                ImGui::Checkbox("Save mesh as VTU", &saveMeshPoint);
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("If checked, save the mesh (moving or fixed) as VTU along with any combination of displacement configurations as shown below.");
                }
                ImGui::Separator();
                /////////////////////////////////////

                ImGui::Checkbox("Only save first frame mesh", &onlySaveFirstFrameMesh);
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("If checked, VTU file will always use the first-frame mesh location for all frames.\nOtherwise the mesh will move in different frames.");
                }

                ImGui::Checkbox("Save accumulative displacement (absolute)", &saveAccumulativeDisplacement);
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("Save the displacement relative to the first frame");
                }
                ImGui::Checkbox("Save accumulative displacement (relative)", &saveAccumulativeDisplacement_relative);
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("Save the displacement relative to the first frame. The mean of all the displacements will be subtracted.");
                }
                ImGui::Checkbox("Save incremental displacement (absolute)", &saveIncrementalDisplacement);
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("Save the displacement relative to the previous frame");
                }
                ImGui::Checkbox("Save incremental displacement (relative)", &saveIncrementalDisplacement_relative);
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("Save the displacement relative to the previous frame. The mean of all the displacements will be subtracted.");
                }
                ImGui::Checkbox("Save jacobian", &saveMeshCell);
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("Save the jacobian of the moving mesh to VTU file (DEBUG ONLY)");
                }

                ImGui::Separator();
                /////////////////////////////////////

                if (ImGui::Button("SaveMeshToVTU (point data)"))
                {
                    meshPointSaveFlag = SaveMeshToVTU_point(onlySaveFirstFrameMesh, saveAccumulativeDisplacement, saveAccumulativeDisplacement_relative, saveIncrementalDisplacement, saveIncrementalDisplacement_relative);
                }
                if (ImGui::Button("SaveMeshToVTU (cell data)"))
                {
                    meshCellSaveFlag = SaveMeshToVTU_cell(onlySaveFirstFrameMesh);
                }
                if (ImGui::Button("SaveImageToTiff"))
                {
                    imageSaveFlag = SaveImageToTIFF(saveMarkerImage, saveCellImage, cellChannel);
                }

                /////////////////////////////////////

                if (ImGui::Button("Export to OBJ"))
                {
                    SaveMeshToOBJ();
                }
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("(DEBUG ONLY)");
                }
                if (ImGui::Button("Export displacement to TXT"))
                {
                    SaveDisplacementToTXT();
                }
                if (showTooltip && ImGui::IsItemHovered())
                {
                    ImGui::SetTooltip("Raw displacement w.r.t. the first frame, without subtracting mean global movement. (DEBUG ONLY)");
                }

                ImGui::TreePop();
                ImGui::Separator();
            }
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    // Optimization

    void GUI::GetMarkersInAvgDispArea(std::vector<bool> &markerInAvgDispArea) {

        const int N = markerArray[0].num;
        markerInAvgDispArea.resize(N);
        int count = 0;

        int x0 = (meanCrop.r0 == -1) ? 0 : meanCrop.r0;
        int x1 = (meanCrop.r1 == -1) ? imgRows : meanCrop.r1;
        int y0 = (meanCrop.c0 == -1) ? 0 : meanCrop.c0;
        int y1 = (meanCrop.c1 == -1) ? imgCols : meanCrop.c1;

        for (int i = 0; i < N; i++) {
            if (markerArray[0].loc(i, 0) > x0 && markerArray[0].loc(i, 0) < x1 &&
                markerArray[0].loc(i, 1) > y0 && markerArray[0].loc(i, 1) < y1) {
                markerInAvgDispArea[i] = true;
                count += 1;
            } else {
                markerInAvgDispArea[i] = false;
            }
        }
        logger().info("#markerInAvgDispArea = {}", count);
    }

    bool GUI::OptimizeAllFrames(bool logEnergy) {

        bool res = true;
        int currentFrame;

        for (currentFrame = 0; currentFrame < currentLoadedFrames - 1; currentFrame++) {

            // OptimizeOneFrame(currentFrame);
            /// Note: why don't we optimize here anymore?
            /// Experiments show that there are rare cases where the optimization may
            /// diverge to somewhere wrong. This is potentially due to wrong depth.

            // apply optical flow result on the next frame
            ApplyOpticalFlow(currentFrame);

            // depth correction for the frame that was just updated
            if (!MarkerRecursiveDepthCorrection(currentFrame + 1, depthCorrectionNum, depthCorrectionGap, logEnergy, true, true)) {
                res = false;
                logger().warn("Depth search unsuccessful on frame {}. See above for detailed reasons.", currentFrame + 1);
            }
        }

        // update visualization variable
        UpdateMarkerPointLocArray();

        return res;
    }

    void GUI::ApplyOpticalFlow(int prevFrameIdx)
    {
        // Apply optical flow results to this frame

        const int N = markerArray[prevFrameIdx].num;

        for (int i = 0; i < N; i++)
        {

            markerArray[prevFrameIdx + 1].loc(i, 0) = markerArray[prevFrameIdx].loc(i, 0) + opticalFlowCorrection[prevFrameIdx](i, 0); // x
            markerArray[prevFrameIdx + 1].loc(i, 1) = markerArray[prevFrameIdx].loc(i, 1) + opticalFlowCorrection[prevFrameIdx](i, 1); // y
            markerArray[prevFrameIdx + 1].loc(i, 2) = markerArray[prevFrameIdx].loc(i, 2) + opticalFlowCorrection[prevFrameIdx](i, 2); // z
        }
    }

    void GUI::OptimizeOneFrame(int prevFrameIdx)
    {
        // Optimize from frame "prevFrameIdx" to frame "prevFrameIdx+1"

        /////////////////////////////////////////////////
        // prepare LBFGS
        LBFGSpp::LBFGSParam<double> param;
        param.epsilon = optimEpsilon;
        param.max_iterations = optimMaxIt;

        // prepare for parallel optimization
        const int N = markerArray[prevFrameIdx].num;
        markerArray[prevFrameIdx + 1].num = N;

        logger().info("Optimization from frame {} to {}", prevFrameIdx, prevFrameIdx + 1);
        logger().info("Optimization #starting points = {}", N);

        // Optimization
        logger().info(">>>>>>>>>> Before optimization >>>>>>>>>>");
        tbb::parallel_for(tbb::blocked_range<int>(0, N),
            //////////////////////////////////////
            // lambda function for parallel_for //
            //////////////////////////////////////
            [this /*.markerArray[?], &bsplineArray[?], .opticalFlowCorrection*/, &param, prevFrameIdx](const tbb::blocked_range<int> &r) {
                // NOTE: LBFGSSolver is NOT thread safe. This must be instantiated for every thread
                LBFGSpp::LBFGSSolver<double> solver(param);

                // NOTE: the "variable count" used by "Autodiff" will be stored in
                //       thread-local memory, so this must be set for every thread
                DiffScalarBase::setVariableCount(3);

                for (int ii = r.begin(); ii != r.end(); ++ii)
                {

                    Eigen::VectorXd vec(3, 1);
                    vec(0) = markerArray[prevFrameIdx].loc(ii, 0) + opticalFlowCorrection[prevFrameIdx](ii, 0); // x
                    vec(1) = markerArray[prevFrameIdx].loc(ii, 1) + opticalFlowCorrection[prevFrameIdx](ii, 1); // y
                    vec(2) = markerArray[prevFrameIdx].loc(ii, 3);                                              // r
                    double res;

                    ///////////////////////////////////
                    // lambda function for optimizer //
                    ///////////////////////////////////
                    auto func = [this /*.markerArray[?], .bsplineArray[?], .opticalFlowCorrection*/, ii, prevFrameIdx](const Eigen::VectorXd &x, Eigen::VectorXd &grad) {
                        DScalar ans;

                        if (!cylinder::IsValid(bsplineArray[prevFrameIdx + 1], x(0), x(1), markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2), x(2), cylinder::H))
                        {
                            grad.setZero();
                            return 1.0;
                        }
                        cylinder::EvaluateCylinder(bsplineArray[prevFrameIdx + 1], DScalar(0, x(0)), DScalar(1, x(1)), markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2), DScalar(2, x(2)), cylinder::H, ans, reverseColor);
                        grad.resize(3, 1);
                        grad = ans.getGradient();
                        return ans.getValue();
                    };
                    // NOTE: the template of "solver.minimize" does not accept a temprary variable (due to non-const argument)
                    //       so we define a "func" and pass it in
                    int it = solver.minimize(func, vec, res);
                    ///////////////////////////////////

                    markerArray[prevFrameIdx + 1].loc(ii, 0) = vec(0);                                                                            // x
                    markerArray[prevFrameIdx + 1].loc(ii, 1) = vec(1);                                                                            // y
                    markerArray[prevFrameIdx + 1].loc(ii, 2) = markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2); // z
                    markerArray[prevFrameIdx + 1].loc(ii, 3) = vec(2);                                                                            // r
                    markerArray[prevFrameIdx + 1].energy(ii) = res;                                                                               // energy
                }
            });
        logger().info("<<<<<<<<<< After optimization <<<<<<<<<<");
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    // Export

    bool GUI::SaveMeshToVTU_point(bool onlySaveFirstFrameMesh, bool saveAccumulativeDisplacement, bool saveAccumulativeDisplacement_relative, bool saveIncrementalDisplacement, bool saveIncrementalDisplacement_relative)
    {
        // return true if successful

        static VTUWriter vtuWriter;
        std::vector<Eigen::MatrixXd> markerPointLocArray_phy;
        GetPhysicalLocation(markerPointLocArray, resolutionX, resolutionY, resolutionZ, markerPointLocArray_phy);

        // determine markers that are used to compute global displacement
        std::vector<bool> markerInAvgDispArea;
        GetMarkersInAvgDispArea(markerInAvgDispArea);

        for (int i = 0; i < currentLoadedFrames; i++)
        {

            // prepare filename
            std::string vtuFileName;
            if (onlySaveFirstFrameMesh)
                vtuFileName = GetFileName(imagePath, i, ".vtu", "FixedMesh-point");
            else
                vtuFileName = GetFileName(imagePath, i, ".vtu", "MovingMesh-point");
            // prepare absolute displacement
            Eigen::MatrixXd accumulativeDisplacement, incrementalDisplacement;
            GetDisplacement(markerPointLocArray_phy, i, false, accumulativeDisplacement);
            GetDisplacement(markerPointLocArray_phy, i, true, incrementalDisplacement);
            if (saveAccumulativeDisplacement)
            {
                vtuWriter.add_field("accumulative displacement (absolute)", accumulativeDisplacement);
            }
            if (saveIncrementalDisplacement)
            {
                vtuWriter.add_field("incremental displacement (absolute)", incrementalDisplacement);
            }
            // prepare relative displacement
            if (saveAccumulativeDisplacement_relative)
            {
                Eigen::RowVectorXd meanDisp(1, 3);
                meanDisp << 0.0f, 0.0f, 0.0f;
                int count = 0;
                for (int j = 0; j < markerInAvgDispArea.size(); j++)
                {
                    if (!markerInAvgDispArea[j])
                        continue;
                    meanDisp += accumulativeDisplacement.row(j);
                    count += 1;
                }
                meanDisp /= double(count);

                accumulativeDisplacement.rowwise() -= meanDisp;
                vtuWriter.add_field("accumulative displacement (relative)", accumulativeDisplacement);
            }
            if (saveIncrementalDisplacement_relative)
            {
                Eigen::RowVectorXd meanDisp(1, 3);
                meanDisp << 0.0f, 0.0f, 0.0f;
                int count = 0;
                for (int j = 0; j < markerInAvgDispArea.size(); j++)
                {
                    if (!markerInAvgDispArea[j])
                        continue;
                    meanDisp += incrementalDisplacement.row(j);
                    count += 1;
                }
                meanDisp /= double(count);

                incrementalDisplacement.rowwise() -= meanDisp;
                vtuWriter.add_field("incremental displacement (relative)", incrementalDisplacement);
            }

            // prepare mesh point array
            Eigen::MatrixXd meshPoint;
            if (onlySaveFirstFrameMesh)
                meshPoint = markerPointLocArray_phy[0];
            else
                meshPoint = markerPointLocArray_phy[i];
            // correct depth
            meshPoint.col(2).array() -= layerBegin;

            // write to VTU
            if (!vtuWriter.write_tet_mesh(vtuFileName, meshPoint, markerMeshArray))
            {
                return false;
            }
        }

        return true;
    }

    bool GUI::SaveMeshToVTU_cell(bool onlySaveFirstFrameMesh)
    {
        // return true if successful

        static VTUWriter vtuWriter;
        std::vector<Eigen::MatrixXd> markerPointLocArray_phy;
        GetPhysicalLocation(markerPointLocArray, resolutionX, resolutionY, resolutionZ, markerPointLocArray_phy);

        for (int i = 0; i < currentLoadedFrames; i++)
        {

            // prepare filename
            std::string vtuFileName;
            if (onlySaveFirstFrameMesh)
                vtuFileName = GetFileName(imagePath, i, ".vtu", "FixedMesh-cell");
            else
                vtuFileName = GetFileName(imagePath, i, ".vtu", "MovingMesh-cell");

            Eigen::MatrixXd incrementalDisplacement;
            GetDisplacement(markerPointLocArray_phy, i, true, incrementalDisplacement);
            // prepare Jacobian
            Eigen::MatrixXd jacobian;
            GetJacobian(markerPointLocArray_phy[i], markerMeshArray, incrementalDisplacement, jacobian);
            vtuWriter.add_field("jacobian", jacobian);
            // prepare mesh point array
            Eigen::MatrixXd meshPoint;
            if (onlySaveFirstFrameMesh)
                meshPoint = markerPointLocArray_phy[0];
            else
                meshPoint = markerPointLocArray_phy[i];
            meshPoint.col(2).array() -= layerBegin;

            // write to VTU
            if (!vtuWriter.write_tet_mesh(vtuFileName, meshPoint, markerMeshArray))
            {
                return false;
            }
        }

        return true;
    }

    bool GUI::SaveImageToTIFF(bool saveMarkerImage, bool saveCellImage, int cellChannel)
    {
        // return true if successful

        imageData_t cellImgData;
        if (saveCellImage)
        {

            // only support loading one channel
            std::vector<bool> channelVec(channelPerSlice, false);
            channelVec[cellChannel] = true;
            // Read all desired frame to "imgData"
            ReadTif(imagePath, layerPerImg, channelVec, desiredFrames, cellImgData, imageCrop.r0, imageCrop.c0, imageCrop.r1, imageCrop.c1);
        }

        for (int i = 0; i < currentLoadedFrames; i++)
        {

            if (saveMarkerImage)
            {

                std::string tifFileName = GetFileName(imagePath, i, ".tif", "marker_channel");
                if (!WriteTif(tifFileName, imgData[i], layerBegin, layerEnd))
                {
                    return false;
                }
            }

            if (saveCellImage)
            {

                std::string tifFileName = GetFileName(imagePath, i, ".tif", "cell_channel");
                if (!WriteTif(tifFileName, cellImgData[i], layerBegin, layerEnd))
                {
                    return false;
                }
            }
        }

        return true;
    }

    void GUI::SaveMeshToOBJ()
    {

        std::string objFileName;
        Eigen::MatrixXd meshPoint;
        std::vector<Eigen::MatrixXd> markerPointLocArray_phy;
        GetPhysicalLocation(markerPointLocArray, resolutionX, resolutionY, resolutionZ, markerPointLocArray_phy);

        for (int i = 0; i < currentLoadedFrames; i++)
        {

            // prepare filename
            objFileName = GetFileName(imagePath, i, ".obj", "marker");

            // prepare mesh point array
            meshPoint = markerPointLocArray_phy[i];
            // correct depth
            meshPoint.col(2).array() -= layerBegin;

            // write to OBJ
            if (!igl::writeOBJ(objFileName, meshPoint, markerMeshArray))
            {
                logger().warn("Write to OBJ failed");
            }
        }
    }

    void GUI::SaveDisplacementToTXT()
    {

        FILE *pFile;
        std::string fileName;
        Eigen::MatrixXd displacement;
        std::vector<Eigen::MatrixXd> markerPointLocArray_phy;
        GetPhysicalLocation(markerPointLocArray, resolutionX, resolutionY, resolutionZ, markerPointLocArray_phy);

        fileName = GetFileName(imagePath, 0, ".txt", "displacement");
        pFile = fopen(fileName.c_str(), "w+");

        for (int i = 0; i < currentLoadedFrames; i++)
        {

            GetDisplacement(markerPointLocArray_phy, i, false, displacement);

            fprintf(pFile, "## Raw displacement with respective to the first frame, without subtracting mean global displacement\n");
            fprintf(pFile, "## Frame %d\n", i);
            for (int r = 0; r < displacement.rows(); r++)
            {
                for (int c = 0; c < displacement.cols(); c++)
                {

                    fprintf(pFile, "%f ", displacement(r, c));
                }
                fprintf(pFile, "\n");
            }

            fprintf(pFile, "\n");
        }

        fclose(pFile);
    }

} // namespace zebrafish
