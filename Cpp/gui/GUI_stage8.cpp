#include <zebrafish/Cylinder.h>
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>
#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/VTUwriter.h>
#include <zebrafish/TiffReader.h>

#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>

#include <LBFGS.h>
#include <string>
#include <ctime>
#include <stdio.h>


namespace zebrafish {

namespace {

using std::string;
string GetFileName(const string &imagePath, int index, const string &ext, const string &extra) {

    string fileName;
    // erase extension (.tif / .tiff)
    size_t lastdot = imagePath.find_last_of(".");
    if (lastdot == string::npos)
        fileName = imagePath;
    else
        fileName = imagePath.substr(0, lastdot);
    // erase path
    size_t lastDelimiter = fileName.find_last_of("/");
    if (lastDelimiter != string::npos) fileName.erase(0, lastDelimiter + 1);
    lastDelimiter = fileName.find_last_of("\\");
    if (lastDelimiter != string::npos) fileName.erase(0, lastDelimiter + 1);
    // add extra
    fileName += "-";
    fileName += extra;
    // add index
    fileName += "-frame";
    fileName += std::to_string(index);
    // add time
    /*
    char buffer[80];
    time_t rawTime;
    struct tm *timeInfo;
    time(&rawTime);
    timeInfo = localtime(&rawTime);
    strftime(buffer, 80, "%Y-%m-%dT%H-%M-%S", timeInfo);
    fileName += "-";
    fileName += buffer;
    */
    // add extension
    fileName += ext;

    return fileName;
}


void GetDisplacement(const std::vector<Eigen::MatrixXd> &markerPointLocArray, int currFrameIdx, bool incremental, Eigen::MatrixXd &disp) {

    disp.resizeLike(markerPointLocArray[currFrameIdx]);

    // loc[currFrame] - loc[lastFrame]
    if (currFrameIdx == 0) {
        disp.setZero();
    } else {
        if (incremental)
            disp = markerPointLocArray[currFrameIdx] - markerPointLocArray[currFrameIdx - 1];
        else  // accumulative
            disp = markerPointLocArray[currFrameIdx] - markerPointLocArray[0];
    }
}


double GetTriArea(const Eigen::MatrixXd &meshPoint, int v1, int v2, int v3) {
// should have used cross product...

    double a2 = (meshPoint(v1, 0)-meshPoint(v2, 0)) * (meshPoint(v1, 0)-meshPoint(v2, 0)) + 
                (meshPoint(v1, 1)-meshPoint(v2, 1)) * (meshPoint(v1, 1)-meshPoint(v2, 1)) + 
                (meshPoint(v1, 2)-meshPoint(v2, 2)) * (meshPoint(v1, 2)-meshPoint(v2, 2));
    double b2 = (meshPoint(v1, 0)-meshPoint(v3, 0)) * (meshPoint(v1, 0)-meshPoint(v3, 0)) + 
                (meshPoint(v1, 1)-meshPoint(v3, 1)) * (meshPoint(v1, 1)-meshPoint(v3, 1)) + 
                (meshPoint(v1, 2)-meshPoint(v3, 2)) * (meshPoint(v1, 2)-meshPoint(v3, 2));
    double c2 = (meshPoint(v3, 0)-meshPoint(v2, 0)) * (meshPoint(v3, 0)-meshPoint(v2, 0)) + 
                (meshPoint(v3, 1)-meshPoint(v2, 1)) * (meshPoint(v3, 1)-meshPoint(v2, 1)) + 
                (meshPoint(v3, 2)-meshPoint(v2, 2)) * (meshPoint(v3, 2)-meshPoint(v2, 2));

    return 0.25 * std::sqrt(4.0 * a2 * b2 - (a2+b2-c2) * (a2+b2-c2));  // Helen's formula
}


void GetPerpendicularVector(const Eigen::MatrixXd &meshPoint, int v1, int v2, int v3, Eigen::Vector3d &V31_perp) {

    using Eigen::Vector3d;

    Vector3d Vik, Vij, Vip, Vpj;
    Vij << meshPoint(v2, 0)-meshPoint(v1, 0), meshPoint(v2, 1)-meshPoint(v1, 1), meshPoint(v2, 2)-meshPoint(v1, 2);
    Vik << meshPoint(v3, 0)-meshPoint(v1, 0), meshPoint(v3, 1)-meshPoint(v1, 1), meshPoint(v3, 2)-meshPoint(v1, 2);
    Vip = Vik.normalized() * (Vij.dot(Vik.normalized()));
    Vpj = Vij - Vip;
    V31_perp = Vpj.normalized() * Vik.norm();
}


void GetJacobian(const Eigen::MatrixXd &meshPoint, const Eigen::MatrixXi &markerMeshArray, const Eigen::MatrixXd &displacement, Eigen::MatrixXd &jacobian) {
// Calculate gradient on a mesh for x, y, z
// Ref: Geometric Modeling - Daniele Panozzo - Single Patch Parametrization

    const int N = markerMeshArray.rows();  // number of triangles
    int v1, v2, v3;
    double area;
    Eigen::Vector3d Vki_perp, Vij_perp;
    jacobian.resize(N, 3);

    for (int i=0; i<N; i++) {  // i-th triangle

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

        for (int dim=0; dim<3; dim++) {  // x, y, z

            jacobian(i, dim) = (displacement(v2, dim) - displacement(v1, dim)) * Vki_perp(dim) + 
                               (displacement(v3, dim) - displacement(v1, dim)) * Vij_perp(dim);
            jacobian(i, dim) /= (2.0 * area);
        }
    }
}

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Stage 8: Dispalcement & Export

void GUI::DrawStage8() {

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
                    markerPointColor
                );
            } else {
                // show markers in manually selected frames
                for (int i=0; i<markerPointStatusArray.rows(); i++) {
                    if (!markerPointStatusArray(i)) continue;
                    viewer.data().add_points(
                        markerPointLocArray[i],
                        markerPointColor
                    );
                }
            }
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

    // Visualize meshes
    DrawMarkerMesh();

    // visualize optical flow points
    // FIXME: DO NOT DO THIS
    if (showOpticalFlow) {

        viewer.data().point_size = pointSize;
        Eigen::MatrixXd pointColor(1, 3);
        pointColor << 0.91, 0.69, 0.01;

        if (!opticalFlowCorrection.empty() && frameToShow < opticalFlowCorrection.size()) {
            // show optical flow corrected location in the frame that is currently focused
            Eigen::MatrixXd temp(markerPointLocArray[frameToShow].rows(), markerPointLocArray[frameToShow].cols());
            temp.col(0) = markerPointLocArray[frameToShow].col(0) + opticalFlowCorrection[frameToShow].col(1);
            temp.col(1) = markerPointLocArray[frameToShow].col(1) - opticalFlowCorrection[frameToShow].col(0);
            temp.col(2) = markerPointLocArray[frameToShow].col(2) + opticalFlowCorrection[frameToShow].col(2);
            viewer.data().add_points(
                temp,
                pointColor
            );
        }
    }

    ImGui::Separator(); /////////////////////////////////////////

    if (ImGui::CollapsingHeader("Displacement", ImGuiTreeNodeFlags_DefaultOpen)) {

        if (ImGui::TreeNode("Advanced depth correction")) {

            const float inputWidth = ImGui::GetWindowWidth() / 3.0;
            ImGui::PushItemWidth(inputWidth);

            ImGui::SliderFloat("DC gap", &depthCorrectionGap, 0, 0.5, "%.3f pixels");
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Depth correction gap in pixels");
            }
            ImGui::SliderInt("DC trial numbers", &depthCorrectionNum, 0, 16, "%d * gap");
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Depth correction trial numbers. A vertical interval of length [range]x[gap] pixels will be searched to determine whether the depth should be corrected.");
            }

            ImGui::PopItemWidth();

            ImGui::TreePop();
            ImGui::Separator();
        }

        static std::string calcDispStr = "";
        if (ImGui::Button("Calculate Displacement")) {
            if (OptimizeAllFrames())
                calcDispStr = "Successful";
            else
                calcDispStr = "Failed";
            logger().debug("   <button> Calculate Displacement");
        }
        ImGui::SameLine();
        ImGui::Text("%s", calcDispStr.c_str());
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Calculate displacement for all loaded frames");
        }
    }

    if (ImGui::TreeNode("Advanced visualization  ")) {

        const float inputWidth = ImGui::GetWindowWidth() / 3.0;
        ImGui::PushItemWidth(inputWidth);
        ImGui::Checkbox("Show background image", &showBackgroundImage);
        ImGui::Checkbox("Show marker centers", &showMarkerPoints);
        ImGui::Checkbox("Show mesh", &showMarkerMesh);
        ImGui::SliderInt("Point size", &pointSize, 1, 30);
        ImGui::SliderFloat("Line width", &lineWidth, 1, 16);
        ImGui::PopItemWidth();

        ImGui::TreePop();
        ImGui::Separator();
    }

    ImGui::Separator(); /////////////////////////////////////////

    static bool onlySaveFirstFrameMesh = true;
    static bool saveAccumulativeDisplacement = true;
    static bool saveAccumulativeDisplacement_relative = true;
    static bool saveIncrementalDisplacement = true;
    static bool saveIncrementalDisplacement_relative = true;
    static bool saveMarkerImage = true;
    static bool saveCellImage = false;
    static int cellChannel = -1;
    if (cellChannel < 0) {
        // make a guess
        if (channelPerSlice < 2) {
            cellChannel = 0;
            saveCellImage = false;
        } else if (channelToLoad == 0) {
            cellChannel = 1;
        } else {
            cellChannel = 0;
        }
    }

    // static bool saveEntireImage = false;
    if (ImGui::CollapsingHeader("Save & Export", ImGuiTreeNodeFlags_DefaultOpen)) {

        static std::string saveStr;
        static bool meshPointSaveFlag, meshCellSaveFlag, imageSaveFlag;
        if (ImGui::Button("Export VTU and TIFF")) {

            meshPointSaveFlag = true;
            meshCellSaveFlag = true;
            imageSaveFlag = true;

            // Save mesh to VTU (point data)
            meshPointSaveFlag = SaveMeshToVTU_point(onlySaveFirstFrameMesh, saveAccumulativeDisplacement, saveAccumulativeDisplacement_relative, saveIncrementalDisplacement, saveIncrementalDisplacement_relative);
            // Save mesh to VTU (cell data)
            meshCellSaveFlag = SaveMeshToVTU_cell(onlySaveFirstFrameMesh);
            // Save image to TIFF
            imageSaveFlag = SaveImageToTIFF(saveMarkerImage, saveCellImage, cellChannel);

            saveStr = (meshPointSaveFlag && meshCellSaveFlag && imageSaveFlag) ? "Data saved" : "Save failed";

            logger().debug("   <button> Export VTU & TIFF");
        }
        if (showTooltip && ImGui::IsItemHovered()) {
            ImGui::SetTooltip("Export the displacements and the gradient of the displacement. The cropped may also be saved as new TIFF images.");
        }
        ImGui::SameLine();
        ImGui::Text("%s", saveStr.c_str());

        if (ImGui::TreeNode("Advanced export")) {
        
            ImGui::Checkbox("Only save first frame mesh", &onlySaveFirstFrameMesh);
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("If checked, the exported VTU will save the first frame mesh location for all frames.\nOtherwise different frames will have their corresponding meshes.");
            }
            ImGui::Separator();
            ImGui::Checkbox("Save accumulative displacement (absolute)", &saveAccumulativeDisplacement);
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Save the displacement relative to the first frame");
            }
            ImGui::Checkbox("Save accumulative displacement (relative)", &saveAccumulativeDisplacement_relative);
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Save the displacement relative to the first frame. The mean of all the displacements in one frame will be subtracted.");
            }
            ImGui::Checkbox("Save incremental displacement (absolute)", &saveIncrementalDisplacement);
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Save the displacement relative to the previous frame");
            }
            ImGui::Checkbox("Save incremental displacement (relative)", &saveIncrementalDisplacement);
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("Save the displacement relative to the previous frame. The mean of all the displacements in one frame will be subtracted.");
            }
            ImGui::Separator();
            ImGui::Checkbox("Save cropped image (marker channel)", &saveMarkerImage);
            ImGui::Checkbox("Save cropped image (cell channel)", &saveCellImage);
            ImGui::SliderInt("Cell channel", &cellChannel, 0, channelPerSlice-1);
            if (showTooltip && ImGui::IsItemHovered()) {
                ImGui::SetTooltip("The channel with cell image.\nNote: the loaded channel is the one with marker information.");
            }

            /////////////////////////////////////

            if (ImGui::Button("SaveMeshToVTU (point data)")) {
                meshPointSaveFlag = SaveMeshToVTU_point(onlySaveFirstFrameMesh, saveAccumulativeDisplacement, saveAccumulativeDisplacement_relative, saveIncrementalDisplacement, saveIncrementalDisplacement_relative);
            }
            if (ImGui::Button("SaveMeshToVTU (cell data)")) {
                meshCellSaveFlag = SaveMeshToVTU_cell(onlySaveFirstFrameMesh);
            }
            if (ImGui::Button("SaveImageToTiff")) {
                imageSaveFlag = SaveImageToTIFF(saveMarkerImage, saveCellImage, cellChannel);
            }

            /////////////////////////////////////

            if (ImGui::Button("Export to OBJ")) {
                SaveMeshToOBJ();
            }
            if (ImGui::Button("Export displacement to TXT")) {
                SaveDisplacementToTXT();
            }

            ImGui::TreePop();
            ImGui::Separator();
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// Optimization


bool GUI::OptimizeAllFrames() {

    bool res = true;
    int currentFrame;

    for (currentFrame=0; currentFrame<currentLoadedFrames-1; currentFrame++) {

        OptimizeOneFrame(currentFrame);
        // depth correction for the frame that was just updated
        if (!MarkerDepthCorrection(currentFrame + 1, depthCorrectionNum, depthCorrectionGap)) {
            res = false;
            logger().warn("Depth correction failure: frame {}", currentFrame + 1);
        }
    }

    // update visualization variable
    UpdateMarkerPointLocArray();

    return res;
}


void GUI::OptimizeOneFrame(int prevFrameIdx) {

    /////////////////////////////////////////////////
    // prepare LBFGS
    LBFGSpp::LBFGSParam<double> param;
    param.epsilon = optimEpsilon;
    param.max_iterations = optimMaxIt;

    // prepare for parallel optimization
    const int N = markerArray[prevFrameIdx].num;
    markerArray[prevFrameIdx+1].num = N;

    logger().info("Optimization from frame {} to {}", prevFrameIdx, prevFrameIdx+1);
    logger().info("Optimization #starting points = {}", N);

    // Optimization
    logger().info(">>>>>>>>>> Before optimization >>>>>>>>>>");
    tbb::parallel_for( tbb::blocked_range<int>(0, N),
        //////////////////////////////////////
        // lambda function for parallel_for //
        //////////////////////////////////////
        [this/*.markerArray[?], &bsplineArray[?], .opticalFlowCorrection*/, &param, prevFrameIdx]
        (const tbb::blocked_range<int> &r) {

        // NOTE: LBFGSSolver is NOT thread safe. This must be instantiated for every thread
        LBFGSpp::LBFGSSolver<double> solver(param);

        // NOTE: the "variable count" used by "Autodiff" will be stored in 
        //       thread-local memory, so this must be set for every thread
        DiffScalarBase::setVariableCount(3);

        for (int ii = r.begin(); ii != r.end(); ++ii) {    

                Eigen::VectorXd vec(3, 1);
                vec(0) = markerArray[prevFrameIdx].loc(ii, 0) + opticalFlowCorrection[prevFrameIdx](ii, 0);  // x
                vec(1) = markerArray[prevFrameIdx].loc(ii, 1) + opticalFlowCorrection[prevFrameIdx](ii, 1);  // y
                vec(2) = markerArray[prevFrameIdx].loc(ii, 3);  // r
                double res;

                ///////////////////////////////////
                // lambda function for optimizer //
                ///////////////////////////////////
                auto func = [this/*.markerArray[?], .bsplineArray[?], .opticalFlowCorrection*/, ii, prevFrameIdx]
                (const Eigen::VectorXd& x, Eigen::VectorXd& grad) {

                        DScalar ans;

                        if (!cylinder::IsValid(bsplineArray[prevFrameIdx+1], x(0), x(1), markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2), x(2), 3)) {
                            grad.setZero();
                            return 1.0;
                        }
                        cylinder::EvaluateCylinder(bsplineArray[prevFrameIdx+1], DScalar(0, x(0)), DScalar(1, x(1)), markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2), DScalar(2, x(2)), 3, ans, reverseColor);
                        grad.resize(3, 1);
                        grad = ans.getGradient();
                        return ans.getValue();
                    };
                // NOTE: the template of "solver.minimize" does not accept a temprary variable (due to non-const argument)
                //       so we define a "func" and pass it in
                int it = solver.minimize(func, vec, res);
                ///////////////////////////////////

                markerArray[prevFrameIdx+1].loc(ii, 0) = vec(0);  // x
                markerArray[prevFrameIdx+1].loc(ii, 1) = vec(1);  // y
                markerArray[prevFrameIdx+1].loc(ii, 2) = markerArray[prevFrameIdx].loc(ii, 2) + opticalFlowCorrection[prevFrameIdx](ii, 2);  // z
                markerArray[prevFrameIdx+1].loc(ii, 3) = vec(2);  // r
                markerArray[prevFrameIdx+1].energy(ii) = res;     // energy
        }
    });
    logger().info("<<<<<<<<<< After optimization <<<<<<<<<<");
}


////////////////////////////////////////////////////////////////////////////////////////
// Export


bool GUI::SaveMeshToVTU_point(bool onlySaveFirstFrameMesh, bool saveAccumulativeDisplacement, bool saveAccumulativeDisplacement_relative, bool saveIncrementalDisplacement, bool saveIncrementalDisplacement_relative) {
// return true if successful

    static VTUWriter vtuWriter;

    for (int i=0; i<currentLoadedFrames; i++) {

        // prepare filename
        std::string vtuFileName;
        if (onlySaveFirstFrameMesh)
            vtuFileName = GetFileName(imagePath, i, ".vtu", "point");
        else
            vtuFileName = GetFileName(imagePath, i, ".vtu", "MovingMesh-point");
        // prepare absolute displacement
        Eigen::MatrixXd accumulativeDisplacement, incrementalDisplacement;
        GetDisplacement(markerPointLocArray, i, false, accumulativeDisplacement);
        GetDisplacement(markerPointLocArray, i, true, incrementalDisplacement);
        if (saveAccumulativeDisplacement) {
            vtuWriter.add_field("accumulative displacement (absolute)", accumulativeDisplacement);
        }
        if (saveIncrementalDisplacement) {
            vtuWriter.add_field("incremental displacement (absolute)", incrementalDisplacement);
        }
        // prepare relative displacement
        if (saveAccumulativeDisplacement_relative) {
            Eigen::RowVectorXd meanAccumulativeDisp = accumulativeDisplacement.colwise().mean();
            accumulativeDisplacement.rowwise() -= meanAccumulativeDisp;
            vtuWriter.add_field("accumulative displacement (relative)", accumulativeDisplacement);
        }
        if (saveIncrementalDisplacement_relative) {
            Eigen::RowVectorXd meanIncrementalDisp = incrementalDisplacement.colwise().mean();
            incrementalDisplacement.rowwise() -= meanIncrementalDisp;
            vtuWriter.add_field("incremental displacement (relative)", incrementalDisplacement);
        }

        // prepare mesh point array
        Eigen::MatrixXd meshPoint;
        if (onlySaveFirstFrameMesh)
            meshPoint = markerPointLocArray[0];
        else
            meshPoint = markerPointLocArray[i];
        meshPoint.col(2).array() -= layerBegin;

        // write to VTU
        if (!vtuWriter.write_tet_mesh(vtuFileName, meshPoint, markerMeshArray)) {
            return false;
        }
    }

    return true;
}


bool GUI::SaveMeshToVTU_cell(bool onlySaveFirstFrameMesh) {
// return true if successful

    static VTUWriter vtuWriter;

    for (int i=0; i<currentLoadedFrames; i++) {

        // prepare filename
        std::string vtuFileName;
        if (onlySaveFirstFrameMesh)
            vtuFileName = GetFileName(imagePath, i, ".vtu", "cell");
        else
            vtuFileName = GetFileName(imagePath, i, ".vtu", "MovingMesh-cell");

        Eigen::MatrixXd incrementalDisplacement;
        GetDisplacement(markerPointLocArray, i, true, incrementalDisplacement);
        // prepare Jacobian
        Eigen::MatrixXd jacobian;
        GetJacobian(markerPointLocArray[i], markerMeshArray, incrementalDisplacement, jacobian);
        vtuWriter.add_field("jacobian", jacobian);
        // prepare mesh point array
        Eigen::MatrixXd meshPoint;
        if (onlySaveFirstFrameMesh)
            meshPoint = markerPointLocArray[0];
        else
            meshPoint = markerPointLocArray[i];
        meshPoint.col(2).array() -= layerBegin;

        // write to VTU
        if (!vtuWriter.write_tet_mesh(vtuFileName, meshPoint, markerMeshArray)) {
            return false;
        }
    }

    return true;
}


bool GUI::SaveImageToTIFF(bool saveMarkerImage, bool saveCellImage, int cellChannel) {
// return true if successful

    imageData_t cellImgData;
    if (saveCellImage) {

        // only support loading one channel
        std::vector<bool> channelVec(channelPerSlice, false);
        channelVec[cellChannel] = true;
        // Read all desired frame to "imgData"
        ReadTif(imagePath, layerPerImg, channelVec, desiredFrames, cellImgData, r0, c0, r1, c1);
    }

    for (int i=0; i<currentLoadedFrames; i++) {

        if (saveMarkerImage) {

            std::string tifFileName = GetFileName(imagePath, i, ".tif", "image-marker");
            if (!WriteTif(tifFileName, imgData[i], layerBegin, layerEnd)) {
                return false;
            }
        }

        if (saveCellImage) {
            
            std::string tifFileName = GetFileName(imagePath, i, ".tif", "image-cell");
            if (!WriteTif(tifFileName, cellImgData[i], layerBegin, layerEnd)) {
                return false;
            }
        }
    }

    return true;
}


void GUI::SaveMeshToOBJ() {

    std::string objFileName;
    Eigen::MatrixXd meshPoint;

    for (int i=0; i<currentLoadedFrames; i++) {

        // prepare filename
        objFileName = GetFileName(imagePath, i, ".obj", "marker");

        // prepare mesh point array
        meshPoint = markerPointLocArray[i];
        meshPoint.col(2).array() -= layerBegin;

        // write to OBJ
        if (!igl::writeOBJ(objFileName, meshPoint, markerMeshArray)) {
            logger().warn("Write to OBJ failed");
        }
    }
}


void GUI::SaveDisplacementToTXT() {

    FILE *pFile;
    std::string fileName;
    Eigen::MatrixXd displacement;

    fileName = GetFileName(imagePath, 0, ".txt", "displacement");
    pFile = fopen(fileName.c_str(), "w+");

    for (int i=0; i<currentLoadedFrames; i++) {

        GetDisplacement(markerPointLocArray, i, false, displacement);

        fprintf(pFile, "## Frame %d\n", i);
        for (int r=0; r<displacement.rows(); r++) {
            for (int c=0; c<displacement.cols(); c++) {

                fprintf(pFile, "%f ", displacement(r, c));
            }
            fprintf(pFile, "\n");
        }

        fprintf(pFile, "\n");
    }

    fclose(pFile);
}

}  // namespace zebrafish
