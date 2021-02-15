#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>
#include <igl/project.h>

#include <sstream>


namespace zebrafish {

namespace {

template <typename T>
std::string ToStringWithPrecision(const T a_value, const int n = 0)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

}  // anonymous namespace

void GUI::DrawText(
    Eigen::Vector3d pos, 
    const std::string &text,
    const Eigen::Vector4f color) {
// DO NOT USE THIS

    Eigen::Vector3f coord = igl::project(Eigen::Vector3f(pos.cast<float>()),
    viewer.core().view, viewer.core().proj, viewer.core().viewport);

    // Draw text labels slightly bigger than normal text
    double pixel_ratio_ = 2;
    ImDrawList* drawList = ImGui::GetWindowDrawList();
    drawList->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 1.2,
        ImVec2(coord[0]/pixel_ratio_, (viewer.core().viewport[3] - coord[1])/pixel_ratio_),
        ImGui::GetColorU32(ImVec4(
            color(0),
            color(1),
            color(2),
            color(3))),
        &text[0], &text[0] + text.size());
}

////////////////////////////////////////////////////////////////////////////////////////
// Special module about visualizing some effects

void GUI::DrawReferenceDots() {

    if (!show_refPoints) return;
    static Eigen::Matrix3d loc = (Eigen::Matrix3d() << 0, 0, 1, imgCols, imgRows, layerPerImg, imgCols-1, imgRows-1, layerPerImg-1).finished();
    static Eigen::MatrixXd referencePointColor = [] {
        Eigen::MatrixXd tmp(1, 3);
        tmp << 0.33, 0.83, 0.33;
        return tmp;
    } ();
    viewer.data(visualID).add_points(loc, referencePointColor);
}


void GUI::DrawAxisDots() {

    if (!show_axisPoints) return;

    const auto Linspace = [](double s, double e, double gap) -> Eigen::VectorXd {
        int n = std::ceil((e - s) / gap) + 1;
        Eigen::VectorXd res(n);
        for (int i=0; i<n; i++)
            res(i) = i * gap;
        res(n-1) = e;
        return res;
    };

    const auto ToMat = [](const Eigen::VectorXd &v1, const Eigen::VectorXd &v2, const Eigen::VectorXd &v3) -> Eigen::MatrixXd {
        Eigen::MatrixXd res(v1.rows(), 3);
        res << v1, v2, v3;
        return res;
    };

    static double gapr = (imgRows > 200) ? double(imgRows)/40.0 : 5;
    static double gapc = (imgCols > 200) ? double(imgCols)/40.0 : 5;
    static double gapz = (layerPerImg > 60) ? double(layerPerImg)/20.0 : 3;
    static int nr = std::ceil(double(imgRows) / gapr)+1;
    static int nc = std::ceil(double(imgCols) / gapc)+1;
    static int nz = std::ceil(double(layerPerImg) / gapz)+1;
    static Eigen::MatrixXd loc_r = ToMat(Eigen::VectorXd::Zero(nr), Linspace(0, imgRows, gapr), Eigen::VectorXd::Zero(nr));
    static Eigen::MatrixXd label_loc_r = ToMat(Eigen::VectorXd::Ones(nr).array() * (-2.0), Linspace(0, imgRows, gapr).array() + 0.5, Eigen::VectorXd::Zero(nr));
    static Eigen::MatrixXd loc_c = ToMat(Linspace(0, imgCols, gapc), Eigen::VectorXd::Ones(nc).array() * imgRows, Eigen::VectorXd::Zero(nc));
    static Eigen::MatrixXd label_loc_c = ToMat(Linspace(0, imgCols, gapc).array() - 0.5, Eigen::VectorXd::Ones(nc).array() * (imgRows+2.0), Eigen::VectorXd::Zero(nc));
    static Eigen::MatrixXd loc_z = ToMat(Eigen::VectorXd::Ones(nz).array() * imgCols, Eigen::VectorXd::Ones(nz).array() * imgRows, Linspace(0, layerPerImg, gapz));
    static Eigen::MatrixXd label_loc_z = ToMat(Eigen::VectorXd::Ones(nz).array() * (imgCols+1.0), Eigen::VectorXd::Ones(nz).array() * imgRows, Linspace(0, layerPerImg, gapz).array() + 1.2);

    static Eigen::MatrixXd referencePointColor = [] {
        Eigen::MatrixXd tmp(1, 3);
        tmp << 0.0, 0.0, 0.3;
        return tmp;
    } ();

    // add points
    viewer.data(visualID).add_points(loc_r, referencePointColor);
    viewer.data(visualID).add_points(loc_c, referencePointColor);
    viewer.data(visualID).add_points(loc_z, referencePointColor);

    // add labels
    for (int i=0; i<nr; i++)
        viewer.data(visualID).add_label(label_loc_r.row(i), ToStringWithPrecision(loc_r(i, 1)));
    for (int i=0; i<nc; i++)
        viewer.data(visualID).add_label(label_loc_c.row(i), ToStringWithPrecision(loc_c(i, 0)));
    for (int i=0; i<nz; i++)
        viewer.data(visualID).add_label(label_loc_z.row(i), ToStringWithPrecision(loc_z(i, 2)));
}


void GUI::ShowAllMarkerIndex() {

    if (markerArray.empty()) return;
    const markerRecord_t &markerRecord = markerArray[frameToShow];
    const int N = markerRecord.num;
    if (N == 0) return;

    viewer.data().show_labels = true;
    Eigen::Vector3d p;
    for (int i=0; i<N; i++) {
        p << 0.5+markerRecord.loc(i, 1), (imgRows-0.5)-markerRecord.loc(i, 0), markerRecord.loc(i, 2) + 0.003;
        viewer.data().add_label(p, std::to_string(i));
    }
}


void GUI::PlotBadDCPoints() {

    if (markerDepthCorrectionSuccess.size()-1 < frameToShow) return;
    const markerRecord_t &markerRecord = markerArray[frameToShow];
    const std::vector<int> &succ = markerDepthCorrectionSuccess[frameToShow];
    const int N = succ.size();
    if (N == 0) return;

    // colors
    Eigen::MatrixXd colors(4, 3);
    colors << 0, 0, 0, 
              0.8, 0.0, 0.1, 
              0.6, 0.7, 0.4, 
              0.6, 0.7, 1.0;

    viewer.data().show_labels = true;
    Eigen::Vector3d p;
    for (int i=0; i<N; i++) {
        if (succ[i] == 0) continue;  // if this point is successful
        p << 0.5+markerRecord.loc(i, 1), (imgRows-0.5)-markerRecord.loc(i, 0), markerRecord.loc(i, 2) + 0.003;
        // label
        viewer.data(visualID).add_label(p, std::to_string(i));
        // point
        viewer.data(visualID).add_points(p.transpose(), colors.row(succ[i]));
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// Draw 3D image & marker mesh


void GUI::Draw3DImage() {

    // do not draw if there is no image or this is not needed
    if (currentLoadedFrames == 0) return;
    if (!showBackgroundImage) return;

    static const Eigen::Vector3d ambient = Eigen::Vector3d(146./255., 172./255., 178./255.);
    static const Eigen::Vector3d diffuse = Eigen::Vector3d(146./255., 172./255., 178./255.);
    static const Eigen::Vector3d specular = Eigen::Vector3d(0., 0., 0.);
    static Eigen::MatrixXd UV = [] {
        Eigen::MatrixXd tmp(4, 2);
        tmp << 0, 1, 1, 1, 1, 0, 0, 0;
        return tmp;
    }();
    static int layerBegin_cache = -1, layerEnd_cache = -1;
    static int imgCols_cache = -1, imgRows_cache = -1;

    if (imageViewerType == 0) {
        // compressed view

        V_texture << 0, 0, 0, imgCols, 0, 0, imgCols, imgRows, 0, 0, imgRows, 0;
        if (layerBegin_cache != layerBegin || layerEnd_cache != layerEnd) {
            // if depth crop has updated
            // Note: only frame 0 in stage 1 could reach here

            const int num = imgData[0].size();
            // the user could provide false metadata
            try {
                ComputeCompressedTextureMax(imgData[0], 0);
            } catch (const std::exception &e) {
                logger().warn("Failed to draw 3D image. This is often due to wrong metadata like incorrect number of channels or slices.");
                std::cerr << "Failed to draw 3D image. This is often due to wrong metadata like incorrect number of channels or slices." << std::endl;
            }
            layerBegin_cache = layerBegin;
            layerEnd_cache = layerEnd;
        }
        texture = compressedImgTextureArray[frameToShow];
    } else {
        // per slice view

        V_texture << 0, 0, sliceToShow, 
             imgCols, 0, sliceToShow, 
             imgCols, imgRows, sliceToShow, 
             0, imgRows, sliceToShow;
        image_t &imgChosen = imgData[frameToShow];
        texture = (imgChosen[sliceToShow].array() * 255).cast<unsigned char>();
        texture.transposeInPlace();
    }

    if (stage <= 2) viewer.core().align_camera_center(V_texture);
    viewer.data().set_mesh(V_texture, F_texture);
    viewer.data().set_uv(UV);
    viewer.data().uniform_colors(ambient, diffuse, specular);
    viewer.data().show_faces = true;
    viewer.data().show_lines = false;
    viewer.data().show_texture = true;
    viewer.data().set_texture(texture, texture, texture);
}


void GUI::DrawMarkerMesh() {

    viewer.data(meshID).clear();

    // do not draw if not needed
    if (!showMarkerMesh) return;
    if (markerMeshArray.size() == 0 || markerPointLocArray.empty()) return;
    if (stage < 6) return;  // do not visualize mesh before stage 6: mesh

    // show mesh
    static Eigen::MatrixXd tempV;
    tempV = markerPointLocArray[frameToShow];
    tempV.col(2).array() -= 0.3;  // colliding with marker points
    viewer.data(meshID).show_faces = false;
    viewer.data(meshID).show_lines = true;
    viewer.data(meshID).show_texture = false;
    // viewer.data(meshID).set_colors(Eigen::RowVector3d(0.78, 0.82, 0.83));
    viewer.data(meshID).line_width = lineWidth;
    viewer.data(meshID).shininess = 0;
    viewer.data(meshID).set_mesh(
        tempV,  // mesh V
        markerMeshArray  // mesh F
    );

    if (!analysisInputPath.empty()) {
        viewer.core().align_camera_center(tempV);
    }
}


////////////////////////////////////////////////////////////////////////////////////////
// Update marker point array


void GUI::UpdateMarkerPointLocArray() {
// Critical marker location visualization array

    markerPointLocArray.resize(currentLoadedFrames);
    // only show the markers in one frame
    markerPointStatusArray.resize(currentLoadedFrames, 1);
    for (int i=0; i<markerPointStatusArray.rows(); i++)
        markerPointStatusArray(i) = false;
    markerPointStatusArray(frameToShow) = true;

    for (int i=0; i<currentLoadedFrames; i++) {

        markerRecord_t &markerRecord = markerArray[i];
        const int N = markerRecord.num;

        markerPointLocArray[i].resize(N, 3);
        markerPointLocArray[i].col(0) = markerRecord.loc.col(1).array() + 0.5;
        markerPointLocArray[i].col(1) = (imgRows-0.5) - markerRecord.loc.col(0).array();
        markerPointLocArray[i].col(2) = markerRecord.loc.col(2);
    }

    logger().trace("   [Vis] Marker point location array updated: all frames");
}


void GUI::UpdateAnalysisPointLocArray() {

    int N = analysisPara.V.size();
    markerPointLocArray.resize(N);
    // only show the markers in one frame
    markerPointStatusArray.resize(N, 1);
    for (int i=0; i<markerPointStatusArray.rows(); i++)
        markerPointStatusArray(i) = false;
    markerPointStatusArray(frameToShow) = true;

    for (int i=0; i<N; i++) {

        const int M = analysisPara.V[i].rows();

        markerPointLocArray[i].resize(M, 3);
        markerPointLocArray[i].col(0) = analysisPara.V[i].col(1);
        markerPointLocArray[i].col(1) = -analysisPara.V[i].col(0);
        markerPointLocArray[i].col(2) = analysisPara.V[i].col(2);
    }

    // to show marker mesh
    markerMeshArray = analysisPara.F;
    showMarkerMesh = true;
}

}  // namespace zebrafish
