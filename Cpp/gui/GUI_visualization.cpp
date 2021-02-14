#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>


namespace zebrafish {

namespace {

}  // anonymous namespace

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
    viewer.data().add_points(loc, referencePointColor);
}


void GUI::ShowAllMarkerIndex() {

    const markerRecord_t &markerRecord = markerArray[frameToShow];
    const int N = markerRecord.num;
    if (N == 0) return;

    Eigen::Vector3d p;
    for (int i=0; i<N; i++) {
        p << markerRecord.loc(i, 0), markerRecord.loc(i, 1), markerRecord.loc(i, 2) + 0.003;
        viewer.data().add_label(p, std::to_string(i));
    }
}

}  // namespace zebrafish
