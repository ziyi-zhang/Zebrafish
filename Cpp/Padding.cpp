#include <zebrafish/Common.h>
#include <zebrafish/Padding.h>
#include <igl/harmonic.h>

#include <unordered_set>
#include <array>
#include <vector>
#include <algorithm>


namespace zebrafish {

typedef std::array<int, 2> rc_t;
typedef std::array<rc_t, 3> tri_t;

namespace {
}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// class padding

void padding::ComputeOneRing(
    const Eigen::MatrixXd &V, 
    const Eigen::MatrixXi &F, 
    const RCMap_t &RCMap, 
    Eigen::MatrixXd &appendV, 
    Eigen::MatrixXi &appendF,
    RCMap_t &appendRCMap) {

    std::unordered_set<int> vids, vids_new;
    RCMap_t RCMap_new;  // map newly created vertices to RC
    double avgEdgeLength = 0.0;

    const auto Smaller = [](const rc_t &x, const rc_t &y) -> bool {
        if (x[0] < y[0]) return true;
        if (x[0] == y[0] && x[1] < y[1]) return true;
        return false;
    };

    const auto FormTriangle = [&Smaller](const rc_t &rc0, const rc_t &rc1, const rc_t &rc2) -> tri_t {
        tri_t tri{rc0, rc1, rc2};
        for (int i=0; i<2; i++)
            for (int j=i; j<3; j++) {
                if (!Smaller(tri[i], tri[j])) {std::swap(tri[i], tri[j]);}
            }
        return tri;
    };

    const auto GetVID = [&vids, &vids_new, &RCMap, &RCMap_new, &V, &appendV, &avgEdgeLength](const rc_t &rc, int central_vid, int dir) -> int {
        // in old RCMap?
        for (int vid : vids) {
            if (RCMap.at(vid) == rc) {
                return vid;
            }
        }
        // in new appendV?
        for (int vid : vids_new) {
            if (RCMap_new.at(vid) == rc) {
                return vid;
            }
        }
        // create one and return the new id
        int N = appendV.rows();
        appendV.conservativeResize(N+1, 3);
        appendV(N, 2) = V(central_vid, 2);  // temporary z-value
        // Update New Vertex Loc
        // y(row)      5    4
        // ^        0  center  3
        // |           1    2
        // | ----> x(col)
        switch (dir) {
        case 0:
            appendV(N, 0) = V(central_vid, 0) - avgEdgeLength;
            appendV(N, 1) = V(central_vid, 1);
            break;
        case 1:
            appendV(N, 0) = V(central_vid, 0) - avgEdgeLength * 0.5;
            appendV(N, 1) = V(central_vid, 1) - avgEdgeLength * std::sqrt(3) / 2.0;
            break;
        case 2:
            appendV(N, 0) = V(central_vid, 0) + avgEdgeLength * 0.5;
            appendV(N, 1) = V(central_vid, 1) - avgEdgeLength * std::sqrt(3) / 2.0;
            break;
        case 3:
            appendV(N, 0) = V(central_vid, 0) + avgEdgeLength;
            appendV(N, 1) = V(central_vid, 1);
            break;
        case 4:
            appendV(N, 0) = V(central_vid, 0) + avgEdgeLength * 0.5;
            appendV(N, 1) = V(central_vid, 1) + avgEdgeLength * std::sqrt(3) / 2.0;
            break;
        case 5:
            appendV(N, 0) = V(central_vid, 0) - avgEdgeLength * 0.5;
            appendV(N, 1) = V(central_vid, 1) + avgEdgeLength * std::sqrt(3) / 2.0;
            break;
        default:
            break;
        }
        // Update vids_new & RCMap_new
        vids_new.insert(V.rows() + N);
        RCMap_new.insert({V.rows() + N, rc});

        return V.rows() + N;
    };


    //////////////////////////////////////////////////////////////////////////////

    const int old_v_cnt = V.rows();
    appendV.resize(0, V.cols());
    appendF.resize(0, 3);

    // get vids in V
    for (int i=0; i<F.rows(); i++)
        for (int j=0; j<F.cols(); j++) {
            vids.insert(F(i, j));
        }

    // estimate avgEdgeLength
    for (int i=0; i<F.rows(); i++)
        for (int j=0; j<F.cols(); j++) {
            avgEdgeLength += (V.row(F(i, j)) - V.row(F(i, (j+1)%3))).norm();
        }
    avgEdgeLength /= double(F.rows() * 3);

    // put known triangles to set
    std::vector<tri_t> tri_set;
    for (int i=0; i<F.rows(); i++) {
        tri_set.push_back(FormTriangle(RCMap.at(F(i, 0)), RCMap.at(F(i, 1)), RCMap.at(F(i, 2))));
    }

    // find ring
    for (int vid : vids) {
        
        rc_t rc0{RCMap.at(vid)[0], RCMap.at(vid)[1]};
        int r = rc0[0];
        int c = rc0[1];

        int col_correction = r & 1;  // for odd rows, the column index of adjacent dots needs to add 1
        int rc_loop[6][2]{
            // anti-clockwise
            {r, c-1},
            {r-1, c-1+col_correction},
            {r-1, c+col_correction},
            {r, c+1}, 
            {r+1, c+col_correction},
            {r+1, c-1+col_correction}
        };

        for (int i=0; i<6; i++) {
            // triangle rc0 x rc1 x rc2
            rc_t rc1{rc_loop[i][0], rc_loop[i][1]};
            rc_t rc2{rc_loop[(i+1)%6][0], rc_loop[(i+1)%6][1]};

            tri_t tri = FormTriangle(rc0, rc1, rc2);
            auto it = std::find(tri_set.begin(), tri_set.end(), tri);
            if (it != tri_set.end()) {
                continue;
            } else {
                // find a new triangle to insert
                int v1 = GetVID(rc1, vid, i);
                int v2 = GetVID(rc2, vid, (i+1)%6);

                // real update
                tri_set.push_back(tri);
                int appendF_cnt = appendF.rows();
                appendF.conservativeResize(appendF_cnt+1, 3);
                appendF.row(appendF_cnt) << vid, v1, v2;  // order matters
            }
        }
    }

    appendRCMap = RCMap_new;
}


template <typename T>
void padding::AddOneRing(const Eigen::MatrixXd &appendV, const Eigen::MatrixXi &appendF, T &V, Eigen::MatrixXi &F) {

    int N = V.rows();
    V.conservativeResize(V.rows() + appendV.rows(), V.cols());
    V.block(N, 0, appendV.rows(), V.cols()) = appendV;

    N = F.rows();
    F.conservativeResize(F.rows() + appendF.rows(), 3);
    F.block(N, 0, appendF.rows(), 3) = appendF;
}
// explicit template instantiation
template void padding::AddOneRing(const Eigen::MatrixXd &appendV, const Eigen::MatrixXi &appendF, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
template void padding::AddOneRing(const Eigen::MatrixXd &appendV, const Eigen::MatrixXi &appendF, Eigen::MatrixX4d &V, Eigen::MatrixXi &F);  // test purpose


void padding::AddOneRingForAll(const Eigen::MatrixXd &appendV, const Eigen::MatrixXi &appendF, const RCMap_t &appendRCMap, std::vector<Eigen::MatrixXd> &V, Eigen::MatrixXi &F, RCMap_t &RCMap) {

    const int N = V.size();
    const int M = F.rows();
    for (int i=0; i<N; i++) {
        padding::AddOneRing<Eigen::MatrixXd>(appendV, appendF, V[i], F);
    }
    RCMap.insert(appendRCMap.begin(), appendRCMap.end());

    // Fix F (we append it for N times)
    F.conservativeResize(M + appendF.rows(), 3);
}


void padding::Harmonic(Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int rawMeshVerts, int lastPadIndex, int order) {

    const int Nv = V.rows();
    const int Nv_rep = lastPadIndex - rawMeshVerts;  // #z-coordinates to be replaced
    Eigen::MatrixXd V_xy = V.block(0, 0, Nv, 2);  // we only need xy location for input
    Eigen::VectorXd Z, bc;
    // boundary
    Eigen::VectorXi b = Eigen::VectorXi::LinSpaced(Nv, 0, Nv-1);
    b.segment(rawMeshVerts, Nv-lastPadIndex) = b.tail(Nv-lastPadIndex);
    b.conservativeResize(Nv-Nv_rep, 1);
    // boundary condition
    bc.resize(b.rows());
    for (int i=0; i<b.rows(); i++) {
        int ii = b(i);
        bc(i) = V(ii, 2);
    }

    // harmonic
    igl::harmonic(V_xy, F, b, bc, order, Z);
    // replace z-coordinate
    V.block(rawMeshVerts, 2, Nv_rep, 1) = Z.segment(rawMeshVerts, Nv_rep);
    // assert
    double m = (V.block(0, 2, rawMeshVerts, 1) - Z.head(rawMeshVerts)).minCoeff();
    if (m != 0)
        std::cerr << "Z-coordinate not preversed " << m << std::endl;
}


void padding::HarmonicForAll(std::vector<Eigen::MatrixXd> &V, const Eigen::MatrixXi &F, int rawMeshVerts, int lastPadIndex, int order) {

    // set outmost ring to have z = mean_z
    double meanZ = V[0].block(0, 2, rawMeshVerts, 1).mean();
    V[0].block(lastPadIndex, 2, V[0].rows()-lastPadIndex, 1) = Eigen::VectorXd::Constant(V[0].rows()-lastPadIndex, meanZ);

    // harmonic
    padding::Harmonic(V[0], F, rawMeshVerts, lastPadIndex, order);

    // replace for all frames
    for (int i=1; i<V.size(); i++) {
        V[i].block(rawMeshVerts, 2, V[i].rows()-rawMeshVerts, 1) = V[0].block(rawMeshVerts, 2, V[0].rows()-rawMeshVerts, 1);
    }
}

}  // namespace zebrafish
