#include <zebrafish/ICP.h>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <limits>
#include <math.h>


namespace zebrafish {


void ICP::MatchPoints(const Eigen::MatrixXd &pt, const Eigen::MatrixXd &q, Eigen::VectorXi &matchIdx, Eigen::VectorXd &distSqVec) {
// Find the closest point in Q for every point in P
// brute force

    int i, j, minIdx;
    const int N = pt.cols(), M = q.cols();
    double minDist, dist;

    for (i=0; i<N; i++) {

        minDist = std::numeric_limits<double>::max();
        for (j=0; j<M; j++) {

            dist = (pt(0, i) - q(0, j))*(pt(0, i) - q(0, j)) + (pt(1, i) - q(1, j))*(pt(1, i) - q(1, j)) + (pt(2, i) - q(2, j))*(pt(2, i) - q(2, j));
            if (dist < minDist) {
                minIdx = j;
                minDist = dist;
            }
        }

        matchIdx(i) = minIdx;
        distSqVec(i) = minDist;
    }
}


void ICP::CalcTransMat(Eigen::MatrixXd pt, Eigen::MatrixXd qt, RMat_t &R, TMat_t &T) {

    Eigen::Vector3d pt_mean = pt.rowwise().mean();
    Eigen::Vector3d qt_mean = qt.rowwise().mean();

    // minus centroid
    pt.colwise() -= pt_mean;
    qt.colwise()  -= qt_mean;

    // SVD
    Eigen::MatrixXd covMat;
    covMat = qt * pt.transpose();
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(covMat, Eigen::ComputeFullV | Eigen::ComputeFullU);
    R = svd.matrixU() * svd.matrixV().transpose();
    T = qt_mean - R * pt_mean;
}


double ICP::RunICP(const Eigen::MatrixXd &p, const Eigen::MatrixXd &q, RMat_t &R_res, TMat_t &T_res) {

    assert(p.rows() == 3 && "P should be of size 3xN");
    assert(q.rows() == 3 && "Q should be of size 3xM");

    int iter, i;
    const int N = p.cols(), M = q.cols();
    double RMSe;
    Eigen::MatrixXd pt = p;  // store transformed p
    Eigen::MatrixXd qt;
    Eigen::VectorXi matchIdx;  // p[i] corresponds to q[matchIdx(i)]
    Eigen::VectorXd distSqVec;
    matchIdx.resize(N, 1);
    distSqVec.resize(N, 1);
    RMat_t R;
    TMat_t T;

    R_res << 1.0, 0.0, 0.0, 
             0.0, 1.0, 0.0, 
             0.0, 0.0, 1.0;
    T_res << 0.0, 0.0, 0.0;

    for (iter=0; iter<maxIt; iter++) {

        // Match closest point
        MatchPoints(pt, q, matchIdx, distSqVec);

        // Form a temporary q
        qt.resize(3, N);
        for (i=0; i<N; i++) {
            qt(0, i) = q(0, matchIdx(i));
            qt(1, i) = q(1, matchIdx(i));
            qt(2, i) = q(2, matchIdx(i));
        }

        // Minimize for R and T
        CalcTransMat(pt, qt, R, T);

        // inplace
        R_res = R * R_res;
        T_res = R * T_res + T;
        pt = R_res * p;
        pt.colwise() += T_res;

        // calculate error
        RMSe = std::sqrt( distSqVec.sum() / distSqVec.size() );
        printf("iter = %d   RMSE = %f\n", iter, RMSe);
    }

    return RMSe;
}


////////////////////////////////////////////////////////////////////////////////////////
// maintenance methods


ICP::ICP() : maxIt(10) {

}


ICP::~ICP() {

}

}  // namespace zebrafish
