#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef Eigen::Matrix<double, 3, 1> TMat_t;  // translation vector
typedef Eigen::Matrix<double, 3, 3> RMat_t;  // rotation matrix

///////////////////////////////////////
// Iterative Closest Point
// Reference: Martin Kjer and Jakob Wilm, Technical University of Denmark, 2012

class ICP {

private:
    int maxIt;  // max iteration

    void MatchPoints(const Eigen::MatrixXd &pt, const Eigen::MatrixXd &q, Eigen::VectorXi &matchIdx, Eigen::VectorXd &distSqVec);
    void CalcTransMat(Eigen::MatrixXd pt, Eigen::MatrixXd qt, RMat_t &R, TMat_t &T);

public:
    double RunICP(const Eigen::MatrixXd &p, const Eigen::MatrixXd &q, RMat_t &R, TMat_t &T);
    /// Run iterative closest point from point set p to q.
    /// Minimize the L2 error of sum(R*p + T - q).
    ///
    /// @param[in]   p     { [3 x #input] input point cloud }
    /// @param[in]   q     { [3 x #ref]   reference point cloud }
    /// @param[out]  R     { [3 x 3] rotation matrix }
    /// @param[out]  T     { [3 x 1] translation matrix }
    /// @return      { RMS error after ICP }
    ///

    // maintenance methods
    ICP();
    ~ICP();
};

}  // namespace zebrafish
