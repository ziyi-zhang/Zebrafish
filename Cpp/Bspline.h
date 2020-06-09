#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/autodiff.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

class bspline {

private:
    Eigen::VectorXd controlPoints;  // #points control points array
    int numX, numY, numZ;  // the dimension of control points (numX * numY * numZ == #points)
    std::vector<double> centersX, centersY, centersZ;  // location of the center of a control point's basis function

public:
    double gapX, gapY, gapZ;  // the interval between two control points along a direction

    void CalcControlPts(const image_t &image, const double xratio, const double yratio, const double zratio);
    /// Use least square to calculate an array of control points and store the result 
    /// in private variables. This function must be called before any evaluation.
    ///
    /// @param[in]   xratio     { the ratio of (#control points) to (#sample points) along x-axis }
    /// @param[in]   yratio     { the ratio of (#control points) to (#sample points) along y-axis }
    /// @param[in]   zratio     { the ratio of (#control points) to (#sample points) along z-axis }

    void Interp3D(const Eigen::MatrixX3d &sample, Eigen::VectorXd &res) const;
    void Interp3D(const Eigen::Matrix<DScalar, Eigen::Dynamic, 3> &sample, Eigen::Matrix<DScalar, Eigen::Dynamic, 1> &res) const;
    /// Calculate the interpolated B-spline result at "sample" points.
    /// Note: this function does not check for input validity
    ///
    /// @param[in]   sample    { #sample x 3 input sample points }
    /// @param[out]  res       { #sample x 1 output interpolated results}
    ///

    // maintenance methods
    bspline();
    ~bspline();
};

}  // namespace zebrafish
