#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/common.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

class bspline {

private:
    Eigen::VectorXd controlPoints;  // #points control points array
    int numX, numY, numZ;  // the dimension of control points (numX * numY * numZ == #points)
    double gapX, gapY, gapZ;  // the interval between two control points along a direction

public:

    void CalcControlPts(const image_t &image, const double xratio, const double yratio, const double zratio);
    /// Use least square to calculate an array of control points and store the result 
    /// in private variables. This function must be called before any evaluation.
    ///
    /// @param[in]   xratio     { the ratio of (#control points) to (#sample points) along x-axis }
    /// @param[in]   yratio     { the ratio of (#control points) to (#sample points) along y-axis }
    /// @param[in]   zratio     { the ratio of (#control points) to (#sample points) along z-axis }

    void Interp3D(const image_t &image, const Eigen::MatrixX3d &sample, Eigen::VectorXd &res);
    /// Calculate the interpolated B-spline result at "sample" points.
    /// Note: this function does not check for input validity
    ///
    /// @param[in]   image     { 3D matrix of the image }
    /// @param[in]   sample    { #sample x 3 input sample points }
    /// @param[out]  res       { #sample x 1 output interpolated results}
    ///

    // maintenance methods
    bspline();
    ~bspline();
};

}  // namespace zebrafish
