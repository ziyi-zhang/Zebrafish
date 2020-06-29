#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/autodiff.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

class bspline {

private:
    Eigen::Matrix<double, 27, Eigen::Dynamic, Eigen::ColMajor> controlPointsCache; //TODO: Optional, probably remove me?

    Eigen::VectorXd controlPoints;  // [#points] control points value
    int Nx, Ny, Nz;         // the dimension of sample points (Nx * Ny * Nz == #pixels)
    int numX, numY, numZ;   // the dimension of control points (numX * numY * numZ == #control points)
    double resolutionX, resolutionY, resolutionZ;  // The distance between two pixels (in micrometers)
    int solverMaxIt;                  // Hypre solver max iterations
    double solverConvTol, solverTol;  // Hypre solver "convergence tolerance" and "tolerance"

    Eigen::Matrix< std::function<double(double)>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> basisXd, basisYd, basisZd;
    Eigen::Matrix<std::function<DScalar(DScalar)>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> basisX, basisY; //basisZ
        // Pre-calculated lambda basis functions (double & DScalar)

    void CalcLeastSquareMat(Eigen::SparseMatrix<double, Eigen::RowMajor> &A);
    template <typename T>
    void CalcBasisFunc(Eigen::Matrix< std::function<T(T)>, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> &basisT, const int& numT, const double& gapT);
    void CreateControlPtsCache(); //Optional

public:
    int degree;  // B-spline degree
                 // "2" for quadratic clamped B-spline
                 // "3" for cubic clamped B-spline
    double gapX, gapY, gapZ;  // the interval between two control points along a direction

    void SetResolution(const double resX, const double resY, const double resZ);
    /// Set the microscope resolution in X, Y and Z direction
    /// Unit in micrometers (um)

    void CalcControlPts(   const image_t &image, const double xratio, const double yratio, const double zratio, const int degree);
    void CalcControlPts_um(const image_t &image, const double distX,  const double distY,  const double distZ,  const int degree);
    /// Use least square to calculate an array of control points and store the result
    /// in private variables. This function must be called before any evaluation.
    ///
    /// @param[in]   xratio     { the ratio of (#control points) to (#sample points) along x-axis }
    /// @param[in]   yratio     { the ratio of (#control points) to (#sample points) along y-axis }
    /// @param[in]   zratio     { the ratio of (#control points) to (#sample points) along z-axis }
    /// @param[in]   distX      { the distance between two control points in X-axis. Unit: micrometer }
    /// @param[in]   distY      { the distance between two control points in Y-axis. Unit: micrometer }
    /// @param[in]   distZ      { the distance between two control points in Z-axis. Unit: micrometer }
    /// @param[in]   degree     { the degree of B-spline. Can be 2 or 3. }
    ///

    template <typename T>
    T Interp3D(const T &x, const T &y, const T &z) const;
    void Interp3D(const Eigen::Matrix<DScalar, Eigen::Dynamic, 2> &sample, const double z, Eigen::Matrix<DScalar, Eigen::Dynamic, 1> &res) const;
    void Interp3D(const Eigen::Matrix<double , Eigen::Dynamic, 2> &sample, const double z, Eigen::Matrix<double , Eigen::Dynamic, 1> &res) const;
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
