#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>
#include <zebrafish/Bspline.h>

#include <Eigen/Dense>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef struct cylinder_t {
// uniquely defines a cylinder
    double x, y, z, r, h;  // bottom x(row), y(col), z + radius + height
} cylinder_t;


typedef struct samplePoints_t {
// The sample points coordinates [x, y, z] and their weights for one cylinder
    Eigen::Matrix<DScalar, Eigen::Dynamic, 2> innerPoints;  // [#points x 2] point positions
    Eigen::Matrix<DScalar, Eigen::Dynamic, 2> outerPoints;
    Eigen::Matrix<DScalar, Eigen::Dynamic, 1> weights;  // [#points] (Gaussian quadrature weight)*(energy function constants)*(subtraction function)
    Eigen::Matrix<double, Eigen::Dynamic, 1> zArray;   // [#depth] depth along z-axis
} samplePoints_t;

///////////////////////////////////////

//Quadrature contains:
// xy for disc + weight
// xy for ring + weight
// z quadarture + weight
//Quadrature quad(10); //10 is the order
//quad.xyrigh quad.xydisc

class cylinder {

private:
    int xmax, ymax, zmax;  // image size
    int numPts;            // size of sample points
    int heightLayers;      // desired height layers used in quadrature
    double minRadius;      // minimal radius of a cylinder

    // hardcoded quadrature weights and locations
    Eigen::Matrix<double, Eigen::Dynamic, 2> xyArray;
    Eigen::Matrix<double, Eigen::Dynamic, 1> weightArray;

    template<typename T>
    bool BoundaryCheck(const T &x, const T &y, const double z, const T &r, const double h) const;
    template<typename T>
    void EnergyHelper(const zebrafish::bspline &bsp, const Eigen::Matrix<double, Eigen::Dynamic, 1> &zArray,
                      const T &r, const T &x, const T &y, T &resT);

public:
    //energy(const Quadrature &quad, const BSpline &image, x, y, z, r, h) -> double
    //energyGrad(const Quadrature &quad, const BSpline &image, x, y, z, r, h) -> Vector3d (or double + Vector3d)
    //is_valid(const BSpline &image, x, y, z, r, h) -> boolean check for all of them, radius, and x,y,z inside image

    //cylinder::energy(...)

    void UpdateBoundary(const zebrafish::image_t &image);
    /// Record the size of the input image
    /// This will only be used to do boundary check

    template <typename T>
    bool EvaluateCylinder(const zebrafish::bspline &bsp, T x, T y, double z, T r, double h, T &res);
    /// Calculate sample points for the given cylinder and evaluate the energy.
    /// Return false if the cylinder is invalid.
    /// Must call "UpdateBoundary" before evaluating any cylinder
    ///
    /// @param[in]   bsp       { a B-spline solver with control points calculated }
    /// @param[in]   x, y, z   { coordinate of the cylinder bottom center }
    /// @param[in]   r         { cylinder radius. Must be positive }
    /// @param[in]   h         { cylinder height. Must be positive }
    /// @param[out]  res       { evaluated energy }
    /// @return      { whether the cylinder is inside the boundary of the image }
    ///

    void SubtractionHelper(const Eigen::MatrixXd &points, const Eigen::VectorXd &weight, Eigen::VectorXd &resWeight);
    /// [Deprecated]
    /// Multiply the weight array with a subtraction array such that the integral result
    /// measures the subtraction of two areas.
    /// The subtraction array can be hard or soft (sigmoid).
    ///
    /// @param[in]   points     { [#points x 2] sample locations used to calculate disk quadrature }
    /// @param[in]   weight     { [#points] weights corresponding to the sample points }
    /// @param[out]  resWeight  { [#points] returned weight array }
    ///

    void LoadQuadParas(int diskQuadMethod);
    /// Select a disk quadrature method and update the inner storage xy & weight array
    /// This can be changed on-the-fly

    // maintenance methods
    cylinder(int diskQuadMethod = 6);
    ~cylinder();
};

}  // namespace zebrafish
