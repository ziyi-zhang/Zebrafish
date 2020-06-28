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

class cylinder {

private:
    int xmax, ymax, zmax;  // image size
    int numPts;            // size of sample points
    int heightLayers;      // desired height layers used in quadrature
    double minRadius;      // minimal radius of a cylinder

    // variables used when evaluating energy
    Eigen::Matrix<double, Eigen::Dynamic, 2> points;  // store inner points / outer points location
    Eigen::Matrix<double, Eigen::Dynamic, 1> interpRes;  // store the results of interpolation
    Eigen::Matrix<double, Eigen::Dynamic, 1> zArray;  // store the array of depths

    // hardcoded quadrature weights and locations
    Eigen::Matrix<double, Eigen::Dynamic, 2> xyArray;
    Eigen::Matrix<double, Eigen::Dynamic, 1> weightArray;

    bool BoundaryCheck(double x, double y, double z, double r, double h) const;
    void EnergyHelper(const zebrafish::bspline &bsp, double r, double x, double y, double &resT);

public:
    void UpdateBoundary(const zebrafish::image_t &image);
    /// Record the size of the input image
    /// This will only be used to do boundary check

    bool EvaluateCylinder(const zebrafish::bspline &bsp, double x, double y, double z, double r, double h, double &res);
/*
    bool SampleCylinder(const zebrafish::image_t &image, const zebrafish::bspline &bspline, 
                        double x = -1, double y = -1, double z = -1, double r = -1, double h = -1);
    /// Calculate the sample points (both interior and exterior point) and 
    /// save them to corresponding private variables
    /// 'image' is only used to do boundary check
    ///
    /// @param[in]   image     { 3D matrix of the image }
    /// @return      { whether the cylinder is inside the boundary of the image }
    ///

    void SubtractionHelper(const Eigen::MatrixXd &points, const Eigen::VectorXd &weight, Eigen::VectorXd &resWeight);
    /// Multiply the weight array with a subtraction array such that the integral result
    /// measures the subtraction of two areas.
    /// The subtraction array can be hard or soft (sigmoid).
    ///
    /// @param[in]   points     { [#points x 2] sample locations used to calculate disk quadrature }
    /// @param[in]   weight     { [#points] weights corresponding to the sample points }
    /// @param[out]  resWeight  { [#points] returned weight array }
    ///

    DScalar EvaluateCylinder(const zebrafish::image_t &image, const zebrafish::bspline &bspline);
    /// Calculate the energy for this cylinder
*/

    void LoadQuadParas(int diskQuadMethod);
    /// Select a disk quadrature method and update the inner storage xy & weight array
    /// This can be changed on-the-fly

    // maintenance methods
    cylinder(int diskQuadMethod = 6);
    ~cylinder();
};

}  // namespace zebrafish
