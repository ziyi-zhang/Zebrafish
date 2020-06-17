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
//
    Eigen::Matrix<DScalar, Eigen::Dynamic, 2> points;  // [#points x 2] point positions
    Eigen::Matrix<DScalar, Eigen::Dynamic, 1> weights;  // [#points] (Gaussian quadrature weight)*(energy function constants)*(subtraction function)
    Eigen::Matrix<double, Eigen::Dynamic, 1> zArray;   // [#depth] depth along z-axis
} samplePoints_t;

///////////////////////////////////////

class cylinder {

private:
    cylinder_t cyl;
    samplePoints_t samplePoints;

public:
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

    // maintenance methods
    cylinder();
    ~cylinder();

private:
    // hardcoded quadrature weights and locations
    // "locations" multiplied by sqrt(2)
    // "weight" multiplied by 1 or -1 depending on distance to center
        // static Eigen::Matrix<double, 57, 3> cools_kim_1;  // 57 samples, degree = 17
        // static Eigen::Matrix<double, 900, 3> lether;  // 900 samples, degree = 59
    Eigen::Matrix<double, Eigen::Dynamic, 2> xyArray;
    Eigen::Matrix<double, Eigen::Dynamic, 1> weightArray;
};

}  // namespace zebrafish
