#include <zebrafish/Common.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/autodiff.h>
#include <cmath>

namespace zebrafish {

namespace {

// Generate a constant with corresponding type
template<typename T>
class GenConst {
public:
    T operator()(const double x) const {
        return T(x);
    }
};

template<>
class GenConst<double> {
public:
    double operator()(const double x) const {
        return x;
    }
};

// Get the value the variable
template<typename T>
class GetVal {
public:
    double operator()(const T x) const {
        return x.getValue();
    }
};

template<>
class GetVal<double> {
public:
    double operator()(const double x) const {
        return x;
    }
};

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
bool cylinder::BoundaryCheck(const T &x_, const T &y_, const double z, const T &r_, const double h) const {

    assert(h > 0);  // we do not optimize h. It should be manually set to be a positive number
    assert(xmax > 0 && ymax > 0 && zmax > 0);  // must update boundary before eval

    static const GetVal<T> getVal;
    static double x = getVal(x_);
    static double y = getVal(y_);
    static double r = getVal(r_);

    // The extended cylinder (union of inner & outer) has radius sqrt(2)*r
    // Also avoid interpolating points lying on the surface of the sample grid
    if (x - 1.5*r < 0 || x + 1.5*r > xmax - 1) return false;  // x-axis
    if (y - 1.5*r < 0 || y + 1.5*r > ymax - 1) return false;  // y-axis
    if (z < 0 || z+h >= zmax - 1) return false;  // depth-axis
    if (r < minRadius) return false;  // radius

    return true;
}
// explicit template instantiation
template bool cylinder::BoundaryCheck(const DScalar &x_, const DScalar &y_, const double z, const DScalar &r_, const double h) const;
template bool cylinder::BoundaryCheck(const double  &x_, const double  &y_, const double z, const double  &r_, const double h) const;

template <typename T>
//xy-point, xy-weight, z-point, z-weight
void cylinder::EnergyHelper(const zebrafish::bspline &bsp, const Eigen::Matrix<double, Eigen::Dynamic, 1> &zArray,
                            const T &r, const T &x, const T &y, T &resT)
{
    // This function calculates an intermediate value used by energy evaluation function
    // Note: "weightScalar" not multiplied here

    static const GenConst<T> genConst;
    static Eigen::Matrix<T, Eigen::Dynamic, 2> points;     // store inner points / outer points location
    static Eigen::Matrix<T, Eigen::Dynamic, 1> interpRes;  // store the results of interpolation
    int i, depth;
    //TODO remove static and move *r and +x/y to spline

    points.resize(numPts, 2);

    resT = T(0.0);
    for (i=0; i<numPts; i++) {
        // * r + [x, y]
        points(i, 0) = xyArray(i, 0) * r + x;
        points(i, 1) = xyArray(i, 1) * r + y;
    }
    for (depth=0; depth<heightLayers; depth++) {

        bsp.Interp3D(points, zArray(depth), interpRes);

        assert(weightArray.size() == interpRes.size());
        for (i=0; i<numPts; i++) {
            resT = resT + interpRes(i) * weightArray(i); //xy.weight * z.weight
        }
    }
}
// explicit template instantiation
template void cylinder::EnergyHelper(const zebrafish::bspline &bsp, const Eigen::Matrix<double, Eigen::Dynamic, 1> &zArray, const DScalar &r, const DScalar &x, const DScalar &y, DScalar &resT);
template void cylinder::EnergyHelper(const zebrafish::bspline &bsp, const Eigen::Matrix<double, Eigen::Dynamic, 1> &zArray, const double  &r, const double  &x, const double  &y, double  &resT);


void cylinder::UpdateBoundary(const zebrafish::image_t &image) {

    assert(image.size() > 0);

    xmax = image[0].rows();
    ymax = image[0].cols();
    zmax = image.size();
}


template <typename T>
bool cylinder::EvaluateCylinder(const zebrafish::bspline &bsp, T x, T y, double z, T r, double h, T &res) {

    static const GenConst<T> genConst;

    //TODO: drop this one and move this in the line search
    if (!BoundaryCheck(x, y, z, r, h)) {
        res = genConst(1);
        return false;
    }
    
    //kill the staic, zArray should be an input
    static Eigen::Matrix<double, Eigen::Dynamic, 1> zArray;  // store the array of depths
    static T resInner, resExt;
    assert(numPts > 0);  // Quadrature method has been chosen

    // depth array
    /// Note: This is simply a very small array with a few doubles
    double pad = h / (heightLayers * 2.0);
    zArray.resize(heightLayers, 1);
    zArray = Eigen::VectorXd::LinSpaced(heightLayers, z+pad, z+h-pad);
    //this is quadrature

    // weight correction term
    //////////////////////////////////////////////////////////////////////////////////
    /// Notes about the quadrature weight correction term "scalar":
    /// scalar = [disk quadrature jacobian] / [cylinder volumn] * [depth layer weight]
    /// For the inner cylinder with equidistant depth layer - 
    ///                r^2                 H
    ///     scalar = --------------- * ---------
    ///                pi * r^2 * H     #layers
    /// For the extended cylinder with equidistant depth layer - 
    ///                    (sqrt(2)*r)^2            H
    ///     scalar = ------------------------ * ---------
    ///               pi * (sqrt(2)*r)^2 * H     #layers
    /// Thus they share the same correction term.
    //////////////////////////////////////////////////////////////////////////////////
    // Note: Theoretically "weightScalar" should be DScalar, but the radius got
    // cancelled in the formula so it is OK to use double
    double weightScalar = 1.0 / (M_PI * double(heightLayers));

    // Inner area
    EnergyHelper<T>(bsp, zArray, r, x, y, resInner);

    // Extended area
    EnergyHelper<T>(bsp, zArray, r * sqrt(2), x, y, resExt);

    ////////////////////////////////////////////////////
    /// Notes about enerygy function evaluation:
    /// Energy = (Inner density) - (Periperal density)
    ///        = S_in / V_in - S_peri / V_peri
    /// [note: we enforced V_in == V_peri]
    ///        = (1 / V_in) * (S_in + S_in - S_in - S_peri)
    ///        = (1 / V_in) * (2*S_in - S_ext)
    ///        = 2 * (Inner density - Extended density)
    ////////////////////////////////////////////////////
    // double quadrature solution of the energy
    res = 2.0*(resInner - resExt) * weightScalar;
    return true;
}
//
template bool cylinder::EvaluateCylinder(const zebrafish::bspline &bsp, DScalar x, DScalar y, double z, DScalar r, double h, DScalar &res);
template bool cylinder::EvaluateCylinder(const zebrafish::bspline &bsp, double  x, double  y, double z, double  r, double h, double  &res);


void cylinder::SubtractionHelper(const Eigen::MatrixXd &points, const Eigen::VectorXd &weight, Eigen::VectorXd &resWeight) {

    // this function can be used to implement sigmoid subtraction function
    // [deprecated]
}


void cylinder::LoadQuadParas(int diskQuadMethod) {

    switch (diskQuadMethod) {
        #include <zebrafish/Quad.ipp>

        default:
            assert(false);
    }

    numPts = xyArray.rows();
}


cylinder::cylinder(int diskQuadMethod) {

    // Autodiff variables
    DiffScalarBase::setVariableCount(3);  // x, y, r
    // Quadrature method (by default use #6)
    LoadQuadParas(diskQuadMethod);
    // Boundary
    xmax = 0;
    ymax = 0;
    zmax = 0;
    // height layer
    heightLayers = 4;
    // minimal radius of a cylinder
    minRadius = 2.0;
}


cylinder::~cylinder() {

}

}  // namespace zebrafish
