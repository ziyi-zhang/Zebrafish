#include <zebrafish/Common.h>
#include <zebrafish/Cylinder.h>
#include <cmath>

namespace zebrafish {

bool cylinder::SampleCylinder(const zebrafish::image_t &image, const zebrafish::bspline &bsp, 
                              double x_, double y_, double z_, double r_, double h_) {

    assert((x_==-1 && y_==-1 && z_==-1 && r_==-1 && h_==-1) || 
           (x_>0 && y_>0 && z_>0 && r_>0 && h_>0));

    // update cylinder
    if (x_ != -1) {
        cyl.x = x_;
        cyl.y = y_;
        cyl.z = z_;
        cyl.r = r_;
        cyl.h = h_;
    }

    // alias
    const double &x = cyl.x;
    const double &y = cyl.y;
    const double &z = cyl.z;
    const double &r = cyl.r;
    const double &h = cyl.h;
    const int &xmax = image[0].rows();
    const int &ymax = image[0].cols();
    const int &zmax = image.size();
    Eigen::MatrixX2d &points = samplePoints.points;
    Eigen::VectorXd &weights = samplePoints.weights;
    Eigen::VectorXd &zArray = samplePoints.zArray;

    // boundary check
    if (x - 1.5*r < 2*bsp.gapX || x + 1.5*r > xmax - 2*bsp.gapX ||
        y - 1.5*r < 2*bsp.gapY || y + 1.5*r > ymax - 2*bsp.gapY ||
        z < 2*bsp.gapZ || z+h > zmax - 2*bsp.gapZ || 
        r < minRadius) {
        return false;
    }

    // depth array
    double pad = h / (heightLayers * 2.0);
    zArray = Eigen::VectorXd::LinSpaced(heightLayers, z+pad, z+h-pad);
    
    // alias of location & weights (this is changeable)
    const Eigen::MatrixXd &xyArray = cools_kim_1.block<57, 2>(0, 0);
    const Eigen::VectorXd &weightArray = cools_kim_1.block<57, 1>(0, 2);

    // points xy (rc) array
    points.resize(xyArray.rows(), 2);  // Note: xyArray already multiplied by sqrt(2)
    points = xyArray.array() * r;  // * r
    Eigen::Vector2d xyVec;
    xyVec << x, y;
    points.rowwise() += xyVec.transpose();  // + [x, y], broadcast operation

    // weight
    weights.resize(weightArray.rows(), 1);
    double scalar;
    scalar = 2 * r * r;  // disk quadrature jacobian
                         // NOTE: this is not density function jacobian
    scalar /= r * r * h * M_PI;  // V_peri
    scalar *= h / double(heightLayers);
    weights = weightArray.array() * scalar;  // Note: weightArray already has subtraction function multiplied

    return true;
}


void cylinder::SubtractionHelper(const Eigen::MatrixXd &points, const Eigen::VectorXd &weight, Eigen::VectorXd &resWeight) {

    // this function can be used to implement sigmoid subtraction function
    // but it seems this is not necessary
}


double cylinder::EvaluateCylinder(const zebrafish::image_t &image, const zebrafish::bspline &bsp) {

    int depth, i;
    double res = 0;
    Eigen::MatrixX3d query;
    Eigen::VectorXd interpRes, temp;
    const int H = samplePoints.zArray.size();
    const int N = samplePoints.points.rows();

    query.resize(N, 3);
    query.block(0, 0, N, 2) = samplePoints.points;
    for (depth=0; depth<H; depth++) {

        temp.setConstant(N, samplePoints.zArray(depth));
        query.block(0, 2, N, 1) = temp;  // set z to current depth
        bsp.Interp3D(query, interpRes);  // interp this layer
        
        assert(samplePoints.weights.size() == interpRes.size());
        res += interpRes.dot(samplePoints.weights);
    }

    return res;
}


cylinder::cylinder() {

}


cylinder::~cylinder() {}


////////////////////////////////////////////////////
// hardcoded quadrature weights and locations

Eigen::Matrix<double, 57, 3> cylinder::cools_kim_1 = []{
/// 57 samples, degree = 17
    Eigen::Matrix<double, 57, 3> tmp;
    tmp << 0.000000, 0.000000, 0.114983, 1.254362, 0.000000, -0.042666, -1.254362, 0.000000, -0.042666, 0.000000, 1.254362, -0.042666, 0.000000, -1.254362, -0.042666, 0.600523, 0.000000, 0.087938, -0.600523, 0.000000, 0.087938, 0.000000, 0.600523, 0.087938, 0.000000, -0.600523, 0.087938, 0.982128, 0.000000, 0.076207, -0.982128, 0.000000, 0.076207, 0.000000, 0.982128, 0.076207, 0.000000, -0.982128, 0.076207, 0.972772, 0.972772, -0.019157, -0.972772, 0.972772, -0.019157, 0.972772, -0.972772, -0.019157, -0.972772, -0.972772, -0.019157, 0.843787, 0.843787, -0.062086, -0.843787, 0.843787, -0.062086, 0.843787, -0.843787, -0.062086, -0.843787, -0.843787, -0.062086, 0.333221, 0.333221, 0.095665, -0.333221, 0.333221, 0.095665, 0.333221, -0.333221, 0.095665, -0.333221, -0.333221, 0.095665, 0.442575, 0.776319, 0.085163, -0.442575, 0.776319, 0.085163, 0.442575, -0.776319, 0.085163, -0.442575, -0.776319, 0.085163, 0.776319, 0.442575, 0.085163, -0.776319, 0.442575, 0.085163, 0.776319, -0.442575, 0.085163, -0.776319, -0.442575, 0.085163, 1.359359, 0.245872, -0.020201, -1.359359, 0.245872, -0.020201, 1.359359, -0.245872, -0.020201, -1.359359, -0.245872, -0.020201, 0.245872, 1.359359, -0.020201, -0.245872, 1.359359, -0.020201, 0.245872, -1.359359, -0.020201, -0.245872, -1.359359, -0.020201, 0.431883, 1.117731, -0.056835, -0.431883, 1.117731, -0.056835, 0.431883, -1.117731, -0.056835, -0.431883, -1.117731, -0.056835, 1.117731, 0.431883, -0.056835, -1.117731, 0.431883, -0.056835, 1.117731, -0.431883, -0.056835, -1.117731, -0.431883, -0.056835, 1.201195, 0.654357, -0.024269, -1.201195, 0.654357, -0.024269, 1.201195, -0.654357, -0.024269, -1.201195, -0.654357, -0.024269, 0.654357, 1.201195, -0.024269, -0.654357, 1.201195, -0.024269, 0.654357, -1.201195, -0.024269, -0.654357, -1.201195, -0.024269;
    return tmp;
}();

}  // namespace zebrafish
