#include <zebrafish/Common.h>
#include <zebrafish/Cylinder.h>
#include <zebrafish/autodiff.h>
#include <cmath>

namespace zebrafish {

bool cylinder::SampleCylinder(const zebrafish::image_t &image, const zebrafish::bspline &bsp, 
                              double x_, double y_, double z_, double r_, double h_) {

    int i;

    DiffScalarBase::setVariableCount(3);  // x, y, r

    // update cylinder
    if (x_ != -1) {
        cyl.x = x_;
        cyl.y = y_;
        cyl.z = z_;
        cyl.r = r_;
        cyl.h = h_;
    }

    // alias
    DScalar xDS = DScalar(0, cyl.x);
    DScalar yDS = DScalar(1, cyl.y);
    DScalar rDS = DScalar(2, cyl.r);
    const double x = cyl.x;
    const double y = cyl.y;
    const double z = cyl.z;
    const double r = cyl.r;
    const double h = cyl.h;
    const int &xmax = image[0].rows();
    const int &ymax = image[0].cols();
    const int &zmax = image.size();
    auto &points = samplePoints.points;
    auto &weights = samplePoints.weights;
    auto &zArray = samplePoints.zArray;

    // boundary check
    if (x - 1.5*r < 2*bsp.gapX || x + 1.5*r > xmax - 2*bsp.gapX ||
        y - 1.5*r < 2*bsp.gapY || y + 1.5*r > ymax - 2*bsp.gapY ||
        z < 2*bsp.gapZ || z+h > zmax - 2*bsp.gapZ || 
        r < minRadius) {
        return false;
    }

    // depth array
    double pad = h / (heightLayers * 2.0);
    zArray.resize(heightLayers, 1);
    zArray = Eigen::VectorXd::LinSpaced(heightLayers, z+pad, z+h-pad);

    // alias of location & weights (this is changeable)
      // const Eigen::MatrixXd &xyArray = cools_kim_1.block<57, 2>(0, 0);
      // const Eigen::VectorXd &weightArray = cools_kim_1.block<57, 1>(0, 2);
    const Eigen::MatrixXd &xyArray = lether.block<900, 2>(0, 0);
    const Eigen::VectorXd &weightArray = lether.block<900, 1>(0, 2).cwiseAbs();

    // points xy (rc) array
    points.resize(xyArray.rows(), 2);  // Note: xyArray already multiplied by sqrt(2)
    /// points = xyArray.array() * r;  // * r
    for (i=0; i<xyArray.rows(); i++) {
        points(i, 0) = xyArray(i, 0) * rDS;
        points(i, 1) = xyArray(i, 1) * rDS;
    }
    Eigen::Matrix<DScalar, 2, 1> xyVec;
    xyVec(0) = xDS;
    xyVec(1) = yDS;
    for (i=0; i<points.rows(); i++) {
        points(i, 0) += xDS;
        points(i, 1) += yDS;
    }
    // points.rowwise() += xyVec.transpose();  // + [x, y], broadcast operation

    // weight
    weights.resize(weightArray.rows(), 1);
    DScalar scalar;
    scalar = 2 * rDS * rDS;  // disk quadrature jacobian
                         // NOTE: this is not density function jacobian
    scalar = scalar / (rDS * rDS * h * M_PI);  // V_peri
    scalar *= h / double(heightLayers);
    /// weights = weightArray.array() * scalar;  // Note: weightArray already has subtraction function multiplied
    for (i=0; i<weights.rows(); i++)
        weights(i) = weightArray(i) * scalar;

    return true;
}


void cylinder::SubtractionHelper(const Eigen::MatrixXd &points, const Eigen::VectorXd &weight, Eigen::VectorXd &resWeight) {

    // this function can be used to implement sigmoid subtraction function
    // but it seems this is not necessary
}


DScalar cylinder::EvaluateCylinder(const zebrafish::image_t &image, const zebrafish::bspline &bsp) {

    int depth, i;
    DScalar res = DScalar(0);
    Eigen::Matrix<DScalar, Eigen::Dynamic, 3> query;
    Eigen::Matrix<DScalar, Eigen::Dynamic, 1> interpRes, temp;
    const int H = samplePoints.zArray.size();
    const int N = samplePoints.points.rows();

    query.resize(N, 3);
    temp.resize(N, 1);
    query.block(0, 0, N, 2) = samplePoints.points;
    for (depth=0; depth<H; depth++) {

        temp.setConstant(N, DScalar(samplePoints.zArray(depth)));
        query.block(0, 2, N, 1) = temp;  // set z to current depth
        bsp.Interp3D(query, interpRes);  // interp this layer

        assert(samplePoints.weights.size() == interpRes.size());
        ////// DEBUG
        // Eigen::VectorXd const_ = Eigen::VectorXd::Ones(900, 1).array() * 14.5;
        // Eigen::VectorXd error = query.col(0) + query.col(1) - interpRes;
        // Eigen::VectorXd error = (query.col(0)-const_)*(query.col(0)-const_) + (query.col(1)-const_)*(query.col(1)-const_) - interpRes;
        //std::cout << "error= " << error.maxCoeff() << std::endl;
        //std::cout << "depth=" << samplePoints.zArray(depth) << " RES=" << interpRes.dot(samplePoints.weights) << std::endl;
        ////// DEBUG

        /// res += interpRes.dot(samplePoints.weights);
        for (i=0; i<interpRes.rows(); i++)
            res = res + interpRes(i) * samplePoints.weights(i);
    }

    return res;
}


cylinder::cylinder() {

}


cylinder::~cylinder() {

}


////////////////////////////////////////////////////
// hardcoded quadrature weights and locations

Eigen::Matrix<double, 57, 3> cylinder::cools_kim_1 = []{
/// 57 samples, degree = 17
    Eigen::Matrix<double, 57, 3> tmp;
    tmp << 0.000000, 0.000000, 0.114983, 1.254362, 0.000000, -0.042666, -1.254362, 0.000000, -0.042666, 0.000000, 1.254362, -0.042666, 0.000000, -1.254362, -0.042666, 0.600523, 0.000000, 0.087938, -0.600523, 0.000000, 0.087938, 0.000000, 0.600523, 0.087938, 0.000000, -0.600523, 0.087938, 0.982128, 0.000000, 0.076207, -0.982128, 0.000000, 0.076207, 0.000000, 0.982128, 0.076207, 0.000000, -0.982128, 0.076207, 0.972772, 0.972772, -0.019157, -0.972772, 0.972772, -0.019157, 0.972772, -0.972772, -0.019157, -0.972772, -0.972772, -0.019157, 0.843787, 0.843787, -0.062086, -0.843787, 0.843787, -0.062086, 0.843787, -0.843787, -0.062086, -0.843787, -0.843787, -0.062086, 0.333221, 0.333221, 0.095665, -0.333221, 0.333221, 0.095665, 0.333221, -0.333221, 0.095665, -0.333221, -0.333221, 0.095665, 0.442575, 0.776319, 0.085163, -0.442575, 0.776319, 0.085163, 0.442575, -0.776319, 0.085163, -0.442575, -0.776319, 0.085163, 0.776319, 0.442575, 0.085163, -0.776319, 0.442575, 0.085163, 0.776319, -0.442575, 0.085163, -0.776319, -0.442575, 0.085163, 1.359359, 0.245872, -0.020201, -1.359359, 0.245872, -0.020201, 1.359359, -0.245872, -0.020201, -1.359359, -0.245872, -0.020201, 0.245872, 1.359359, -0.020201, -0.245872, 1.359359, -0.020201, 0.245872, -1.359359, -0.020201, -0.245872, -1.359359, -0.020201, 0.431883, 1.117731, -0.056835, -0.431883, 1.117731, -0.056835, 0.431883, -1.117731, -0.056835, -0.431883, -1.117731, -0.056835, 1.117731, 0.431883, -0.056835, -1.117731, 0.431883, -0.056835, 1.117731, -0.431883, -0.056835, -1.117731, -0.431883, -0.056835, 1.201195, 0.654357, -0.024269, -1.201195, 0.654357, -0.024269, 1.201195, -0.654357, -0.024269, -1.201195, -0.654357, -0.024269, 0.654357, 1.201195, -0.024269, -0.654357, 1.201195, -0.024269, 0.654357, -1.201195, -0.024269, -0.654357, -1.201195, -0.024269;
    return tmp;
}();

Eigen::Matrix<double, 900, 3> cylinder::lether = []{
/// 900 samples, degree = 59
    Eigen::Matrix<double, 900, 3> tmp;
    tmp << 1.406958, -0.142629, -0.000008, 1.385265, -0.283795, -0.000033, 1.349357, -0.422048, -0.000072, 1.299603, -0.555971, -0.000126, 1.236513, -0.684189, -0.000190, 1.160735, -0.805386, -0.000264, 1.073046, -0.918318, -0.000343, 0.974346, -1.021828, -0.000424, 0.865649, -1.114852, -0.000505, 0.748068, -1.196436, -0.000582, 0.622811, -1.265743, -0.000651, 0.491164, -1.322062, -0.000710, 0.354476, -1.364815, -0.000757, 0.214151, -1.393563, -0.000789, 0.071629, -1.408011, -0.000805, -0.071629, -1.408011, -0.000805, -0.214151, -1.393563, -0.000789, -0.354476, -1.364815, -0.000757, -0.491164, -1.322062, -0.000710, -0.622811, -1.265743, -0.000651, -0.748068, -1.196436, -0.000582, -0.865649, -1.114852, -0.000505, -0.974346, -1.021828, -0.000424, -1.073046, -0.918318, -0.000343, -1.160735, -0.805386, -0.000264, -1.236513, -0.684189, -0.000190, -1.299603, -0.555971, -0.000126, -1.349357, -0.422048, -0.000072, -1.385265, -0.283795, -0.000033, -1.406958, -0.142629, -0.000008, 1.406958, -0.140737, -0.000019, 1.385265, -0.280030, -0.000076, 1.349357, -0.416449, -0.000168, 1.299603, -0.548595, -0.000291, 1.236513, -0.675112, -0.000441, 1.160735, -0.794701, -0.000611, 1.073046, -0.906135, -0.000794, 0.974346, -1.008271, -0.000983, 0.865649, -1.100061, -0.001170, 0.748068, -1.180563, -0.001348, 0.622811, -1.248951, -0.001508, 0.491164, -1.304523, -0.001646, 0.354476, -1.346708, -0.001754, 0.214151, -1.375075, -0.001829, 0.071629, -1.389331, -0.001867, -0.071629, -1.389331, -0.001867, -0.214151, -1.375075, -0.001829, -0.354476, -1.346708, -0.001754, -0.491164, -1.304523, -0.001646, -0.622811, -1.248951, -0.001508, -0.748068, -1.180563, -0.001348, -0.865649, -1.100061, -0.001170, -0.974346, -1.008271, -0.000983, -1.073046, -0.906135, -0.000794, -1.160735, -0.794701, -0.000611, -1.236513, -0.675112, -0.000441, -1.299603, -0.548595, -0.000291, -1.349357, -0.416449, -0.000168, -1.385265, -0.280030, -0.000076, -1.406958, -0.140737, -0.000019, 1.406958, -0.137354, -0.000030, 1.385265, -0.273298, -0.000118, 1.349357, -0.406438, -0.000261, 1.299603, -0.535407, -0.000454, 1.236513, -0.658883, -0.000687, 1.160735, -0.775597, -0.000952, 1.073046, -0.884353, -0.001238, 0.974346, -0.984034, -0.001532, 0.865649, -1.073617, -0.001824, 0.748068, -1.152184, -0.002101, 0.622811, -1.218928, -0.002351, 0.491164, -1.273164, -0.002565, 0.354476, -1.314335, -0.002734, 0.214151, -1.342020, -0.002850, 0.071629, -1.355933, -0.002910, -0.071629, -1.355933, -0.002910, -0.214151, -1.342020, -0.002850, -0.354476, -1.314335, -0.002734, -0.491164, -1.273164, -0.002565, -0.622811, -1.218928, -0.002351, -0.748068, -1.152184, -0.002101, -0.865649, -1.073617, -0.001824, -0.974346, -0.984034, -0.001532, -1.073046, -0.884353, -0.001238, -1.160735, -0.775597, -0.000952, -1.236513, -0.658883, -0.000687, -1.299603, -0.535407, -0.000454, -1.349357, -0.406438, -0.000261, -1.385265, -0.273298, -0.000118, -1.406958, -0.137354, -0.000030, 1.406958, -0.132515, -0.000040, 1.385265, -0.263670, -0.000159, 1.349357, -0.392119, -0.000352, 1.299603, -0.516545, -0.000611, 1.236513, -0.635670, -0.000926, 1.160735, -0.748273, -0.001283, 1.073046, -0.853197, -0.001668, 0.974346, -0.949366, -0.002066, 0.865649, -1.035793, -0.002459, 0.748068, -1.111592, -0.002832, 0.622811, -1.175984, -0.003169, 0.491164, -1.228310, -0.003458, 0.354476, -1.268031, -0.003685, 0.214151, -1.294740, -0.003842, 0.071629, -1.308163, -0.003922, -0.071629, -1.308163, -0.003922, -0.214151, -1.294740, -0.003842, -0.354476, -1.268031, -0.003685, -0.491164, -1.228310, -0.003458, -0.622811, -1.175984, -0.003169, -0.748068, -1.111592, -0.002832, -0.865649, -1.035793, -0.002459, -0.974346, -0.949366, -0.002066, -1.073046, -0.853197, -0.001668, -1.160735, -0.748273, -0.001283, -1.236513, -0.635670, -0.000926, -1.299603, -0.516545, -0.000611, -1.349357, -0.392119, -0.000352, -1.385265, -0.263670, -0.000159, -1.406958, -0.132515, -0.000040, 1.406958, -0.126271, -0.000050, 1.385265, -0.251247, -0.000199, 1.349357, -0.373644, -0.000440, 1.299603, -0.492207, -0.000763, 1.236513, -0.605719, -0.001155, 1.160735, -0.713016, -0.001601, 1.073046, -0.812997, -0.002081, 0.974346, -0.904635, -0.002577, 0.865649, -0.986990, -0.003067, 0.748068, -1.059218, -0.003533, 0.622811, -1.120576, -0.003954, 0.491164, -1.170436, -0.004314, 0.354476, -1.208285, -0.004597, 0.214151, -1.233736, -0.004793, 0.071629, -1.246527, -0.004893, -0.071629, -1.246527, -0.004893, -0.214151, -1.233736, -0.004793, -0.354476, -1.208285, -0.004597, -0.491164, -1.170436, -0.004314, -0.622811, -1.120576, -0.003954, -0.748068, -1.059218, -0.003533, -0.865649, -0.986990, -0.003067, -0.974346, -0.904635, -0.002577, -1.073046, -0.812997, -0.002081, -1.160735, -0.713016, -0.001601, -1.236513, -0.605719, -0.001155, -1.299603, -0.492207, -0.000763, -1.349357, -0.373644, -0.000440, -1.385265, -0.251247, -0.000199, -1.406958, -0.126271, -0.000050, 1.406958, -0.118689, -0.000060, 1.385265, -0.236160, -0.000236, 1.349357, -0.351208, -0.000522, 1.299603, -0.462652, -0.000906, 1.236513, -0.569348, -0.001372, 1.160735, -0.670202, -0.001901, 1.073046, -0.764179, -0.002472, 0.974346, -0.850315, -0.003061, 0.865649, -0.927725, -0.003643, 0.748068, -0.995615, -0.004196, 0.622811, -1.053289, -0.004696, 0.491164, -1.100155, -0.005124, 0.354476, -1.135732, -0.005460, 0.214151, -1.159654, -0.005693, 0.071629, -1.171677, -0.005812, -0.071629, -1.171677, -0.005812, -0.214151, -1.159654, 
    -0.005693, -0.354476, -1.135732, -0.005460, -0.491164, -1.100155, -0.005124, -0.622811, -1.053289, -0.004696, -0.748068, -0.995615, -0.004196, -0.865649, -0.927725, -0.003643, -0.974346, -0.850315, -0.003061, -1.073046, -0.764179, -0.002472, -1.160735, -0.670202, -0.001901, -1.236513, -0.569348, -0.001372, -1.299603, -0.462652, -0.000906, -1.349357, -0.351208, -0.000522, -1.385265, -0.236160, -0.000236, -1.406958, -0.118689, -0.000060, 1.406958, -0.109849, -0.000068, 1.385265, -0.218570, -0.000271, 1.349357, -0.325049, -0.000599, 1.299603, -0.428192, -0.001040, 1.236513, -0.526941, -0.001575, 1.160735, -0.620284, -0.002182, 1.073046, -0.707261, -0.002837, 0.974346, -0.786981, -0.003512, 0.865649, -0.858625, -0.004181, 0.748068, -0.921459, -0.004815, 0.622811, -0.974837, -0.005389, 0.491164, -1.018212, -0.005879, 0.354476, -1.051139, -0.006266, 0.214151, -1.073280, -0.006533, 0.071629, -1.084408, -0.006669, -0.071629, -1.084408, -0.006669, -0.214151, -1.073280, -0.006533, -0.354476, -1.051139, -0.006266, -0.491164, -1.018212, -0.005879, -0.622811, -0.974837, -0.005389, -0.748068, -0.921459, -0.004815, -0.865649, -0.858625, -0.004181, -0.974346, -0.786981, -0.003512, -1.073046, -0.707261, -0.002837, -1.160735, -0.620284, -0.002182, -1.236513, -0.526941, -0.001575, -1.299603, -0.428192, -0.001040, -1.349357, -0.325049, -0.000599, -1.385265, -0.218570, -0.000271, -1.406958, -0.109849, -0.000068, 1.406958, -0.099844, -0.000077, 1.385265, -0.198663, -0.000303, 1.349357, -0.295444, -0.000670, 1.299603, -0.389194, -0.001162, 1.236513, -0.478949, -0.001760, 1.160735, -0.563790, -0.002439, 1.073046, -0.642846, -0.003171, 0.974346, -0.715305, -0.003927, 0.865649, -0.780424, -0.004674, 0.748068, -0.837535, -0.005383, 0.622811, -0.886052, -0.006025, 0.491164, -0.925477, -0.006573, 0.354476, -0.955405, -0.007005, 0.214151, -0.975529, 0.007303, 0.071629, -0.985643, 0.007455, -0.071629, -0.985643, 0.007455, -0.214151, -0.975529, 0.007303, -0.354476, -0.955405, -0.007005, -0.491164, -0.925477, -0.006573, -0.622811, -0.886052, -0.006025, -0.748068, -0.837535, -0.005383, -0.865649, -0.780424, -0.004674, -0.974346, -0.715305, -0.003927, -1.073046, -0.642846, -0.003171, -1.160735, -0.563790, -0.002439, -1.236513, -0.478949, -0.001760, -1.299603, -0.389194, -0.001162, -1.349357, -0.295444, -0.000670, -1.385265, -0.198663, -0.000303, -1.406958, -0.099844, -0.000077, 1.406958, -0.088781, -0.000084, 1.385265, -0.176651, -0.000332, 1.349357, -0.262708, -0.000733, 1.299603, -0.346070, -0.001273, 1.236513, -0.425880, -0.001927, 1.160735, -0.501320, -0.002671, 1.073046, -0.571616, -0.003472, 0.974346, -0.636047, -0.004299, 0.865649, -0.693950, -0.005118, 0.748068, -0.744733, -0.005894, 0.622811, -0.787874, -0.006597, 0.491164, -0.822931, 0.007197, 0.354476, -0.849542, 0.007670, 0.214151, -0.867437, 0.007996, 0.071629, -0.876430, 0.008163, -0.071629, -0.876430, 0.008163, -0.214151, -0.867437, 0.007996, -0.354476, -0.849542, 0.007670, -0.491164, -0.822931, 0.007197, -0.622811, -0.787874, -0.006597, -0.748068, -0.744733, -0.005894, -0.865649, -0.693950, -0.005118, -0.974346, -0.636047, -0.004299, -1.073046, -0.571616, -0.003472, -1.160735, -0.501320, -0.002671, -1.236513, -0.425880, -0.001927, -1.299603, -0.346070, -0.001273, -1.349357, -0.262708, -0.000733, -1.385265, -0.176651, -0.000332, -1.406958, -0.088781, -0.000084, 1.406958, -0.076777, -0.000090, 1.385265, -0.152766, -0.000357, 1.349357, -0.227187, -0.000789, 1.299603, -0.299277, -0.001370, 1.236513, -0.368296, -0.002074, 1.160735, -0.433536, -0.002874, 1.073046, -0.494327, -0.003737, 0.974346, -0.550046, -0.004626, 0.865649, -0.600121, -0.005507, 0.748068, -0.644037, 0.006342, 0.622811, -0.681345, 0.007099, 0.491164, -0.711661, 0.007744, 0.354476, -0.734675, 0.008253, 0.214151, -0.750150, 0.008605, 0.071629, -0.757927, 0.008784, -0.071629, -0.757927, 0.008784, -0.214151, -0.750150, 0.008605, -0.354476, -0.734675, 0.008253, -0.491164, -0.711661, 0.007744, -0.622811, -0.681345, 0.007099, -0.748068, -0.644037, 0.006342, -0.865649, -0.600121, -0.005507, -0.974346, -0.550046, -0.004626, -1.073046, -0.494327, -0.003737, -1.160735, -0.433536, -0.002874, -1.236513, -0.368296, -0.002074, -1.299603, -0.299277, -0.001370, -1.349357, -0.227187, -0.000789, -1.385265, -0.152766, -0.000357, -1.406958, -0.076777, -0.000090, 1.406958, -0.063959, -0.000096, 1.385265, -0.127261, -0.000378, 1.349357, -0.189258, -0.000837, 1.299603, -0.249312, -0.001452, 1.236513, -0.306808, -0.002199, 1.160735, -0.361156, -0.003047, 1.073046, -0.411798, -0.003961, 0.974346, -0.458215, -0.004904, 0.865649, -0.499929, 0.005838, 0.748068, -0.536514, 0.006724, 0.622811, -0.567593, 0.007525, 0.491164, -0.592848, 0.008210, 0.354476, -0.612020, 0.008749, 0.214151, -0.624911, 0.009122, 0.071629, -0.631390, 0.009312, -0.071629, -0.631390, 0.009312, -0.214151, -0.624911, 0.009122, -0.354476, -0.612020, 0.008749, -0.491164, -0.592848, 0.008210, -0.622811, -0.567593, 0.007525, -0.748068, -0.536514, 0.006724, -0.865649, -0.499929, 0.005838, -0.974346, -0.458215, -0.004904, -1.073046, -0.411798, -0.003961, -1.160735, -0.361156, -0.003047, -1.236513, -0.306808, -0.002199, -1.299603, -0.249312, -0.001452, -1.349357, -0.189258, -0.000837, -1.385265, -0.127261, -0.000378, -1.406958, -0.063959, -0.000096, 1.406958, -0.050463, -0.000100, 1.385265, -0.100408, -0.000396, 1.349357, -0.149322, -0.000875, 1.299603, 
    -0.196705, -0.001519, 1.236513, -0.242069, -0.002300, 1.160735, -0.284948, -0.003187, 1.073046, -0.324904, -0.004144, 0.974346, -0.361527, -0.005130, 0.865649, -0.394439, 0.006107, 0.748068, -0.423304, 0.007034, 0.622811, -0.447825, 0.007872, 0.491164, -0.467751, 0.008588, 0.354476, -0.482877, 0.009153, 0.214151, -0.493048, 0.009542, 0.071629, -0.498160, 0.009741, -0.071629, -0.498160, 0.009741, -0.214151, -0.493048, 0.009542, -0.354476, -0.482877, 0.009153, -0.491164, -0.467751, 0.008588, -0.622811, -0.447825, 0.007872, -0.748068, -0.423304, 0.007034, -0.865649, -0.394439, 0.006107, -0.974346, -0.361527, -0.005130, -1.073046, -0.324904, -0.004144, -1.160735, -0.284948, -0.003187, -1.236513, -0.242069, -0.002300, -1.299603, -0.196705, -0.001519, -1.349357, -0.149322, -0.000875, -1.385265, -0.100408, -0.000396, -1.406958, -0.050463, -0.000100, 1.406958, -0.036432, -0.000103, 1.385265, -0.072490, -0.000409, 1.349357, -0.107804, -0.000905, 1.299603, -0.142012, -0.001570, 1.236513, -0.174763, -0.002377, 1.160735, -0.205720, -0.003294, 1.073046, -0.234566, -0.004282, 0.974346, -0.261006, -0.005302, 0.865649, -0.284767, 0.006311, 0.748068, -0.305606, 0.007269, 0.622811, -0.323309, 0.008135, 0.491164, -0.337695, 0.008876, 0.354476, -0.348615, 0.009459, 0.214151, -0.355958, 0.009862, 0.071629, -0.359649, 0.010067, -0.071629, -0.359649, 0.010067, -0.214151, -0.355958, 0.009862, -0.354476, -0.348615, 0.009459, -0.491164, -0.337695, 0.008876, -0.622811, -0.323309, 0.008135, -0.748068, -0.305606, 0.007269, -0.865649, -0.284767, 0.006311, -0.974346, -0.261006, -0.005302, -1.073046, -0.234566, -0.004282, -1.160735, -0.205720, -0.003294, -1.236513, -0.174763, -0.002377, -1.299603, -0.142012, -0.001570, -1.349357, -0.107804, -0.000905, -1.385265, -0.072490, -0.000409, -1.406958, -0.036432, -0.000103, 1.406958, -0.022015, -0.000106, 1.385265, -0.043804, -0.000418, 1.349357, -0.065143, -0.000924, 1.299603, -0.085814, -0.001604, 1.236513, -0.105604, -0.002429, 1.160735, -0.124311, -0.003366, 1.073046, -0.141742, -0.004376, 0.974346, -0.157718, 0.005418, 0.865649, -0.172077, 0.006449, 0.748068, -0.184669, 0.007427, 0.622811, -0.195367, 0.008313, 0.491164, -0.204059, 0.009069, 0.354476, -0.210658, 0.009665, 0.214151, -0.215096, 0.010076, 0.071629, -0.217326, 0.010286, -0.071629, -0.217326, 0.010286, -0.214151, -0.215096, 0.010076, -0.354476, -0.210658, 0.009665, -0.491164, -0.204059, 0.009069, -0.622811, -0.195367, 0.008313, -0.748068, -0.184669, 0.007427, -0.865649, -0.172077, 0.006449, -0.974346, -0.157718, 0.005418, -1.073046, -0.141742, -0.004376, -1.160735, -0.124311, -0.003366, -1.236513, -0.105604, -0.002429, -1.299603, -0.085814, -0.001604, -1.349357, -0.065143, -0.000924, -1.385265, -0.043804, -0.000418, -1.406958, -0.022015, -0.000106, 1.406958, -0.007364, -0.000107, 1.385265, -0.014653, -0.000422, 1.349357, -0.021791, -0.000934, 1.299603, -0.028706, -0.001621, 1.236513, -0.035326, -0.002455, 1.160735, -0.041584, -0.003402, 1.073046, -0.047415, -0.004422, 0.974346, -0.052759, 0.005476, 0.865649, -0.057562, 0.006518, 0.748068, -0.061775, 0.007507, 0.622811, -0.065353, 0.008402, 0.491164, -0.068261, 0.009166, 0.354476, -0.070468, 0.009768, 0.214151, -0.071953, 0.010184, 0.071629, -0.072699, 0.010397, -0.071629, -0.072699, 0.010397, -0.214151, -0.071953, 0.010184, -0.354476, -0.070468, 0.009768, -0.491164, -0.068261, 0.009166, -0.622811, -0.065353, 0.008402, -0.748068, -0.061775, 0.007507, -0.865649, -0.057562, 0.006518, -0.974346, -0.052759, 0.005476, -1.073046, -0.047415, -0.004422, -1.160735, -0.041584, -0.003402, -1.236513, -0.035326, -0.002455, -1.299603, -0.028706, -0.001621, -1.349357, -0.021791, -0.000934, -1.385265, -0.014653, -0.000422, -1.406958, -0.007364, -0.000107, 1.406958, 0.007364, -0.000107, 1.385265, 0.014653, -0.000422, 1.349357, 0.021791, -0.000934, 1.299603, 0.028706, -0.001621, 1.236513, 0.035326, -0.002455, 1.160735, 0.041584, -0.003402, 1.073046, 0.047415, -0.004422, 0.974346, 0.052759, 0.005476, 0.865649, 0.057562, 0.006518, 0.748068, 0.061775, 0.007507, 0.622811, 0.065353, 0.008402, 0.491164, 0.068261, 0.009166, 0.354476, 0.070468, 0.009768, 0.214151, 0.071953, 0.010184, 0.071629, 0.072699, 0.010397, -0.071629, 0.072699, 0.010397, -0.214151, 0.071953, 0.010184, -0.354476, 0.070468, 0.009768, -0.491164, 0.068261, 0.009166, -0.622811, 0.065353, 0.008402, -0.748068, 0.061775, 0.007507, -0.865649, 0.057562, 0.006518, -0.974346, 0.052759, 0.005476, -1.073046, 0.047415, -0.004422, -1.160735, 0.041584, -0.003402, -1.236513, 0.035326, -0.002455, -1.299603, 0.028706, -0.001621, -1.349357, 0.021791, -0.000934, -1.385265, 0.014653, -0.000422, -1.406958, 0.007364, -0.000107, 1.406958, 0.022015, -0.000106, 1.385265, 0.043804, -0.000418, 1.349357, 0.065143, -0.000924, 1.299603, 0.085814, -0.001604, 1.236513, 0.105604, -0.002429, 1.160735, 0.124311, -0.003366, 1.073046, 0.141742, -0.004376, 0.974346, 0.157718, 0.005418, 0.865649, 0.172077, 0.006449, 0.748068, 0.184669, 0.007427, 0.622811, 0.195367, 0.008313, 0.491164, 0.204059, 0.009069, 0.354476, 0.210658, 0.009665, 0.214151, 0.215096, 0.010076, 0.071629, 0.217326, 0.010286, -0.071629, 0.217326, 0.010286, -0.214151, 0.215096, 0.010076, -0.354476, 0.210658, 0.009665, -0.491164, 0.204059, 0.009069, -0.622811, 0.195367, 0.008313, 
    -0.748068, 0.184669, 0.007427, -0.865649, 0.172077, 0.006449, -0.974346, 0.157718, 0.005418, -1.073046, 0.141742, -0.004376, -1.160735, 0.124311, -0.003366, -1.236513, 0.105604, -0.002429, -1.299603, 0.085814, -0.001604, -1.349357, 0.065143, -0.000924, -1.385265, 0.043804, -0.000418, -1.406958, 0.022015, -0.000106, 1.406958, 0.036432, -0.000103, 1.385265, 0.072490, -0.000409, 1.349357, 0.107804, -0.000905, 1.299603, 0.142012, -0.001570, 1.236513, 0.174763, -0.002377, 1.160735, 0.205720, -0.003294, 1.073046, 0.234566, -0.004282, 0.974346, 0.261006, -0.005302, 0.865649, 0.284767, 0.006311, 0.748068, 0.305606, 0.007269, 0.622811, 0.323309, 0.008135, 0.491164, 0.337695, 0.008876, 0.354476, 0.348615, 0.009459, 0.214151, 0.355958, 0.009862, 0.071629, 0.359649, 0.010067, -0.071629, 0.359649, 0.010067, -0.214151, 0.355958, 0.009862, -0.354476, 0.348615, 0.009459, -0.491164, 0.337695, 0.008876, -0.622811, 0.323309, 0.008135, -0.748068, 0.305606, 0.007269, -0.865649, 0.284767, 0.006311, -0.974346, 0.261006, -0.005302, -1.073046, 0.234566, -0.004282, -1.160735, 0.205720, -0.003294, -1.236513, 0.174763, -0.002377, -1.299603, 0.142012, -0.001570, -1.349357, 0.107804, -0.000905, -1.385265, 0.072490, -0.000409, -1.406958, 0.036432, -0.000103, 1.406958, 0.050463, -0.000100, 1.385265, 0.100408, -0.000396, 1.349357, 0.149322, -0.000875, 1.299603, 0.196705, -0.001519, 1.236513, 0.242069, -0.002300, 1.160735, 0.284948, -0.003187, 1.073046, 0.324904, -0.004144, 0.974346, 0.361527, -0.005130, 0.865649, 0.394439, 0.006107, 0.748068, 0.423304, 0.007034, 0.622811, 0.447825, 0.007872, 0.491164, 0.467751, 0.008588, 0.354476, 0.482877, 0.009153, 0.214151, 0.493048, 0.009542, 0.071629, 0.498160, 0.009741, -0.071629, 0.498160, 0.009741, -0.214151, 0.493048, 0.009542, -0.354476, 0.482877, 0.009153, -0.491164, 0.467751, 0.008588, -0.622811, 0.447825, 0.007872, -0.748068, 0.423304, 0.007034, -0.865649, 0.394439, 0.006107, -0.974346, 0.361527, -0.005130, -1.073046, 0.324904, -0.004144, -1.160735, 0.284948, -0.003187, -1.236513, 0.242069, -0.002300, -1.299603, 0.196705, -0.001519, -1.349357, 0.149322, -0.000875, -1.385265, 0.100408, -0.000396, -1.406958, 0.050463, -0.000100, 1.406958, 0.063959, -0.000096, 1.385265, 0.127261, -0.000378, 1.349357, 0.189258, -0.000837, 1.299603, 0.249312, -0.001452, 1.236513, 0.306808, -0.002199, 1.160735, 0.361156, -0.003047, 1.073046, 0.411798, -0.003961, 0.974346, 0.458215, -0.004904, 0.865649, 0.499929, 0.005838, 0.748068, 0.536514, 0.006724, 0.622811, 0.567593, 0.007525, 0.491164, 0.592848, 0.008210, 0.354476, 0.612020, 0.008749, 0.214151, 0.624911, 0.009122, 0.071629, 0.631390, 0.009312, -0.071629, 0.631390, 0.009312, -0.214151, 0.624911, 0.009122, -0.354476, 0.612020, 0.008749, -0.491164, 0.592848, 0.008210, -0.622811, 0.567593, 0.007525, -0.748068, 0.536514, 0.006724, -0.865649, 0.499929, 0.005838, -0.974346, 0.458215, -0.004904, -1.073046, 0.411798, -0.003961, -1.160735, 0.361156, -0.003047, -1.236513, 0.306808, -0.002199, -1.299603, 0.249312, -0.001452, -1.349357, 0.189258, -0.000837, -1.385265, 0.127261, -0.000378, -1.406958, 0.063959, -0.000096, 1.406958, 0.076777, -0.000090, 1.385265, 0.152766, -0.000357, 1.349357, 0.227187, -0.000789, 1.299603, 0.299277, -0.001370, 1.236513, 0.368296, -0.002074, 1.160735, 0.433536, -0.002874, 1.073046, 0.494327, -0.003737, 0.974346, 0.550046, -0.004626, 0.865649, 0.600121, -0.005507, 0.748068, 0.644037, 0.006342, 0.622811, 0.681345, 0.007099, 0.491164, 0.711661, 0.007744, 0.354476, 0.734675, 0.008253, 0.214151, 0.750150, 0.008605, 0.071629, 0.757927, 0.008784, -0.071629, 0.757927, 0.008784, -0.214151, 0.750150, 0.008605, -0.354476, 0.734675, 0.008253, -0.491164, 0.711661, 0.007744, -0.622811, 0.681345, 0.007099, -0.748068, 0.644037, 0.006342, -0.865649, 0.600121, -0.005507, -0.974346, 0.550046, -0.004626, -1.073046, 0.494327, -0.003737, -1.160735, 0.433536, -0.002874, -1.236513, 0.368296, -0.002074, -1.299603, 0.299277, -0.001370, -1.349357, 0.227187, -0.000789, -1.385265, 0.152766, -0.000357, -1.406958, 0.076777, -0.000090, 1.406958, 0.088781, -0.000084, 1.385265, 0.176651, -0.000332, 1.349357, 0.262708, -0.000733, 1.299603, 0.346070, -0.001273, 1.236513, 0.425880, -0.001927, 1.160735, 0.501320, -0.002671, 1.073046, 0.571616, -0.003472, 0.974346, 0.636047, -0.004299, 0.865649, 0.693950, -0.005118, 0.748068, 0.744733, -0.005894, 0.622811, 0.787874, -0.006597, 0.491164, 0.822931, 0.007197, 0.354476, 0.849542, 0.007670, 0.214151, 0.867437, 0.007996, 0.071629, 0.876430, 0.008163, -0.071629, 0.876430, 0.008163, -0.214151, 0.867437, 0.007996, -0.354476, 0.849542, 0.007670, -0.491164, 0.822931, 0.007197, -0.622811, 0.787874, -0.006597, -0.748068, 0.744733, -0.005894, -0.865649, 0.693950, -0.005118, -0.974346, 0.636047, -0.004299, -1.073046, 0.571616, -0.003472, -1.160735, 0.501320, -0.002671, -1.236513, 0.425880, -0.001927, -1.299603, 0.346070, -0.001273, -1.349357, 0.262708, -0.000733, -1.385265, 0.176651, -0.000332, -1.406958, 0.088781, -0.000084, 1.406958, 0.099844, -0.000077, 1.385265, 0.198663, -0.000303, 1.349357, 0.295444, -0.000670, 1.299603, 0.389194, -0.001162, 1.236513, 0.478949, -0.001760, 1.160735, 0.563790, -0.002439, 1.073046, 0.642846, 
    -0.003171, 0.974346, 0.715305, -0.003927, 0.865649, 0.780424, -0.004674, 0.748068, 0.837535, -0.005383, 0.622811, 0.886052, -0.006025, 0.491164, 0.925477, -0.006573, 0.354476, 0.955405, -0.007005, 0.214151, 0.975529, 0.007303, 0.071629, 0.985643, 0.007455, -0.071629, 0.985643, 0.007455, -0.214151, 0.975529, 0.007303, -0.354476, 0.955405, -0.007005, -0.491164, 0.925477, -0.006573, -0.622811, 0.886052, -0.006025, -0.748068, 0.837535, -0.005383, -0.865649, 0.780424, -0.004674, -0.974346, 0.715305, -0.003927, -1.073046, 0.642846, -0.003171, -1.160735, 0.563790, -0.002439, -1.236513, 0.478949, -0.001760, -1.299603, 0.389194, -0.001162, -1.349357, 0.295444, -0.000670, -1.385265, 0.198663, -0.000303, -1.406958, 0.099844, -0.000077, 1.406958, 0.109849, -0.000068, 1.385265, 0.218570, -0.000271, 1.349357, 0.325049, -0.000599, 1.299603, 0.428192, -0.001040, 1.236513, 0.526941, -0.001575, 1.160735, 0.620284, -0.002182, 1.073046, 0.707261, -0.002837, 0.974346, 0.786981, -0.003512, 0.865649, 0.858625, -0.004181, 0.748068, 0.921459, -0.004815, 0.622811, 0.974837, -0.005389, 0.491164, 1.018212, -0.005879, 0.354476, 1.051139, -0.006266, 0.214151, 1.073280, -0.006533, 0.071629, 1.084408, -0.006669, -0.071629, 1.084408, -0.006669, -0.214151, 1.073280, -0.006533, -0.354476, 1.051139, -0.006266, -0.491164, 1.018212, -0.005879, -0.622811, 0.974837, -0.005389, -0.748068, 0.921459, -0.004815, -0.865649, 0.858625, -0.004181, -0.974346, 0.786981, -0.003512, -1.073046, 0.707261, -0.002837, -1.160735, 0.620284, -0.002182, -1.236513, 0.526941, -0.001575, -1.299603, 0.428192, -0.001040, -1.349357, 0.325049, -0.000599, -1.385265, 0.218570, -0.000271, -1.406958, 0.109849, -0.000068, 1.406958, 0.118689, -0.000060, 1.385265, 0.236160, -0.000236, 1.349357, 0.351208, -0.000522, 1.299603, 0.462652, -0.000906, 1.236513, 0.569348, -0.001372, 1.160735, 0.670202, -0.001901, 1.073046, 0.764179, -0.002472, 0.974346, 0.850315, -0.003061, 0.865649, 0.927725, -0.003643, 0.748068, 0.995615, -0.004196, 0.622811, 1.053289, -0.004696, 0.491164, 1.100155, -0.005124, 0.354476, 1.135732, -0.005460, 0.214151, 1.159654, -0.005693, 0.071629, 1.171677, -0.005812, -0.071629, 1.171677, -0.005812, -0.214151, 1.159654, -0.005693, -0.354476, 1.135732, -0.005460, -0.491164, 1.100155, -0.005124, -0.622811, 1.053289, -0.004696, -0.748068, 0.995615, -0.004196, -0.865649, 0.927725, -0.003643, -0.974346, 0.850315, -0.003061, -1.073046, 0.764179, -0.002472, -1.160735, 0.670202, -0.001901, -1.236513, 0.569348, -0.001372, -1.299603, 0.462652, -0.000906, -1.349357, 0.351208, -0.000522, -1.385265, 0.236160, -0.000236, -1.406958, 0.118689, -0.000060, 1.406958, 0.126271, -0.000050, 1.385265, 0.251247, -0.000199, 1.349357, 0.373644, -0.000440, 1.299603, 0.492207, -0.000763, 1.236513, 0.605719, -0.001155, 1.160735, 0.713016, -0.001601, 1.073046, 0.812997, -0.002081, 0.974346, 0.904635, -0.002577, 0.865649, 0.986990, -0.003067, 0.748068, 1.059218, -0.003533, 0.622811, 1.120576, -0.003954, 0.491164, 1.170436, -0.004314, 0.354476, 1.208285, -0.004597, 0.214151, 1.233736, -0.004793, 0.071629, 1.246527, -0.004893, -0.071629, 1.246527, -0.004893, -0.214151, 1.233736, -0.004793, -0.354476, 1.208285, -0.004597, -0.491164, 1.170436, -0.004314, -0.622811, 1.120576, -0.003954, -0.748068, 1.059218, -0.003533, -0.865649, 0.986990, -0.003067, -0.974346, 0.904635, -0.002577, -1.073046, 0.812997, -0.002081, -1.160735, 0.713016, -0.001601, -1.236513, 0.605719, -0.001155, -1.299603, 0.492207, -0.000763, -1.349357, 0.373644, -0.000440, -1.385265, 0.251247, -0.000199, -1.406958, 0.126271, -0.000050, 1.406958, 0.132515, -0.000040, 1.385265, 0.263670, -0.000159, 1.349357, 0.392119, -0.000352, 1.299603, 0.516545, -0.000611, 1.236513, 0.635670, -0.000926, 1.160735, 0.748273, -0.001283, 1.073046, 0.853197, -0.001668, 0.974346, 0.949366, -0.002066, 0.865649, 1.035793, -0.002459, 0.748068, 1.111592, -0.002832, 0.622811, 1.175984, -0.003169, 0.491164, 1.228310, -0.003458, 0.354476, 1.268031, -0.003685, 0.214151, 1.294740, -0.003842, 0.071629, 1.308163, -0.003922, -0.071629, 1.308163, -0.003922, -0.214151, 1.294740, -0.003842, -0.354476, 1.268031, -0.003685, -0.491164, 1.228310, -0.003458, -0.622811, 1.175984, -0.003169, -0.748068, 1.111592, -0.002832, -0.865649, 1.035793, -0.002459, -0.974346, 0.949366, -0.002066, -1.073046, 0.853197, -0.001668, -1.160735, 0.748273, -0.001283, -1.236513, 0.635670, -0.000926, -1.299603, 0.516545, -0.000611, -1.349357, 0.392119, -0.000352, -1.385265, 0.263670, -0.000159, -1.406958, 0.132515, -0.000040, 1.406958, 0.137354, -0.000030, 1.385265, 0.273298, -0.000118, 1.349357, 0.406438, -0.000261, 1.299603, 0.535407, -0.000454, 1.236513, 0.658883, -0.000687, 1.160735, 0.775597, -0.000952, 1.073046, 0.884353, -0.001238, 0.974346, 0.984034, -0.001532, 0.865649, 1.073617, -0.001824, 0.748068, 1.152184, -0.002101, 0.622811, 1.218928, -0.002351, 0.491164, 1.273164, -0.002565, 0.354476, 1.314335, -0.002734, 0.214151, 1.342020, -0.002850, 0.071629, 1.355933, -0.002910, -0.071629, 1.355933, -0.002910, -0.214151, 1.342020, -0.002850, -0.354476, 1.314335, -0.002734, -0.491164, 1.273164, -0.002565, -0.622811, 1.218928, -0.002351, -0.748068, 1.152184, -0.002101, -0.865649, 1.073617, -0.001824, -0.974346, 0.984034, -0.001532, -1.073046, 
    0.884353, -0.001238, -1.160735, 0.775597, -0.000952, -1.236513, 0.658883, -0.000687, -1.299603, 0.535407, -0.000454, -1.349357, 0.406438, -0.000261, -1.385265, 0.273298, -0.000118, -1.406958, 0.137354, -0.000030, 1.406958, 0.140737, -0.000019, 1.385265, 0.280030, -0.000076, 1.349357, 0.416449, -0.000168, 1.299603, 0.548595, -0.000291, 1.236513, 0.675112, -0.000441, 1.160735, 0.794701, -0.000611, 1.073046, 0.906135, -0.000794, 0.974346, 1.008271, -0.000983, 0.865649, 1.100061, -0.001170, 0.748068, 1.180563, -0.001348, 0.622811, 1.248951, -0.001508, 0.491164, 1.304523, -0.001646, 0.354476, 1.346708, -0.001754, 0.214151, 1.375075, -0.001829, 0.071629, 1.389331, -0.001867, -0.071629, 1.389331, -0.001867, -0.214151, 1.375075, -0.001829, -0.354476, 1.346708, -0.001754, -0.491164, 1.304523, -0.001646, -0.622811, 1.248951, -0.001508, -0.748068, 1.180563, -0.001348, -0.865649, 1.100061, -0.001170, -0.974346, 1.008271, -0.000983, -1.073046, 0.906135, -0.000794, -1.160735, 0.794701, -0.000611, -1.236513, 0.675112, -0.000441, -1.299603, 0.548595, -0.000291, -1.349357, 0.416449, -0.000168, -1.385265, 0.280030, -0.000076, -1.406958, 0.140737, -0.000019, 1.406958, 0.142629, -0.000008, 1.385265, 0.283795, -0.000033, 1.349357, 0.422048, -0.000072, 1.299603, 0.555971, -0.000126, 1.236513, 0.684189, -0.000190, 1.160735, 0.805386, -0.000264, 1.073046, 0.918318, -0.000343, 0.974346, 1.021828, -0.000424, 0.865649, 1.114852, -0.000505, 0.748068, 1.196436, -0.000582, 0.622811, 1.265743, -0.000651, 0.491164, 1.322062, -0.000710, 0.354476, 1.364815, -0.000757, 0.214151, 1.393563, -0.000789, 0.071629, 1.408011, -0.000805, -0.071629, 1.408011, -0.000805, -0.214151, 1.393563, -0.000789, -0.354476, 1.364815, -0.000757, -0.491164, 1.322062, -0.000710, -0.622811, 1.265743, -0.000651, -0.748068, 1.196436, -0.000582, -0.865649, 1.114852, -0.000505, -0.974346, 1.021828, -0.000424, -1.073046, 0.918318, -0.000343, -1.160735, 0.805386, -0.000264, -1.236513, 0.684189, -0.000190, -1.299603, 0.555971, -0.000126, -1.349357, 0.422048, -0.000072, -1.385265, 0.283795, -0.000033, -1.406958, 0.142629, -0.000008;

    return tmp;
}();

}  // namespace zebrafish
