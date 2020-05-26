#include <zebrafish/Bspline.h>
#include <zebrafish/common.h>

#include <Eigen/Dense>
#include <Eigen/Core>


namespace zebrafish {

void Interp3D(const image_t &image, const Eigen::MatrixX3d &sample, Eigen::VectorXd &res) {

    double t1, t2, t3, t4;
    Eigen::MatrixX3d truncSample = sample.array().floor();
    Eigen::MatrixX3d baseIdx = truncSample.array() - 1.0;

    Eigen::Array<double, 64, 1> yArray, BxArray, ByArray, BzArray;

    // yArray
    const Eigen::Matrix4d &z1layer = image[baseIdx(2)  ].block<4, 4>(baseIdx(0), baseIdx(1));
    const Eigen::Matrix4d &z2layer = image[baseIdx(2)+1].block<4, 4>(baseIdx(0), baseIdx(1));
    const Eigen::Matrix4d &z3layer = image[baseIdx(2)+2].block<4, 4>(baseIdx(0), baseIdx(1));
    const Eigen::Matrix4d &z4layer = image[baseIdx(2)+3].block<4, 4>(baseIdx(0), baseIdx(1));
    yArray << Eigen::Map<const Eigen::VectorXd>(z1layer.data(), 16), Eigen::Map<const Eigen::VectorXd>(z2layer.data(), 16),
              Eigen::Map<const Eigen::VectorXd>(z3layer.data(), 16), Eigen::Map<const Eigen::VectorXd>(z4layer.data(), 16);  // advanced initialization feature

    // Bx
    t1 = sample(0) - truncSample(0) + 1;
    t2 = sample(0) - truncSample(0);
    t3 = sample(0) - truncSample(0) - 1;
    t4 = sample(0) - truncSample(0) - 2;
    t1 = - 1.0/6.0 * t1*t1*t1 + t1*t1 - 2 * t1 + 4.0/3.0;
    t2 = 0.5 * t2*t2*t2 - t2*t2 + 2.0/3.0;
    t3 = - 0.5 * t3*t3*t3 - t3*t3 + 2.0/3.0;
    t4 = 1.0/6.0 * t4*t4*t4 + t4*t4 + 2 * t4 + 4.0/3.0;
    BxArray << t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, 
               t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4, t1, t2, t3, t4;
    
    // By
    t1 = sample(1) - truncSample(1) + 1;
    t2 = sample(1) - truncSample(1);
    t3 = sample(1) - truncSample(1) - 1;
    t4 = sample(1) - truncSample(1) - 2;
    t1 = - 1.0/6.0 * t1*t1*t1 + t1*t1 - 2 * t1 + 4.0/3.0;
    t2 = 0.5 * t2*t2*t2 - t2*t2 + 2.0/3.0;
    t3 = - 0.5 * t3*t3*t3 - t3*t3 + 2.0/3.0;
    t4 = 1.0/6.0 * t4*t4*t4 + t4*t4 + 2 * t4 + 4.0/3.0;
    ByArray << t1, t1, t1, t1, t2, t2, t2, t2, t3, t3, t3, t3, t4, t4, t4, t4, t1, t1, t1, t1, t2, t2, t2, t2, t3, t3, t3, t3, t4, t4, t4, t4,
               t1, t1, t1, t1, t2, t2, t2, t2, t3, t3, t3, t3, t4, t4, t4, t4, t1, t1, t1, t1, t2, t2, t2, t2, t3, t3, t3, t3, t4, t4, t4, t4;

    // Bz
    t1 = sample(2) - truncSample(2) + 1;
    t2 = sample(2) - truncSample(2);
    t3 = sample(2) - truncSample(2) - 1;
    t4 = sample(2) - truncSample(2) - 2;
    t1 = - 1.0/6.0 * t1*t1*t1 + t1*t1 - 2 * t1 + 4.0/3.0;
    t2 = 0.5 * t2*t2*t2 - t2*t2 + 2.0/3.0;
    t3 = - 0.5 * t3*t3*t3 - t3*t3 + 2.0/3.0;
    t4 = 1.0/6.0 * t4*t4*t4 + t4*t4 + 2 * t4 + 4.0/3.0;
    BzArray << t1, t1, t1, t1, t1, t1, t1, t1, t1, t1, t1, t1, t1, t1, t1, t1, t2, t2, t2, t2, t2, t2, t2, t2, t2, t2, t2, t2, t2, t2, t2, t2, 
               t3, t3, t3, t3, t3, t3, t3, t3, t3, t3, t3, t3, t3, t3, t3, t3, t4, t4, t4, t4, t4, t4, t4, t4, t4, t4, t4, t4, t4, t4, t4, t4;

    // res = \sum{ yArray * BxArray * ByArray * BzArray }
    res(0) = (yArray * BxArray * ByArray * BzArray).sum();
}

}  // namespace zebrafish
