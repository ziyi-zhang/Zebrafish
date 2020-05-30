#include <zebrafish/Bspline.h>
#include <zebrafish/common.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>


namespace zebrafish {

namespace {

    inline double CubicBasis(double t) {
    /// Defines the cubic basis function centered at t=0 for uniform knot vector B-spline
    /// If t>2 or t<-2, the function will return zero due to local control property

        if (t>2 || t<-2)
            return 0;
        else if (t > 1)
            return -1/6 * t*t*t + t*t - 2 * t + 4/3;
        else if (t > 0)
            return 0.5 * t*t*t - t*t + 2/3;
        else if (t > -1)
            return -0.5* t*t*t - t*t + 2/3;
        else // t > -2
            return 1/6 * t*t*t + t*t + 2*t + 4/3;
    }
}  // anonymous namespace


void bspline::CalcControlPts(const image_t &image, const double xratio, const double yratio, const double zratio) {
// Note: this function will only be called once

    assert(xratio <= 1 && yratio <= 1 && zratio <= 1);

    int Nz = image.size(), Nx = image[0].rows(), Ny = image[0].cols();  // dimension of sample points
    int num, N = Nx*Ny*Nz;
    int i, j, ix, iy, iz, jx, jy, jz;
    std::vector<double> centersX, centersY, centersZ;  // location of the center of a control point's basis function
    double centerX, centerY, centerZ, t;
    Eigen::VectorXd inputPts, vectorY;
    std::vector<Eigen::Triplet<double> > tripletArray;  // to fill in sparse matrix

    numX = ceil(Nx * xratio);  // dimension of control points
    numY = ceil(Ny * yratio);
    numZ = ceil(Nz * zratio);
    num = numX * numY * numZ;
    // Eigen::MatrixXd A;
    Eigen::SparseMatrix<double> A(N, num);  // A[i, j] is the evaluation of basis j for the i-th sample point
    Eigen::SparseMatrix<double> AtransposeA(num, num);

    gapX = Nx / (numX - 1);  // the interval between two control points along a direction
    gapY = Ny / (numY - 1);
    gapZ = Nz / (numZ - 1);
    for (i=0; i<numX; i++)
        centersX.push_back(i * gapX);
    for (i=0; i<numY; i++)
        centersY.push_back(i * gapY);
    for (i=0; i<numZ; i++)
        centersZ.push_back(i * gapZ);

    // map 3D "image" to 1D "inputPts"
    // order matters! z -> x -> y
    inputPts.resize(N);
    for (auto it=image.begin(); it!=image.end(); it++) {
        inputPts << Eigen::Map<const Eigen::VectorXd>(it->data(), Nx*Ny);
    }

    // calculate matrix A
    for (jz=0; jz<numZ; jz++)
        for (jx=0; jx<numX; jx++)
            for (jy=0; jy<numY; jy++) {

                j = jz * numX * numY + jx * numY + jy;  // iterating control points' basis function
                centerX = centersX[jx];
                centerY = centersY[jy];
                centerZ = centersZ[jz];
                for (iz=0; iz<Nz; iz++)
                    for (ix=0; ix<Nx; ix++)
                        for (iy=0; iy<Ny; iy++) {

                            i = iz * Nx * Ny + ix * Ny + iy;  // iterating sample points
                            t = CubicBasis( (ix-centerX) / gapX ) * 
                                CubicBasis( (iy-centerY) / gapY ) * 
                                CubicBasis( (iz-centerZ) / gapZ );
                            if (fabs(t) > 1e-3)
                                tripletArray.push_back(Eigen::Triplet<double>(i, j, t));
                        }
            }
    A.setFromTriplets(tripletArray.begin(), tripletArray.end());

    // calculate control points based on least square
    std::cout << "Starting solving linear system..." << std::endl;
    AtransposeA = (A.transpose() * A).pruned();  // A' * A
    vectorY.resize(N);
    vectorY = A.transpose() * inputPts;  // A' * y
    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(AtransposeA);  // performs Cholesky factorization
    controlPoints = chol.solve(vectorY);
    std::cout << "Control points calculated..." << std::endl;
}


void bspline::Interp3D(const image_t &image, const Eigen::MatrixX3d &sample, Eigen::VectorXd &res) {

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


bspline::bspline() {

    controlPoints.setZero();
}


bspline::~bspline() {

}

}  // namespace zebrafish
