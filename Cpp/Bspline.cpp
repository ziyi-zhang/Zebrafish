#include <zebrafish/Bspline.h>
#include <zebrafish/common.h>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <polysolve/LinearSolver.hpp>
#include <string>
#include <cmath>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <ctime>


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
    int ixMin, ixMax, iyMin, iyMax, izMin, izMax;
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
    t = 0;
    for (auto it=image.begin(); it!=image.end(); it++) {
        inputPts.segment(t*Nx*Ny, Nx*Ny) << Eigen::Map<const Eigen::VectorXd>(it->data(), Nx*Ny);
        t++;
    }

    // calculate matrix A
    for (jz=0; jz<numZ; jz++)
        for (jx=0; jx<numX; jx++)
            for (jy=0; jy<numY; jy++) {

                j = jz * numX * numY + jx * numY + jy;  // iterating control points' basis function
                if (j % 100000 == 0)
                    std::cout << j << "/" << numX*numY*numZ << std::endl;
                centerX = centersX[jx];
                centerY = centersY[jy];
                centerZ = centersZ[jz];

                ixMin = std::max(0, int(floor(centerX-2*gapX)));
                ixMax = std::min(Nx, int(ceil(centerX+2*gapX)));
                iyMin = std::max(0, int(floor(centerY-2*gapY)));
                iyMax = std::min(Ny, int(ceil(centerY+2*gapY)));
                izMin = std::max(0, int(floor(centerZ-2*gapZ)));
                izMax = std::min(Nz, int(ceil(centerZ+2*gapZ)));
                for (iz=izMin; iz<izMax; iz++)
                    for (ix=ixMin; ix<ixMax; ix++)
                        for (iy=iyMin; iy<iyMax; iy++) {

                            i = iz * Nx * Ny + ix * Ny + iy;  // iterating sample points
                            t = CubicBasis( (ix-centerX) / gapX ) * 
                                CubicBasis( (iy-centerY) / gapY ) * 
                                CubicBasis( (iz-centerZ) / gapZ );
                            if (fabs(t) > 1e-3)
                                tripletArray.push_back(Eigen::Triplet<double>(i, j, t));
                        }
            }
    std::cout << "# of non-zero elements = " << tripletArray.size() << std::endl;
    std::cout << "Matrix size = " << A.rows() << " * " << A.cols() << std::endl;
    A.setFromTriplets(tripletArray.begin(), tripletArray.end());
    tripletArray.clear();
    tripletArray.shrink_to_fit();  // give my memory back to me!

    // calculate control points based on least square
    std::cout << "Starting solving linear system..." << std::endl;
    AtransposeA = (A.transpose() * A).pruned();  // A' * A
    vectorY.resize(N);
    vectorY = A.transpose() * inputPts;  // A' * y
    // Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(AtransposeA);  // performs Cholesky factorization
    // controlPoints = chol.solve(vectorY);
    
    ///////////////// TEST ONLY
    Eigen::SparseMatrix<double> testMatrix(3, 3);
    testMatrix.insert(0, 0) = 1;
    testMatrix.insert(1, 1) = 2;
    testMatrix.insert(2, 2) = 3;
    Eigen::Vector3d testVector;
    testVector << 4, 5, 6;
    AtransposeA = testMatrix;
    vectorY = testVector;

    std::cout << std::endl;
    std::cout << AtransposeA << std::endl;
    std::cout << vectorY << std::endl;
    std::cout << std::endl;
    ///////////////// TEST ONLY

    const std::string solverName = "Hypre";
    auto solver = polysolve::LinearSolver::create(solverName, "");
    const nlohmann::json params = {
        {"max_iter", 500}, 
        {"tolerance", 1e-5}
    };
    solver->setParameters(params);

        auto timeStamp = std::chrono::system_clock::now();
        std::time_t time = std::chrono::system_clock::to_time_t(timeStamp);
        std::cout << std::endl << "Start analyzing pattern..." << std::ctime(&time);
    solver->analyzePattern(AtransposeA, AtransposeA.rows());
    std::cout << "Matrix pattern analyzed" << std::endl;
        timeStamp = std::chrono::system_clock::now();
        time = std::chrono::system_clock::to_time_t(timeStamp);
        std::cout << std::endl << "Start factorizing matrix..." << std::ctime(&time);
    solver->factorize(AtransposeA);
    std::cout << "Matrix factorized" << std::endl;
        timeStamp = std::chrono::system_clock::now();
        time = std::chrono::system_clock::to_time_t(timeStamp);
        std::cout << std::endl << "Start solving linear system..." << std::ctime(&time);
    solver->solve(vectorY, controlPoints);

    std::cout << "Control points calculated..." << std::endl;
        timeStamp = std::chrono::system_clock::now();
        time = std::chrono::system_clock::to_time_t(timeStamp);
        std::cout << std::endl << "Control points generated..." << std::ctime(&time);
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
