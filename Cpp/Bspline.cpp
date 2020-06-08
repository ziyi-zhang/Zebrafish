#include <zebrafish/Bspline.h>
#include <zebrafish/Common.h>

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
#include <array>
#include <chrono>
#include <ctime>


namespace zebrafish {

namespace {

    inline double CubicBasis(double t) {
    /// Defines the cubic basis function centered at t=0 for uniform knot vector B-spline
    /// If t>2 or t<-2, the function will return zero due to local control property

        if (t>=2 || t<=-2)
            return 0;
        else if (t > 1)
            return -1.0/6.0 * t*t*t + t*t - 2.0 * t + 4.0/3.0;
        else if (t > 0)
            return 0.5 * t*t*t - t*t + 2.0/3.0;
        else if (t > -1)
            return -0.5 * t*t*t - t*t + 2.0/3.0;
        else // t > -2
            return 1.0/6.0 * t*t*t + t*t + 2.0*t + 4.0/3.0;
    }


    void PrintTime(std::string str) {
    // print human readable time

        auto timeStamp = std::chrono::system_clock::now();
        std::time_t time = std::chrono::system_clock::to_time_t(timeStamp);
        std::cout << str << " " << std::ctime(&time);
    }
}  // anonymous namespace


void bspline::CalcControlPts(const image_t &image, const double xratio, const double yratio, const double zratio) {
// Note: this function will only be called once

    assert(xratio <= 1 && yratio <= 1 && zratio <= 1);
    std::cout << "====================================================" << std::endl;
    PrintTime("Start calculating control points...");

    int Nz = image.size(), Nx = image[0].rows(), Ny = image[0].cols();  // dimension of sample points
    int num, N = Nx*Ny*Nz;
    int i, j, ix, iy, iz, jx, jy, jz, refIdx_x, refIdx_y, refIdx_z;
    int jxMin, jxMax, jyMin, jyMax, jzMin, jzMax;
    double centerX, centerY, centerZ, t;
    Eigen::VectorXd inputPts, vectorY;
    numX = ceil(Nx * xratio);  // dimension of control points
    numY = ceil(Ny * yratio);
    numZ = ceil(Nz * zratio);
    num = numX * numY * numZ;
    // Eigen::MatrixXd A;
    Eigen::SparseMatrix<double> A(N, num);  // A[i, j] is the evaluation of basis j for the i-th sample point
    Eigen::SparseMatrix<double> AtransposeA(num, num);
    std::cout << "numX= " << numX << ", numY= " << numY << ", numZ= " << numZ << std::endl;

    gapX = double(Nx-1) / double(numX-1);  // the interval between two control points along a direction
    gapY = double(Ny-1) / double(numY-1);
    gapZ = double(Nz-1) / double(numZ-1);
    centersX.clear();  // location of the center of a control point's basis function
    centersY.clear();
    centersZ.clear();
    for (i=0; i<numX; i++)
        centersX.push_back(i * gapX);
    for (i=0; i<numY; i++)
        centersY.push_back(i * gapY);
    for (i=0; i<numZ; i++)
        centersZ.push_back(i * gapZ);

    // map 3D "image" to 1D "inputPts"
    // order matters! z -> y -> x
    inputPts.resize(N);
    t = 0;
    for (auto it=image.begin(); it!=image.end(); it++) {
        inputPts.segment(t*Nx*Ny, Nx*Ny) << Eigen::Map<const Eigen::VectorXd>(it->data(), Nx*Ny);
        t++;
    }

    // calculate matrix A
    // order matters!
    // for column (j): z -> x -> y
    // for row (i): z -> y -> x
    int count;
    double rowSum;
    std::array<int, 64> cache_i, cache_j;
    std::array<double, 64> cache_t;
    PrintTime("Start filling least square matrix...");
    for (iz=0; iz<=Nz-1; iz++)
        for (iy=0; iy<=Ny-1; iy++)
            for (ix=0; ix<=Nx-1; ix++) {

                i = iz * Nx * Ny + iy * Nx + ix;  // iterating sample points
                if (i % 5000 == 0) std::cout << i << " / " << Nx*Ny*Nz << std::endl;
                
                refIdx_x = floor(ix / gapX) - 1;
                refIdx_y = floor(iy / gapY) - 1;
                refIdx_z = floor(iz / gapZ) - 1;
                jxMin = std::max(0, refIdx_x);
                jxMax = std::min(numX-1, refIdx_x+3);
                jyMin = std::max(0, refIdx_y);
                jyMax = std::min(numY-1, refIdx_y+3);
                jzMin = std::max(0, refIdx_z);
                jzMax = std::min(numZ-1, refIdx_z+3);
                count = 0;
                rowSum = 0.0;
                for (jz=jzMin; jz<=jzMax; jz++)
                    for (jx=jxMin; jx<=jxMax; jx++)
                        for (jy=jyMin; jy<=jyMax; jy++) {

                            j = jz * numX * numY + jx * numY + jy;  // iterating control points' basis function
                            t = CubicBasis( (ix-centersX[jx]) / gapX) * 
                                CubicBasis( (iy-centersY[jy]) / gapY) * 
                                CubicBasis( (iz-centersZ[jz]) / gapZ);
 
                            if (fabs(t) > 0.000000) {
                                cache_i[count] = i;
                                cache_j[count] = j;
                                cache_t[count] = t;
                                count++;
                                rowSum += t;
                                // A.insert(i, j) = t;  // direct insertion
                                // printf("%d %d : %.5f\n", i, j, t);  // DEBUG only
                            }
                        }
                // Fix a minor error when calculating matrix A
                for (i=0; i<count; i++) {
                    // A.insert(cache_i[i], cache_j[i]) = cache_t[i] / rowSum;
                    // if (rowSum != 1.0) std::cout << rowSum << " ";
                    A.insert(cache_i[i], cache_j[i]) = cache_t[i];
                }
            }

    std::cout << "Control points " << numX << " * " << numY << " * " << numZ << std::endl;
    std::cout << "Matrix A size = " << A.rows() << " * " << A.cols() << std::endl;
    std::cout << "A: # non-zero elements = " << A.nonZeros() << std::endl;

    // calculate control points based on least square
    PrintTime("Calculating A'*A and A'*y...");
    AtransposeA = (A.transpose() * A).pruned();  // A' * A
    std::cout << "A'A: # non-zero elements = " << AtransposeA.nonZeros() << std::endl;
    vectorY.resize(N);
    vectorY = A.transpose() * inputPts;  // A' * y
    // Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(AtransposeA);  // performs Cholesky factorization
    // controlPoints = chol.solve(vectorY);

    ///////////////// TEST ONLY
    /*
    Eigen::SparseMatrix<double> testMatrix(3, 3);
    testMatrix.insert(0, 0) = 1;
    testMatrix.insert(1, 1) = 2;
    testMatrix.insert(2, 2) = 3;
    Eigen::Vector3d testVector;
    testVector << 4, 5, 6;
    AtransposeA = testMatrix;
    vectorY = testVector;
    controlPoints.resize(3, 1);

    std::cout << std::endl;
    std::cout << AtransposeA << std::endl;
    std::cout << vectorY << std::endl;
    std::cout << std::endl;
    */
    ///////////////// TEST ONLY

    const std::string solverName = "Hypre";
    auto solver = polysolve::LinearSolver::create(solverName, "");
    const nlohmann::json params = {
        {"max_iter", 3000}, 
        {"tolerance", 1e-15}
    };
    solver->setParameters(params);
        PrintTime("Start analyzing matrix pattern...");
    solver->analyzePattern(AtransposeA, AtransposeA.rows());
        PrintTime("Start factorizing matrix...");
    solver->factorize(AtransposeA);
        PrintTime("Start solving linear system...");
    controlPoints.resize(num, 1);
    solver->solve(vectorY, controlPoints);
    /////////////// TEST ONLY ///////////
    // std::cout << ">>>>>" << std::endl;
    // std::cout << controlPoints << std::endl;
    /////////////// TEST ONLY ///////////
        PrintTime("Control points calculated...");
        std::cout << "====================================================" << std::endl;
}


void bspline::Interp3D(const Eigen::MatrixX3d &sample, Eigen::VectorXd &res) const {

    assert(numX != 0 && numY != 0 && numZ != 0);

    int i, ix, iy, iz, idx_x, idx_y, idx_z;
    double xcoef, ycoef, zcoef;
    int refIdx_x, refIdx_y, refIdx_z;
    double evalX1, evalX2, evalX3, evalX4, evalY1, evalY2, evalY3, evalY4, evalZ1, evalZ2, evalZ3, evalZ4;
    Eigen::MatrixX3d truncSample = sample.array().floor();
    Eigen::MatrixX3d baseIdx = truncSample.array() - 1.0;

    res.resize(sample.rows(), 1);
    for (i=0; i<sample.rows(); i++) {

        // Get reference index
        refIdx_x = floor(sample(i, 0) / gapX) - 1;
        refIdx_y = floor(sample(i, 1) / gapY) - 1;
        refIdx_z = floor(sample(i, 2) / gapZ) - 1;
        assert(refIdx_x >= 0 && refIdx_y >= 0 && refIdx_z >= 0);
        assert(refIdx_x <= numX-4 && refIdx_y <= numY-4 && refIdx_z <= numZ-4);
        // 
        res(i) = 0;
        for (ix=0; ix<=3; ix++)
            for (iy=0; iy<=3; iy++)
                for (iz=0; iz<=3; iz++) {

                    idx_x = refIdx_x + ix;
                    idx_y = refIdx_y + iy;
                    idx_z = refIdx_z + iz;

                    xcoef = CubicBasis( (sample(i, 0)-centersX[idx_x]) / gapX );
                    ycoef = CubicBasis( (sample(i, 1)-centersY[idx_y]) / gapY );
                    zcoef = CubicBasis( (sample(i, 2)-centersZ[idx_z]) / gapZ );

                    res(i) += controlPoints(numX*numY*idx_z + numY*idx_x + idx_y) * xcoef * ycoef * zcoef;
                }
    }

    /*
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
    */
}


bspline::bspline() {

    controlPoints.setZero();
    numX = 0;
    numY = 0;
    numZ = 0;
}


bspline::~bspline() {

}

}  // namespace zebrafish
