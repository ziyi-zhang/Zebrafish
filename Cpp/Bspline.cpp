#include <zebrafish/Bspline.h>
#include <zebrafish/Common.h>
#include <zebrafish/autodiff.h>

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


    inline DScalar CubicBasis(DScalar t) {
    /// Defines the cubic basis function centered at t=0 for uniform knot vector B-spline
    /// If t>2 or t<-2, the function will return zero due to local control property

        if (t>=2 || t<=-2)
            return DScalar(0);
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
        std::cout << str << "       " << std::ctime(&time);
    }
}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void bspline::CalcControlPts_um(const image_t &image, const double distX, const double distY, const double distZ) {

    double xratio, yratio, zratio;

    xratio = resolutionX / distX;  // if resolution is 1um, want 2um, the ratio should be 0.5
    yratio = resolutionY / distY;
    zratio = resolutionZ / distZ;
    assert(xratio > 0 && yratio > 0 && zratio > 0);
    assert(xratio <= 1 && yratio <= 1 && zratio <= 1);

    CalcControlPts(image, xratio, yratio, zratio);
}


void bspline::CalcControlPts(const image_t &image, const double xratio, const double yratio, const double zratio) {
// Note: this function will only be called once for one 3D image

    assert(xratio <= 1 && yratio <= 1 && zratio <= 1);
        std::cout << "====================================================" << std::endl;
        PrintTime("Start calculating control points...");

    int Nx, Ny, Nz, N, num;
    double centerX, centerY, centerZ, t;
    Eigen::VectorXd inputPts, vectorY;
    
    // dimension of sample points
    Nz = image.size();
    Nx = image[0].rows();
    Ny = image[0].cols();
    N = Nx*Ny*Nz;
        std::cout << "Sample N= " << N << " Nx= " << Nx << " Ny= " << Ny << " Nz= " << Nz << std::endl;
    
    // dimension of control points
    numX = ceil(Nx * xratio);
    numY = ceil(Ny * yratio);
    numZ = ceil(Nz * zratio);
    num = numX * numY * numZ;
        std::cout << "Control num= " << num << " numX= " << numX << " numY= " << numY << " numZ= " << numZ << std::endl;

    // A[i, j]: the evaluation of basis j for the i-th sample point
    Eigen::SparseMatrix<double> A(N, num);
    Eigen::SparseMatrix<double> AtransposeA(num, num);

    // the interval between two control points along a direction
    gapX = double(Nx-1) / double(numX-1);
    gapY = double(Ny-1) / double(numY-1);
    gapZ = double(Nz-1) / double(numZ-1);
    // location of the center of a control point's basis function
    centersX = Eigen::VectorXd::LinSpaced(numX, 0, Nx-1);
    centersY = Eigen::VectorXd::LinSpaced(numY, 0, Ny-1);
    centersZ = Eigen::VectorXd::LinSpaced(numZ, 0, Nz-1);

    // map 3D "image" to 1D "inputPts"
        // order matters! z -> y -> x
    inputPts.resize(N);
    t = 0;
    for (auto it=image.begin(); it!=image.end(); it++) {
        inputPts.segment(t*Nx*Ny, Nx*Ny) << Eigen::Map<const Eigen::VectorXd>(it->data(), Nx*Ny);
        t++;
    }

    // Calculate & fill the least square matrix A
        PrintTime("Calculating & filling least square matrix...");
    CalcLeastSquareMat(A, Nx, Ny, Nz);
        std::cout << "Least square matrix A size = " << A.rows() << " * " << A.cols() << std::endl;
        std::cout << "A: # non-zero elements = " << A.nonZeros() << std::endl;

    // calculate control points based on least square
        PrintTime("Calculating A'*A and A'*y...");
    AtransposeA = (A.transpose() * A).pruned();  // A' * A
        std::cout << "A' * A size = " << AtransposeA.rows() << " * " << AtransposeA.cols() << std::endl;
        std::cout << "A'A: # non-zero elements = " << AtransposeA.nonZeros() << std::endl;
    vectorY.resize(N);
    vectorY = A.transpose() * inputPts;  // A' * y

    // solve linear system for control points
        // Eigen::SimplicialCholesky<Eigen::SparseMatrix<double> > chol(AtransposeA);  // performs Cholesky factorization
        // controlPoints = chol.solve(vectorY);
    const std::string solverName = "Hypre";
    auto solver = polysolve::LinearSolver::create(solverName, "");
    const nlohmann::json params = {
        // {"max_iter", 500}, 
        // {"tolerance", 1e-6}
    };
    solver->setParameters(params);
        PrintTime("Analyzing matrix pattern...");
    solver->analyzePattern(AtransposeA, AtransposeA.rows());
        PrintTime("Factorizing matrix...");
    solver->factorize(AtransposeA);
        PrintTime("Solving linear system...");
    controlPoints.resize(num, 1);
    solver->solve(vectorY, controlPoints);
    /////////////// DEBUG ONLY ///////////
    // std::cout << ">>>>> control points >>>>>" << std::endl;
    // std::cout << controlPoints << std::endl;
    /////////////// DEBUG ONLY ///////////
        PrintTime("Control points calculated...");
        std::cout << "====================================================" << std::endl;
}


void bspline::CalcLeastSquareMat(Eigen::SparseMatrix<double> &A, int Nx, int Ny, int Nz) {

    int count, i, j, ix, iy, iz, jx, jy, jz, refIdx_x, refIdx_y, refIdx_z;
    int jxMin, jxMax, jyMin, jyMax, jzMin, jzMax;
    double rowSum, t;
    std::array<int, 64> cache_i, cache_j;
    std::array<double, 64> cache_t;

    // order matters!
    // for column (j): z -> x -> y
    // for row (i): z -> y -> x
    for (iz=0; iz<=Nz-1; iz++)
        for (iy=0; iy<=Ny-1; iy++)
            for (ix=0; ix<=Nx-1; ix++) {

                i = iz * Nx * Ny + iy * Nx + ix;  // iterating sample points
                if (i % 3000 == 0) std::cout << " " << i << " / " << Nx*Ny*Nz << std::endl;
                
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
                            t = CubicBasis( (ix-centersX(jx)) / gapX) * 
                                CubicBasis( (iy-centersY(jy)) / gapY) * 
                                CubicBasis( (iz-centersZ(jz)) / gapZ);

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
                    A.insert(cache_i[i], cache_j[i]) = cache_t[i] / rowSum;
                }
            }
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

                    xcoef = CubicBasis( (sample(i, 0)-centersX(idx_x)) / gapX );
                    ycoef = CubicBasis( (sample(i, 1)-centersY(idx_y)) / gapY );
                    zcoef = CubicBasis( (sample(i, 2)-centersZ(idx_z)) / gapZ );

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


void bspline::Interp3D(const Eigen::Matrix<DScalar, Eigen::Dynamic, 3> &sampleDS, Eigen::Matrix<DScalar, Eigen::Dynamic, 1> &res) const {

    assert(numX != 0 && numY != 0 && numZ != 0);

    int i, ix, iy, iz, idx_x, idx_y, idx_z;
    DScalar xcoef, ycoef, zcoef;
    int refIdx_x, refIdx_y, refIdx_z;
    double evalX1, evalX2, evalX3, evalX4, evalY1, evalY2, evalY3, evalY4, evalZ1, evalZ2, evalZ3, evalZ4;
    Eigen::MatrixXd sample(sampleDS.rows(), 3);
    for (i=0; i<sampleDS.rows(); i++) {
        sample(i, 0) = sampleDS(i, 0).getValue();
        sample(i, 1) = sampleDS(i, 1).getValue();
        sample(i, 2) = sampleDS(i, 2).getValue();
    }
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
        res(i) = DScalar(0);
        for (ix=0; ix<=3; ix++)
            for (iy=0; iy<=3; iy++)
                for (iz=0; iz<=3; iz++) {

                    idx_x = refIdx_x + ix;
                    idx_y = refIdx_y + iy;
                    idx_z = refIdx_z + iz;

                    xcoef = CubicBasis( (sampleDS(i, 0)-centersX(idx_x)) / gapX );
                    ycoef = CubicBasis( (sampleDS(i, 1)-centersY(idx_y)) / gapY );
                    zcoef = CubicBasis( (sampleDS(i, 2)-centersZ(idx_z)) / gapZ );

                    res(i) += controlPoints(numX*numY*idx_z + numY*idx_x + idx_y) * xcoef * ycoef * zcoef;
                }
    }
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
