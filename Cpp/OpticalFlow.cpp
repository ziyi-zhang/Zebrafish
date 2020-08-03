#include <zebrafish/OpticalFlow.h>

#include <iostream>
#include <algorithm>
#include <vector>


namespace zebrafish {

double GetMatValue(const std::vector<Eigen::MatrixXd> &mat, int x, int y, int z) {

    const int layer = mat.size();
    if (z == -1 || z == layer) return 0.0;

    const int rows = mat[0].rows();
    const int cols = mat[0].cols();
    if (x == -1 || x == rows) return 0.0;
    if (y == -1 || y == cols) return 0.0;

    return mat[z](x, y);
}


void Initialize3DMatrix(std::vector<Eigen::MatrixXd> &mat, int layer, int rows, int cols) {
// Initialize "mat" to be zero with size [rows x cols x layer]

    mat.clear();
    for (int i=0; i<layer; i++)
        mat.push_back(Eigen::MatrixXd::Zero(rows, cols));
}


std::vector<Eigen::MatrixXd> operator*(const std::vector<Eigen::MatrixXd> &matA, const std::vector<Eigen::MatrixXd> &matB) {
// Elementwise product

    int i;
    Eigen::MatrixXd temp;
    std::vector<Eigen::MatrixXd> res;

    // matA and matB must have the same dimension
    assert(matA.size() == matB.size());
    assert(matA[0].rows() == matB[0].rows());
    assert(matA[0].cols() == matB[0].cols());

    res.clear();
    for (i=0; i<matA.size(); i++) {

        temp = matA[i].cwiseProduct(matB[i]);
        res.push_back(temp);
    }

    return res;
}


std::vector<Eigen::MatrixXd> operator/(const std::vector<Eigen::MatrixXd> &matA, const std::vector<Eigen::MatrixXd> &matB) {
// elementwise division

    int i;
    Eigen::MatrixXd temp;
    std::vector<Eigen::MatrixXd> res;

    // matA and matB must have the same dimension
    assert(matA.size() == matB.size());
    assert(matA[0].rows() == matB[0].rows());
    assert(matA[0].cols() == matB[0].cols());

    res.clear();
    for (i=0; i<matA.size(); i++) {

        temp = matA[i].cwiseProduct(matB[i].cwiseInverse());
        res.push_back(temp);
    }

    return res;
}


std::vector<Eigen::MatrixXd> operator+(const std::vector<Eigen::MatrixXd> &matA, const std::vector<Eigen::MatrixXd> &matB) {
// elementwise addition

    int i;
    Eigen::MatrixXd temp;
    std::vector<Eigen::MatrixXd> res;

    // matA and matB must have the same dimension
    assert(matA.size() == matB.size());
    assert(matA[0].rows() == matB[0].rows());
    assert(matA[0].cols() == matB[0].cols());

    res.clear();
    for (i=0; i<matA.size(); i++) {

        temp = matA[i] + matB[i];
        res.push_back(temp);
    }

    return res;
}


std::vector<Eigen::MatrixXd> operator-(const std::vector<Eigen::MatrixXd> &matA, const std::vector<Eigen::MatrixXd> &matB) {
// elementwise subtraction

    int i;
    Eigen::MatrixXd temp;
    std::vector<Eigen::MatrixXd> res;

    // matA and matB must have the same dimension
    assert(matA.size() == matB.size());
    assert(matA[0].rows() == matB[0].rows());
    assert(matA[0].cols() == matB[0].cols());

    res.clear();
    for (i=0; i<matA.size(); i++) {

        temp = matA[i] - matB[i];
        res.push_back(temp);
    }

    return res;
}


void OpticalFlow::ComputePartialDerivative(const image_t &img1, const image_t &img2, PartialDerivative_t &dx, PartialDerivative_t &dy, PartialDerivative_t &dz, PartialDerivative_t &dt) {

    assert(img1.size() > 0 && img2.size() > 0);
    const int layer = img1.size();
    const int rows = img1[0].rows();
    const int cols = img1[0].cols();
    int x, y, z;

    Initialize3DMatrix(dx, layer, rows, cols);
    Initialize3DMatrix(dy, layer, rows, cols);
    Initialize3DMatrix(dz, layer, rows, cols);
    Initialize3DMatrix(dt, layer, rows, cols);

    for (x=0; x<rows; x++)
        for (y=0; y<cols; y++)
            for (z=0; z<layer; z++) {

                // consider an 8-element cube [x-1:x] * [y-1:y] * [z-1:z] for both img1 & img2
                dx[z](x, y) = (1.0 / 8.0) * (
                    GetMatValue(img1, x  , y  , z  ) - GetMatValue(img1, x-1, y  , z  ) +
                    GetMatValue(img1, x  , y-1, z  ) - GetMatValue(img1, x-1, y-1, z  ) +
                    GetMatValue(img1, x  , y  , z-1) - GetMatValue(img1, x-1, y  , z-1) +
                    GetMatValue(img1, x  , y-1, z-1) - GetMatValue(img1, x-1, y-1, z-1) +
                    GetMatValue(img2, x  , y  , z  ) - GetMatValue(img2, x-1, y  , z  ) +
                    GetMatValue(img2, x  , y-1, z  ) - GetMatValue(img2, x-1, y-1, z  ) +
                    GetMatValue(img2, x  , y  , z-1) - GetMatValue(img2, x-1, y  , z-1) +
                    GetMatValue(img2, x  , y-1, z-1) - GetMatValue(img2, x-1, y-1, z-1)
                );

                dy[z](x, y) = (1.0 / 8.0) * (
                    GetMatValue(img1, x  , y  , z  ) - GetMatValue(img1, x  , y-1, z  ) + 
                    GetMatValue(img1, x-1, y  , z  ) - GetMatValue(img1, x-1, y-1, z  ) + 
                    GetMatValue(img1, x  , y  , z-1) - GetMatValue(img1, x  , y-1, z-1) + 
                    GetMatValue(img1, x-1, y  , z-1) - GetMatValue(img1, x  , y-1, z-1) + 
                    GetMatValue(img2, x  , y  , z  ) - GetMatValue(img2, x  , y-1, z  ) + 
                    GetMatValue(img2, x-1, y  , z  ) - GetMatValue(img2, x-1, y-1, z  ) + 
                    GetMatValue(img2, x  , y  , z-1) - GetMatValue(img2, x  , y-1, z-1) + 
                    GetMatValue(img2, x-1, y  , z-1) - GetMatValue(img2, x  , y-1, z-1)
                );

                dz[z](x, y) = (1.0 / 8.0) * (
                    GetMatValue(img1, x  , y  , z  ) - GetMatValue(img1, x  , y  , z-1) +
                    GetMatValue(img1, x-1, y  , z  ) - GetMatValue(img1, x-1, y  , z-1) +
                    GetMatValue(img1, x  , y-1, z  ) - GetMatValue(img1, x  , y-1, z-1) +
                    GetMatValue(img1, x-1, y-1, z  ) - GetMatValue(img1, x-1, y-1, z-1) +
                    GetMatValue(img2, x  , y  , z  ) - GetMatValue(img2, x  , y  , z-1) +
                    GetMatValue(img2, x-1, y  , z  ) - GetMatValue(img2, x-1, y  , z-1) +
                    GetMatValue(img2, x  , y-1, z  ) - GetMatValue(img2, x  , y-1, z-1) +
                    GetMatValue(img2, x-1, y-1, z  ) - GetMatValue(img2, x-1, y-1, z-1)
                );

                dt[z](x, y) = (1.0 / 8.0) * (
                    GetMatValue(img2, x  , y  , z  ) - GetMatValue(img1, x  , y  , z  ) + 
                    GetMatValue(img2, x-1, y  , z  ) - GetMatValue(img1, x-1, y  , z  ) + 
                    GetMatValue(img2, x  , y-1, z  ) - GetMatValue(img1, x  , y-1, z  ) + 
                    GetMatValue(img2, x  , y  , z-1) - GetMatValue(img1, x  , y  , z-1) + 
                    GetMatValue(img2, x-1, y-1, z  ) - GetMatValue(img1, x-1, y-1, z  ) + 
                    GetMatValue(img2, x-1, y  , z-1) - GetMatValue(img1, x-1, y  , z-1) + 
                    GetMatValue(img2, x  , y-1, z-1) - GetMatValue(img1, x  , y-1, z-1) + 
                    GetMatValue(img2, x-1, y-1, z-1) - GetMatValue(img1, x-1, y-1, z-1)
                );

                // multiply by -1
                dx[z](x, y) *= -1.0;
                dy[z](x, y) *= -1.0;
                dz[z](x, y) *= -1.0;
                dt[z](x, y) *= -1.0;
            }
}


void OpticalFlow::ComputeAverageFlow(const FlowVelocity_t &u, FlowVelocity_t &u_bar) {

    assert(u.size() > 0);
    const int layer = u.size();
    const int rows = u[0].rows();
    const int cols = u[0].cols();
    int x, y, z;

    Initialize3DMatrix(u_bar, layer, rows, cols);

    for (x=0; x<rows; x++)
        for (y=0; y<cols; y++)
            for (z=0; z<layer; z++) {
            
                // consider a 27-element cube [x-1:x+1] * [y-1:y+1] * [z-1:z+1]
                u_bar[z](x, y) = (1.0 / 6.0) * (
                    GetMatValue(u, x-1, y  , z  ) +
                    GetMatValue(u, x+1, y  , z  ) +
                    GetMatValue(u, x  , y-1, z  ) +
                    GetMatValue(u, x  , y+1, z  ) +
                    GetMatValue(u, x  , y  , z-1) +
                    GetMatValue(u, x  , y  , z+1)
                );
            }
}


void OpticalFlow::RunOpticalFlow(const image_t &img1, const image_t &img2, double alpha, int iter, FlowVelocity_t &ux, FlowVelocity_t &uy, FlowVelocity_t &uz) {

    PartialDerivative_t dx, dy, dz, dt, dx_square, dy_square, dz_square;
    FlowVelocity_t ux_bar, uy_bar, uz_bar;
    std::vector<Eigen::MatrixXd> constMat;
    const int layer = img1.size();
    const int rows = img1[0].rows();
    const int cols = img1[0].cols();
    int i;

    Initialize3DMatrix(ux, layer, rows, cols);
    Initialize3DMatrix(uy, layer, rows, cols);
    Initialize3DMatrix(uz, layer, rows, cols);
    Initialize3DMatrix(constMat, layer, rows, cols);
    for (i=0; i<layer; i++)
        constMat[i] = Eigen::MatrixXd::Constant(rows, cols, alpha * alpha);  // alpha^2 constant matrix

    // partial derivative
    ComputePartialDerivative(img1, img2, dx, dy, dz, dt);

    // temporory variable
    dx_square = dx * dx;
    dy_square = dy * dy;
    dz_square = dz * dz;

    // iterative solution
    std::cout << "iter res" << std::endl;
    std::vector<Eigen::MatrixXd> numerator, denominator;
    for (i=0; i<iter; i++) {

        ComputeAverageFlow(ux, ux_bar);
        ComputeAverageFlow(uy, uy_bar);
        ComputeAverageFlow(uz, uz_bar);

        numerator = dt + (dx * ux_bar) + (dy * uy_bar) + (dz * uz_bar);
        denominator = dx_square + dy_square + dz_square + constMat;
        ux = ux_bar - dx * (numerator / denominator);
        uy = uy_bar - dy * (numerator / denominator);
        uz = uz_bar - dz * (numerator / denominator);

        std::cout << ux[40](58, 25) << " " << uy[40](58, 25) << " " << uz[40](58, 25) << std::endl;
    }
}

////////////////////////////////////////////////////////////////////////////////////////
// maintenance methods


OpticalFlow::OpticalFlow() {

}


OpticalFlow::~OpticalFlow() {

}

}  // namespace zebrafish
