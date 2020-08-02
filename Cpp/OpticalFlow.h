#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
#include <string>
////////////////////////////////////////////////////////////////////////////////


namespace zebrafish {

typedef std::vector<Eigen::MatrixXd> FlowVelocity_t;  // a 3D double matrix
typedef std::vector<Eigen::MatrixXd> PartialDerivative_t;  // a 3D double matrix

///////////////////////////////////////
// 3D Optical Flow (Horn & Schunck)


class OpticalFlow {

private:
    static void ComputePartialDerivative(const image_t &img1, const image_t &img2, PartialDerivative_t &ix, PartialDerivative_t &iy, PartialDerivative_t &iz, PartialDerivative_t &it);
    /// Compute partial derivative with two 3D images for four dimensions
    static void ComputeAverageFlow(const FlowVelocity_t &u, FlowVelocity_t &u_bar);
    /// Compute the decentered mean for each pixel in u.
    /// This is used when estimating the Laplacian.

public:
    static void RunOpticalFlow(const image_t &img1, const image_t &img2, double alpha, int iter, FlowVelocity_t &ux, FlowVelocity_t &uy, FlowVelocity_t &uz);
    /// Run optical flow from 3D image img1 to img2.
    /// Minimize the Horn Schunck error with weighting factor alpha^2.
    ///
    /// @param[in]   img1        { first frame 3d image }
    /// @param[in]   img2        { second frame 3d image }
    /// @param[in]   alpha       { double weighting factor }
    /// @param[in]   iter        { iterations to calculate the flow velocity }
    /// @param[out]  ux, uy, uz  { 3d flow velocity in x(row), y(col), z(depth) direction }
    ///

    // maintenance methods
    OpticalFlow();
    ~OpticalFlow();
};

}  // namespace zebrafish
