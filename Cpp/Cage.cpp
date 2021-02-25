#include <zebrafish/Common.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Cage.h>
#include <igl/boundary_loop.h>

#include <unordered_set>
#include <array>
#include <vector>
#include <algorithm>

namespace zebrafish
{

    namespace
    {
    } // anonymous namespace

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // class cage

    void cage::ComputeCage(
        const Eigen::MatrixXd &V,
        const Eigen::MatrixXi &F,
        Eigen::MatrixXd &sideV,
        Eigen::MatrixXi &sideF,
        bool above)
    {

        // reserve space for vertices in V
        const int N = V.rows();
        sideV.resize(N * 2, 3);
        sideV << V, V;
        // reserve space for sideF
        Eigen::VectorXi boundary;
        igl::boundary_loop(F, boundary);
        sideF.resize(F.rows() + boundary.size() * 2, 3);
        // flat bottom/top
        Eigen::MatrixXd box_min = V.colwise().minCoeff();
        Eigen::MatrixXd box_max = V.colwise().maxCoeff();
        const double diag = (box_max - box_min).norm();
        const double meanZ = V.col(2).mean();
        const double targetZ = meanZ + (above ? 0.5 : -0.5) * diag;
        sideV.block(N, 2, N, 1) = Eigen::ArrayXd::Constant(N, targetZ);
        sideF.topRows(F.rows()) = F.array() + N;
        // side surface
        int cnt = F.rows();
        for (int i = 0; i < boundary.size(); i++)
        {
            int base0 = boundary(i);
            int base1 = boundary((i + 1) % boundary.size());
            sideF.row(cnt++) << base0 + N, base0, base1;
            sideF.row(cnt++) << base1, base1 + N, base0 + N;
        }
    }

    void cage::AddCage(
        const Eigen::MatrixXd &sideV,
        const Eigen::MatrixXi &sideF,
        Eigen::MatrixXd &V,
        Eigen::MatrixXi &F,
        int baseV,
        bool updateF)
    {

        const int Nv = V.rows();
        // int baseV = std::min(Nv, int(sideV.rows()));
        /*
    for (int i=0; i<baseV; i++) {
        if (!V.row(i).isApprox(sideV.row(i), 0.0)) {
            baseV = i;
            break;
        }
    }
    */
        if (baseV == 0)
        {
            logger().error("AddCage: vertCnt not properly registered");
            return;
        }

        // V
        V.conservativeResize(Nv + sideV.rows() - baseV, 3);
        V.bottomRows(sideV.rows() - baseV) = sideV.bottomRows(sideV.rows() - baseV);
        // F
        if (updateF)
        {
            F.conservativeResize(F.rows() + sideF.rows(), 3);
            F.bottomRows(sideF.rows()) = sideF;
            for (int i = F.rows() - sideF.rows(); i < F.rows(); i++)
                for (int j = 0; j < 3; j++)
                {
                    if (F(i, j) >= baseV)
                        F(i, j) += Nv - baseV;
                }
        }
    }

    bool cage::AddCageForAll(
        const Eigen::MatrixXd &sideV,
        const Eigen::MatrixXi &sideF,
        std::vector<Eigen::MatrixXd> &V,
        Eigen::MatrixXi &F,
        int baseV)
    {

        if (sideV.size() == 0 || sideF.size() == 0)
            return false;
        int frames = V.size();

        for (int i = 0; i < frames; i++)
        {
            cage::AddCage(sideV, sideF, V[i], F, baseV, i == 0);
        }
        return true;
    }

} // namespace zebrafish
