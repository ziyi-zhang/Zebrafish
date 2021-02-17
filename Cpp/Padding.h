#pragma once

////////////////////////////////////////////////////////////////////////////////
#include <zebrafish/Common.h>

#include <Eigen/Dense>
#include <map>
#include <array>
////////////////////////////////////////////////////////////////////////////////

namespace zebrafish {

typedef std::map<int, std::array<int, 2> > RCMap_t;

////////////////////////////////////////////////
// padding [static class]

class padding {

private:

public:
    static void ComputeOneRing(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const RCMap_t &RCMap, Eigen::MatrixXd &appendV, Eigen::MatrixXi &appendF, RCMap_t &appendRCMap);

    template <typename T>
    static void AddOneRing(const Eigen::MatrixXd &appendV, const Eigen::MatrixXi &appendF, T &V, Eigen::MatrixXi &F);
    static void AddOneRingForAll(const Eigen::MatrixXd &appendV, const Eigen::MatrixXi &appendF, const RCMap_t &appendRCMap, std::vector<Eigen::MatrixXd> &V, Eigen::MatrixXi &F, RCMap_t &RCMap);

    static void Harmonic(Eigen::MatrixXd &V, const Eigen::MatrixXi &F, int rawMeshVerts, int lastPadIndex, int order);
    static void HarmonicForAll(std::vector<Eigen::MatrixXd> &V, const Eigen::MatrixXi &F, int rawMeshVerts, int lastPadIndex, int order = 2);
};

}  // namespace zebrafish
