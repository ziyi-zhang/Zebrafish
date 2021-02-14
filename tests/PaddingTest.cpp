// cylinder test
#include <zebrafish/Common.h>
#include <zebrafish/Logger.hpp>
#include <zebrafish/Padding.h>

#include <catch.hpp>

#include <math.h>
#include <stdlib.h>
#include <random>

using namespace std;
using namespace Eigen;
using namespace zebrafish;


const auto Initialize = [](Eigen::MatrixXd &V, Eigen::MatrixXi &F, RCMap_t &RCMap) {
    F.resize(8, 3);
    F << 8, 6, 5, 
         6, 7, 3, 
         4, 5, 2, 
         5, 3, 1, 
         4, 8, 5, 
         5, 6, 3, 
         2, 5, 1, 
         1, 3, 0;
    double s3 = std::sqrt(3);
    double s32 = s3 * 2.0;
    V.resize(9, 3);
    V << 5, 0, 1, 
         3, 0, 1, 
         1, 0, 1, 
         4, s3, 1, 
         0, s3, 1, 
         2, s3, 1, 
         3, s32, 1, 
         5, s32, 1, 
         1, s32, 1;
    RCMap.insert({0, {2, 3}});
    RCMap.insert({1, {2, 2}});
    RCMap.insert({2, {2, 1}});
    RCMap.insert({3, {1, 2}});
    RCMap.insert({4, {1, 0}});
    RCMap.insert({5, {1, 1}});
    RCMap.insert({6, {0, 2}});
    RCMap.insert({7, {7, 3}});
    RCMap.insert({8, {0, 1}});
};


////////////////////////////////////////////////////////////////////////////

TEST_CASE("Padding_9", "[PaddingTest]") {

    Eigen::MatrixXd V, appendV;
    Eigen::MatrixXi F, appendF;
    RCMap_t RCMap;

    Initialize(V, F, RCMap);

    padding::AddOneRing(V, F, RCMap, appendV, appendF);
    cout << appendV << endl;
    cout << "appendF" << endl;
    cout << appendF << endl;
}
