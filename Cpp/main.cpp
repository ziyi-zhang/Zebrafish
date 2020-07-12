#include <zebrafish/ICP.h>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
using namespace Eigen;
using namespace zebrafish;

int main() {

    MatrixXd p, q;
    RMat_t R, R0;
    TMat_t T, T0;

    p = MatrixXd::Random(3, 10);

    const double theta = M_PI / 10.0;
    R0 << cos(theta), -sin(theta), 0, 
          sin(theta),  cos(theta), 0, 
          0, 0, 1;

    T0 << 0.02, 0.03, 0.05;

    q = R0 * p;
    q.colwise() += T0;

    ICP icpSolver;
    icpSolver.RunICP(p, q, R, T);

    cout << ">>>> R0 = " << endl << R0 << endl;
    cout << ">>>> T0 = " << endl << T0 << endl;
    cout << ">>>> R = " << endl << R << endl;
    cout << ">>>> T = " << endl << T << endl;

    cout << "R max err = " << (R0 - R).maxCoeff() << endl;
    cout << "T max err = " << (T0 - T).maxCoeff() << endl;

    // cout << p << endl;
    // cout << q << endl;
}
