#include <zebrafish/Quantile.h>

#include <queue>
#include <math.h>

namespace zebrafish {

double QuantileImage(const zebrafish::image_t &image, double q) {

    int i, j, depth;
    int N = image.size();
    int M = image[0].size();  // rows() * cols()
    int num = floor((1.0-q) * N * M);  // target number for quantile(q)
    std::priority_queue<double, std::vector<double>, std::greater<double> > heap;

    for (depth=0; depth<N; depth++) {

        const Eigen::MatrixXd &layer = image[depth];  // alias
        for (i=0; i<layer.size(); i++) {

            if (heap.size() < num)
                heap.push(layer(i));
            else if (layer(i) > heap.top()) {
                heap.pop();
                heap.push(layer(i));
            }
        }
    }

    return heap.top();
}

}  // namespace zebrafish
