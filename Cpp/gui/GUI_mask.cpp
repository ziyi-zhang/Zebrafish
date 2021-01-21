#include <zebrafish/GUI.h>
#include <zebrafish/Logger.hpp>

#include <cmath>


namespace zebrafish {

namespace {

}  // anonymous namespace

////////////////////////////////////////////////////////////////////////////////////////
// Special module to apply the membrane mask on either points or clusters

bool GUI::PointInMaskArea(double x, double y, double z) {

    y += 1;  // bottom center -> roughly cylinder center

    int xc = std::round(x);
    int yc = std::round(y);
    int zc = std::round(z);
    int validCnt = 0, ttlCnt = 0;
    for (int xi=xc-1; xi<=xc+1; xi++)
        for (int yi=yc-1; yi<=yc+1; yi++)
            for (int zi=zc-1; zi<=zc+1; zi++) {
                // whether in image cropped area
                if (xi>=0 && yi>=0 && zi>=layerBegin &&
                    xi<imgCols && yi<imgRows && zi<=layerEnd)
                    ttlCnt++;
                else
                    continue;
                // whether in mask
                if (membraneMask[zi](xi, yi) == 1)
                    validCnt++;
            }
    if (validCnt > 0)
        logger().debug("{}/{} ", validCnt, ttlCnt);
    if (double(validCnt) / double(ttlCnt) > 1.0/3.0)
        return true;
    else
        return false;
}

}  // namespace zebrafish
