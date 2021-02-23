import numpy as np
import igl
import Geometry
from shapely.geometry import Polygon, Point


def InsertTri(V, F, x0, y0, x1, y1, x2, y2):

    n = V.shape[0]
    V = np.append(V, np.array([[x0, y0], [x1, y1], [x2, y2]]), axis=0)
    F = np.append(F, np.array([[n, n+1, n+2]]), axis=0)
    return V, F


def BackgroundMesh(Vout, maxArea = 0.05):

    xLen = max(Vout[:, 0]) - min(Vout[:, 0])
    yLen = max(Vout[:, 1]) - min(Vout[:, 1])
    L = max(xLen, yLen)
    t = L / np.sqrt(3)
    x0 = min(Vout[:, 0])
    y0 = min(Vout[:, 1])

    # prepare the mesh array
    V = np.empty((0, 2), float)
    F = np.empty((0, 3), int)

    # insert the first bounding triangle
    V, F = InsertTri(V, F, x0-t, y0, x0+L+t, y0, x0+L/2.0, y0+np.sqrt(3) * (t+L/2.0))

    # upsample
    area = np.square(L/2.0 + t) * np.sqrt(3);
    upsampleTime = 0
    while (area > maxArea):
        area /= 4.0
        upsampleTime += 1
    V, F = igl.upsample(V, F, upsampleTime)

    return V, F


def BoundaryCheck(V, F, V_shell_out):

    N = F.shape[0]
    mask = np.full((N), True)
    poly = Polygon(V_shell_out)

    for i in range(0, N):
        # i-th triangle
        pt1 = Point(V[F[i, 0], 0], V[F[i, 0], 1])
        pt2 = Point(V[F[i, 1], 0], V[F[i, 1], 1])
        pt3 = Point(V[F[i, 2], 0], V[F[i, 2], 1])
        if not (pt1.within(poly) or pt2.within(poly) or pt3.within(poly)):
            mask[i] = False  # erase this triangle
    F = F[mask, :]

    return V, F


def UpsampleTriangle(pt1, pt2, pt3):
    
    pt12 = (pt1 + pt2) / 2.0;
    pt13 = (pt1 + pt3) / 2.0;
    pt23 = (pt2 + pt3) / 2.0;
    V = np.array([pt1, pt2, pt3, pt12, pt13, pt23])
    F = np.array([
        [0, 3, 4],
        [3, 1, 5],
        [3, 4, 5],
        [2, 4, 5]
    ])

    return V, F


def Carve(V, F, V_shell_out):
    
    # Discard triangles outside the outer shell
    V, F = BoundaryCheck(V, F, V_shell_out)
    
    N = F.shape[0]
    mask = np.full((N), True)  # whether this triangle will remain
    poly = Polygon(V_shell_out)
    
    newV = np.empty((0, 2), float)  # new triangles
    newF = np.empty((0, 3), int)
    for i in range(0, N):
        # i-th triangle
        pt1 = Point(V[F[i, 0], 0], V[F[i, 0], 1])
        pt2 = Point(V[F[i, 1], 0], V[F[i, 1], 1])
        pt3 = Point(V[F[i, 2], 0], V[F[i, 2], 1])
        if not (pt1.within(poly) or pt2.within(poly) or pt3.within(poly)):
            # Discard this triangle
            mask[i] = False
            triangle = Polygon(pt1, pt2, pt3)
            if triangle.area < 0.005:
                # This is too small to be upsampled
                pass
            else:
                # Upsample & carve it
                Vt, Ft = UpsampleTriangle(V[F[i, 0], :], V[F[i, 1], :], V[F[i, 2], :])
                Vt, Ft = Carve(Vt, Ft, V_shell_out)
                t = newV.shape[0]
                newV = np.concatenate((newV, Vt), axis=0)
                newF = np.concatenate((newF, Ft + t), axis=0)
            
    # Merge old triangles and new triangles
    F = F[mask, :]
    t = V.shape[0]
    V = np.concatenate((V, newV), axis=0)
    F = np.concatenate((F, newF + t), axis=0)
    
    return V, F