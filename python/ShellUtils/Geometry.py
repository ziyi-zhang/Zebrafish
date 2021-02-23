import numpy as np
import igl


def GenDualShell(Vin, Fin):

    Nvert = Vin.shape[0]
    Nvert_per_surface = int(Nvert / 4)

    cntSingular = 0
    cntClose = 0
    for i in range(Fin.shape[0]):
        # for the i-th triangle
        flag = False
        for j in range(3):
            vi = Fin[i][j]
            thres = 1e-15
            if (np.linalg.norm(Vin[vi, :] - Vin[vi+Nvert_per_surface, :]) < thres):
                flag = True
                if ((Vin[vi, :] == Vin[vi+Nvert_per_surface, :]).all()):
                    cntSingular += 1
                else:
                    print(np.linalg.norm(Vin[vi, :] - Vin[vi+Nvert_per_surface, :]))
                    # Vin[vi, :] = Vin[vi+Nvert_per_surface, :]
            if (np.linalg.norm(Vin[vi+Nvert_per_surface*2, :] - Vin[vi+Nvert_per_surface*3, :]) < thres):
                flag = True
                if ((Vin[vi+Nvert_per_surface*2, :] == Vin[vi+Nvert_per_surface*3, :]).all()):
                    cntSingular += 1
                else:
                    print(np.linalg.norm(Vin[vi+Nvert_per_surface*2, :] - Vin[vi+Nvert_per_surface*3, :]))
                    # Vin[vi+Nvert_per_surface*2, :] = Vin[vi+Nvert_per_surface*3, :]
        if (flag):
            cntClose += 1
            print(">>>> Triangle #%d" % (i))
            print("  %.16f  %.16f  %.16f" % (Vin[vi, 0], Vin[vi, 1], Vin[vi, 2]))
            print("  %.16f  %.16f  %.16f" % (Vin[vi+Nvert_per_surface, 0], Vin[vi+Nvert_per_surface, 1], Vin[vi+Nvert_per_surface, 2]))
            print("  %.16f  %.16f  %.16f" % (Vin[vi+Nvert_per_surface*2, 0], Vin[vi+Nvert_per_surface*2, 1], Vin[vi+Nvert_per_surface*2, 2]))
            print("  %.16f  %.16f  %.16f" % (Vin[vi+Nvert_per_surface*3, 0], Vin[vi+Nvert_per_surface*3, 1], Vin[vi+Nvert_per_surface*3, 2]))
            print("<<<<")
    print("#Singularities = %d" % cntSingular)
    print("#Close vertices = %d" % cntClose)

    F = np.concatenate((Fin, Fin+Nvert_per_surface), axis=0)
    F = np.concatenate((F,   Fin+Nvert_per_surface*2), axis=0)
    F = np.concatenate((F,   Fin+Nvert_per_surface*3), axis=0)

    return Vin, F


def PointInTet(pt0, pt1, pt2, pt3, pt_test):

    flip = igl.orient3d(pt0, pt1, pt2, pt3)
    assert(flip == -1)

    return (igl.orient3d(pt0, pt3, pt1, pt_test) <= 0) and (igl.orient3d(pt1, pt3, pt2, pt_test) <= 0) and (igl.orient3d(pt0, pt1, pt2, pt_test) <= 0) and (igl.orient3d(pt0, pt2, pt3, pt_test) <= 0)


def f2e(F):
    return np.vstack([F[:, [0, 2]], F[:, [0, 1]], F[:, [1, 2]]])


def t2e(T):
    return np.vstack([T[:, [0, 1]], T[:, [0, 2]], T[:, [0, 3]], T[:, [1, 2]], T[:, [1, 3]], T[:, [2, 3]]])
