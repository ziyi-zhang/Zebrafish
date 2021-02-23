import numpy as np
import Geometry
import igl


def GenHexagonShell(factor=1.2):

    Vin = np.array([
        [0.5, np.sqrt(3)/2],
        [1, 0],
        [0.5, -np.sqrt(3)/2],
        [-0.5, -np.sqrt(3)/2],
        [-1, 0],
        [-0.5, np.sqrt(3)/2]
    ])
    Vout = Vin * factor
    # np.concatenate((V0, V0*factor), axis=0)
    Fin = np.array([
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 5],
        [5, 0]
    ])
    Fout = Fin
    # np.concatenate((F0, F0+6), axis=0)

    return Vin, Fin, Vout, Fout


def GenDistortedShell():

    Vin = np.array([
        [1, 3],
        [-3, 3],
        [-1, -0.5],
        [-2, -np.sqrt(3)*3],
        [4, -2],
        [2.5, 1]
    ])
    Vout = np.array([
        [1.2, 4],
        [-3.5, 3.3],
        [-2, -0.5],
        [-2.15, -np.sqrt(3)*3-0.25],
        [4.8, -3],
        [3.5, 1.2]
    ])
    Vin = Vin * 0.2
    Vout = Vout * 0.2
    # np.concatenate((V0, V0*factor), axis=0)
    Fin = np.array([
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 4],
        [4, 5],
        [5, 0]
    ])
    Fout = Fin
    # np.concatenate((F0, F0+6), axis=0)

    return Vin, Fin, Vout, Fout


def GenCubeObj():

    Vraw, Fraw = igl.read_triangle_mesh('./data/cube.obj')
    randTransMat = np.array([[-0.700390490817757,  -0.150045418195264,  -0.697810527901858],
                            [-0.114407726497277,  -0.941413431486692,   0.317256399674106],
                            [-0.704531072763853,   0.302038281505938,   0.642190660174339]])
    # Vraw = randTransMat.dot(Vraw.transpose()).transpose()  # rotate
    V_shell = np.concatenate((Vraw, Vraw*1.2), axis=0)
    V_shell = np.concatenate((V_shell, Vraw*1.6), axis=0)
    V_shell = np.concatenate((V_shell, Vraw*1.8), axis=0)
    # F_shell = np.concatenate((Fraw, Fraw+Vraw.shape[0]), axis=0)
    # F_shell = np.concatenate((F_shell, Fraw+2*Vraw.shape[0]), axis=0)
    # F_shell = np.concatenate((F_shell, Fraw+3*Vraw.shape[0]), axis=0)
    Nv = Vraw.shape[0]
    for i in range(Fraw.shape[0]):
        # 1 (split_b)
        Fraw[i, :] = np.sort(Fraw[i, :])
        pts = np.array([
            V_shell[Fraw[i, 0], :],
            V_shell[Fraw[i, 1], :],
            V_shell[Fraw[i, 2], :],
            V_shell[Fraw[i, 0]+Nv, :],
            V_shell[Fraw[i, 1]+Nv, :],
            V_shell[Fraw[i, 2]+Nv, :]])
        if ((igl.orient3d(pts[0, :], pts[3, :], pts[4, :], pts[5, :]) == 1) or
            (igl.orient3d(pts[0, :], pts[5, :], pts[4, :], pts[1, :]) == 1) or
            (igl.orient3d(pts[0, :], pts[1, :], pts[2, :], pts[5, :]) == 1)):
            pass
        else:
            print("{} = 1 split_b".format(i))
            continue
        # 2 (split_a)
        t = Fraw[i, 1]
        Fraw[i, 1] = Fraw[i, 2]
        Fraw[i, 2] = t
        pts = np.array([
            V_shell[Fraw[i, 0], :],
            V_shell[Fraw[i, 1], :],
            V_shell[Fraw[i, 2], :],
            V_shell[Fraw[i, 0]+Nv, :],
            V_shell[Fraw[i, 1]+Nv, :],
            V_shell[Fraw[i, 2]+Nv, :]])
        if ((igl.orient3d(pts[0, :], pts[3, :], pts[4, :], pts[5, :]) == 1) or
            (igl.orient3d(pts[0, :], pts[5, :], pts[4, :], pts[2, :]) == 1) or
            (igl.orient3d(pts[0, :], pts[1, :], pts[2, :], pts[4, :]) == 1)):
            print("split a and b both wrong: i={}" % (i))
        else:
            print("{} = 2 split_a".format(i))
            continue
    F_shell = Fraw
    print(F_shell)

    # igl.write_obj("./data/dualShell_cube.obj", V_shell, F_shell)
    return V_shell, F_shell


def GenSphereObj():

    Vraw, Fraw = igl.read_triangle_mesh('./data/sphere.obj')
    V_shell = np.concatenate((Vraw, Vraw*1.2), axis=0)
    V_shell = np.concatenate((V_shell, Vraw*1.6), axis=0)
    V_shell = np.concatenate((V_shell, Vraw*1.8), axis=0)
    F_shell = np.concatenate((Fraw, Fraw+Vraw.shape[0]), axis=0)
    F_shell = np.concatenate((F_shell, Fraw+2*Vraw.shape[0]), axis=0)
    F_shell = np.concatenate((F_shell, Fraw+3*Vraw.shape[0]), axis=0)

    # igl.write_obj("./data/dualShell_sphere.obj", V_shell, F_shell)
    return V_shell, F_shell


def ply2obj(filename):
    # Load Zhongshi's h5 model and save as std ply input
    Vin, Fin = igl.read_triangle_mesh('./dataset/' + filename + '.stl.h5.ply')  # './dataset/100070.stl.h5.ply'

    V, F = Geometry.GenDualShell(Vin, Fin)
    # igl.write_triangle_mesh("./dataset/dualShell_" + filename + ".ply", V, F, False)
    print("F = %dx%d" % (F.shape[0], F.shape[1]))
    # print("./dataset/dualShell_" + filename + ".ply")
    return V, F
