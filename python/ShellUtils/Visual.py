import numpy as np
import meshplot as mp


def PlotShell(Vin, Fin, Vout, Fout):

    V = np.concatenate((Vin, Vout))
    F = np.concatenate((Fin, Fout + Vin.shape[0]))
    plt = mp.plot(V, F, return_plot=True, shading={
        "line_color": "orange",
        "line_width": 3.0
    })
    return plt


def PlotMesh(plt, V, F):

    plt.add_mesh(V, F, shading={"wireframe": True})
    return


def PlotMeshColor(plt, V, F, cat):

    plt.add_mesh(V, F, c=cat, shading={"wireframe": True})
    return


def Explode(tetV, tetT, alpha=0.5):

    VT = tetV[tetT]
    mean = VT.mean(axis=1, keepdims=True)
    V = (VT - mean)*alpha + mean
    mp.plot(V.reshape(-1, 3), np.arange(4*len(tetT)).reshape(-1, 4), shading={"wireframe": True})
