import numpy as np


def GenHexagonShell(factor = 1.2):

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
