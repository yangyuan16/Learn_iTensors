import os, sys
os.chdir(sys.path[0])
#print('sys.path:', str(sys.path))
import h5py
import numpy as np
import matplotlib.pyplot as plt
#
def plot_entropy(x,y):
    plt.plot(x,y,"-", marker = "o",markersize=5 )
    plt.title("entropy")
    plt.xlabel("sites")
    plt.ylabel("Entropy")
    plt.show()
    return



if __name__ == "__main__":
    Lz = 2
    Ly = 3
    Lx = 4
    N = Lz * Ly * Lx
    maxdim = 200
    bc = "obc"
    #
    # load entropy
    #
    loadpath =  "./output/obs/entropy/"
    filename_entropy_1 = "Lz%d_Ly%d_Lx%d_N%d_S%g_maxdim%d_" % (Lz, Ly, Lx, N, 0.5,maxdim)
    filename_entropy_2 = bc
    filename_entropy_3 = "_entropy.h5"
    filename_entropy = loadpath + filename_entropy_1 + filename_entropy_2 + filename_entropy_3
    hdfFile = h5py.File(filename_entropy, 'r')
    data = hdfFile.get('entropy')
    entropy = list(data)
    hdfFile.close()
    print("entropy:", entropy)
    plot_entropy(x=entropy[0],y=entropy[1])
    

    