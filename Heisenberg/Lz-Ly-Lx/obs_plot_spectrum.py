import os, sys
os.chdir(sys.path[0])
#print('sys.path:', str(sys.path))
import h5py
import numpy as np
import matplotlib.pyplot as plt
#
def plot_e_spectrum(x,y):
    plt.plot(x,y,"-", marker = "o",markersize=5 )
    plt.title("spectrum")
    plt.xlabel("chi")
    plt.ylabel("e_spectrum")
    plt.show()
    return
#
if __name__ == "__main__":
    Lz = 2
    Ly = 3
    Lx = 4
    N = Lz * Ly * Lx
    maxdim = 200
    bc = "obc"
    #
    # load entanglement spectrum
    loadpath = "./output/obs/e_spectrum/"
    filename_spectrum_1 = "Lz%d_Ly%d_Lx%d_N%d_S%g_maxdim%d_" % (Lz, Ly, Lx, N, 0.5,maxdim)
    filename_spectrum_2 = bc
    filename_spectrum_3 = "_e_spectrum_site%d.h5"%(int(N/2))
    filename_spectrum = loadpath + filename_spectrum_1 + \
        filename_spectrum_2 + filename_spectrum_3
    hdfFile = h5py.File(filename_spectrum, 'r')
    data = hdfFile.get('spectrum')
    spectrum = list(data)
    #print("spectrum:", spectrum)
    hdfFile.close()
    plot_e_spectrum(x=range(len(spectrum[0])),y=spectrum[0]) 