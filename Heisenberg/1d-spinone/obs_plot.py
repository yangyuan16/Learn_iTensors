import os, sys
os.chdir(sys.path[0])
#print('sys.path:', str(sys.path))
import h5py
import numpy
import matplotlib.pyplot as plt

def plot_magz(x,y):
    plt.plot(x,y,"-", marker = "o",markersize=5)
    plt.title("magz")
    plt.xlabel("sites")
    plt.ylabel("Sz")
    plt.show()
    return
#
def plot_magx(x,y):
    plt.plot(x,y,"-", marker = "o",markersize=5 )
    plt.title("magx")
    plt.xlabel("sites")
    plt.ylabel("Sx")
    plt.show()
    return
#
def plot_corr_z(x,y, matrix):
    plt.plot(x,y,"-", marker = "o",markersize=5)
    plt.title("corr_z")
    plt.xlabel("sites")
    plt.ylabel("<SzSz>")
    plt.show()
    plt.imshow(matrix)
    plt.title("corr_z")
    plt.colorbar()
    plt.show()
    return

def plot_corr_spsm(x,y, matrix):
    plt.plot(x,y,"-", marker = "o",markersize=5 )
    plt.title("corr_spsm")
    plt.xlabel("sites")
    plt.ylabel("<SpSm>")
    plt.show()
    plt.imshow(matrix)
    plt.title("corr_spsm")
    plt.colorbar()
    plt.show()
    return

def plot_corr_smsp(x,y,matrix):
    plt.plot(x,y,"-", marker = "o",markersize=5)
    plt.title("corr_smsp")
    plt.xlabel("sites")
    plt.ylabel("<SmSp>")
    plt.show()
    plt.imshow(matrix)
    plt.title("corr_smsp")
    plt.colorbar()
    plt.show()
    return

def plot_entropy(x,y):
    plt.plot(x,y,"-", marker = "o",markersize=5 )
    plt.title("entropy")
    plt.xlabel("sites")
    plt.ylabel("Entropy")
    plt.show()
    return

def plot_e_spectrum(x,y):
    plt.plot(x,y,"-", marker = "o",markersize=5 )
    plt.title("spectrum")
    plt.xlabel("chi")
    plt.ylabel("e_spectrum")
    plt.show()
    return 

if __name__ == "__main__":
    print()
    N = 100
    maxdim = 200
    bc = "ob"
    #
    # load magz:
    #
    loadpath =  "./output/obs/magz/"
    filename_magz_1 = "1dchain_N%d_S1_maxdim%d_" % (N, maxdim)
    filename_magz_2 = bc
    filename_magz_3 = "_magz.h5"
    filename_magz = loadpath + filename_magz_1 + filename_magz_2 + filename_magz_3
    hdfFile = h5py.File(filename_magz, 'r')
    data = hdfFile.get('magz')
    magz = list(data)
    hdfFile.close()
    plot_magz(x=range(len(magz)),y=magz)
    #
    # load magx:
    #
    loadpath =  "./output/obs/magx/"
    filename_magx_1 = "1dchain_N%d_S1_maxdim%d_" % (N, maxdim)
    filename_magx_2 = bc
    filename_magx_3 = "_magx.h5"
    filename_magx = loadpath + filename_magx_1 + filename_magx_2 + filename_magx_3
    hdfFile = h5py.File(filename_magx, 'r')
    data = hdfFile.get('magx')
    magx = list(data)
    hdfFile.close()
    plot_magx(x=range(len(magx)),y=magx)
    #
    # load corr_z
    #     
    loadpath = "./output/obs/corr_z/"
    filename_corr_z_1 = "1dchain_N%d_S1_maxdim%d_" % (N, maxdim)
    filename_corr_z_2 = bc
    filename_corr_z_3 = "_corr_z.h5"
    filename_corr_z = loadpath + filename_corr_z_1 + filename_corr_z_2 + filename_corr_z_3
    hdfFile = h5py.File(filename_corr_z, "r")
    data = hdfFile.get("corr_z")
    corr_z = list(data)
    hdfFile.close()
    s1 = 20
    s2 = 80
    plot_corr_z(x=range(20,80),y=corr_z[s1][s1:s2],matrix=corr_z)
    #
    # load corr_spsm
    #
    loadpath = "./output/obs/corr_spsm/"
    filename_corr_spsm_1 = "1dchain_N%d_S1_maxdim%d_" % (N, maxdim)
    filename_corr_spsm_2 = bc
    filename_corr_spsm_3 = "_corr_spsm.h5"
    filename_corr_spsm = loadpath + filename_corr_spsm_1 + \
        filename_corr_spsm_2 + filename_corr_spsm_3
    hdfFile = h5py.File(filename_corr_spsm, "r")
    data = hdfFile.get("corr_spsm")
    corr_spsm = list(data)
    hdfFile.close()
    s1 = 20
    s2 = 80
    plot_corr_spsm(x=range(20,80),y=corr_spsm[s1][s1:s2],matrix=corr_spsm)
    #
    # load corr_smsp
    #
    loadpath = "./output/obs/corr_smsp/"
    filename_corr_smsp_1 = "1dchain_N%d_S1_maxdim%d_" % (N, maxdim)
    filename_corr_smsp_2 = bc
    filename_corr_smsp_3 = "_corr_smsp.h5"
    filename_corr_smsp = loadpath + filename_corr_smsp_1 + \
        filename_corr_smsp_2 + filename_corr_smsp_3
    hdfFile = h5py.File(filename_corr_smsp, "r")
    data = hdfFile.get("corr_smsp")
    corr_smsp = list(data)
    hdfFile.close()
    s1 = 20
    s2 = 80
    plot_corr_smsp(x=range(20,80),y=corr_smsp[s1][s1:s2],matrix=corr_smsp)
    #
    # load entropy
    # 
    loadpath =  "./output/obs/entropy/"
    filename_entropy_1 = "1dchain_N%d_S1_maxdim%d_" % (N, maxdim)
    filename_entropy_2 = bc
    filename_entropy_3 = "_entropy.h5"
    filename_entropy = loadpath + filename_entropy_1 + filename_entropy_2 + filename_entropy_3
    hdfFile = h5py.File(filename_entropy, 'r')
    data = hdfFile.get('entropy')
    entropy = list(data)
    hdfFile.close()
    plot_entropy(x=entropy[0],y=entropy[1])    
    #
    # load spectrum
    #
    loadpath = "./output/obs/e_spectrum/"
    filename_spectrum_1 = "1dchain_N%d_S1_maxdim%d_" % (N, maxdim)
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



