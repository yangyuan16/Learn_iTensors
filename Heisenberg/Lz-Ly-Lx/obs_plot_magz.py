# plot the magz of layered squared lattice
#
import os, sys
os.chdir(sys.path[0])
#print('sys.path:', str(sys.path))
import h5py
import numpy as np
import matplotlib.pyplot as plt
#
def plot_dots(coordinates, magz_value, title, magz_max, magz_min, magz_mean):
    fig, ax = plt.subplots(figsize=(9,9))
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    plt.axis("off")
    for i1 in range(len(coordinates)):
        for i2 in range(len(coordinates)):
            if np.sqrt((coordinates[i1][0] - coordinates[i2][0])**2 + (coordinates[i1][1]-coordinates[i2][1])**2) < 1.1:
                ax.plot([coordinates[i1][0], coordinates[i2][0]], [coordinates[i1][1], coordinates[i2][1]], "-k", linewidth=1)
    #
    for i in range(len(coordinates)):
        #ax.plot(coordinates[i][0], coordinates[i][1], "ro", markersize=10)
        if magz_value[i] > 0:
            marker = r"${\uparrow}$"
            color = "red"
            mz_norm = (magz_value[i] - magz_min) / (magz_max - magz_min)
            markersize = mz_norm * 40
            ax.plot(coordinates[i][0], coordinates[i][1], marker = marker, 
                    markersize = markersize, color= color )
        elif magz_value[i] < 0:
            marker = r"${\downarrow}$"
            color = "green"
            mz_norm = (magz_value[i] - magz_min) / (magz_max - magz_min)
            markersize = mz_norm * 40
            ax.plot(coordinates[i][0], coordinates[i][1], marker = marker, 
                    markersize = markersize, color= color )
        else:
            raise Exception("magzations is zero")
    plt.text(0.5, 0.5, "magz_max=%d"%(magz_max))
    plt.text(0.5, 0.75, "magz_min=%d"%(magz_min))
    plt.text(0.5, 1, "magz_mean=%d"%(magz_mean))
    plt.title(title)
    #plt.savefig("square.eps")
    plt.show()
    return
#
def magz_numbered(magz):  
    c2 = magz
    c1 = []
    for i in range(len(c2)):
        c1.append(i+1)
    magz_n = [c1,c2]
    return magz_n
#
def magz_layerd(magz_n):
    magz_layer1 = []
    it_layer1 = []
    magz_layer2 = []
    it_layer2 = []
    #
    for it in range(len(magz_n[0])):
        if it % 2 == 0:
            it_layer1.append(magz_n[0][it])
            magz_layer1.append(magz_n[1][it])
        elif it % 2 == 1:
            it_layer2.append(magz_n[0][it])
            magz_layer2.append(magz_n[1][it])
        else:
            raise Exception("wrong input of magz_n")        
    magz_1 = [it_layer1, magz_layer1]
    magz_2 = [it_layer2, magz_layer2]
    return magz_1, magz_2
#
def magz_coord(magz_1, magz_2, Ly, Lx):
    #magz_1 = []
    coor = []
    for i in range(Lx):
        for j in range(Ly):
            coor.append((i,j))
    #
    magz_1_coor = [magz_1, coor]
    magz_2_coor = [magz_2, coor]
    #
    return magz_1_coor, magz_2_coor
#
def plot_layer1(magz):
    print(magz[1])
    coordinates = magz[1]
    title = "first layer"
    magz_value = magz[0][1]
    magz_max = np.array(magz_value).max()
    magz_min = np.array(magz_value).min()
    magz_mean = np.array(magz_value).mean()
    plot_dots(coordinates=coordinates, magz_value=magz_value, title=title,
                magz_max=magz_max, magz_min=magz_min, magz_mean=magz_mean)
    return
#
def plot_layer2(magz):
    print(magz[1])
    coordinates = magz[1]
    title = "second layer"
    magz_value = magz[0][1]
    magz_max = np.array(magz_value).max()
    magz_min = np.array(magz_value).min()
    magz_mean = np.array(magz_value).mean()
    plot_dots(coordinates=coordinates, magz_value=magz_value, title=title,
                magz_max=magz_max, magz_min=magz_min, magz_mean=magz_mean)
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
    # load magz:
    #
    loadpath =  "./output/obs/magz/"
    filename_magz_1 = "Lz%d_Ly%d_Lx%d_N%d_S%g_maxdim%d_" % (Lz, Ly, Lx, N, 0.5,maxdim)
    filename_magz_2 = bc
    filename_magz_3 = "_magz.h5"
    filename_magz = loadpath + filename_magz_1 + filename_magz_2 + filename_magz_3
    hdfFile = h5py.File(filename_magz, 'r')
    data = hdfFile.get('magz')
    magz = list(data)
    hdfFile.close()
    #
    # group magz into layer1 and layer2  
    #
    magz_n = magz_numbered(magz=magz) # number the magz
    magz_1, magz_2 = magz_layerd(magz_n=magz_n) # two layers of magz
    m1, m2 = magz_coord(magz_1=magz_1, magz_2=magz_2, Ly=Ly, Lx=Lx) # add coordinates 
    print("m1:\n", m1)
    print("m2:\n", m2)
    #
    plot_layer1(magz=m1)
    plot_layer2(magz=m2)