import os, sys
os.chdir(sys.path[0])
#print('sys.path:', str(sys.path))
import h5py
import numpy as np
import matplotlib.pyplot as plt
#
def plot_bonds(corrz_layer, title):
    #
    corrz = corrz_layer["corrz"]
    site_pair = corrz_layer["site_pair"]
    coord = corrz_layer["coord"]
    site = corrz_layer["site"]
    #
    fig, ax = plt.subplots(figsize=(9,9))
    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.95, top=0.95)
    plt.axis("off")
    coordinates = coord
    for i1 in range(len(coordinates)):
        for i2 in range(len(coordinates)):
            if np.sqrt((coordinates[i1][0] - coordinates[i2][0])**2 + (coordinates[i1][1]-coordinates[i2][1])**2) < 1.1:
                coord1 = (coordinates[i1][0], coordinates[i1][1])   # 得到相邻 bond 的第一个点的坐标
                coord2 = (coordinates[i2][0], coordinates[i2][1])   # 得到相邻 bond 的第二个点的坐标
                #print("(coord1, coord2):", (coord1, coord2))
                index1 = coordinates.index(coord1) # 得到相邻 bond 的第一个点的坐标索引
                index2 = coordinates.index(coord2) # 得到相邻 bond 的第二个点的坐标索引
                s1 = site[index1]  # 根据坐标索引得到相邻 bond 的第一个点的晶格格点序号
                s2 = site[index2]  # 根据坐标索引得到相邻 bond 的第二个点的晶格格点序号
                #print("(s1, s2):", (s1, s2))
                index_corrz = site_pair.index((s1,s2))  # 根据格点序号得到关联函数索引
                #print("index_corrz:", index_corrz)
                nn_corrz = corrz[index_corrz] # 根据索引得到满足条件的最近邻关联
                if nn_corrz > 0:
                    color = "red"
                else:
                    color = "green"
                ax.plot([coordinates[i1][0], coordinates[i2][0]], [coordinates[i1][1], coordinates[i2][1]], color=color, linewidth=3)
    for i in range(len(coordinates)):
        ax.plot(coordinates[i][0], coordinates[i][1], 's', color = "w", markersize=10)
    # plt.savefig('square.eps') 
    plt.title(title)
    #plt.savefig("square.eps")
    plt.show()
    return
#
def plot_corrz(x,y, title):
    plt.plot(x,y,"-", marker = "o",markersize=5)
    plt.title("corrz")
    plt.xlabel("sites")
    plt.ylabel("correlation")
    plt.title(title)
    plt.show()
    return 
#
def corr_z_layered(corr_z, Lz, Ly, Lx):
    corrz_1 = []  # layer 1
    corrz_2 = []  # layer 2
    corrz_12= [] # between two layers
    sites_layer1 = [] 
    sites_layer2 = []
    sites_layer12 = []
    #
    for it1 in range(Lz * Ly * Lx):
        for it2 in range(Lz * Ly *  Lx):
            if (it1 % 2 == 0) and (it2 % 2 ==0): 
                corrz_1.append(corr_z[it1, it2]) # layer 1
                sites_layer1.append((it1+1, it2+1))
            elif (it1 % 2 == 1) and (it2 % 2 ==1):
                corrz_2.append(corr_z[it1, it2]) # layer 2
                sites_layer2.append((it1+1, it2+1))
            else:
                corrz_12.append(corr_z[it1, it2]) # between layer 1 and layer 2
                sites_layer12.append((it1+1, it2+1))
    #
    assert(len(corrz_1)==len(sites_layer1))
    assert(len(corrz_2)==len(sites_layer2))
    assert(len(corrz_12)==len(sites_layer12))
    #
    coord = []
    for it1 in range(Lx):
        for it2 in range(Ly):
            coord.append((it1, it2))
    #
    si_layer1 = []
    si_layer2 = []
    for it in range(Ly * Lx):
        si_layer1.append(2 * it + 1)
        si_layer2.append(2 * it + 2)
    
    assert(len(si_layer1) == len(coord))
    assert(len(si_layer2) == len(coord)) 
    
    #
    corr_z_lay1 = {"corrz": corrz_1, "site_pair":sites_layer1, "coord": coord, "site": si_layer1}
    corr_z_lay2 = {"corrz": corrz_2, "site_pair":sites_layer2, "coord": coord, "site": si_layer2}
    corr_z_lay12 = {"corrz": corrz_12, "site_pair":sites_layer12}
    #
    return corr_z_lay1, corr_z_lay2, corr_z_lay12
#
def get_corrz_all_bonds(corrz_1, corrz_2, corrz_12):
    corrz_all_bonds = {}
    corrz_all_bonds["corrz"] = corrz_1["corrz"] + corrz_2["corrz"] + corrz_12["corrz"]
    corrz_all_bonds["site_pair"] = corrz_1["site_pair"] + corrz_2["site_pair"] + corrz_12["site_pair"]
    assert(len(corrz_all_bonds["corrz"]) == len(corrz_all_bonds["site_pair"]) )
    return corrz_all_bonds

def get_corrz_along_path(corrz, path):
    corrz_path = []
    for it in path:
        index = corrz["site_pair"].index(it)
        corrz_path.append(corrz["corrz"][index])
    return corrz_path 

#
if __name__ == "__main__":
    Lz = 2
    Ly = 3
    Lx = 4
    N = Lz * Ly * Lx
    maxdim = 200
    bc = "obc"
    #
    # load corrz:
    #
    loadpath =  "./output/obs/corr_z/"
    filename_corr_z_1 = "Lz%d_Ly%d_Lx%d_N%d_S%g_maxdim%d_" % (Lz, Ly, Lx, N, 0.5,maxdim)
    filename_corr_z_2 = bc
    filename_corr_z_3 = "_corr_z.h5"
    filename_corr_z = loadpath + filename_corr_z_1 + filename_corr_z_2 + filename_corr_z_3
    hdfFile = h5py.File(filename_corr_z, 'r')
    data = hdfFile.get('corr_z')
    corr_z = np.array(list(data))
    hdfFile.close()
    #
    print("corr_z:\n", corr_z)
    corr_z_lay1, corr_z_lay2, corr_z_lay12 = corr_z_layered(corr_z=corr_z, Lz=Lz, Ly=Ly, Lx=Lx)
    #
    print("corr_z_lay1:\n", corr_z_lay1["corrz"])
    print("len(corr_z_lay1):", len(corr_z_lay1["corrz"]))
    print("corr_z_lay1[site_pair]:", corr_z_lay1["site_pair"])
    print("corr_z_lay1[coord]:", corr_z_lay1["coord"])
    print("corr_z_lay1[site]", corr_z_lay1["site"])
    #
    print("corr_z_lay2:\n", corr_z_lay2["corrz"])
    print("len(corr_z_lay2):", len(corr_z_lay2["corrz"]))
    print("corr_z_lay2[site_pair]:", corr_z_lay2["site_pair"])
    print("corr_z_lay2[coord]:", corr_z_lay2["coord"])
    print("corr_z_lay2[site]:", corr_z_lay2["site"])
    #
    print("len(corr_z_lay12):", len(corr_z_lay12["corrz"]))
    print("corr_z_lay12[sites]:", corr_z_lay12["site_pair"])
    #
    # plot the figure
    #
    #plot_bonds(corrz_layer=corr_z_lay1, title="layer 1")
    #plot_bonds(corrz_layer=corr_z_lay2, title="layer 2")
    #
    # plot corrz along path1 
    corrz_all_bonds = get_corrz_all_bonds(corrz_1=corr_z_lay1, corrz_2=corr_z_lay2, corrz_12=corr_z_lay12)
    path1 = []
    fix_site = 1
    x = []
    for i in range(1,24,2):
        path1.append((fix_site, i))
        x.append(i)
    print("path1:", path1)
    corrz_along_path = get_corrz_along_path(corrz=corrz_all_bonds,path=path1)
    title = "fix_site = 1"
    plot_corrz(x=x, y=corrz_along_path, title=title)
    #
    # plot corrz along path2
    path2 = []
    fix_site = 2
    x = []
    for i in range(2,25,2):
        path2.append((fix_site, i))
        x.append(i)
    print("path2:", path2)
    corrz_along_path = get_corrz_along_path(corrz=corrz_all_bonds,path=path2)
    title = "fix_site = 2"
    plot_corrz(x=x, y=corrz_along_path, title=title)
    #
    # plot corrz along path3
    path3 = []
    fix_site = 3
    x = []
    for i in range(fix_site,25,1):
        path3.append((fix_site, i))
        x.append(i)
    print("path3", path3)
    corrz_along_path =  get_corrz_along_path(corrz=corrz_all_bonds,path=path3)
    title = "fix_site = 3"
    plot_corrz(x=x, y=corrz_along_path, title=title)
    #