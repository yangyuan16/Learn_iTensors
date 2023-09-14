# 画出观测量
#
using ITensors
using ITensors.HDF5
using GR
#=
plot the magz
=#
let 
    N = 100
    maxdim = [200]
    bc = "ob"

    # plot magz
    loadpath = "./chain-Heisenberg/output/obs/magz/"
    filename_magz_1 = "1dchain_N$(N)_S1_maxdim$(last(maxdim))_"
    filename_magz_2 = bc
    filename_magz_3 = "_magz.h5"
    filename_magz = join([loadpath, filename_magz_1, filename_magz_2, filename_magz_3])

    f = h5open(filename_magz, "r")
    magz = read(f, "magz")
    close(f)
    #@show magz
    sites = [i for i = 1:size(magz)[1]]
    
    #
    x = sites
    y = magz
    legend("magz")
    xlabel("sites")
    ylabel("Sz")
    #xlim([0,N])
    #ylim([-1,1])
    title("Heisenberg chain 1d")
    plot((x,y))
    
end

#=
plot the correlation 
=#

