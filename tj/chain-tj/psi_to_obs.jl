using ITensors
using ITensors.HDF5
# 通过基态波函数计算观测值
#
let 
    N = 100
    maxdim = [400]
    bc = "obc"   # boundary condition 
    loadpath = "./tj/chain-tj/output/psi/"
    filename_psi_1 = "1d_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([loadpath, filename_psi_1, filename_psi_2, filename_psi_3])

    f = h5open(filename_psi, "r")
    psi = read(f, "psi", MPS)
    close(f)

    # calculate the magz
    magz = expect(psi, "Sz")     # 测量 magz 
    @show magz

    #
    savepath = "./tj/chain-tj/output/obs/magz/"
    filename_magz_1 = "1d_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_magz_2 = bc
    filename_magz_3 = "_magz.h5"
    filename_magz = join([savepath, filename_magz_1, filename_magz_2, filename_magz_3])

    f = h5open(filename_magz, "w")
    write(f, "magz", magz)
    close(f)

    # calculate the Nup
    Nup = expect(psi, "Nup")
    @show Nup
    savepath = "./tj/chain-tj/output/obs/Nup/"
    filename_Nup_1 = "1d_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_Nup_2 = bc
    filename_Nup_3 = "_Nup.h5"
    filename_Nup = join([savepath, filename_Nup_1, filename_Nup_2, filename_Nup_3])

    f = h5open(filename_Nup, "w")
    write(f, "Nup", Nup)
    close(f)

    # calculate the Ndn
    Ndn = expect(psi, "Ndn")
    @show Ndn
    savepath = "./tj/chain-tj/output/obs/Ndn/"
    filename_Ndn_1 = "1d_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_Ndn_2 = bc
    filename_Ndn_3 = "filename_Ndn_1.h5"
    filename_Ndn = join([savepath, filename_Ndn_1, filename_Ndn_2, filename_Ndn_3])

    f = h5open(filename_Ndn, "w")
    write(f, "Ndn", Ndn)
    close(f)

    # calculate the Ntot
    Ntot = expect(psi, "Ntot")
    @show Ntot
    savepath = "./tj/chain-tj/output/obs/Ntot/"
    filename_Ntot_1 = "1d_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_Ntot_2 = bc
    filename_Ntot_3 = "filename_Ntot.h5"
    filename_Ntot = join([savepath, filename_Ntot_1, filename_Ntot_2, filename_Ntot_3])

    f = h5open(filename_Ntot, "w")
    write(f, "Ntot", Ntot)
    close(f)   
end
