using ITensors
using ITensors.HDF5
# 通过基态波函数计算观测值

let 
    # load psi
    Lz = 2
    Ly = 3
    Lx = 4
    N = Lz * Ly * Lx
    maxdim = [200]
    bc = "obc"   # boundary condition 
    loadpath = "./Heisenberg/Lz-Ly-Lx/output/psi/"
    filename_psi_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([loadpath, filename_psi_1, filename_psi_2, filename_psi_3])

    f = h5open(filename_psi, "r")
    psi = read(f, "psi", MPS)
    close(f)
    #@show psi

    # calculate the magz
    magz = expect(psi, "Sz")     # 测量 magz 
    @show magz
    
    savepath = "./Heisenberg/Lz-Ly-Lx/output/obs/magz/"
    filename_magz_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_magz_2 = bc
    filename_magz_3 = "_magz.h5"
    filename_magz = join([savepath, filename_magz_1, filename_magz_2, filename_magz_3])

    f = h5open(filename_magz, "w")
    write(f, "magz", magz)
    close(f)

    # calculate the magx
    #=    Sz = 0 子空间得到的 state 不能测量 magx
    magx = expect(psi, "Sx")     # 测量 magx
    @show magx

    savepath = "./Heisenberg/Lz-Ly-Lx/output/obs/magx/"
    filename_magx_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_magx_2 = bc
    filename_magx_3 = "_magx.h5"
    filename_magx = join([savepath, filename_magx_1, filename_magx_2, filename_magx_3])

    f = h5open(filename_magx, "w")
    write(f, "magx", magx)
    close(f)
    =#
    
    # calculate the correlation 
    corr_z = correlation_matrix(psi,"Sz", "Sz")
    corr_spsm = correlation_matrix(psi, "S+", "S-")
    corr_smsp = correlation_matrix(psi, "S-", "S+")
    #@show corr_z
    #@show corr_spsm
    #@show corr_smsp

    savepath = "./Heisenberg/Lz-Ly-Lx/output/obs/corr_z/"
    filename_corr_z_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_corr_z_2 = bc
    filename_corr_z_3 = "_corr_z.h5"
    filename_corr_z = join([savepath, filename_corr_z_1, filename_corr_z_2, filename_corr_z_3])
    f = h5open(filename_corr_z, "w")
    write(f, "corr_z", corr_z)
    close(f)

    savepath = "./Heisenberg/Lz-Ly-Lx/output/obs/corr_spsm/"
    filename_corr_spsm_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_corr_spsm_2 = bc
    filename_corr_spsm_3 = "_corr_spsm.h5"
    filename_corr_spsm = join([savepath, filename_corr_spsm_1, filename_corr_spsm_2, 
                                filename_corr_spsm_3])
    f = h5open(filename_corr_spsm, "w")
    write(f, "corr_spsm", corr_spsm)
    close(f)

    savepath = "./Heisenberg/Lz-Ly-Lx/output/obs/corr_smsp/"
    filename_corr_smsp_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_corr_smsp_2 = bc
    filename_corr_smsp_3 = "_corr_smsp.h5"
    filename_corr_smsp = join([savepath, filename_corr_smsp_1, filename_corr_smsp_2, 
                                filename_corr_smsp_3])
    f = h5open(filename_corr_smsp, "w")
    write(f, "corr_smsp", corr_smsp)
    close(f)

    # calculate the entanglement_spectrum and  entanglement_entropy
    
    entropy = fill(0.0,(N-2,2))
    for b = 2:N-1
        orthogonalize!(psi, b) # 将格点 b 作为正交中心
        U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
        SvN = 0.0
        spectrum = fill(0.0, (dim(S,1),1))
        for n = 1:dim(S,1)
            p = S[n,n]^2
            SvN -= p * log(p)
            spectrum[n] = S[n,n]
        end
        # save entanglement spectrum
        savepath = "./Heisenberg/Lz-Ly-Lx/output/obs/e_spectrum/"
        filename_e_spectrum_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
        filename_e_spectrum_2 = bc
        filename_e_spectrum_3 = "_e_spectrum_site$(b).h5"
        filename_e_spectrum = join([savepath, filename_e_spectrum_1, filename_e_spectrum_2, 
                                    filename_e_spectrum_3])
        f = h5open(filename_e_spectrum, "w")
        write(f, "spectrum", spectrum)
        close(f)
        #
        entropy[b-1,1] = b
        entropy[b-1,2] = SvN
    end

    # save entropy
    savepath = "./Heisenberg/Lz-Ly-Lx/output/obs/entropy/"
    filename_entropy_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_entropy_2 = bc
    filename_entropy_3 = "_entropy.h5"
    filename_entropy = join([savepath, filename_entropy_1, filename_entropy_2, 
                                filename_entropy_3])
    f = h5open(filename_entropy, "w")
    write(f, "entropy", entropy)
    close(f)

    return
end