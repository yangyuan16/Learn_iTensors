using ITensors
using ITensors.HDF5
let
    N = 100
    bc = "ob" # open boundary condition
    sites = siteinds("S=1",N) # S=1 spins
    @show sites[1]

    os = OpSum()
    for j=1:N-1
        os += "Sz",j,"Sz",j+1
        os += 1/2,"S+",j,"S-",j+1
        os += 1/2,"S-",j,"S+",j+1
    end
    H = MPO(os,sites)

    psi0 = randomMPS(sites,10) # 10 is the linkdims
    nsweeps = 5
    maxdim = [10,20,100,100,200]
    cutoff = [1E-10]

    energy, psi = dmrg(H,psi0; nsweeps, maxdim, cutoff)
    @show energy


    # 将 psi 保留下来
    @show psi
    workpath = "./chain-Heisenberg/output/psi/"
    filename_psi_1 = "1dchain_N$(N)_S1_maxdim$(last(maxdim))_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    f = h5open(filename_psi, "w")
    write(f, "psi", psi)
    close(f)
    
    return
    
end