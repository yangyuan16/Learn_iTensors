# Lz layerd squared lattice Heisenberg model
#
using ITensors
using ITensors.HDF5
# 
function lattice_obc(Lz, Ly, Lx)    
    b1 = fill(0, (Ly * Lx, 2))
    for i =1:Ly*Lx
        b1[i,1] = (2 * i - 1)
        b1[i,2] = 2 * i
    end
    @show b1
    #
    b2 = fill(0, ((Lz*(Ly-1))*Lx, 2))
    for i = 1:Lx
        s0 = (i-1) * Lz * Ly + 1
        for j = 1:(Lz * (Ly - 1))
            b2[(i-1) * (Lz * (Ly-1)) + j, 1] = s0 + (j-1)
            b2[(i-1) * (Lz * (Ly-1)) + j, 2] = s0 + (j+1)
        end
    end
    @show b2
    #
    b3 = fill(0, ((Lz * Ly) * (Lx-1),2))
    for i = 1:(Lz*Ly)*(Lx-1)
        b3[i, 1] = i 
        b3[i, 2] = i + (Lz * Ly)
    end
    @show b3
    
    b12 = cat(b1,b2,dims=1)
    b123 = cat(b12,b3, dims=1)
    
    @show b123

    return b123
end 

function lattice_pbc(Lz, Ly, Lx)
    b1 = fill(0, (Ly * Lx, 2))
    for i =1:Ly*Lx
        b1[i,1] = (2 * i - 1)
        b1[i,2] = 2 * i
    end
    @show b1
    #
    b2 = fill(0, ((Lz*(Ly-1))*Lx, 2))
    for i = 1:Lx
        s0 = (i-1) * Lz * Ly + 1
        for j = 1:(Lz * (Ly - 1))
            b2[(i-1) * (Lz * (Ly-1)) + j, 1] = s0 + (j-1)
            b2[(i-1) * (Lz * (Ly-1)) + j, 2] = s0 + (j+1)
        end
    end
    @show b2
    #
    b3 = fill(0, ((Lz * Ly) * (Lx-1),2))
    for i = 1:(Lz*Ly)*(Lx-1)
        b3[i, 1] = i 
        b3[i, 2] = i + (Lz * Ly)
    end
    @show b3
    #
    b4 = fill(0, (Lx * Lz, 2))  # consider the periodic boundary
    for i = 1:Lx
        s0 = (i - 1)*(Lz * Ly) + 1
        for j = 1:Lz
            b4[(i-1) * Lz + j, 1] = s0 + j - 1
            b4[(i-1) * Lz + j, 2] = s0 + j - 1 + Lz * (Ly - 1) 
        end 
    end
    @show b4

    b12 = cat(b1,b2,dims=1)
    b123 = cat(b12,b3, dims=1)
    b1234 = cat(b123, b4, dims=1)

    return b1234

end
#
function rundmrg() 
    bc = "obc" # ob: open boundary condition; pb: periodic boundary condition
    Lz = 2
    Ly = 3
    Lx = 4
    N = Lz * Ly * Lx
    
    sites = siteinds("S=1/2", N; conserve_qns=true) # S=1/2 spins

    if bc == "obc"
        b = lattice_obc(Lz,Ly,Lx)
        #@show b
    else
        b = lattice_pbc(Lz, Ly, Lx)
    end

    @show size(b)

    os = OpSum()
    for j = 1:size(b)[1]
        os += "Sz", b[j,1], "Sz", b[j,2]
        os += 0.5, "S+", b[j,1], "S-", b[j,2]
        os += 0.5, "S-", b[j,1], "S+", b[j,2]
    end
    H = MPO(os, sites)

    state = [isodd(n) ? "Up" : "Dn" for n = 1:N]
    psi0 = MPS(sites, state)
    @show flux(psi0)

    nsweeps = 5
    maxdim = [10,20,100,100,200]
    cutoff = [1E-10]

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)

    @show energy

    workpath = "./Heisenberg/Lz-Ly-Lx/output/psi/"
    filename_psi_1 = "Lz$(Lz)_Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_maxdim$(last(maxdim))_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    f = h5open(filename_psi, "w")
    write(f, "psi", psi)
    close(f)
    return
end

let 
    @time rundmrg()
end