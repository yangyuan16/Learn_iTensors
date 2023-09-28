# tJ Square lattice on cylinder with symmetry condition
# 控制粒子数 和 total spin
# 只考虑最近邻的正方格子
using ITensors
using ITensors.HDF5
#
function tune_dopping(N, dop_level)
    L_cell = 2 * dop_level[2]
    L_emp = 2 * dop_level[1]
    N_cell = div(N, L_cell)
    @show N_cell
    if N % L_cell > 0
        println("N / L_cell must be an integer")
        sqrt(-1)
    else
        state = fill("a", N)
        L_spin = L_cell - L_emp
        @show L_spin
        if L_spin % 2 > 0
            print("L_spin must be an even to keep total Sz=0")
            sqrt(-1)
        else
            for it = 1:N_cell 
                for j = 1:L_spin
                    if j % 2 ==1
                        state[(it-1) * L_cell + j] = "Up"
                        #@show (it-1) * L_cell + j
                    else
                        state[(it-1) * L_cell + j] = "Dn"
                        #@show (it-1) * L_cell + j
                    end
                end
                for jj = L_spin+1 : L_cell
                    state[(it-1) * L_cell + jj] = "Emp"
                    #@show (it-1) * L_cell + jj
                end
            end
        end
        return state
    end
end
#
function rundmrg()
    bc = "cylinder"
    Lx = 4
    Ly = 4
    N = Lx * Ly
    t = 1.0
    J = 1.0 / 3
    dop_level = (1,4)
    #
    sites = siteinds("tJ",N; conserve_qns = true) # 加上 conserve_qns 条件
    # 2D sqaure wrapped on a cylinder
    lattice = square_lattice(Lx,Ly, yperiodic = true)
    
    os = OpSum()
    for b in lattice
        @show b
        os .+= -t, "Cdagup", b.s1, "Cup", b.s2
        os .+= -t, "Cdagup", b.s2, "Cup", b.s1
        os .+= -t, "Cdagdn", b.s1, "Cdn", b.s2
        os .+= -t, "Cdagdn", b.s2, "Cdn", b.s1
        os .+= J, "Sz", b.s1, "Sz", b.s2
        os .+= 0.5 * J, "S+", b.s1, "S-", b.s2
        os .+= 0.5 * J, "S-", b.s1, "S+", b.s2 
        os .+= -0.25 * J, "Ntot", b.s1, "Ntot", b.s2  
    end
    
    H = MPO(os, sites)
    
    state = tune_dopping(N, dop_level)
    @show state
    #psi0 = randomMPS(sites, state, L)
    psi0 = MPS(sites, state)

    nsweeps = 10
    maxdim = [20, 60, 100, 100, 200, 400]
    cutoff = [1E-4]

    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
    @show energy
    #
    per_energy = energy / N
    @show per_energy
    #
    workpath = "./tj/lattice-tj/output/psi/"
    filename_psi_1 = "Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_t$(t)_J$(round(J;digits=4))_maxdim$(last(maxdim))_dop$(dop_level)_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    f = h5open(filename_psi, "w")
    write(f, "psi", psi)
    close(f)
    return
end
#
let 
    @time rundmrg()
end
#
