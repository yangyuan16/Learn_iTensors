# tJ Square lattice on cylinder with symmetry condition
# 控制粒子数 和 total spin
# 考虑具有最近邻和次近邻相互作用的格子
#
using ITensors
using ITensors.HDF5
#
function tune_dopping(N, dop_level)
    # 设置一个加了对称性的初始 state
    #state = [isodd(n) ? "Up" : "Dn" for n = 1:N]
    #=
    调整掺杂浓度，可以这样设置，这里我们
    将粒子放置在每4个站点中的第一个和第二个站点上，
    这样，每4个格点中有两个格点有粒子。相当于 1/4 粒子浓度。
    可以更改该这里的代码逻辑以放置更多或更少的粒子。
    =#

    if dop_level  == 0
        state = [isodd(n) ? "Up" : "Dn" for n = 1:N] # 半满的状态， 一个格点上一个spin
        return state
    elseif N % dop_level > 0
        println("N / dop_level must be an integer")
        sqrt(-1)
    else
        state = fill("a", N)
        for j = 1:N
            if j % dop_level == 1
                state[j] = "Up"
            elseif j % dop_level == 2
                state[j] = "Dn"
            else
                state[j] = "Emp"
            end
        end
        return state
    end  
end
#
function add_NNN(Ly, Lx)  # add NNN bonds
    #
    bonds1 = fill(0, ((Ly-1) * (Lx-1), 2))
    bonds2 = fill(0, ((Ly-1) * (Lx-1), 2))
    for i = 1 : (Lx-1)
        s0 = (i - 1) * Ly + 1
        for j = 1:(Ly-1)
            bonds1[(i-1) * (Ly-1) + j, 1] = s0 + j - 1
            bonds1[(i-1) * (Ly-1) + j, 2] = s0 + j - 1 + (Ly + 1)

            bonds2[(i-1) * (Ly-1) + j, 1] = s0 + j 
            bonds2[(i-1) * (Ly-1) + j, 2] = s0 + j + (Ly-1)
        end
    end
    @show bonds1
    @show bonds2

    bonds3 = fill(0, (Lx-1, 2))
    bonds4 = fill(0, (Lx-1, 2))
    for i = 1 : Lx - 1
        bonds3[i,1] = (i-1) * Ly + 1
        bonds3[i,2] = (i-1) * Ly + Ly + Ly 

        bonds4[i,1] = i * Ly
        bonds4[i,2] = i * Ly + 1
    end
    @show bonds3
    @show bonds4

    b12 = cat(bonds1,bonds2,dims=1)
    b123 = cat(b12,bonds3, dims=1)
    bonds = cat(b123, bonds4, dims=1)
    return bonds
end
#
function rundmrg()
    bc = "cylinder"
    Ly = 4
    Lx = 5
    N = Lx * Ly
    t = 1.0
    J = 1.0 / 3
    dop_level = 4
    #
    sites = siteinds("tJ",N; conserve_qns = true) # 加上 conserve_qns 条件
    # 2D sqaure wrapped on a cylinder
    lattice = square_lattice(Lx,Ly, yperiodic = true)
    NNN_bonds = add_NNN(Ly, Lx)
    @show NNN_bonds
    #
    os = OpSum()
    for b in lattice
        os .+= -t, "Cdagup", b.s1, "Cup", b.s2
        os .+= -t, "Cdagup", b.s2, "Cup", b.s1
        os .+= -t, "Cdagdn", b.s1, "Cdn", b.s2
        os .+= -t, "Cdagdn", b.s2, "Cdn", b.s1
        os .+= J, "Sz", b.s1, "Sz", b.s2
        os .+= 0.5 * J, "S+", b.s1, "S-", b.s2
        os .+= 0.5 * J, "S-", b.s1, "S+", b.s2 
        os .+= -0.25 * J, "Ntot", b.s1, "Ntot", b.s2  
    end
    for j = 1:size(NNN_bonds)[1]
        os .+= -t, "Cdagup", NNN_bonds[j,1], "Cup", NNN_bonds[j,2]
        os .+= -t, "Cdagup", NNN_bonds[j,2], "Cup", NNN_bonds[j,1]
        
        os .+= -t, "Cdagdn", NNN_bonds[j,1], "Cdn", NNN_bonds[j,2]
        os .+= -t, "Cdagdn", NNN_bonds[j,2], "Cdn", NNN_bonds[j,1]

        os .+= J, "Sz", NNN_bonds[j,1], "Sz", NNN_bonds[j,2]
        os .+= 0.5 * J, "S+", NNN_bonds[j,1], "S-", NNN_bonds[j,2]
        os .+= 0.5 * J, "S-", NNN_bonds[j,1], "S+", NNN_bonds[j,2] 
        os .+= -0.25 * J, "Ntot", NNN_bonds[j,1], "Ntot", NNN_bonds[j,2]
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
    
    @show (Ly, Lx, t, J, dop_level)
    
    per_energy = energy / N
    @show per_energy
    #
    workpath = "./tj/Ly-Lx-t-t-j/output/psi/"
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
