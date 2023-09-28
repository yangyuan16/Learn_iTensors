# 一维 t-J 模型
#
using ITensors
using ITensors.HDF5
#
function chain_obc(N)
    b = fill(0, (N-1, 2))
    for i =1:N-1
        b[i,1] = i
        b[i,2] = i+1
    end
    @show b
end
#
function chain_pbc(N)
    b = fill(0, (N, 2))
    for i =1:N-1
        b[i,1] = i
        b[i,2] = i + 1
    end
    b[N, 1] = 1
    b[N, 2] = N
    @show b 
end
#
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
    bc = "pbc" 
    t = 1 
    J = 1
    N = 100
    dop_level = (1,5)
    sites = siteinds("tJ",N; conserve_qns = true) # 加上 conserve_qns 条件
    #    
    if bc == "obc"
        b = chain_obc(N)
    else
        b = chain_pbc(N)
    end
    #
    os = OpSum()
    for j = 1:size(b)[1]
        os += -t, "Cdagup", b[j,1], "Cup", b[j,2]
        os += -t, "Cdagup", b[j,2], "Cup", b[j,1]
        os += -t, "Cdagdn", b[j,1], "Cdn", b[j,2]
        os += -t, "Cdagdn", b[j,2], "Cdn", b[j,1]
        os += J, "Sz", b[j,1], "Sz", b[j,2]
        os += 0.5 * J, "S+", b[j,1], "S-", b[j,2]
        os += 0.5 * J, "S-", b[j,1], "S+", b[j,2] 
        os += -0.25 * J, "Ntot", b[j,1], "Ntot", b[j,2]   
    end
    H = MPO(os, sites)
    # 设置加了对称性的初始态
    #state = [isodd(n) ? "Up" : "Dn" for n = 1:N] # 半满的状态， 一个格点上一个spin
    state = tune_dopping(N, dop_level)  # 调整 dopping 浓度
    @show state
    # 设置初始态
    #psi0 = randomMPS(sites, state, N)
    psi0 = MPS(sites, state)
    #
    nsweeps = 10
    maxdim = [20, 60, 100, 100, 200, 400]
    cutoff = [1E-4]
    #
    energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff)
    @show energy
    #
    per_energy = energy / N
    @show per_energy
    #
    workpath = "./tj/chain-tj/output/psi/"
    filename_psi_1 = "1d_N$(N)_S$(0.5)_maxdim$(last(maxdim))_dop$(dop_level)_"
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
# maxdim: 400; cutoff: 1E-4; energy: -68.8577438