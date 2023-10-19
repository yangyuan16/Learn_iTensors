# tJ Square lattice on cylinder with symmetry condition
# 控制粒子数 和 total spin
#
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
function build_sites(N,Is_conserve_qns)   
    if Is_conserve_qns == true
        sites = siteinds("tJ",N; conserve_qns = true, conserve_nf=true, conserve_sz=true, conserve_nfparity=true) # 加上 conserve_qns 条件
    elseif Is_conserve_qns == false
        sites = siteinds("tJ", N)
    else
        println("wrong input of Is_conserve_qns")
        sqrt(-1)
    end
    return sites
end
#
function get_psi0(N, dop_level, sites, Is_conserve_qns,way_of_psi0,loadfile)
    if way_of_psi0 == "default"
        println("get psi0 by the way of default")
        # get the initial state psi0
        if Is_conserve_qns == true
            state = tune_dopping(N, dop_level)
            @show state
            println("initial state psi0 is random but with conserved sites")
            psi0 = randomMPS(sites, state, N)
            println("setting initial psi0")
            @show flux(psi0)
        elseif Is_conserve_qns == false
            state = tune_dopping(N, dop_level)
            @show state
            println("initial state psi0 is random !!!<without>!!! conserved sites")
            psi0 = randomMPS(sites, state, N)
            println("setting initial psi0")
            @show flux(psi0)
        else
            println("wrong input of Is_conserve_qns")
            sqrt(-1)
        end
        return psi0
    elseif way_of_psi0 == "readfile"
        println("get psi0 by the way of readfile")
        println("filename is:")
        println(loadfile)
        f = h5open(loadfile, "r")
        psi0 = read(f, "psi", MPS)
        close(f)
        return psi0
    else
        println("wrong input of way_of_psi0")
        sqrt(-1)
    end
end
#
#
function rundmrg()
    #-----------------------------------------------------------------------
    # control parameter
    bc = "cylinder"
    Ly = 4
    Lx = 8
    N = Lx * Ly
    t = 1.0
    J = 1.0 
    dop_level = (0,8) # delta
    #
    Is_conserve_qns = true
    way_of_psi0_list = ["default", "readfile"] # {readfile or default}
    nepoch = 3
    #
    nsweeps_list = [5]
    maxdim_preset = [[100, 200, 400],[400,500,500]] # pre setting cutoff dimemsion
    write_disk_dim = 500
    cutoff = [1E-4]
    noise = [1E-6, 1E-7, 1E-8, 1E-9, 1E-10, 1E-11, 0.0]
    #
    workpath = "./tj/Ly-Lx-tj/output/psi/"
    dim_lable = 500
    #
    filename_psi_1 = "Ly$(Ly)_Lx$(Lx)_N$(N)_S$(0.5)_t$(t)_J$(round(J;digits=4))_dimlabel$(dim_lable)_dop$(dop_level)_"
    filename_psi_2 = bc
    filename_psi_3 = "_psi.h5"
    filename_psi = join([workpath, filename_psi_1, filename_psi_2, filename_psi_3])
    #------------------------------------------------------------------------------
    # get lattice bonds
    # 2D sqaure wrapped on a cylinder
    lattice = square_lattice(Lx,Ly, yperiodic = true)
    #
    os = OpSum()
    for b in lattice
        #@show (b.s1, b.s2)
        os += -t, "Cdagup", b.s1, "Cup", b.s2
        os += -t, "Cdagup", b.s2, "Cup", b.s1
        os += -t, "Cdagdn", b.s1, "Cdn", b.s2
        os += -t, "Cdagdn", b.s2, "Cdn", b.s1
        os += J, "Sz", b.s1, "Sz", b.s2
        os += 0.5 * J, "S+", b.s1, "S-", b.s2
        os += 0.5 * J, "S-", b.s1, "S+", b.s2 
        os += -0.25 * J, "Ntot", b.s1, "Ntot", b.s2  
    end
    #-------------------------------------------------------------------------------
    # begin dmrg
    for in = 1:nepoch
        println("**************begin $(in) sweeps epoch************************")
        if in == 1
            way_of_psi0 = way_of_psi0_list[1]
        else
            way_of_psi0 = way_of_psi0_list[2]
        end
        # get the hamiltonian MPO
        sites = build_sites(N, Is_conserve_qns)
        psi0 = get_psi0(N, dop_level,sites,Is_conserve_qns,way_of_psi0,filename_psi) # 
        sites = siteinds(psi0)    
        H = MPO(os, sites)
        #@show psi0
        for it = 1:size(nsweeps_list)[1]
            println("begin $(it) sweeps period")
            nsweeps = nsweeps_list[it]
            maxdim = maxdim_preset[it]
            #
            energy, psi = dmrg(H, psi0; nsweeps, maxdim, cutoff, noise,write_when_maxdim_exceeds=write_disk_dim)
            #
            println("-----------------------------------------------")
            @show energy
            @show flux(psi)
            @show maxlinkdim(psi)
            #
            @show (Ly, Lx, t, J, dop_level)
            per_energy = energy / N
            @show per_energy
            #Compute the energy variance of an MPS to check whether it is an eigenstate.
            H2 = inner(H,psi,H,psi)
            E = inner(psi',H,psi)
            var = H2-E^2
            @show var
            println("-----------------------------------------------")
            psi0 = psi            
            #
            #-----------------------------------------------------
            # save psi
            f = h5open(filename_psi, "w")
            write(f, "psi", psi)
            close(f)
        end
    end
    return
end
#
let 
    @time rundmrg()
end
#
