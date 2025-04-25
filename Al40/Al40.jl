using ASEconvert
using DFTK
using JLD2
using LinearMaps
using LinearAlgebra
using Printf
using Dates
using Random
using ForwardDiff

disable_threading()
save_to_dir = "Al40/data_logs/"

function setup_model(; repeat=10, rattle_intensity=0.05, Ecut=40, kgrid=(1, 3, 3), temperature=0.001, tol=1e-14, debug=false)
    model_specs = (; repeat=repeat,
        rattle_intensity=rattle_intensity,
        Ecut=Ecut,
        kgrid=kgrid,
        temperature=temperature)

    println("------ Setting up model ... ------")
    mixing = KerkerMixing()
    system = ase.build.bulk("Al", cubic=true) * pytuple((repeat, 1, 1))
    system = pyconvert(AbstractSystem, system)
    system = attach_psp(system; Al="hgh/pbe/Al-q3")
    model = model_PBE(system, temperature=temperature, symmetries=false)
    basis = PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)
    println(show(stdout, MIME("text/plain"), basis))
    println("------ Running SCF ... ------")
    DFTK.reset_timer!(DFTK.timer)
    scfres = self_consistent_field(basis; tol=tol, mixing=mixing)
    println(DFTK.timer)

    println("------ Computing rhs ... ------")

    ρ, ψ, ham, basis, occupation, εF, eigenvalues = scfres.ρ, scfres.ψ, scfres.ham, scfres.basis, scfres.occupation, scfres.εF, scfres.eigenvalues
    positions = model.positions
    lattice = model.lattice
    atoms = model.atoms
    R = [zeros(3) for pos in positions]
    Random.seed!(1234)
    for iR in 1:length(R)
        R[iR] = -ones(3) + 2 * rand(3)
    end
    function V(ε)
        T = typeof(ε)
        pos = positions + ε * R
        modelV = Model(Matrix{T}(lattice), atoms, pos; model_name="potential",
            terms=[DFTK.AtomicLocal(), DFTK.AtomicNonlocal()], symmetries=false)
        basisV = PlaneWaveBasis(modelV; Ecut, kgrid)
        jambon = Hamiltonian(basisV)
        DFTK.total_local_potential(jambon)
    end
    δV = ForwardDiff.derivative(V, 0.0)
    println("||δV|| = ", norm(δV))
    flush(stdout)
    DFTK.reset_timer!(DFTK.timer)
    δρ0 = apply_χ0(ham, ψ, occupation, εF, eigenvalues, δV; tol=1e-16)
    println(DFTK.timer)
    println("||δρ0|| = ", norm(δρ0))

    println("------ Saving model ... ------")
    save_to_file = save_to_dir * "repeat" * string(repeat) * (debug ? "/00debug_" : "/") * "scf.jld2"
    save(save_to_file, "model_specs", model_specs,
        "ρ", ρ, "ψ", ψ, "δρ0", δρ0)

    println("------ Computing CG_tol_scale used in experiments ... ------")
    scfres = self_consistent_field(basis; tol=1e-12, ρ=ρ, ψ=ψ, mixing=mixing, maxiter=2)
    num_kpoints = length(basis.kpoints)
    apply_χ0_info = DFTK.get_apply_χ0_info(scfres.ham, ψ, occupation, scfres.εF, scfres.eigenvalues)
    CG_tol_scale = apply_χ0_info.CG_tol_scale
    Nocc_ks = [length(CG_tol_scale[ik]) for ik in 1:num_kpoints]
    Nocc = sum(Nocc_ks)
    fn_occ = [occupation[ik][maskk] for (ik, maskk) in enumerate(apply_χ0_info.mask_occ)]
    CG_tol_scale = [fn_occ[ik] * basis.kweights[ik] for ik in 1:num_kpoints] * Nocc * sqrt(prod(basis.fft_size)) / basis.model.unit_cell_volume
    
    kcoef = zeros(num_kpoints)
    for k in 1:num_kpoints
    accum = zeros(basis.fft_size)
    for n in 1:Nocc_ks[k]
        accum += (abs2.(real.(ifft(basis, basis.kpoints[k], ψ[k][:, n]))))
    end
    kcoef[k] = sqrt(maximum(accum)) * basis.kweights[k]
    end

    CG_tol_scale_grt = [fn_occ[ik] * kcoef[ik] for ik in 1:num_kpoints] * sqrt(Nocc) * sqrt(prod(basis.fft_size)) / sqrt(basis.model.unit_cell_volume)

    jldopen(save_to_file, "r+") do file
        if haskey(file, "CG_tol_scale")
            delete!(file, "CG_tol_scale")
        end
        if haskey(file, "CG_tol_scale_grt")
            delete!(file, "CG_tol_scale_grt")
        end
        write(file, "CG_tol_scale", CG_tol_scale)
        write(file, "CG_tol_scale_grt", CG_tol_scale_grt)
    end

    println("------ Model setup done ------")
end

function load_model(; tol=1e-12, debug=false, repeat=10)
    load_from_file = save_to_dir * "repeat" * string(repeat) * (debug ? "/00debug_" : "/") * "scf.jld2"
    model_specs, ρ, ψ, δρ0 = load(load_from_file, "model_specs", "ρ", "ψ", "δρ0")
    mixing = KerkerMixing()
    system = ase.build.bulk("Al", cubic=true) * pytuple((model_specs.repeat, 1, 1))
    system = pyconvert(AbstractSystem, system)
    system = attach_psp(system; Al="hgh/pbe/Al-q3")
    model = model_PBE(system, temperature=model_specs.temperature, symmetries=false)
    basis = PlaneWaveBasis(model; Ecut=model_specs.Ecut, kgrid=model_specs.kgrid)
    scfres = self_consistent_field(basis; tol=tol, ρ=ρ, ψ=ψ, mixing=mixing, maxiter=2)

    return (; ρ=scfres.ρ, ψ=scfres.ψ, ham=scfres.ham, basis=scfres.ham.basis,
        occupation=scfres.occupation, εF=scfres.εF, eigenvalues=scfres.eigenvalues, δρ0=δρ0)
end

function setup_model(debug::String; repeat=10)
    if lowercase(debug) == "debug"
        println("Setting up model in debug mode")
        # use a small model for debugging
        rattle_intensity = 0.3
        Ecut = 10
        kgrid = (repeat == 1 ? (3, 3, 3) : (1, 3, 3))
        # check if save_to_dir * "repeat" * string(repeat) * "/" exists and create if not
        if !isdir(save_to_dir * "repeat" * string(repeat) * "/")
            mkdir(save_to_dir * "repeat" * string(repeat) * "/")
        end
        open(save_to_dir * "repeat" * string(repeat) * "/00debug_log_scf.log", "w") do io
            redirect_stdout(io) do
                setup_model(;repeat=repeat, rattle_intensity=rattle_intensity, Ecut=Ecut, kgrid=kgrid, debug=true)
            end
        end
    elseif lowercase(debug) == "full"
        println("Setting up full model with repeat = ", repeat)
        # otherwise use the full model
        open(save_to_dir * "repeat" * string(repeat) * "/log_scf.log", "w") do io
            redirect_stdout(io) do
                if repeat == 1
                    setup_model(;repeat=repeat,kgrid=(3, 3, 3))
                else
                    setup_model(;repeat=repeat)
                end
            end
        end
    else
        error("Invalid model setup choice")
    end
end

function run_gmres(; debug=false, restart=20, tol=1e-12, adaptive=true, CG_tol_scale_choice="agr", repeat=10, precon=false)
    save_to_dir_new = save_to_dir * "repeat" * string(repeat) * (debug ? "/00debug_" : "/")
    save_to_file = ((precon == true ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                * string(restart) * "_" * string(Int(log10(tol))))

    ρ, ψ, ham, basis, occupation, εF, eigenvalues, δρ0 = load_model(debug=debug, repeat=repeat)
    mixing = KerkerMixing()
    # Define some helper functions and variables
    pack(δρ) = vec(δρ)
    unpack(δρ) = reshape(δρ, size(ρ))

    num_kpoints = length(basis.kpoints)
    apply_χ0_info = DFTK.get_apply_χ0_info(ham, ψ, occupation, εF, eigenvalues)
    CG_tol_scale = apply_χ0_info.CG_tol_scale
    Nocc_ks = [length(CG_tol_scale[ik]) for ik in 1:num_kpoints]
    Nocc = sum(Nocc_ks)

    if adaptive == true
        if CG_tol_scale_choice == "agr"
            CG_tol_scale = [[1.0 for _ in 1:Nocc_ks[ik]] for ik in 1:num_kpoints]
        elseif CG_tol_scale_choice == "hdmd" 
            CG_tol_scale = load(save_to_dir_new * "scf.jld2", "CG_tol_scale")
        elseif CG_tol_scale_choice == "grt"
            CG_tol_scale = load(save_to_dir_new * "scf.jld2", "CG_tol_scale_grt")
        else
            error("Invalid CG tolerance scaling choice")
        end
    elseif adaptive == "D10" || adaptive == "D100" || adaptive == "D10_n"
        CG_tol_scale = [[1.0 for _ in 1:Nocc_ks[ik]] for ik in 1:num_kpoints]
    else
        error("Invalid adaptive choice")
    end

    # k_to_k_minus_q, mask_occ, ψ_occ, ψ_extra, ε_occ, ε_minus_q_occ, CG_tol_scale)
    apply_χ0_info = (; k_to_k_minus_q=apply_χ0_info.k_to_k_minus_q,
        mask_occ=apply_χ0_info.mask_occ,
        ψ_occ=apply_χ0_info.ψ_occ,
        ψ_extra=apply_χ0_info.ψ_extra,
        ε_occ=apply_χ0_info.ε_occ,
        ε_minus_q_occ=apply_χ0_info.ε_minus_q_occ,
        CG_tol_scale=CG_tol_scale)

    normδV_all = Float64[]
    tol_sternheimer_all = Float64[]
    CG_niters_all = Vector{Vector{Int64}}[]
    CG_niters_all_exact = Vector{Vector{Int64}}[]

    inds = [0, 1, 1] # i, ik, n
    inds_exact = [0, 1, 1] # i, ik, n
    function sternheimer_callback(CG_niters)
        function callback(info)
            if inds[3] > Nocc_ks[inds[2]]
                inds[3] = 1
                inds[2] += 1
            end
            CG_niters[inds[2]][inds[3]] = info.res.n_iter
            #println("| ", inds[1], " | ", tol_sternheimer, " | ", inds[2], " | ", inds[3], " | ", tol_sternheimer / CG_tol_scale[inds[2]][inds[3]], " | ", info.res.n_iter, " |") Δtime
            #@printf("| %-7d | %-10.3e | %-7d | %-11d | %-10.3e | %-8d | %-7.3f |\n", inds[1], tol_sternheimer, inds[2], inds[3], tol_sternheimer / CG_tol_scale[inds[2]][inds[3]], info.res.n_iter, (runtime_ns[2] - runtime_ns[1]) / 1e9)
            #flush(stdout)
            inds[3] += 1
        end
    end

    function sternheimer_callback_exact(CG_niters)
        function callback(info)
            if inds_exact[3] > Nocc_ks[inds_exact[2]]
                inds_exact[3] = 1
                inds_exact[2] += 1
            end
            CG_niters[inds_exact[2]][inds_exact[3]] = info.res.n_iter
            #println("| ", inds[1], " | ", tol_sternheimer, " | ", inds[2], " | ", inds[3], " | ", tol_sternheimer / CG_tol_scale[inds[2]][inds[3]], " | ", info.res.n_iter, " |")
            #@printf("| %-7d | %-7d | %-11d | %-8d | %-7.3f |\n", inds_exact[1], inds_exact[2], inds_exact[3], info.res.n_iter, (runtime_ns[2] - runtime_ns[1]) / 1e9)
            #flush(stdout)
            inds_exact[3] += 1
        end
    end

    function operators_a(tol_sternheimer)
        function eps_fun(δρ)
            δρ = unpack(δρ)
            δV = apply_kernel(basis, δρ; ρ)

            #println("| GMRES i | τ_i | k-point | Sternheimer | τ_CG | CG iters |")
            #@printf("| %-7s | %-10s | %-7s | %-11s | %-10s | %-8s | %-7s |\n", "GMRES i", "τ_i", "k-point", "Sternheimer", "τ_CG", "CG iters", "Δtime")
            inds[1] += 1
            inds[2:3] = [1, 1]
            push!(normδV_all, norm(DFTK.symmetrize_ρ(basis, δV)))
            if adaptive == true && CG_tol_scale_choice == "grt"
                tol_sternheimer = tol_sternheimer ./ (2*normδV_all[end])
            end
            push!(tol_sternheimer_all, tol_sternheimer)
            CG_niters = [zeros(Int64, Nocc_ks[i]) for i in 1:num_kpoints]

            println("---- τ_CG's used for each Sternheimer equation (row) of each k-point (column) ----")
            τ_CG_table = [max.(0.5*eps(Float64), tol_sternheimer ./ CG_tol_scale[ik]) for ik in 1:num_kpoints]
            @printf("| %-7s ", "k-point")
            for n in 1:maximum(Nocc_ks)
                @printf("| %-8d ", n)
            end
            @printf("|\n")
            for (k, row) in enumerate(τ_CG_table)
                @printf("| %-7d ", k)
                for τ in row[1:end]
                    @printf("| %-8.2e ", τ)
                end
                @printf("|\n")
            end
            @printf("| %-10s | %-10s | %-10s | %-10s |\n", "τ_i", "min τ_CG", "mean τ_CG", "max τ_CG")
            @printf("| %-10.3e | %-10.3e | %-10.3e | %-10.3e |\n\n", tol_sternheimer, minimum(reduce(vcat, τ_CG_table)), exp(sum(log.(reduce(vcat, τ_CG_table)))/Nocc), maximum(reduce(vcat, τ_CG_table)))
            flush(stdout)

            t1 = Dates.now()
            χ0δV = apply_χ0(ham, ψ, occupation, εF, eigenvalues, δV; tol=tol_sternheimer, callback=sternheimer_callback(CG_niters), apply_χ0_info=apply_χ0_info)
            t2 = Dates.now()

            push!(CG_niters_all, CG_niters)
            println("no. CG iters for each Sternheimer equation (row) of each k-point (column):")
            @printf("| %-7s ", "k-point")
            for n in 1:maximum(Nocc_ks)
                @printf("| %-3d ", n)
            end
            @printf("|\n")
            for (k, row) in enumerate(CG_niters)
                @printf("| %-7d ", k)
                for niters in row[1:end]
                    @printf("| %-3d ", niters)
                end
                @printf("|\n")
            end
            @printf("| %-10s | %-10s | %-10s | %-10s | %-10s |\n", "min", "mean", "max", "sum", "total")
            @printf("| %-10d | %-10.3f | %-10d | %-10d | %-10d |\n\n", minimum(reduce(vcat, CG_niters)), sum(reduce(vcat, CG_niters)) / Nocc, maximum(reduce(vcat, CG_niters)), sum(reduce(vcat, CG_niters)), sum(sum.(sum(CG_niters_all))))
            println("χ0Time = ", canonicalize(t2 - t1), ", time now: ", Dates.format(t2, "yyyy-mm-dd HH:MM:SS"))
            flush(stdout)
            #println("ave CG iters = ", sum(reduce(vcat, CG_niters)) / Nocc)

            if precon
                return pack(DFTK.mix_density(mixing, basis, δρ - χ0δV))
            else
                return pack(δρ - χ0δV)
            end
        end
        return LinearMap(eps_fun, prod(size(δρ0)))
    end

    function eps_fun_exact(δρ)
        normδρ = norm(δρ)
        δρ = δρ ./ normδρ
        δρ = unpack(δρ)
        δV = apply_kernel(basis, δρ; ρ)

        #println("| GMRES i | τ_i | k-point | Sternheimer | τ_CG | CG iters |")
        #@printf("| %-7s | %-7s | %-11s | %-8s | %-7s |\n", "GMRES i", "k-point", "Sternheimer", "CG iters", "Δtime")
        inds_exact[1] += 1
        inds_exact[2:3] = [1, 1]
        CG_niters = [zeros(Int64, Nocc_ks[i]) for i in 1:num_kpoints]

        t1 = Dates.now()
        χ0δV = apply_χ0(ham, ψ, occupation, εF, eigenvalues, δV; tol=1e-16, callback=sternheimer_callback_exact(CG_niters))
        t2 = Dates.now()

        push!(CG_niters_all_exact, CG_niters)

        println("no. CG iters:")
        @printf("| %-10s | %-10s | %-10s | %-10s | %-10s |\n", "min", "mean", "max", "sum", "total")
        @printf("| %-10d | %-10.3f | %-10d | %-10d | %-10d |\n\n", minimum(reduce(vcat, CG_niters)), sum(reduce(vcat, CG_niters)) / Nocc, maximum(reduce(vcat, CG_niters)), sum(reduce(vcat, CG_niters)), sum(sum.(sum(CG_niters_all_exact))))
        println("χ0Time = ", canonicalize(t2 - t1), ", time now: ", Dates.format(t2, "yyyy-mm-dd HH:MM:SS"))
        pack(δρ - χ0δV) .* normδρ
    end

    open(save_to_dir_new * "log_" * save_to_file * ".log", "w") do io
        redirect_stdout(io) do
            println("----- running GMRES: tol=", tol, ", restart=", restart, ", adaptive=", adaptive, ", CG_tol_scale_choice=", CG_tol_scale_choice, " -----")
            Pδρ0 = δρ0
            if precon
                Pδρ0 = DFTK.mix_density(mixing, basis, δρ0)
            end
            println("||Pδρ0|| = ", norm(Pδρ0))
            DFTK.reset_timer!(DFTK.timer)
            if adaptive == "D10"
                results_a = DFTK.gmres(operators_a(tol / 10), pack(Pδρ0); restart=restart, tol=tol, verbose=1, debug=true)
            elseif adaptive == "D10_n"
                results_a = DFTK.gmres(operators_a(tol / 10 / norm(Pδρ0)), pack(Pδρ0); restart=restart, tol=tol, verbose=1, debug=true)
            elseif adaptive == "D100"
                results_a = DFTK.gmres(operators_a(tol / 100), pack(Pδρ0); restart=restart, tol=tol, verbose=1, debug=true)
            elseif adaptive == true
                results_a = DFTK.gmres(operators_a, pack(Pδρ0); restart=restart, tol=tol, verbose=1, debug=true)
            else
                error("Invalid adaptive choice")
            end
            println(DFTK.timer)

            if haskey(results_a, :minsvdvals)
                save(save_to_dir_new * "results_" * save_to_file * ".jld2", "X", results_a.x, "res_tilde", results_a.residuals,
                    "MvTime", results_a.MvTime, "restart_inds", results_a.restart_inds, "minsvdvals", results_a.minsvdvals,
                    "normδV_all", normδV_all, "tol_sternheimer_all", tol_sternheimer_all, "CG_niters_all", CG_niters_all)
            else
                save(save_to_dir_new * "results_" * save_to_file * ".jld2", "X", results_a.x, "res_tilde", results_a.residuals,
                    "MvTime", results_a.MvTime, "restart_inds", results_a.restart_inds,
                    "normδV_all", normδV_all, "tol_sternheimer_all", tol_sternheimer_all, "CG_niters_all", CG_niters_all)
            end
            if length(results_a.restart_inds) > 0
                computeToIndex = max(1,results_a.restart_inds[end] - max(ceil(Int,1.3*restart),10))
            else
                computeToIndex = 1
            end
            println("----- running exact applications of χ0 from i = ", size(results_a.x, 2), " to ", computeToIndex, " -----")
            flush(stdout)
            res = zeros(Float64, size(results_a.x, 2))
            DFTK.reset_timer!(DFTK.timer)
            for i in size(results_a.x, 2):-1:computeToIndex
                res[i] = norm(eps_fun_exact(results_a.x[:, i]) - pack(δρ0))
                println("i = ", i, ", true res = ", res[i], "\n")
                flush(stdout)
                jldopen(save_to_dir_new * "results_" * save_to_file * ".jld2", "r+") do file
                    if haskey(file, "res")
                        delete!(file, "res")
                    end
                    write(file, "res", res)
                end
            end
            println(DFTK.timer)
        end
    end
end

function run_test_repeat(; repeat=10, rattle_intensity=0.05)
    open(save_to_dir * "log_" * string(repeat) * ".log", "w") do io
        redirect_stdout(io) do
            Ecut = 40
            kgrid = (1, 3, 3)
            temperature = 0.001
            tol_scf = 1e-12
            println("------ Setting up model ... ------")
            mixing = KerkerMixing()
            system = ase.build.bulk("Al", cubic=true) * pytuple((repeat, 1, 1))
            system = pyconvert(AbstractSystem, system)
            system = attach_psp(system; Al="hgh/pbe/Al-q3")
            model = model_PBE(system, temperature=temperature, symmetries=false)
            basis = PlaneWaveBasis(model; Ecut=Ecut, kgrid=kgrid)
            println(show(stdout, MIME("text/plain"), basis))
            println("------ Running SCF ... ------")
            DFTK.reset_timer!(DFTK.timer)
            flush(stdout)
            scfres = self_consistent_field(basis; tol=tol_scf, mixing=mixing)
            println(DFTK.timer)

            println("------ Setting up rattle model ... ------")
            system_r = ase.build.bulk("Al", cubic=true) * pytuple((repeat, 1, 1))
            system_r.rattle(stdev=rattle_intensity, seed=42)
            system_r = pyconvert(AbstractSystem, system_r)
            system_r = attach_psp(system_r; Al="hgh/pbe/Al-q3")
            model_r = model_PBE(system_r, temperature=temperature, symmetries=false)
            basis_r = PlaneWaveBasis(model_r; Ecut=Ecut, kgrid=kgrid)
            println("------ Running SCF for rattle model ... ------")
            DFTK.reset_timer!(DFTK.timer)
            scfres_r = self_consistent_field(basis_r; tol=tol_scf, mixing=mixing)
            println(DFTK.timer)

            ρ = scfres.ρ
            ψ = scfres.ψ
            ham = scfres.ham
            basis = ham.basis
            occupation = scfres.occupation
            εF = scfres.εF
            eigenvalues = scfres.eigenvalues
            δρ0 = scfres_r.ρ - ρ

            println("------ Computing CG_tol_scale used in experiments ... ------")
            num_kpoints = length(basis.kpoints)
            apply_χ0_info = DFTK.get_apply_χ0_info(scfres.ham, ψ, occupation, scfres.εF, scfres.eigenvalues)
            CG_tol_scale = apply_χ0_info.CG_tol_scale
            Nocc_ks = [length(CG_tol_scale[ik]) for ik in 1:num_kpoints]
            Nocc = sum(Nocc_ks)
            fn_occ = [occupation[ik][maskk] for (ik, maskk) in enumerate(apply_χ0_info.mask_occ)]
            CG_tol_scale = [fn_occ[ik] * basis.kweights[ik] for ik in 1:num_kpoints] * Nocc * sqrt(prod(basis.fft_size)) / basis.model.unit_cell_volume

            # gmres
            restart = 20
            tol = 1e-12
            # Define some helper functions and variables
            pack(δρ) = vec(δρ)
            unpack(δρ) = reshape(δρ, size(ρ))

            num_kpoints = length(basis.kpoints)
            apply_χ0_info = DFTK.get_apply_χ0_info(ham, ψ, occupation, εF, eigenvalues)

            apply_χ0_info = (; k_to_k_minus_q=apply_χ0_info.k_to_k_minus_q,
                mask_occ=apply_χ0_info.mask_occ,
                ψ_occ=apply_χ0_info.ψ_occ,
                ψ_extra=apply_χ0_info.ψ_extra,
                ε_occ=apply_χ0_info.ε_occ,
                ε_minus_q_occ=apply_χ0_info.ε_minus_q_occ,
                CG_tol_scale=CG_tol_scale)

            normδV_all = Float64[]
            tol_sternheimer_all = Float64[]
            CG_niters_all = Vector{Vector{Int64}}[]
            CG_niters_all_exact = Vector{Vector{Int64}}[]

            inds = [0, 1, 1] # i, ik, n
            inds_exact = [0, 1, 1] # i, ik, n
            runtime_ns = UInt64[time_ns(), time_ns()]
            function sternheimer_callback(CG_niters, tol_sternheimer)
                function callback(info)
                    if inds[3] > Nocc_ks[inds[2]]
                        inds[3] = 1
                        inds[2] += 1
                    end
                    CG_niters[inds[2]][inds[3]] = info.res.n_iter
                    runtime_ns[1] = runtime_ns[2]
                    runtime_ns[2] = time_ns()
                    #println("| ", inds[1], " | ", tol_sternheimer, " | ", inds[2], " | ", inds[3], " | ", tol_sternheimer / CG_tol_scale[inds[2]][inds[3]], " | ", info.res.n_iter, " |") Δtime
                    @printf("| %-7d | %-10.3e | %-7d | %-11d | %-10.3e | %-8d | %-7.3f |\n", inds[1], tol_sternheimer, inds[2], inds[3], tol_sternheimer / CG_tol_scale[inds[2]][inds[3]], info.res.n_iter, (runtime_ns[2] - runtime_ns[1]) / 1e9)
                    flush(stdout)
                    inds[3] += 1
                end
            end

            function sternheimer_callback_exact(CG_niters)
                function callback(info)
                    if inds_exact[3] > Nocc_ks[inds_exact[2]]
                        inds_exact[3] = 1
                        inds_exact[2] += 1
                    end
                    CG_niters[inds_exact[2]][inds_exact[3]] = info.res.n_iter
                    runtime_ns[1] = runtime_ns[2]
                    runtime_ns[2] = time_ns()
                    #println("| ", inds[1], " | ", tol_sternheimer, " | ", inds[2], " | ", inds[3], " | ", tol_sternheimer / CG_tol_scale[inds[2]][inds[3]], " | ", info.res.n_iter, " |")
                    @printf("| %-7d | %-7d | %-11d | %-8d | %-7.3f |\n", inds_exact[1], inds_exact[2], inds_exact[3], info.res.n_iter, (runtime_ns[2] - runtime_ns[1]) / 1e9)
                    flush(stdout)
                    inds_exact[3] += 1
                end
            end

            function operators_a(tol_sternheimer)
                function eps_fun(δρ)
                    δρ = unpack(δρ)
                    δV = apply_kernel(basis, δρ; ρ)

                    #println("| GMRES i | τ_i | k-point | Sternheimer | τ_CG | CG iters |")
                    @printf("| %-7s | %-10s | %-7s | %-11s | %-10s | %-8s | %-7s |\n", "GMRES i", "τ_i", "k-point", "Sternheimer", "τ_CG", "CG iters", "Δtime")
                    inds[1] += 1
                    inds[2:3] = [1, 1]
                    push!(normδV_all, norm(DFTK.symmetrize_ρ(basis, δV)))
                    push!(tol_sternheimer_all, tol_sternheimer)
                    CG_niters = [zeros(Int64, Nocc_ks[i]) for i in 1:num_kpoints]

                    χ0δV = apply_χ0(ham, ψ, occupation, εF, eigenvalues, δV; tol=tol_sternheimer, callback=sternheimer_callback(CG_niters, tol_sternheimer), apply_χ0_info=apply_χ0_info)

                    push!(CG_niters_all, CG_niters)
                    println("ave CG iters = ", sum(reduce(vcat, CG_niters)) / Nocc)

                    pack(δρ - χ0δV)
                end
                return LinearMap(eps_fun, prod(size(δρ0)))
            end

            function eps_fun_exact(δρ)
                δρ = unpack(δρ)
                δV = apply_kernel(basis, δρ; ρ)

                #println("| GMRES i | τ_i | k-point | Sternheimer | τ_CG | CG iters |")
                @printf("| %-7s | %-7s | %-11s | %-8s | %-7s |\n", "GMRES i", "k-point", "Sternheimer", "CG iters", "Δtime")
                inds_exact[1] += 1
                inds_exact[2:3] = [1, 1]
                CG_niters = [zeros(Int64, Nocc_ks[i]) for i in 1:num_kpoints]

                χ0δV = apply_χ0(ham, ψ, occupation, εF, eigenvalues, δV; tol=1e-16, callback=sternheimer_callback_exact(CG_niters))

                push!(CG_niters_all_exact, CG_niters)
                println("ave CG iters = ", sum(reduce(vcat, CG_niters)) / Nocc)

                pack(δρ - χ0δV)
            end

            println("----- running GMRES: tol=", tol, ", restart=", restart, " -----")
            DFTK.reset_timer!(DFTK.timer)
            results_a = DFTK.gmres(operators_a, pack(δρ0); restart=restart, tol=tol, verbose=1, debug=true)
            println(DFTK.timer)
            save(save_to_dir * "results_" * string(repeat) * ".jld2", "X", results_a.x, "res_tilde", results_a.residuals,
                "MvTime", results_a.MvTime, "restart_inds", results_a.restart_inds, "minsvdvals", results_a.minsvdvals,
                "normδV_all", normδV_all, "tol_sternheimer_all", tol_sternheimer_all, "CG_niters_all", CG_niters_all)
            println("----- running exact applications of χ0 -----")
            res = zeros(Float64, size(results_a.x, 2))
            DFTK.reset_timer!(DFTK.timer)
            for i in size(results_a.x, 2):-1:min(1, size(results_a.x, 2) - 10)
                res[i] = norm(eps_fun_exact(results_a.x[:, i]) - pack(δρ0))
                jldopen(save_to_dir * "results_" * string(repeat) * ".jld2", "r+") do file
                    if haskey(file, "res")
                        delete!(file, "res")
                    end
                    write(file, "res", res)
                end
            end
            println(DFTK.timer)
        end
    end
end

if false
    begin
        restart_list = [5,10,20]
        tol_list = [1e-3,1e-6,1e-9,1e-12]

        restart = restart_list[1]
        tol = tol_list[1]

        repeat_list = [1,3,5,8,10]
        configs = [Dict("adaptive" => true, "CG_tol_scale_choice" => "agr"),
            Dict("adaptive" => true, "CG_tol_scale_choice" => "hdmd"),
            Dict("adaptive" => "D10", "CG_tol_scale_choice" => nothing),
            Dict("adaptive" => "D100", "CG_tol_scale_choice" => nothing)]

        for repeat in repeat_list
            for (i, config) in enumerate(configs)
                adaptive = config["adaptive"]
                CG_tol_scale_choice = config["CG_tol_scale_choice"]
                println("-----currently running restart=", restart, ", tol=", tol, ", repeat=", repeat, ", adaptive=", adaptive, ", CG_tol_scale_choice=", CG_tol_scale_choice, "-----")
                flush(stdout)
                run_gmres(; debug=false, restart=restart, tol=tol, adaptive=adaptive, CG_tol_scale_choice=CG_tol_scale_choice, repeat=repeat)
            end
        end
    end
end

if false
    begin
        restart_list = [5,10,20]

        restart = restart_list[1]
        tol = 1e-12

        repeat_list = [1,3,5,8]
        configs = [Dict("adaptive" => true, "CG_tol_scale_choice" => "agr"),
            Dict("adaptive" => true, "CG_tol_scale_choice" => "hdmd"),
            Dict("adaptive" => "D10", "CG_tol_scale_choice" => nothing),
            Dict("adaptive" => "D100", "CG_tol_scale_choice" => nothing)]

        for repeat in repeat_list
            for (i, config) in enumerate(configs)
                adaptive = config["adaptive"]
                CG_tol_scale_choice = config["CG_tol_scale_choice"]
                println("-----currently running restart=", restart, ", tol=", tol, ", repeat=", repeat, ", adaptive=", adaptive, ", CG_tol_scale_choice=", CG_tol_scale_choice, "-----")
                flush(stdout)
                run_gmres(; debug=false, restart=restart, tol=tol, adaptive=adaptive, CG_tol_scale_choice=CG_tol_scale_choice, repeat=repeat)
            end
        end
    end
end
if false
    begin
        restart_list = [5,10,20]
        configs = [Dict("adaptive" => true, "CG_tol_scale_choice" => "agr"),
            Dict("adaptive" => true, "CG_tol_scale_choice" => "hdmd"),
            Dict("adaptive" => "D10", "CG_tol_scale_choice" => nothing),
            Dict("adaptive" => "D100", "CG_tol_scale_choice" => nothing)]

        restart = restart_list[1]
        tol = 1e-12
        repeat = 10
        config = configs[1]
        adaptive = config["adaptive"]
        CG_tol_scale_choice = config["CG_tol_scale_choice"]
        
        println("-----currently running restart=", restart, ", tol=", tol, ", repeat=", repeat, ", adaptive=", adaptive, ", CG_tol_scale_choice=", CG_tol_scale_choice, "-----")
        flush(stdout)
        run_gmres(; debug=false, restart=restart, tol=tol, adaptive=adaptive, CG_tol_scale_choice=CG_tol_scale_choice, repeat=repeat)
    end
end