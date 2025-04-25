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
save_to_dir = "silicon/data_logs/"

function setup_model(;a = 10.0, Ecut=40, kgrid=[8, 8, 8], tol=1e-12, debug=false)
    model_specs = (;
        a=a,
        Ecut=Ecut,
        kgrid=kgrid)

    # set up model
    println("------ Setting up model ... ------")
    lattice = a / 2 * [[0 1 1.];
                        [1 0 1.];
                        [1 1 0.]]
    Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))
    atoms     = [Si, Si]
    positions = [ones(3)/8, -ones(3)/8]
    
    model = model_PBE(lattice, atoms, positions, symmetries=false)
    basis = PlaneWaveBasis(model; Ecut, kgrid)
    
    println(show(stdout, MIME("text/plain"), basis))
    println("------ Running SCF ... ------")
    DFTK.reset_timer!(DFTK.timer)
    scfres = self_consistent_field(basis; tol=tol)
    println(DFTK.timer)

    println("------ Computing RHS ... ------")
    
    R = [zeros(3), ones(3)]
    function V(ε)
        T = typeof(ε)
        pos = positions + ε*R
        modelV = Model(Matrix{T}(lattice), atoms, pos; model_name="potential",
                        terms=[DFTK.AtomicLocal(),
                                DFTK.AtomicNonlocal()], symmetries=false)
        basisV = PlaneWaveBasis(modelV; Ecut, kgrid)
        jambon = Hamiltonian(basisV)
        DFTK.total_local_potential(jambon)
    end
    δV = ForwardDiff.derivative(V, 0.0)
    println("||δV|| = ", norm(δV))
    flush(stdout)
    DFTK.reset_timer!(DFTK.timer)
    δρ0 = apply_χ0(scfres, δV; tol=1e-16)
    println(DFTK.timer)
    normδρ0 = norm(δρ0)
    println("||δρ0|| = ", normδρ0)

    println("------ Saving model ... ------")
    save_to_file = save_to_dir * (debug ? "00debug_" : "") * "scf.jld2"
    ψ = scfres.ψ
    ρ = scfres.ρ
    save(save_to_file, "model_specs", model_specs,
        "ρ", ρ, "ψ", ψ, "δρ0", δρ0)

    println("------ Computing CG_tol_scale used in experiments ... ------")
    scfres = self_consistent_field(basis; tol=tol, ρ=ρ, ψ=ψ,maxiter=2)
    occupation = scfres.occupation
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

function load_model(; tol=1e-12, debug=false)
    load_from_file = save_to_dir * (debug ? "00debug_" : "") * "scf.jld2"
    model_specs, ρ, ψ, δρ0 = load(load_from_file, "model_specs", "ρ", "ψ", "δρ0")

    lattice = model_specs.a / 2 * [[0 1 1.];
                                    [1 0 1.];
                                    [1 1 0.]]
    Si = ElementPsp(:Si, psp=load_psp("hgh/lda/Si-q4"))
    atoms     = [Si, Si]
    positions = [ones(3)/8, -ones(3)/8]
    
    model_new = model_PBE(lattice, atoms, positions, symmetries=false)
    basis = PlaneWaveBasis(model_new; model_specs.Ecut, model_specs.kgrid)
    scfres = self_consistent_field(basis; tol=tol, ρ=ρ, ψ=ψ,maxiter=2)

    return (; ρ=scfres.ρ, ψ=scfres.ψ, ham=scfres.ham, basis=scfres.ham.basis,
        occupation=scfres.occupation, εF=scfres.εF, eigenvalues=scfres.eigenvalues, δρ0=δρ0)
end


function setup_model(debug::String)
    if lowercase(debug) == "debug"
        println("Setting up model in debug mode")
        a = 10.0
        Ecut = 5
        kgrid = (8, 8, 8)
        open(save_to_dir * "00debug_log_scf.log", "w") do io
            redirect_stdout(io) do
                setup_model(a=a, Ecut=Ecut, kgrid=kgrid, debug=true)
            end
        end
    elseif lowercase(debug) == "full"
        println("Setting up full model")
        # otherwise use the full model
        open(save_to_dir * "log_scf.log", "w") do io
            redirect_stdout(io) do
                setup_model()
            end
        end
    else
        error("Invalid debug choice")
    end
end



function run_gmres(; debug=false, restart=20, tol=1e-12, adaptive=true, CG_tol_scale_choice="agr", precon=false)
    save_to_dir_new = save_to_dir * (debug ? "00debug_" : "")
    save_to_file = ((precon == true ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                * string(restart) * "_" * string(Int(log10(tol))))

    ρ, ψ, ham, basis, occupation, εF, eigenvalues, δρ0 = load_model(debug=debug)
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
            
            println("----- running exact applications of χ0 of the last iteration -----")
            flush(stdout)
            DFTK.reset_timer!(DFTK.timer)
            res = norm(eps_fun_exact(results_a.x[:, end]) - pack(δρ0))
            println("i = ", size(results_a.x, 2), ", true res = ", res, "\n")
            flush(stdout)
            jldopen(save_to_dir_new * "results_" * save_to_file * ".jld2", "r+") do file
                if haskey(file, "res")
                    delete!(file, "res")
                end
                write(file, "res", res)
            end
            println(DFTK.timer)
        end
    end
end