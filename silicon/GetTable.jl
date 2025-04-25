using JLD2
using LinearAlgebra
strategies = [Dict("adaptive" => true, "CG_tol_scale_choice" => "agr"),
    Dict("adaptive" => true, "CG_tol_scale_choice" => "hdmd"),
    Dict("adaptive" => true, "CG_tol_scale_choice" => "grt"),
    Dict("adaptive" => "D10", "CG_tol_scale_choice" => nothing),
    Dict("adaptive" => "D100", "CG_tol_scale_choice" => nothing),
    Dict("adaptive" => "D10_n", "CG_tol_scale_choice" => nothing)]

debug = true
result_table = zeros(length(strategies), 4)
println("res_end CG_iter efficiency")
for (i, config) in enumerate(strategies)
    adaptive = config["adaptive"]
    CG_tol_scale_choice = config["CG_tol_scale_choice"]
    load_from_file = ("silicon/data_logs/" * (debug ? "00debug_" : "") * "results_" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                      * string(20) * "_" * string(Int(log10(1e-9))))
    results = load(load_from_file * ".jld2")

    # 77.47952930755548 or 327.37936457757695
    normδρ0 = norm(load("silicon/data_logs/"* (debug ? "00debug_" : "") *"scf.jld2")["δρ0"])
    
    debug ? 77.47952930755548 : 327.37936457757695
    res_end = results["res"][end]
    CG_niters_all = results["CG_niters_all"]
    gmres_iters = length(CG_niters_all)
    CG_iter_all = sum(sum([sum.(CG_niters_all[i]) for i in 1:gmres_iters]))
    efficiency = (log10(normδρ0) - log10(res_end)) / CG_iter_all * 10000
    # print res_end, CG_iter_all, efficiency
    println(res_end, " ", CG_iter_all, " ", efficiency)
    result_table[i, 1] = res_end
    result_table[i, 2] = CG_iter_all
    result_table[i, 3] = efficiency
end
result_table[:, 4] = result_table[:, 3] ./ result_table[4, 3]