using FileIO
using Plots
using LaTeXStrings
using Statistics
using DataFrames
using LinearAlgebra
using Printf
using ASEconvert
using PrettyTables
results = load("Al40/data_logs/Al_extracted.jld2")

load_from_dir = "metals/Al40/data_logs/"
repeat_list = [1, 3, 5, 10]
restart_list = [5, 10, 20]
tol_list = [1e-6, 1e-9, 1e-12]

strategies = [Dict("adaptive" => true, "CG_tol_scale_choice" => "agr"),
    Dict("adaptive" => true, "CG_tol_scale_choice" => "hdmd"),
    Dict("adaptive" => true, "CG_tol_scale_choice" => "grt"),
    Dict("adaptive" => "D10", "CG_tol_scale_choice" => nothing),
    Dict("adaptive" => "D100", "CG_tol_scale_choice" => nothing),
    Dict("adaptive" => "D10_n", "CG_tol_scale_choice" => nothing)]

#columns = [:agr, :hdmd, :D10, :D100, :Pagr, :Phdmd, :Pgrt, :PD10, :PD100, :PD10_n]
columns = [:Pagr, :Phdmd, :Pgrt, :PD10, :PD100, :PD10_n, :agr, :hdmd, :D10, :D100]
df_PRE = DataFrame(; vcat([:config => String[]], [col => Float64[] for col in columns])...)
df_RE = DataFrame(; vcat([:config => String[]], [col => Float64[] for col in columns])...)
df_E = DataFrame(; vcat([:config => String[]], [col => Float64[] for col in columns])...)
df_CG_iters = DataFrame(; vcat([:config => String[]], [col => Float64[] for col in columns])...)
df_res = DataFrame(; vcat([:config => String[]], [col => Float64[] for col in columns])...)

normδρ0_list = Dict(1 => 97.66346537186564, 3 => 539.1145933154517, 5 => 1307.9139596204495, 10 => 6617.248515007278)
for repeat in repeat_list
    normδρ0 = normδρ0_list[repeat]
    for restart in restart_list
        for tol in tol_list
            current_row = ["ℓ=" * string(repeat) * ", m=" * string(restart) * ", τ=10^-" * string(-Int(log10(tol))); fill(NaN, length(columns))]
            push!(df_PRE, current_row)
            push!(df_RE, current_row)
            push!(df_E, current_row)
            push!(df_CG_iters, current_row)
            push!(df_res, current_row)
            filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * "D10" * "_"
                        * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
            filenameP = (load_from_dir * "repeat" * string(repeat) * "/results_P" * "D10" * "_"
                            * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
            filenameifP = haskey(results, filenameP) ? filenameP : filename
            base_efficiencyP = (log10(normδρ0) - log10(results[filenameifP]["res"][end])) / results[filenameifP]["CG_niters_all"]
            base_efficiency = (log10(normδρ0) - log10(results[filename]["res"][end])) / results[filename]["CG_niters_all"]
            for strategy in strategies
                adaptive = strategy["adaptive"]
                CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
                filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                            * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
                filenameP = (load_from_dir * "repeat" * string(repeat) * "/results_P" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
                if haskey(results, filename)
                    res = results[filename]["res"]
                    efficiency = (log10(normδρ0) - log10(res[end])) / results[filename]["CG_niters_all"]
                    col = Symbol((adaptive == true ? CG_tol_scale_choice : adaptive))
                    df_PRE[end, col] = efficiency / base_efficiencyP
                    df_RE[end, col] = efficiency / base_efficiency
                    df_E[end, col] = efficiency * 10000
                    df_CG_iters[end, col] = results[filename]["CG_niters_all"]
                    df_res[end, col] = res[end]
                end
                if haskey(results, filenameP)
                    res = results[filenameP]["res"]
                    efficiency = (log10(normδρ0) - log10(res[end])) / results[filenameP]["CG_niters_all"]
                    col = Symbol("P" * (adaptive == true ? CG_tol_scale_choice : adaptive))
                    df_PRE[end, col] = efficiency / base_efficiencyP
                    df_RE[end, col] = efficiency / base_efficiency
                    df_E[end, col] = efficiency * 10000
                    df_CG_iters[end, col] = results[filenameP]["CG_niters_all"]
                    df_res[end, col] = res[end]
                end
            end
        end
    end
end

columns_speedup = [Symbol("ℓ=" * string(repeat)) for repeat in repeat_list]
df_speedup = DataFrame(; vcat([:config => String[]], [col => Float64[] for col in columns_speedup])... )
for restart in restart_list[2:end]
    for tol in tol_list
        for strategy in [strategies[1:2];strategies[4:5]]
            adaptive = strategy["adaptive"]
            CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
            push!(df_speedup, ["m=" * string(restart) * ", τ=10^-" * string(-Int(log10(tol))) * ", " * (adaptive == true ? CG_tol_scale_choice : adaptive);
                    fill(NaN, length(columns_speedup))])
            for (i,repeat) in enumerate(repeat_list)
                i_row = findfirst(df_E.config .== "ℓ=" * string(repeat) * ", m=" * string(restart) * ", τ=10^-" * string(-Int(log10(tol))))
                col = Symbol((adaptive == true ? CG_tol_scale_choice : adaptive))
                colP = Symbol("P" * (adaptive == true ? CG_tol_scale_choice : adaptive))
                df_speedup[end, Symbol("ℓ=" * string(repeat))] = df_E[i_row, colP] / df_E[i_row, col]
            end
        end
    end
end

ft_nonan(v::Float64, i::Int, j::Int) = isnan(v) ? "" : v
ft_nonan(v::String, i::Int, j::Int) = v == "NaN" ? "" : v

output_file = "Al40/data_logs/summary_tables.md"
open(output_file, "w") do io
    write(io, "# Summary tables for aluminium systems\n")
    write(io, "See Section 4.2 & Supplementary Material SM7 of the paper for detailed setups for computations.\n")
    write(io, "## Relative efficiency referencing `D10`\n")
    pretty_table(io, df_RE; backend = Val(:markdown), formatters = (ft_nonan,ft_printf("%1.2f")),)
    write(io, "\n")
    write(io, "## Number of Hamiltonian applications (CG iterations)\n")
    pretty_table(io, df_CG_iters; backend = Val(:markdown), formatters = (ft_nonan,ft_printf("%7.0f")))
    write(io, "\n")
    write(io, "## Accuracy of the final solution\n")
    pretty_table(io, df_res; backend = Val(:markdown), formatters = (ft_nonan,ft_printf("%1.2e")))
    write(io, "\n")
end