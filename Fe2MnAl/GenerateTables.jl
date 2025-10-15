using FileIO
using Plots
using LaTeXStrings
using Statistics
using DataFrames
using LinearAlgebra
using Printf
using ASEconvert
using PrettyTables
results = load("Fe2MnAl/data_logs/Fe2MnAl_extracted.jld2")

load_from_dir = "Fe2MnAl/data_logs/"
restart_list = [5, 10, 20]
tol_list = [1e-6, 1e-9, 1e-12]
normδρ0 = 4709.452671852897
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

for restart in restart_list
    for tol in tol_list
        current_row = ["m=" * string(restart) * ", τ=10^-" * string(-Int(log10(tol))); fill(NaN, length(columns))]
        push!(df_PRE, current_row)
        push!(df_RE, current_row)
        push!(df_E, current_row)
        push!(df_CG_iters, current_row)
        push!(df_res, current_row)
        filename = (load_from_dir * "results_" * "D10" * "_"
                    * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
        filenameP = (load_from_dir * "results_P" * "D10" * "_"
                     * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
        filenameifP = haskey(results, filenameP) ? filenameP : filename
        base_efficiencyP = (log10(normδρ0) - log10(results[filenameifP]["res"][end])) / results[filenameifP]["CG_niters_all"]
        base_efficiency = (log10(normδρ0) - log10(results[filename]["res"][end])) / results[filename]["CG_niters_all"]
        for strategy in strategies
            adaptive = strategy["adaptive"]
            CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
            filename = (load_from_dir * "results_" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                        * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
            filenameP = (load_from_dir * "results_P" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                         * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
            if haskey(results, filename)
                res = results[filename]["res"]
                efficiency = (log10(normδρ0) - log10(res[end])) / results[filename]["CG_niters_all"]
                df_PRE[end, Symbol((adaptive == true ? CG_tol_scale_choice : adaptive))] = efficiency / base_efficiencyP
                df_RE[end, Symbol((adaptive == true ? CG_tol_scale_choice : adaptive))] = efficiency / base_efficiency
                df_E[end, Symbol((adaptive == true ? CG_tol_scale_choice : adaptive))] = efficiency * 10000
                df_CG_iters[end, Symbol((adaptive == true ? CG_tol_scale_choice : adaptive))] = results[filename]["CG_niters_all"]
                df_res[end, Symbol((adaptive == true ? CG_tol_scale_choice : adaptive))] = res[end]
            end
            if haskey(results, filenameP)
                res = results[filenameP]["res"]
                efficiency = (log10(normδρ0) - log10(res[end])) / results[filenameP]["CG_niters_all"]
                df_PRE[end, Symbol("P" * (adaptive == true ? CG_tol_scale_choice : adaptive))] = efficiency / base_efficiencyP
                df_RE[end, Symbol("P" * (adaptive == true ? CG_tol_scale_choice : adaptive))] = efficiency / base_efficiency
                df_E[end, Symbol("P" * (adaptive == true ? CG_tol_scale_choice : adaptive))] = efficiency * 10000
                df_CG_iters[end, Symbol("P" * (adaptive == true ? CG_tol_scale_choice : adaptive))] = results[filenameP]["CG_niters_all"]
                df_res[end, Symbol("P" * (adaptive == true ? CG_tol_scale_choice : adaptive))] = res[end]
            end
        end
    end
end

ft_nonan(v::Float64, i::Int, j::Int) = isnan(v) ? "" : v
ft_nonan(v::String, i::Int, j::Int) = v == "NaN" ? "" : v

output_file = "Fe2MnAl/data_logs/summary_tables.md"
open(output_file, "w") do io
    write(io, "# Summary tables for Heusler systems\n")
    write(io, "See Section 4.3 & Supplementary Material SM7 of the paper for detailed setups for computations.\n")
    write(io, "## Relative efficiency referencing `D10`\n")
    pretty_table(io, df_RE; backend = Val(:markdown), formatters = (ft_nonan,ft_printf("%1.2f")),)
    write(io, "\n")
    write(io, "## Number of Hamiltonian applications (CG iterations)\n")
    pretty_table(io, df_CG_iters; backend = Val(:markdown), formatters = (ft_nonan,ft_printf("%8.0f")))
    write(io, "\n")
    write(io, "## Accuracy of the final solution\n")
    pretty_table(io, df_res; backend = Val(:markdown), formatters = (ft_nonan,ft_printf("%1.2e")))
    write(io, "\n")
end