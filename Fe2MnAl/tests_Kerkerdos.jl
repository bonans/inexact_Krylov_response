using FileIO
using Plots
using LaTeXStrings

include("Fe2MnAlPDos.jl")
debug = false # "debug" or "full"

setup_model((debug ? "debug" : "full"))
run_gmres(; debug=debug, restart=10, tol=1e-9, adaptive=true, CG_tol_scale_choice="hdmd", precon=true)
run_gmres(; debug=debug, restart=10, tol=1e-9, adaptive="D100", CG_tol_scale_choice="agr", precon=true)

# Do the plots and summary table

lw_true = 4
legend_labels_P = [latexstring("\\texttt{Pagr}"), latexstring("\\texttt{Pbal}"), 
latexstring("\\texttt{Pgrt}"), latexstring("\\texttt{PD10}"), 
latexstring("\\texttt{PD100}"), latexstring("\\texttt{PD10n}")]
colors = [palette(:tab10)[1], palette(:tab10)[3], palette(:tab10)[10], palette(:tab10)[4], palette(:tab10)[7], palette(:tab10)[2]]

save_to_dir = "Fe2MnAl/figures/"
if !isdir(save_to_dir)
    mkdir(save_to_dir)
end
ρ, ψ, ham, basis, occupation, εF, eigenvalues, δρ0 = load_model(debug=debug)
normδρ0 = norm(δρ0)
results_all = load("Fe2MnAl/data_logs/" * (debug ? "00debug_" : "") * "results_Phdmd_10_-9_PDos.jld2")
first_nonzero = findfirst(x -> x > 0, results_all["res"])
results_all["res"][1:first_nonzero-1] = results_all["res_tilde"][1:first_nonzero-1]
res_bal = copy(results_all["res"])
CG_iter_bal = sum([sum(x) for x in sum(results_all["CG_niters_all"])])
efficiency_bal = (log10(normδρ0) - log10(res_bal[end])) / CG_iter_bal

results_all = load("Fe2MnAl/data_logs/" * (debug ? "00debug_" : "") * "results_PD100_10_-9_PDos.jld2")
first_nonzero = findfirst(x -> x > 0, results_all["res"])
results_all["res"][1:first_nonzero-1] = results_all["res_tilde"][1:first_nonzero-1]
res_D100 = copy(results_all["res"])
CG_iter_D100 = sum([sum(x) for x in sum(results_all["CG_niters_all"])])
efficiency_D100 = (log10(normδρ0) - log10(res_D100[end])) / CG_iter_D100

max_xaxis = max(length(res_bal), length(res_D100))
fig1 = plot()
hline!(fig1, [1e-9], label=latexstring("\\tau = 10^{" * string(Int(log10(1e-9))) * "}"), lw=lw_true, color=:black, ls=:dash)
plot!(fig1, x=1:1, 1e-9 * ones(1), label=legend_labels_P[2], lw=lw_true, color=colors[2])
plot!(fig1, x=1:1, 1e-9 * ones(1), label=legend_labels_P[5], lw=lw_true, color=colors[5], ls=:dash)
plot!(fig1, res_bal,lw=lw_true, color=colors[2], label = "")
plot!(fig1, res_D100,lw=lw_true, color=colors[5], ls=:dash, label = "")

plot!(fig1, yscale=:log10, yticks=10.0 .^ (Vector(-14:2:ceil(log10(normδρ0)))), 
    xtickfont=font(14, "times"), ytickfont=font(14, "times"),
    size=(945, 405), dpi=300, margin=0Plots.mm, legend=:topright, legendfont=font(12, "times"),
    xlims=(1, max_xaxis+2),bottom_margin=1Plots.mm,ylims=(1e-10,1e3))
savefig(fig1, save_to_dir * "Hconvergence_10_-9_PDos.svg")

# Summary table: res_end, CG_iter, Rel efficiency (D100 baseline = 1.0)
rel_eff_D100 = 1.0
rel_eff_bal = efficiency_bal / efficiency_D100

println("Summary:")
header = @sprintf("%15s %12s %12s","", "Bal", "D100")
println(header)
println(""*repeat('-', length(header)))
println(@sprintf("%15s %12.3e %12.3e", "res_end", res_bal[end], res_D100[end]))
println(@sprintf("%15s %12d %12d", "CG_iter", Int(CG_iter_bal), Int(CG_iter_D100)))
println(@sprintf("%15s %12.3f %12.3f", "Rel efficiency", rel_eff_bal, rel_eff_D100))
