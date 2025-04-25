include("Al40.jl")
using FileIO
using Plots
using LaTeXStrings
using Statistics

repeat = 3
debug = true
tol = 1e-12
restart = 20
precon = true

setup_model("debug"; repeat=repeat)

run_gmres(; debug=debug, restart=restart, tol=tol, adaptive=true, CG_tol_scale_choice="hdmd", repeat=repeat, precon=precon)
run_gmres(; debug=debug, restart=restart, tol=tol, adaptive="D10_n", CG_tol_scale_choice=nothing, repeat=repeat, precon=precon)

function get_CG_niters(results, mode::String)
    gmres_iters = length(results["CG_niters_all"])
    if mode == "all"
        return sum(sum([sum.(results["CG_niters_all"][i]) for i in 1:gmres_iters]))
    end
    if mode == "i_mean"
        return [mean(reduce(vcat, results["CG_niters_all"][i])) for i in 1:gmres_iters]
    end
end

load_from_dir = "Al40/data_logs/" * "repeat" * string(repeat) * (debug ? "/00debug_" : "/")
save_to_dir = "Al40/figures/" * "repeat" * string(repeat) * (debug ? "/00debug_" : "/")
save_to_file = ("_" * string(restart) * "_" * string(Int(log10(tol))))
if !isdir("Al40/figures/")
    mkdir("Al40/figures/")
end
if !isdir("Al40/figures/repeat" * string(repeat) * "/")
    mkdir("Al40/figures/repeat" * string(repeat) * "/")
end
configs = [Dict("adaptive" => true, "CG_tol_scale_choice" => "Phdmd"),
    Dict("adaptive" => "PD10_n", "CG_tol_scale_choice" => nothing)]

lw_true = 2
lw_tilde = 1.5
op_tilde = 0.5
m_restart = :circle
m_sz_restart = 4

max_xaxis = 0
for (i, config) in enumerate(configs)
    adaptive = config["adaptive"]
    CG_tol_scale_choice = config["CG_tol_scale_choice"]
    load_from_file = (load_from_dir * "results_" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                      * string(restart) * "_" * string(Int(log10(tol))))
    results = load(load_from_file * ".jld2")
    global max_xaxis = max(max_xaxis, length(results["res"]))
end

β₀ = 1.0
fig = plot()
tau_label = "τ = 1e-" * string(Int(log10(tol)))
plot!(fig, x=1:max_xaxis, tol * ones(max_xaxis), label=tau_label, lw=lw_true, color=:black, ls=:dash)
for (i, config) in enumerate(configs)
    adaptive = config["adaptive"]
    CG_tol_scale_choice = config["CG_tol_scale_choice"]
    load_from_file = (load_from_dir * "results_" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                      * string(restart) * "_" * string(Int(log10(tol))))
    results = load(load_from_file * ".jld2")
    global β₀ = max(results["res_tilde"][1], β₀)
    legend_label = (adaptive == true ? CG_tol_scale_choice : adaptive) * ":"
    scatter!([1], [1e-2], label=legend_label, ms=0, mc=:white, msc=:white)
    scatter!([1], [1e-2], label="N=" * string(get_CG_niters(results, "all")), ms=0, mc=:white, msc=:white)
    # find the first non-zero for results["res"]
    first_nonzero = findfirst(x -> x > 0, results["res"])
    plot!(fig, first_nonzero:length(results["res"]),results["res"][first_nonzero:end], label="true residual",
        lw=lw_true, color=i)
    plot!(fig, results["res_tilde"], label="estimated residual",
        lw=lw_tilde, color=i, ls=:dash, opacity=op_tilde)
    # marker for restarts
    scatter!(fig, results["restart_inds"], results["res_tilde"][results["restart_inds"]],
        label="", m=m_restart, color=i, markersize=m_sz_restart)
end
plot!(fig, yscale=:log10, xlims=(0, floor(1.05 * max_xaxis)),
    yticks=10.0.^(Vector(-14:2:ceil(log10(β₀)))),
    xtickfont=font(14, "times"), ytickfont=font(14, "times"))
plot!(fig, size=(945, 405), dpi=300, margin=5Plots.mm)
xlabel!(fig, "iteration", xguidefont=font(16, "times"))
ylabel!(fig, "residual norm ", yguidefont=font(16, "times"))
plot!(fig, legendfont=font(12, "times"), legend=:outertopright)
# save to svg
savefig(fig, save_to_dir * "convergence" * save_to_file * ".svg")
