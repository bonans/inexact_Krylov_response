using FileIO
using Plots
using LaTeXStrings
using Statistics
using DataFrames
using LinearAlgebra
using Printf
#31.11396874956486
results_all = load("Fe2MnAl/data_logs/Fe2MnAl_extracted.jld2")

restart = 10
tol = 1e-9
save_to_dir = "Fe2MnAl/figures/"
if !isdir(save_to_dir)
    mkdir(save_to_dir)
end
save_to_file = ("_" * string(restart) * "_" * string(Int(log10(tol))))
lw_true = 4
lw_tilde = 3
op_tilde = 1.0
m_restart = :circle
m_sz_restart = 8
m_sz = 4
markers = [:utriangle,:utriangle,:utriangle,:dtriangle,:dtriangle,:dtriangle]
lagend_labels = [latexstring("\\texttt{agr} / \\texttt{Pagr}"), latexstring("\\texttt{bal} / \\texttt{Pbal}"), 
    latexstring("\\texttt{Pgrt}"), latexstring("\\texttt{D10} / \\texttt{PD10}"), 
    latexstring("\\texttt{D100} / \\texttt{PD100}"), latexstring("\\texttt{PD10n}")]
lagend_labels_noP = [latexstring("\\texttt{agr}"), latexstring("\\texttt{bal}"), 
    latexstring("\\texttt{grt}"), latexstring("\\texttt{D10}"), 
    latexstring("\\texttt{D100}"), latexstring("\\texttt{D10n}")]
lagend_labels_P = [latexstring("\\texttt{Pagr}"), latexstring("\\texttt{Pbal}"), 
latexstring("\\texttt{Pgrt}"), latexstring("\\texttt{PD10}"), 
latexstring("\\texttt{PD100}"), latexstring("\\texttt{PD10n}")]
colors = palette(cgrad(:roma), 6)#["#1f77b4", "#2ca02c", "#17becf", "#ff7f0e", "#d62728", "#ffbf00"]

colors = [palette(:tab10)[1], palette(:tab10)[3], palette(:tab10)[10], palette(:tab10)[4], palette(:tab10)[7], palette(:tab10)[2]]

CG_err_res_Al = load("Fe2MnAl/data_logs/CG_err_res_2_58.jld2")
CG_err_res_H = load("Fe2MnAl/data_logs/CG_err_res_732_26_H.jld2")

fig1 = plot()
plot!(fig1, CG_err_res_Al["residuals"],lw=lw_true, color=colors[1])
plot!(fig1, CG_err_res_Al["errors"],lw=lw_true, color=colors[4])
plot!(fig1, CG_err_res_Al["residuals_bound"],lw=lw_true, color=colors[2])
plot!(fig1, yscale=:log10, yticks=10.0 .^ (Vector(-16:2:0)), 
xtickfont=font(14, "times"), ytickfont=font(14, "times"),dpi=300, margin=0Plots.mm,legend=false)
annotate!(fig1, 23, 5e-2, Plots.text("(A)", "times", 16), :black)
fig2 = plot()
plot!(fig2, CG_err_res_H["residuals"], label=L"\| \textbf{r}_k \|",lw=lw_true, color=colors[1])
plot!(fig2, CG_err_res_H["errors"], label=L"\| \textbf{x}_k - \textbf{x}_\star \|",lw=lw_true, color=colors[4])
plot!(fig2, CG_err_res_H["residuals_bound"], label=L"\| \textbf{r}_k \|/(\epsilon_{N_{\mathrm{occ}+1}} - \epsilon_{N_{\mathrm{occ}}})",lw=lw_true, color=colors[2])
plot!(fig2, yscale=:log10, yticks=10.0 .^ (Vector(-16:2:0)),
xtickfont=font(14, "times"), ytickfont=font(14, "times"),dpi=300, 
margin=0Plots.mm,legend=:topright, legendfont=font(12, "times"))
annotate!(fig2, 23.5, 2e1, Plots.text("(B)", "times", 16), :black)
l = @layout [a b]
fig = plot(fig1, fig2, layout=l, size=(945, 405), dpi=300, 
left_margin=6Plots.mm, bottom_margin=7Plots.mm, top_margin=5Plots.mm)
savefig(fig, save_to_dir * "CG_err_res.svg")


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

########### plot convergence
max_xaxis = 0
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    filename = (load_from_dir * "results_P" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        global max_xaxis = max(max_xaxis, length(results["res"]))
    end
end


fig1 = plot()
hline!(fig1, [tol], label=latexstring("\\tau = 10^{" * string(Int(log10(tol))) * "}"), lw=lw_true, color=:black, ls=:dash)
for (i, legend_label) in enumerate(lagend_labels_P)
    if i <= 3
        plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i])
    else
        plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i], ls=:dash)
    end
end
for Pindex in 2:2
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    strategy_name = (Pindex == 2 ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive)
    filename = (load_from_dir * "results_" * strategy_name * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if adaptive != true && contains(adaptive, "_")
        adaptive = replace(adaptive, "_" => "\\_")
    end
    if haskey(results_all, filename)
        results = results_all[filename]
        first_nonzero = findfirst(x -> x > 0, results["res"])
        results["res"][1:first_nonzero-1] = results["res_tilde"][1:first_nonzero-1]
        if adaptive == true
            plot!(fig1, results["res"], 
            lw=lw_true, color=colors[i], label = "")#, label=lagend_label)# markershape=markers[i], markersize=m_sz)
        else
            plot!(fig1, results["res"], 
            lw=lw_true, color=colors[i], label = "", ls=:dash)
        end
    end
end
end
plot!(fig1, yscale=:log10, yticks=10.0 .^ (Vector(-14:2:ceil(log10(normδρ0)))), 
    xtickfont=font(14, "times"), ytickfont=font(14, "times"),
    size=(945, 405), dpi=300, margin=0Plots.mm, legend=:topright, legendfont=font(12, "times"),
    xlims=(1, max_xaxis+2),bottom_margin=1Plots.mm,ylims=(1e-10,1e3))
savefig(fig1, save_to_dir * "Hconvergence" * save_to_file * ".svg")


########### plot CG convergence
max_xaxis = 0
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    filename = (load_from_dir * "results_P" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        global max_xaxis = max(max_xaxis, results["CG_niters_all"])
    end
end


fig1 = plot()
hline!(fig1, [tol], label=latexstring("\\tau = 10^{" * string(Int(log10(tol))) * "}"), lw=lw_true, color=:black, ls=:dash)
for (i, legend_label) in enumerate(lagend_labels_P)
    if i <= 3
        plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i])
    else
        plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i], ls=:dash)
    end
end
for Pindex in 2:2
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    strategy_name = (Pindex == 2 ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive)
    filename = (load_from_dir * "results_" * strategy_name * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if adaptive != true && contains(adaptive, "_")
        adaptive = replace(adaptive, "_" => "\\_")
    end
    if haskey(results_all, filename)
        #lagend_label = latexstring("\\texttt{" * (Pindex == 2 ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive) * "}")
        results = results_all[filename]
        num_Mv = length(results["tol_sternheimer_all"])
        iter_inds = setdiff(1:num_Mv, results["restart_inds"] .+ (1:length(results["restart_inds"])) )
        CG_i_accu = results["CG_niters_accu"]
        first_nonzero = findfirst(x -> x > 0, results["res"])
        if adaptive == true
        plot!(fig1, CG_i_accu[iter_inds][first_nonzero:end], results["res"][first_nonzero:end], lw=lw_true, 
        color=colors[i], label="")
        else
        plot!(fig1, CG_i_accu[iter_inds][first_nonzero:end], results["res"][first_nonzero:end], lw=lw_true,
        color=colors[i], label="", ls=:dash)
        end
    end
end
end
plot!(fig1, yscale=:log10, yticks=10.0 .^ (Vector(-14:2:ceil(log10(normδρ0)))), 
    xticks=5e6:5e6:2e7, xaxis=(formatter = x -> string(round(Int, x/1e6)) * " M"),
    xtickfont=font(14, "times"), ytickfont=font(14, "times"),
    size=(945, 405), dpi=300, legend=false, legendfont=font(12, "times"),
    xlims=(1, 2.1e7),margin=0Plots.mm,bottom_margin=2Plots.mm,ylims=(1e-10,1e3))
savefig(fig1, save_to_dir * "HCG_convergence" * save_to_file * ".svg")



########### plot CG iterations and tolerance
max_xaxis = 0
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    filename = (load_from_dir * "results_P" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        global max_xaxis = max(max_xaxis, length(results["res"]))
    end
end

fig1 = plot()
axis2 = twinx()
for (i, legend_label) in enumerate(lagend_labels_P)
    plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i])
end
for Pindex in 2:2
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    if adaptive != true
        continue
    end
    lagend_label = ""#lagend_labels_P[i]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    strategy_name = (Pindex == 2 ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive)
    filename = (load_from_dir * "results_" * strategy_name * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        num_Mv = length(results["tol_sternheimer_all"])
        iter_inds = setdiff(1:num_Mv, results["restart_inds"] .+ (1:length(results["restart_inds"])) )
        CG_niters_i_mean = results["CG_niters_i"]
        plot!([iter_inds], CG_niters_i_mean[iter_inds], lw=lw_true, color=colors[i], label=lagend_label)
        if CG_tol_scale_choice == "hdmd"
            plot!(axis2, [iter_inds], max.(0.5*eps(), results["tol_sternheimer_all"][iter_inds]./31.11396874956486), lw=lw_true, color=colors[i], ls=:dash, label="")
        else
            plot!(axis2, [iter_inds], max.(0.5*eps(), results["tol_sternheimer_all"][iter_inds]), lw=lw_true, color=colors[i], ls=:dash, label="")
        end
    end
end
end
plot!(xtickfont=font(14, "times"), ytickfont=font(14, "times"),
    size=(945, 405), dpi=300,xlims=(0, max_xaxis+1), margin=0Plots.mm,ylims=(-1,31))
plot!(axis2,yscale=:log10,yticks=10.0 .^ (Vector( -16:2:0)),ylims=(1e-16,1e-1),yformatter=_->"")
plot!(legendfont=font(12, "times"),legend=:right)

fig2 = plot()
axis2 = twinx();
for Pindex in 2:2
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    if adaptive == true
        continue
    end
    lagend_label = ""#lagend_labels_P[i]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    strategy_name = (Pindex == 2 ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive)
    filename = (load_from_dir * "results_" * strategy_name * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        num_Mv = length(results["tol_sternheimer_all"])
        iter_inds = setdiff(1:num_Mv, results["restart_inds"] .+ (1:length(results["restart_inds"])) )
        CG_niters_i_mean = results["CG_niters_i"]
        plot!([iter_inds], CG_niters_i_mean[iter_inds], lw=lw_true, color=colors[i], label=lagend_label)
        if CG_tol_scale_choice == "hdmd"
            plot!(axis2, [iter_inds], max.(0.5*eps(), results["tol_sternheimer_all"][iter_inds]./31.11396874956486), lw=lw_true, color=colors[i], ls=:dash, label="")
        else
            plot!(axis2, [iter_inds], max.(0.5*eps(), results["tol_sternheimer_all"][iter_inds]), lw=lw_true, color=colors[i], ls=:dash, label="")
        end
    end
end
end
plot!(xtickfont=font(14, "times"), ytickfont=font(1, "times"),ytickfontcolor=:white,
    size=(945, 405), dpi=300,xlims=(0, max_xaxis+4), margin=0Plots.mm,ylims=(-1,31))
plot!(axis2,yscale=:log10,yticks=10.0 .^ (Vector( -16:2:0)),ylims=(1e-16,1e-1),ytickfont=font(12, "times"))
plot!(legendfont=font(12, "times"),legend=:right)

fig = plot(fig1, fig2, layout=grid(1,2, widths=(2/3,1/3)), size=(945, 405), dpi=300, 
left_margin=0Plots.mm, bottom_margin=1Plots.mm, top_margin=0Plots.mm)
savefig(fig, save_to_dir * "HCG_niters" * save_to_file * ".svg")