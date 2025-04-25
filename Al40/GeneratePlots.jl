using FileIO
using Plots
using LaTeXStrings
using Statistics
using DataFrames
using LinearAlgebra
using Printf
using ASEconvert

results_all = load("Al40/data_logs/Al_extracted.jld2")

load_from_dir = "metals/Al40/data_logs/"
repeat_list = [1, 3, 5, 10]

strategies = [Dict("adaptive" => true, "CG_tol_scale_choice" => "agr"),
    Dict("adaptive" => true, "CG_tol_scale_choice" => "hdmd"),
    Dict("adaptive" => true, "CG_tol_scale_choice" => "grt"),
    Dict("adaptive" => "D10", "CG_tol_scale_choice" => nothing),
    Dict("adaptive" => "D100", "CG_tol_scale_choice" => nothing),
    Dict("adaptive" => "D10_n", "CG_tol_scale_choice" => nothing)]

repeat = 10
restart = 20
tol = 1e-9
save_to_dir = "Al40/figures/" * "repeat" * string(repeat) * "/"
if !isdir("Al40/figures/")
    mkdir("Al40/figures/")
end
if !isdir(save_to_dir)
    mkdir(save_to_dir)
end
save_to_file = ("_" * string(restart) * "_" * string(Int(log10(tol))))
lw_true = 4
lw_tilde = 3
op_tilde = 1.0
m_restart = :circle
m_sz_restart = 8
m_sz = 8
markers = [:utriangle,:utriangle,:utriangle,:dtriangle,:circle,:dtriangle]
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

fig1 = plot()
for (i,repeat) in enumerate(repeat_list)
    normδV_all = results_all[load_from_dir * "repeat" * string(repeat) * "/results_hdmd_20_-9.jld2"]["normδV_all"]
    scatter!(fig1, repeat * ones(length(normδV_all)), normδV_all, m=:circle, ms=m_sz,color=(i==1 ? colors[5] : colors[i]))
end
plot!(fig1, 1:0.1:10.1, 7.65339^2/π .* (1:0.1:10.1).^2, lw=lw_true, color=:black, ls=:dash)
plot!(fig1,legend=false, yticks=0:250:2000, xticks=repeat_list,
    xtickfont=font(14, "times"), ytickfont=font(14, "times"),margin=0Plots.mm)
fig2 = plot()
for (i,repeat) in enumerate(repeat_list)
    normδV_all = results_all[load_from_dir * "repeat" * string(repeat) * "/results_Phdmd_20_-9.jld2"]["normδV_all"]
    scatter!(fig2, repeat * ones(length(normδV_all)), normδV_all, m=:circle, ms=m_sz,color=(i==1 ? colors[5] : colors[i]) )
end
plot!(fig2, 1:0.1:10.1, 7.65339^2/π .* (1:0.1:10.1).^2, lw=lw_true, color=:black, ls=:dash,
    label="bound of the Hartree kernel")
plot!(fig2,legend=false, yticks=0:250:2000, xticks=repeat_list,
    xtickfont=font(14, "times"), ytickfont=font(14, "times"),margin=0Plots.mm)
l = @layout [a b]
fig = plot(fig1, fig2, layout=l, size=(945, 405), dpi=300, left_margin=5Plots.mm, bottom_margin=6Plots.mm)
savefig(fig, save_to_dir * "normdeltaV.svg")

# plot supercell
system = ase.build.bulk("Al", cubic=true) * pytuple((repeat, 1, 1))
ase.io.write(save_to_dir * "Al.eps",system, rotation="10x,-20y,0z")

########### plot convergence v.s. i for D10 and D100 with est res
normδρ0_list = Dict(1 => 97.66346537186564, 3 => 539.1145933154517, 5 => 1307.9139596204495, 10 => 6617.248515007278)
max_xaxis = 0
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        global max_xaxis = max(max_xaxis, length(results["res"]))
    end
end

fig1 = plot()
hline!(fig1, [tol], label=latexstring("\\tau = 10^{" * string(Int(log10(tol))) * "}"), lw=lw_true, color=:black, ls=:dash)
jj = 1

plot!(fig1, x=1:1, 1e-20 * ones(1), label=latexstring("\\texttt{D10}\\ \\Vert \\mathbf{r}_i \\Vert"), 
    lw=lw_true, color=1, markershape=markers[4], markersize=6)
plot!(fig1, x=1:1, 1e-20 * ones(1), label=latexstring("\\texttt{D10}\\ \\Vert \\tilde{\\mathbf{r}}_i \\Vert"),
    lw=4, color=1, ls=:dash, opacity=op_tilde)
plot!(fig1, x=1:1, 1e-20 * ones(1), label=latexstring("\\texttt{D100}\\ \\Vert \\mathbf{r}_i \\Vert"), 
    lw=lw_true, color=2, markershape=markers[5], markersize=6)
plot!(fig1, x=1:1, 1e-20 * ones(1), label=latexstring("\\texttt{D100}\\ \\Vert \\tilde{\\mathbf{r}}_i \\Vert"),
    lw=4, color=2, ls=:dash, opacity=op_tilde)

jj = 1
for Pindex in 1:1
    for (i, strategy) in enumerate(strategies)
        # if i < 4 then continue end
        if i <= 3 || i >= 6
            continue
        end
        adaptive = strategy["adaptive"]
        CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
        strategy_name = (Pindex == 2 ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive)
        filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * strategy_name * "_"
                    * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
        if adaptive != true && contains(adaptive, "_")
            adaptive = replace(adaptive, "_" => "\\_")
        end
        if haskey(results_all, filename)
            results = results_all[filename]
            first_nonzero = findfirst(x -> x > 0, results["res"])
            results["res"][1:first_nonzero-1] = results["res_tilde"][1:first_nonzero-1]
            plot!(fig1, results["res"], 
            lw=lw_true, color=jj, label = "", markershape=markers[i], markersize=6)
            plot!(fig1, results["res_tilde"], lw=4, color=jj, ls=:dash, opacity=op_tilde, label="")
            jj += 1
        end
    end
end
plot!(fig1, yscale=:log10, yticks=10.0 .^ (Vector(-14:2:ceil(log10(normδρ0_list[repeat])))), 
    xtickfont=font(14, "times"), ytickfont=font(14, "times"),
    size=(945, 405), dpi=300, margin=0Plots.mm,bottom_margin=1Plots.mm,
    legend=:topright, legendfont=font(12, "times"),
    xlims=(1, max_xaxis+1),ylims=(1e-10,0.99*1e4))
savefig(fig1, save_to_dir * "convergence_naive" * save_to_file * ".svg")


########### plot convergence v.s. i for all (no est res)
max_xaxis = 0
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        global max_xaxis = max(max_xaxis, length(results["res"]))
    end
end
fig1 = plot()
hline!(fig1, [tol], label=latexstring("\\tau = 10^{" * string(Int(log10(tol))) * "}"), lw=lw_true, color=:black, ls=:dash)
for (i, legend_label) in enumerate(lagend_labels)
    if i <= 3
        plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i])
    else
        plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i], ls=:dash)
    end
end
for Pindex in 1:2
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    strategy_name = (Pindex == 2 ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive)
    filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * strategy_name * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
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
plot!(fig1, yscale=:log10, yticks=10.0 .^ (Vector(-14:2:ceil(log10(normδρ0_list[repeat])))), 
    xtickfont=font(14, "times"), ytickfont=font(14, "times"),
    size=(945, 405), dpi=300, margin=0Plots.mm, bottom_margin=1Plots.mm,
    legend=:topright, legendfont=font(12, "times"),
    xlims=(1, max_xaxis+1),ylims=(1e-10,0.99*1e4))

x = collect(range(0, 25, length= 100))
y1 = 5e-10 * ones(100)
y2 = 1e4*exp.(-1.05.*x)
plot!(fig1,x,y1, fillrange = y2, fillalpha = 0.35, c = :grey, label="w/ Kerker")
annotate!(fig1,8,1e-3, Plots.text("w/ Kerker", 16, rotation = -50, font = "times"))
annotate!(fig1,25,1e-1, Plots.text("w/o Kerker", 16, rotation = -25, font = "times"))
savefig(fig1, save_to_dir * "convergence" * save_to_file * ".svg")


########### plot CG iter vs residual
max_xaxis = 0
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        global max_xaxis = max(max_xaxis, results["CG_niters_all"])
    end
end
fig1 = plot()
hline!(fig1, [tol], label=latexstring("\\tau = 10^{" * string(Int(log10(tol))) * "}"), lw=lw_true, color=:black, ls=:dash)
for (i, legend_label) in enumerate(lagend_labels)
    if i <= 3
        plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i])
    else
        plot!(fig1, x=1:1, tol * ones(1), label=legend_label, lw=lw_true, color=colors[i], ls=:dash)
    end
end
for Pindex in 1:2
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    strategy_name = (Pindex == 2 ? "P" : "") * (adaptive == true ? CG_tol_scale_choice : adaptive)
    filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * strategy_name * "_"
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
        #,markershape=markers[i], markersize=m_sz)
        # scatter!(fig, CG_i_accu[results["restart_inds"]], results["res_tilde"][results["restart_inds"]], 
        # m=m_restart, color=(i <= 3) ? i : i+2, markersize=m_sz_restart, label="")
    end
end
end
plot!(fig1, yscale=:log10, yticks=10.0 .^ (Vector(-14:2:ceil(log10(normδρ0_list[repeat])))), 
    xtickfont=font(14, "times"), ytickfont=font(14, "times"),
    size=(945, 405), dpi=300, legend=false, legendfont=font(12, "times"),
    xlims=(1, max_xaxis+1),ylims=(1e-10,0.99*1e4),margin=0Plots.mm,bottom_margin=1Plots.mm)
plot!(fig1, xticks=1e5:1e5:5e5,xaxis=(formatter = x -> string(round(Int, x/1e3)) * " k"))
x = collect(range(0, 3e5, length= 100))
y1 = 5e-10 * ones(100)
y2 = 2e4*exp.(-0.0001 .* x.^0.99)
plot!(fig1,x,y1, fillrange = y2, fillalpha = 0.35, c = :grey, label="w/ Kerker")
annotate!(fig1,7.5e4,1e-3, Plots.text("w/ Kerker", 16, rotation = -50, font = "times"))
annotate!(fig1,3e5,1e-1, Plots.text("w/o Kerker", 16, rotation = -25, font = "times"))
savefig(fig1, save_to_dir * "CG_convergence" *  save_to_file * ".svg")


########### plot CG iterations and tolerance
max_xaxis = 0
for (i, strategy) in enumerate(strategies)
    adaptive = strategy["adaptive"]
    CG_tol_scale_choice = strategy["CG_tol_scale_choice"]
    filename = (load_from_dir * "repeat" * string(repeat) * "/results_P" * (adaptive == true ? CG_tol_scale_choice : adaptive) * "_"
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
    filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * strategy_name * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        num_Mv = length(results["tol_sternheimer_all"])
        iter_inds = setdiff(1:num_Mv, results["restart_inds"] .+ (1:length(results["restart_inds"])) )
        CG_niters_i_mean = results["CG_niters_i"]
        plot!([iter_inds], CG_niters_i_mean[iter_inds], lw=lw_true, color=colors[i], label=lagend_label)
        if CG_tol_scale_choice == "hdmd"
            plot!(axis2, [iter_inds], max.(0.5*eps(), results["tol_sternheimer_all"][iter_inds]./17.663871318548775), lw=lw_true, color=colors[i], ls=:dash, label="")
        else
            plot!(axis2, [iter_inds], max.(0.5*eps(), results["tol_sternheimer_all"][iter_inds]), lw=lw_true, color=colors[i], ls=:dash, label="")
        end
        # scatter!(axis4, results["restart_inds"], max.(0.5*eps(), results["tol_sternheimer_all"][results["restart_inds"].+(1:length(results["restart_inds"]))]),
        # label="", m=m_restart, color=i, markersize=m_sz_restart)
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
    filename = (load_from_dir * "repeat" * string(repeat) * "/results_" * strategy_name * "_"
                * string(restart) * "_" * string(Int(log10(tol))) * ".jld2")
    if haskey(results_all, filename)
        results = results_all[filename]
        num_Mv = length(results["tol_sternheimer_all"])
        iter_inds = setdiff(1:num_Mv, results["restart_inds"] .+ (1:length(results["restart_inds"])) )
        CG_niters_i_mean = results["CG_niters_i"]
        plot!([iter_inds], CG_niters_i_mean[iter_inds], lw=lw_true, color=colors[i], label=lagend_label)
        if CG_tol_scale_choice == "hdmd"
            plot!(axis2, [iter_inds], max.(0.5*eps(), results["tol_sternheimer_all"][iter_inds]./17.663871318548775), lw=lw_true, color=colors[i], ls=:dash, label="")
        else
            plot!(axis2, [iter_inds], max.(0.5*eps(), results["tol_sternheimer_all"][iter_inds]), lw=lw_true, color=colors[i], ls=:dash, label="")
        end
        # scatter!(axis4, results["restart_inds"], max.(0.5*eps(), results["tol_sternheimer_all"][results["restart_inds"].+(1:length(results["restart_inds"]))]),
        # label="", m=m_restart, color=i, markersize=m_sz_restart)
    end
end
end
plot!(xtickfont=font(14, "times"), ytickfont=font(1, "times"),ytickfontcolor=:white,
    size=(945, 405), dpi=300,xlims=(0, max_xaxis-2), margin=0Plots.mm,ylims=(-1,31))
plot!(axis2,yscale=:log10,yticks=10.0 .^ (Vector( -16:2:0)),ylims=(1e-16,1e-1),ytickfont=font(12, "times"))
plot!(legendfont=font(12, "times"),legend=:right)

fig = plot(fig1, fig2, layout=grid(1,2, widths=(2/3,1/3)), size=(945, 405), dpi=300, 
left_margin=0Plots.mm, bottom_margin=1Plots.mm, top_margin=0Plots.mm)
savefig(fig, save_to_dir * "CG_niters" * save_to_file * ".svg")
