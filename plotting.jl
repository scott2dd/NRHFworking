## SECTIONS
#1: Stacked plots for  Node label, 2d 3d euc lattice....
#2: LB comparison
#3: CPLEX
#4: Table print

#rest is misc?

# Color Key:
#     NODE: red
#     LABEL: grey    
#     Else:       dashed (LB)
using CairoMakie
using Gadfly
using HybridUAVPlanning
CairoMakie.activate!()
using HybridUAVPlanning
using JLD2, SparseArrays
import Statistics.mean
println("Packages loaded")
################################################################################
## AVG PLOT WITH BOX PLOTS PER PROBLEM
# EUC 3D
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("euc", "euc_probs_disc", algo = "_node")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("euc", "euc_probs_disc", algo = "_label")
times_NODE_LB, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("euc", "euc_probs_disc", algo = "_node", heur = "_eucLB")
times_LABEL_LB, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("euc", "euc_probs_disc", algo = "_label",  heur = "_eucLB")
#get max and main_across both times
max_time = max(maximum(times_NODE), maximum(times_LABEL), maximum(times_NODE_LB), maximum(times_LABEL_LB))
min_time = min(minimum(times_NODE), minimum(times_LABEL), minimum(times_NODE_LB), minimum(times_LABEL_LB))
max_time = 10^ceil(log10(max_time))
min_time = 10^floor(log10(min_time))

f = Figure(resolution = (480,800))
ax1= Axis(f[1,1],
    ylabel="Mean Time to Solve (s)",
    yscale = log10,xticklabelsvisible = false, xticksvisible=false)
ax2 = Axis(f[2,1],
    ylabel="Time to Solve (s)",
    yscale = log10,xticklabelsvisible = false, xticksvisible=false)
ax3 = Axis(f[3,1],
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values])
nodeline=lines!(ax1,Nvec1, avg_times_NODE, color = :red)
lblline=lines!(ax1,Nvec2, avg_times_LABEL, color = :grey)
nodelbline = lines!(ax1,Nvec3, avg_times_NODE_euc_LB, linestyle = :dot, color = "tomato", linewidth = 2)
lbllbline = lines!(ax1,Nvec4, avg_times_LABEL_euc_LB, linestyle = :dot, color = "gray", linewidth = 2.5)
nodebox = boxplot!(ax2, repeat(Nvec1,inner=10), reshape(times_NODE',(length(times_NODE),1))[:], color = :red, width = 800)
labelbox = boxplot!(ax3, repeat(Nvec2,inner=10), reshape(times_LABEL',(length(times_LABEL),1))[:], color = :grey, width = 800)
axislegend(ax1, [nodeline,lblline, nodelbline, lbllbline],["NODE-SUP", "LABEL-SUP", "NODE-SLD", "LABEL-SLD"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
axislegend(ax2, [nodebox],["NODE-SUP"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
axislegend(ax3, [labelbox],["LABEL-SUP"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
[ylims!(ax, min_time, max_time) for ax in [ax1, ax2, ax3]]

save("Figs\\main_euc3D.pdf",f)
f

## EUC 2D
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("euc", "euc_probs2D", algo = "_node", conn = "")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("euc", "euc_probs2D", algo = "_label", conn = "")
times_NODE_LB, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("euc", "euc_probs2D", algo = "_node", heur = "_eucLB",conn = "")
times_LABEL_LB, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("euc", "euc_probs2D", algo = "_label",  heur = "_eucLB", conn = "")
#get max and main_across both times
max_time = max(maximum(times_NODE), maximum(times_LABEL), maximum(times_NODE_LB), maximum(times_LABEL_LB))
min_time = min(minimum(times_NODE), minimum(times_LABEL), minimum(times_NODE_LB), minimum(times_LABEL_LB))
max_time = 10^ceil(log10(max_time))
min_time = 10^floor(log10(min_time))

f = Figure(resolution = (480,800))
ax1= Axis(f[1,1],
    ylabel="Mean Time to Solve (s)",
    yscale = log10,xticklabelsvisible = false, xticksvisible=false)
ax2 = Axis(f[2,1],
    ylabel="Time to Solve (s)",
    yscale = log10,xticklabelsvisible = false, xticksvisible=false)
ax3 = Axis(f[3,1],
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values])
linkxaxes!(ax1,ax2,ax3)
nodeline=lines!(ax1,Nvec1, avg_times_NODE, color = :red)
lblline=lines!(ax1,Nvec2, avg_times_LABEL, color = :grey)
nodelbline = lines!(ax1,Nvec3, avg_times_NODE_euc_LB, linestyle = :dot, color = "tomato", linewidth = 2)
lbllbline = lines!(ax1,Nvec4, avg_times_LABEL_euc_LB, linestyle = :dot, color = "gray", linewidth = 2.5)
nodebox = boxplot!(ax2, repeat(Nvec1,inner=10), reshape(times_NODE',(length(times_NODE),1))[:], color = :red, width = 800)
labelbox = boxplot!(ax3, repeat(Nvec2,inner=10), reshape(times_LABEL',(length(times_LABEL),1))[:], color = :grey, width = 800)
axislegend(ax1, [nodeline,lblline, nodelbline, lbllbline],["NODE-SUP", "LABEL-SUP", "NODE-SLD", "LABEL-SLD"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
axislegend(ax2, [nodebox],["NODE-SUP"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
axislegend(ax3, [labelbox],["LABEL-SUP"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
[ylims!(ax, min_time, max_time) for ax in [ax1, ax2, ax3]]
save("Figs\\main_euc2D.pdf",f)
f

## lattice 3D
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node", conn = "")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label", conn = "")
times_NODE_LB, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node", heur = "_eucLB",conn = "")
times_LABEL_LB, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label",  heur = "_eucLB",conn = "")
#get max and main_across both times
max_time = max(maximum(times_NODE), maximum(times_LABEL), maximum(times_NODE_LB), maximum(times_LABEL_LB))
min_time = min(minimum(times_NODE), minimum(times_LABEL), minimum(times_NODE_LB), minimum(times_LABEL_LB))
max_time = 10^ceil(log10(max_time))
min_time = 10^floor(log10(min_time))

maxlat3D = max_time
minlat3D = min_time

f = Figure(resolution = (480,800))
ax1= Axis(f[1,1],
    ylabel="Mean Time to Solve (s)",
    yscale = log10,xticklabelsvisible = false, xticksvisible=false)
ax2 = Axis(f[2,1],
    ylabel="Time to Solve (s)",
    yscale = log10,xticklabelsvisible = false, xticksvisible=false)
ax3 = Axis(f[3,1],
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values])
nodeline=lines!(ax1,5Nvec1.^2, avg_times_NODE, color = :red)
lblline=lines!(ax1,5Nvec2.^2, avg_times_LABEL, color = :grey)
nodelbline = lines!(ax1,5Nvec3.^2, avg_times_NODE_euc_LB, linestyle = :dot, color = "tomato", linewidth = 2)
lbllbline = lines!(ax1,5Nvec4.^2, avg_times_LABEL_euc_LB, linestyle = :dot, color = "gray", linewidth = 2.5)
linkaxes!(ax1,ax2,ax3)
nodebox = boxplot!(ax2, repeat(5Nvec1.^2,inner=10), reshape(times_NODE',(length(times_NODE),1))[:], color = :red, width = 200)
labelbox = boxplot!(ax3, repeat(5Nvec2.^2,inner=10), reshape(times_LABEL',(length(times_LABEL),1))[:], color = :grey, width = 200)
axislegend(ax1, [nodeline,lblline, nodelbline, lbllbline],["NODE-SUP", "LABEL-SUP", "NODE-SLD", "LABEL-SLD"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
axislegend(ax2, [nodebox],["NODE-SUP"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
axislegend(ax3, [labelbox],["LABEL-SUP"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
[ylims!(ax, min_time, max_time) for ax in [ax1, ax2, ax3]]

save("Figs\\main_lattice3D.pdf",f)
f
 
## lattice 2D
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node", conn = "")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label", conn = "")
times_NODE_LB, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node", heur = "_eucLB",conn = "")
times_LABEL_LB, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label",  heur = "_eucLB",conn = "")
#get max and main_across both times
max_time = max(maximum(times_NODE), maximum(times_LABEL), maximum(times_NODE_LB), maximum(times_LABEL_LB))
min_time = min(minimum(times_NODE), minimum(times_LABEL), minimum(times_NODE_LB), minimum(times_LABEL_LB))
max_time = 10^ceil(log10(max_time))
min_time = 10^floor(log10(min_time))

maxlat2D = max_time
minlat2D = min_time

max_max = max(maxlat2D, maxlat3D)
min_min = min(minlat2D, minlat3D)

f = Figure(resolution = (480,800))
ax1= Axis(f[1,1],
    ylabel="Mean Time to Solve (s)",
    yscale = log10,xticklabelsvisible = false, xticksvisible=false)
ax2 = Axis(f[2,1],
    ylabel="Time to Solve (s)",
    yscale = log10,xticklabelsvisible = false, xticksvisible=false)
ax3 = Axis(f[3,1],
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values])
linkaxes!(ax1,ax2,ax3)
nodeline=lines!(ax1,Nvec1.^2, avg_times_NODE, color = :red)
lblline=lines!(ax1,Nvec2.^2, avg_times_LABEL, color = :grey)
nodelbline = lines!(ax1,Nvec3.^2, avg_times_NODE_euc_LB, linestyle = :dot, color = "tomato", linewidth = 2)
lbllbline = lines!(ax1,Nvec4.^2, avg_times_LABEL_euc_LB, linestyle = :dot, color = "gray", linewidth = 2.5)
nodebox = boxplot!(ax2, repeat(Nvec1.^2,inner=10), reshape(times_NODE',(length(times_NODE),1))[:], color = :red, width = 50)
labelbox = boxplot!(ax3, repeat(Nvec2.^2,inner=10), reshape(times_LABEL',(length(times_LABEL),1))[:], color = :grey, width = 50)
axislegend(ax1, [nodeline,lblline, nodelbline, lbllbline],["NODE-SUP", "LABEL-SUP", "NODE-SLD", "LABEL-SLD"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
axislegend(ax2, [nodebox],["NODE-SUP"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
axislegend(ax3, [labelbox],["LABEL-SUP"], position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )
[ylims!(ax, min_time, max_time) for ax in [ax1, ax2, ax3]]
save("Figs\\main_lattice2D.pdf",f)
f

#################################################################################################
## Get connectivity_expr results

connvec = [4,5,8,10,12]
node_times_mat, node_avg_times = [], []
label_times_mat, label_avg_times = [], []
labelNvec, nodeNvec = Vector[],Vector[]
for conn in connvec
    node_times, mean_node_times, nodeNvec_conni= get_sol_vec_conn_expr(conn, "_node")
    label_times,mean_label_times,labelNvec_conni = get_sol_vec_conn_expr(conn, "_label")
    push!(node_times_mat, node_times)
    push!(node_avg_times, mean_node_times)
    push!(label_times_mat, label_times)
    push!(label_avg_times, mean_label_times)
    push!(labelNvec, labelNvec_conni)
    push!(nodeNvec, nodeNvec_conni)
end
#now plot avgs for each connectivity all on one plot
f = Figure(size = (480, 400))
ax = Axis(f[1,1],
xlabel="Number of Nodes",
ylabel="Mean Time to Solve (s)",
yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values])
line_styles = [:dash, :solid, :dot, :solid, :dash]
label_colors = [:paleturquoise4, :black, :paleturquoise4, :black, :paleturquoise4]
linewidths = [2.5, 1, 2.5, 1, 2.5]
for i in 1:length(connvec)
    Nveclabel = labelNvec[i]
    Nvecnode = nodeNvec[i]
    lines!(ax, Nvecnode, node_avg_times[i], color = :red, linestyle = line_styles[i], label = "NODE $(connvec[i])conn", linewidth = 2)
    lines!(ax, Nveclabel, label_avg_times[i], color = label_colors[i], linestyle = line_styles[i], label = "LABEL $(connvec[i])conn", linewidth = linewidths[i])
end
axislegend(ax, position = :rb, labelsize=10, patchsize = (30.0f0, 10.0f0), markersize = 20 )

save("Figs\\connectivity_expr.pdf",f)
f

##################################################################################################
## lower bound plots.....
##euc 3D LB....
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("euc", "euc_probs_disc", algo = "_node")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("euc", "euc_probs_disc", algo = "_label")
times_NODE, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("euc", "euc_probs_disc", algo = "_node", heur = "_eucLB")
times_LABEL, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("euc", "euc_probs_disc", algo = "_label",  heur = "_eucLB")

f = Figure(resolution = (800,600))
ax = Axis(f[1,1],
xlabel="Number of Nodes",
ylabel="Mean Time to Solve (s)",
yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values],)
lin1 = lines!(Nvec1, avg_times_NODE, color = :red, linewidth = 2)
lin2 = lines!(Nvec2, avg_times_LABEL, color = :grey, linewidth = 2.5)
lin3 = lines!(Nvec3, avg_times_NODE_euc_LB, linestyle = :dashdot, color = "tomato", linewidth = 3)
lin4 = lines!(Nvec4, avg_times_LABEL_euc_LB, linestyle = :dashdot, color = "gray", linewidth = 3)
axislegend(ax, [lin1,lin2,lin3,lin4],["NODE  ", "LABEL", "NODE-EucLB", "LABEL-EucLB"], position = :rb)
f
save("Figs\\LB_euc3D.pdf",f)


##euc 2D LB....
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("euc", "euc_probs2D", algo = "_node", conn = "")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("euc", "euc_probs2D", algo = "_label",conn = "")
times_NODE, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("euc", "euc_probs2D", algo = "_node", heur = "_eucLB",conn = "")
times_LABEL, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("euc", "euc_probs2D", algo = "_label",  heur = "_eucLB", conn = "")

f = Figure(resolution = (800,600))
ax = Axis(f[1,1],
xlabel="Number of Nodes",
ylabel="Mean Time to Solve (s)",
yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values])
lin1 = lines!(Nvec1, avg_times_NODE, color = :red, linewidth = 2)
lin2 = lines!(Nvec2, avg_times_LABEL, color = :grey, linewidth = 2.5)
lin3 = lines!(Nvec3, avg_times_NODE_euc_LB, linestyle = :dashdot, color = "tomato", linewidth = 3)
lin4 = lines!(Nvec4, avg_times_LABEL_euc_LB, linestyle = :dashdot, color = "gray", linewidth = 3)
axislegend(ax, [lin1,lin2,lin3,lin4],["NODE  ", "LABEL", "NODE-EucLB", "LABEL-EucLB"], position = :rb)
f
save("Figs\\LB_euc2D.pdf",f)


##Lattice 3D LB....
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node",conn = "")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label",conn = "")
times_NODE, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node", heur = "_eucLB",conn = "")
times_LABEL, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label",  heur = "_eucLB",conn = "")
f = Figure(resolution = (800,600))
ax = Axis(f[1,1],
xlabel="Number of Nodes",
ylabel="Mean Time to Solve (s)",
yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values])
lin1 = lines!(5Nvec1.^2, avg_times_NODE, color = :red, linewidth = 2)
lin2 = lines!(5Nvec2.^2, avg_times_LABEL, color = :grey, linewidth = 2.5)
lin3 = lines!(5Nvec3.^2, avg_times_NODE_euc_LB, linestyle = :dashdot, color = "tomato", linewidth = 3)
lin4 = lines!(5Nvec4.^2, avg_times_LABEL_euc_LB, linestyle = :dashdot, color = "gray", linewidth = 3)
axislegend(ax, [lin1,lin2,lin3,lin4],["NODE  ", "LABEL", "NODE-EucLB", "LABEL-EucLB"], position = :rb)
f
save("Figs\\LB_lattice3D.pdf",f)



##Lattice 2D LB....
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node",conn = "")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label",conn = "")
times_NODE, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node", heur = "_eucLB",conn = "")
times_LABEL, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label",  heur = "_eucLB",conn = "")
f = Figure(resolution = (800,600))
ax = Axis(f[1,1],
xlabel="Number of Nodes",
ylabel="Mean Time to Solve (s)",
yscale = log10,xtickformat = values -> ["$(Int(value))" for value in values])
lin1 = lines!(Nvec1.^2, avg_times_NODE, color = :red, linewidth = 2)
lin2 = lines!(Nvec2.^2, avg_times_LABEL, color = :grey, linewidth = 2.5)
lin3 = lines!(Nvec3.^2, avg_times_NODE_euc_LB, linestyle = :dashdot, color = "tomato", linewidth = 3)
lin4 = lines!(Nvec4.^2, avg_times_LABEL_euc_LB, linestyle = :dashdot, color = "gray", linewidth = 3)
axislegend(ax, [lin1,lin2,lin3,lin4],["NODE  ", "LABEL", "NODE-EucLB", "LABEL-EucLB"], position = :rb)
f
save("Figs\\LB_lattice2D.pdf",f)


#############################################################################
## CPLEX Plots
_, avg_times1, Nvec1 = HybridUAVPlanning.get_sol_vec_old("euc", "euc_probs_disc", algo = "_CPLEX", prob = "MILP")
_, avg_times2, Nvec2 = HybridUAVPlanning.get_sol_vec_old("euc", "euc_probs2D", algo="_CPLEX", prob="MILP", conn="")
_, avg_times3, Nvec3 = HybridUAVPlanning.get_sol_vec_old("lattice", "lattice_probs_disc_old_10_23_23", algo="_CPLEX", prob="MILP", conn="")
_, avg_times4, Nvec4 = HybridUAVPlanning.get_sol_vec_old("lattice", "lattice_probs_2D", algo="_CPLEX", prob="MILP", conn="")



f = Figure(resolution = (500,300))
ax = Axis(f[1,1],
xlabel="Number of Nodes",
ylabel="Mean Time to Solve (s)",
yscale = log10,)
lin1 = lines!(Nvec1, avg_times1, color = :red)
lin2 = lines!(Nvec2, avg_times2, color = :grey)
lin3 = lines!(Nvec3, avg_times3, color = :cyan)
lin4 = lines!(Nvec4, avg_times4, color = :magenta)
axislegend(ax, [lin1,lin2,lin3,lin4],["Euc3D", "Euc2D", "Lattice3D", "Lattice2D"], position = :rb)
save("Figs\\CPLEX_TTS.pdf",f)
f





#######################################################################################
##Plot a solution
N = 550
k  = 5
string1 = "euc_probs2D"
@load "Problems\\$(string1)\\$(N)_$(k)" euc_inst
@load "Solutions\\$(string1)\\$(N)_$(k)_node" tdp cost pathL gen
plt_sol = plot_euc_graph(euc_inst, path = pathL, gen = gen, color_ends = false);
draw(PDF("Figs/TASE_path_gen_ex.pdf", 10cm, 10cm), plt)




## Plot euc and lattice graph examples
@load "Problems\\euc_probs2D\\550_5" euc_inst
@load "Problems\\lattice_probs_2D\\25_1" lattice_inst

plt_euc = plot_euc_graph(euc_inst)
plt_lattice = plot_euc_graph(lattice_inst)


draw(PDF("Figs/TASE_euc_graph.pdf", 10cm, 10cm), plt_euc)  
draw(PDF("Figs/TASE_lattice_graph.pdf", 10cm, 10cm), plt_lattice)  


#############################################################################
## Big boy table for EUclidean results
#first load all our data....
Nvec = [50:500:2000; 2000:1000:20000]
_, NODE_3D, Nvec1 = get_sol_vec("euc", "euc_probs_disc", algo = "_node", conn = "_4conn")
_, NODE_3D_LB, Nvcec2 = get_sol_vec("euc", "euc_probs_disc", algo = "_node", conn = "_4conn", heur = "_eucLB")
_, LABEL_3D, Nvec3 = get_sol_vec("euc", "euc_probs_disc", algo = "_label", conn = "_4conn")
_, LABEL_3D_LB, Nvec4 = get_sol_vec("euc", "euc_probs_disc", algo = "_label", conn = "_4conn", heur = 
"_eucLB")

_, NODE_2D, Nvec5 = get_sol_vec("euc", "euc_probs2D", algo = "_node", conn = "")
_, NODE_2D_LB, Nvec6 = get_sol_vec("euc", "euc_probs2D", algo = "_node", conn = "", heur = "_eucLB")
_, LABEL_2D, Nvec7 = get_sol_vec("euc", "euc_probs2D", algo = "_label", conn = "")
_, LABEL_2D_LB, Nvec8 = get_sol_vec("euc", "euc_probs2D", algo = "_label", conn = "", heur = 
"_eucLB")

_, CPLEX_3D, Nvec9 = get_sol_vec("euc", "euc_probs_disc", algo = "_CPLEX", prob = "MILP")
_, CPLEX_2D, Nvec10 = get_sol_vec("euc", "euc_probs2D", algo = "_CPLEX", prob = "MILP", conn = "")

CPLEX_3D = vcat(CPLEX_3D, fill("--", length(Nvec)-length(CPLEX_3D)))
CPLEX_2D = vcat(CPLEX_2D, fill("--", length(Nvec)-length(CPLEX_2D)))

NODE_3D_LB = vcat(NODE_3D_LB, fill("--", length(Nvec)-length(NODE_3D_LB)))
NODE_2D_LB = vcat(NODE_2D_LB, fill("--", length(Nvec)-length(NODE_2D_LB)))

NODE_2D = vcat(NODE_2D, fill("--", length(Nvec)-length(NODE_2D)))

LABEL_3D_LB = vcat(LABEL_3D_LB, fill("--", length(Nvec)-length(LABEL_3D_LB)))
LABEL_2D_LB = vcat(LABEL_2D_LB, fill("--", length(Nvec)-length(LABEL_2D_LB)))

data_table = hcat(Nvec, CPLEX_3D, LABEL_3D, LABEL_3D_LB, NODE_3D, NODE_3D_LB, Nvec, CPLEX_2D, LABEL_2D, LABEL_2D_LB, NODE_2D, NODE_2D_LB)


#now print table......
println("----------------TABLE DATA  EUC Problems   ---------- \n")
idx = 0
for N in Nvec #each row is a size... print avg for each case... with " &" after
    idx += 1
    print(" $(data_table[idx,1]) ")
    for data in data_table[idx, 2:end]
        if data isa Int
            print(" & $(data)")
        elseif data isa Number
            print("&   $(round(data,digits=3))")
        else #if a string (--)
            print("&   $(data)")
        end
    end
    println(" \\\\ \\hline")

end

############################################################################
## Big boy table for LATTICE results
#first load all our data....
Nvec = 5:50

_, NODE_3D, Nvec1 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node", conn = "")
_, NODE_3D_LB, Nvec2 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node", conn = "", heur = "_eucLB")
_, LABEL_3D, Nvec3 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label", conn = "")
_, LABEL_3D_LB, Nvec4 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label", conn = "", heur = "_eucLB")

_, NODE_2D, Nvec5 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node", conn = "")
_, NODE_2D_LB, Nvec6 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node", conn = "", heur = "_eucLB")
_, LABEL_2D, Nvec7 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label", conn = "")
_, LABEL_2D_LB, Nvec8 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label", conn = "", heur = "_eucLB")

_, CPLEX_3D, Nvec9 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_CPLEX", prob = "MILP", conn="")
_, CPLEX_2D, Nvec10 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_CPLEX", prob = "MILP", conn = "")

CPLEX_3D = vcat(CPLEX_3D, fill("--", length(Nvec)-length(CPLEX_3D)))
CPLEX_2D = vcat(CPLEX_2D, fill("--", length(Nvec)-length(CPLEX_2D)))

NODE_3D_LB = vcat(NODE_3D_LB, fill("--", length(Nvec)-length(NODE_3D_LB)))
NODE_2D_LB = vcat(NODE_2D_LB, fill("--", length(Nvec)-length(NODE_2D_LB)))

LABEL_3D_LB = vcat(LABEL_3D_LB, fill("--", length(Nvec)-length(LABEL_3D_LB)))
LABEL_2D_LB = vcat(LABEL_2D_LB, fill("--", length(Nvec)-length(LABEL_2D_LB)))

LABEL_3D = vcat(LABEL_3D, fill("--", length(Nvec)-length(LABEL_3D)))
NODE_3D = vcat(NODE_3D, fill("--", length(Nvec)-length(NODE_3D)))


data_table = hcat(5*(5:50).^2, CPLEX_3D, LABEL_3D, LABEL_3D_LB, NODE_3D, NODE_3D_LB,
(5:50).^2, CPLEX_2D, LABEL_2D, LABEL_2D_LB, NODE_2D, NODE_2D_LB)

#now print table......
println("----------------TABLE DATA  LATTICE Problems   ---------- \n")
idx = 0
for N in Nvec #each row is a size... print avg for each case... with " &" after
    idx += 1
    print(" $(data_table[idx,1]) ")
    for data in data_table[idx, 2:end]
        if data isa Int
            print(" & $(data)")
        elseif data isa Number
            print(" & $(round(data,digits=3))")
        else #if a string (--)
            print(" & $(data)")
        end
    end
    println(" \\\\ \\hline")

end


###############################################################################################################################################################################################################
## Get average % NR edges
perc_NR = []
for N in [50:500:2000; 2000:1000:20000]
    println("N = $N")
    for k in 1:10
        @load "Problems\\euc_probs_disc\\$(N)_4conn_$(k)" euc_inst
        total_edges = sum(length.(euc_inst.Alist))
        nr_edges = length(nonzeros(euc_inst.GFlipped))
        perc_NR_k = nr_edges/total_edges
        push!(perc_NR, perc_NR_k)
        @load "Problems\\euc_probs2D\\$(N)_$(k)" euc_inst
        total_edges = sum(length.(euc_inst.Alist))
        nr_edges = length(nonzeros(euc_inst.GFlipped))
        perc_NR_k = nr_edges/total_edges
        push!(perc_NR, perc_NR_k)
    end
end

for N in 5:50
    println("N = $N")
    for k in 1:10
        @load "Problems\\lattice_probs_disc\\$(N)_$(k)" lattice_inst
        total_edges = sum(length.(lattice_inst.Alist))
        nr_edges = length(nonzeros(lattice_inst.GFlipped))
        perc_NR_k = nr_edges/total_edges
        push!(perc_NR, perc_NR_k)
        @load "Problems\\lattice_probs_2D\\$(N)_$(k)" lattice_inst
        total_edges = sum(length.(lattice_inst.Alist))
        nr_edges = length(nonzeros(lattice_inst.GFlipped))
        perc_NR_k = nr_edges/total_edges
        push!(perc_NR, perc_NR_k)
    end
end

mean_NR = mean(perc_NR)

## Old boxplots....
# EUC3D
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("euc", "euc_probs_disc", algo = "_node")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("euc", "euc_probs_disc", algo = "_label")
stacked_plt = get_stacked_plots([Nvec1, Nvec2], [times_NODE, times_LABEL], [avg_times_NODE, avg_times_LABEL])
draw(PDF("Figs/main_results_EUC3D.pdf", 12cm, 20cm), stacked_plt)  
# EUC2D
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("euc", "euc_probs2D", algo = "_node", conn = "")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("euc", "euc_probs2D", algo = "_label", conn = "")
stacked_plt = get_stacked_plots([Nvec1, Nvec2], [times_NODE, times_LABEL], [avg_times_NODE, avg_times_LABEL])
draw(PDF("Figs/main_results_EUC2D.pdf", 12cm, 20cm), stacked_plt)  


# LATTICE 3D
times1, avg_times1, Nvec1 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node", conn = "")
times2, avg_times2, Nvec2 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label", conn = "")
stacked_plt = get_stacked_plots([Vector(Nvec1), Vector(Nvec2)], [times1, times2], [avg_times1, avg_times2])
draw(PDF("Figs/main_results_lattice3D.pdf", 12cm, 20cm), stacked_plt)  

# LATTICE 2D    
times1, avg_times1, Nvec1 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node", conn = "")
times2, avg_times2, Nvec2 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label", conn = "")
stacked_plt = get_stacked_plots([Vector(Nvec1), Vector(Nvec2)], [times1, times2], [avg_times1, avg_times2])
draw(PDF("Figs/main_results_lattice2D.pdf", 12cm, 20cm), stacked_plt)  

## OLD LB PLOTS
# EUC3D-LB
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("euc", "euc_probs_disc", algo = "_node", heur = "_eucLB")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("euc", "euc_probs_disc", algo = "_label",  heur = "_eucLB")
stacked_plt = get_stacked_plots([Nvec1, Nvec2], [times_NODE, times_LABEL], [avg_times_NODE, avg_times_LABEL])
draw(PDF("Figs/main_results_EUC3D_LB.pdf", 12cm, 20cm), stacked_plt)  
# EUC2D
times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("euc", "euc_probs2D", algo = "_node", conn = "",  heur = "_eucLB")
times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("euc", "euc_probs2D", algo = "_label", conn = "",  heur = "_eucLB")
stacked_plt = get_stacked_plots([Nvec1, Nvec2], [times_NODE, times_LABEL], [avg_times_NODE, avg_times_LABEL])
draw(PDF("Figs/main_results_EUC2D_LB.pdf", 12cm, 20cm), stacked_plt)  


# LATTICE 3D
times1, avg_times1, Nvec1 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node", conn = "",  heur = "_eucLB")
times2, avg_times2, Nvec2 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label", conn = "",  heur = "_eucLB")
stacked_plt = get_stacked_plots([Vector(Nvec1), Vector(Nvec2)], [times1, times2], [avg_times1, avg_times2])
draw(PDF("Figs/main_results_lattice3D_LB.pdf", 12cm, 20cm), stacked_plt)  

# LATTICE 2D
times1, avg_times1, Nvec1 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node", conn = "",  heur = "_eucLB")
times2, avg_times2, Nvec2 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label", conn = "",  heur = "_eucLB")
stacked_plt = get_stacked_plots([Vector(Nvec1), Vector(Nvec2)], [times1, times2], [avg_times1, avg_times2])
draw(PDF("Figs/main_results_lattice2D_LB.pdf", 12cm, 20cm), stacked_plt)  

# # plt = plot(
#     # layer(x = Nvec1, y = avg_times_NODE, Geom.line, Geom.point, Theme(default_color = "red"), shape = [:circle]),
#     # layer(x = Nvec2, y = avg_times_LABEL, Geom.line, Geom.point, Theme(default_color = "grey"), shape = [:circle]),
#     # layer(x = Nvec3, y = avg_times_NODE_euc_LB, Geom.line, Geom.point, Theme(default_color = "tomato", line_style = [:dashdot]), shape = [:dtriangle]),
#     # layer(x = Nvec4, y = avg_times_LABEL_euc_LB, Geom.line, Geom.point, Theme(default_color = "lightgrey", line_style = [:dashdot]), shape= [:dtriangle]),
#     # Scale.y_log10,
#     # Scale.x_continuous(format= :plain),
#     # Guide.xlabel(" Number of Nodes"),
#     # Guide.ylabel(" Mean Time to Solve (s)"),
#     # Guide.title("EUCLIDEAN GRAPH 3D - Straight Line Distance LB"),
#     # # Theme(key_position = :top),
#     # Guide.manual_color_key("", ["NODE ", "LABEL ", "NODE-EucLB   ", "LABEL-EucLB "], ["red","grey", "tomato", "lightgrey"]),
#     # Guide.shapekey(title = "", labels = ["NODE ", "LABEL ", "NODE-EucLB   ", "LABEL-EucLB "])
# # )
# # set_default_plot_size(10cm, 10cm)
# # draw(PDF("Figs/LB_euc3D.pdf", 10cm, 10cm), plt)  
# ## lattice 3D LB....
# times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node",conn = "")
# times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label",conn = "")
# times_NODE, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_node", heur = "_eucLB",conn = "")
# times_LABEL, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("lattice", "lattice_probs_disc", algo = "_label",  heur = "_eucLB",conn = "")
# plt = plot(
#     layer(x = Nvec1, y = avg_times_NODE, Geom.line, Geom.point, Theme(default_color = "red")),
#     layer(x = Nvec2, y = avg_times_LABEL, Geom.line, Geom.point, Theme(default_color = "grey")),
#     layer(x = Nvec3, y = avg_times_NODE_euc_LB, Geom.line, Geom.point, Theme(default_color = "tomato"), shape=[Shape.vline]),
#     layer(x = Nvec4, y = avg_times_LABEL_euc_LB, Geom.line, Geom.point, Theme(default_color = "lightgrey"), shape=[Shape.vline]),
#     Scale.y_log10,
#     Scale.x_continuous(format= :plain),
#     Guide.xlabel(" Number of Nodes"),
#     Guide.ylabel(" Mean Time to Solve (s)"),
#     Guide.title("LATTICE GRAPH 3D - Straight Line Distance LB"),
#     Theme(key_position = :top),
#     Guide.manual_color_key("", ["NODE  ", "LABEL  ", "NODE-EucLB  ", "LABEL-EucLB  "], ["red","grey", "tomato", "lightgrey"]),
# )
# draw(PDF("Figs/LB_euc2D.pdf", 10cm, 10cm), plt)  


# # lattice 2D LB....
# times_NODE, avg_times_NODE, Nvec1 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node",conn = "")
# times_LABEL, avg_times_LABEL, Nvec2 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label",conn = "")
# times_NODE, avg_times_NODE_euc_LB, Nvec3 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_node", heur = "_eucLB",conn = "")
# times_LABEL, avg_times_LABEL_euc_LB, Nvec4 = get_sol_vec("lattice", "lattice_probs_2D", algo = "_label",  heur = "_eucLB",conn = "")
# plt = plot(
#     layer(x = Nvec1, y = avg_times_NODE, Geom.line, Geom.point, Theme(default_color = "red")),
#     layer(x = Nvec2, y = avg_times_LABEL, Geom.line, Geom.point, Theme(default_color = "grey")),
#     layer(x = Nvec3, y = avg_times_NODE_euc_LB, Geom.line, Geom.point, Theme(default_color = "tomato"), shape=[Shape.vline]),
#     layer(x = Nvec4, y = avg_times_LABEL_euc_LB, Geom.line, Geom.point, Theme(default_color = "lightgrey"), shape=[Shape.vline]),
#     Scale.y_log10,
#     Scale.x_continuous(format= :plain),
#     Guide.xlabel(" Number of Nodes"),
#     Guide.ylabel(" Mean Time to Solve (s)"),
#     Guide.title("LATTICE GRAPH 3D - Straight Line Distance LB"),
#     Theme(key_position = :top),
#     Guide.manual_color_key("", ["NODE  ", "LABEL  ", "NODE-EucLB  ", "LABEL-EucLB  "], ["red","grey", "tomato", "lightgrey"]),
# )



## Old CPLEX
set_default_plot_size(15cm, 15cm)
plt = plot(layer(x=Nvec1, y=avg_times1, Geom.line, Geom.point, Theme(default_color="red")),
    layer(x=Nvec2, y=avg_times2, Geom.line, Geom.point, Theme(default_color="grey")),
    layer(x=5 * Nvec3 .^ 2, y=avg_times3, Geom.line, Geom.point, Theme(default_color="magenta")),
    layer(x=Nvec4 .^ 2, y=avg_times4, Geom.line, Geom.point, Theme(default_color="blue")),
    Scale.y_log10,
    Guide.xlabel("Number of Nodes "),
    Guide.ylabel("TIme to Solve (s)"),
    Guide.manual_color_key("", ["Euc3D", "Euc2D", "lattice3D", "lattice2D"], ["red", "grey", "magenta", "blue"]),
    Theme(key_position=:top),
    Guide.title("Time to Solve - CPLEX")
)




draw(PDF("Figs/CPLEX_time_to_solve.pdf", 15cm, 15cm), stacked_plt)


## Plot average path length vs graph size....
Nvec = [50:500:2000; 2000:1000:20000]
lengths = zeros(length(Nvec), 10) 
avg_lengths = zeros(length(Nvec))

for N in Nvec
    for k in 1:10
        @load 
    end
end





##plot1: euc3D label labeleuc node solve_euc2D_with_node_euc
Nvec =  [50:500:2000; 2000:1000:20000]
Nvec1 = [50:500:2000; 2000:1000:4000]
Nvec2 = [50;550]
times_NODE, avg_NODE = get_sol_vec(Nvec, "euc_probs_disc", algo = "_node", conn = "_4conn")
times_NODE_LB, avg_NODE_LB = get_sol_vec(Nvec2, "euc_probs_disc", algo = "_node", conn = "_4conn", heur = "_eucLB")

times_LABEL, avg_LABEL = get_sol_vec(Nvec, "euc_probs_disc", algo = "_label", conn = "_4conn")
times_LABEL_LB, avg_LABEL_LB = get_sol_vec(Nvec1, "euc_probs_disc", algo = "_label", conn = "_4conn", heur = "_eucLB")


plt = plot(
    layer(x = Nvec, y = avg_NODE, Geom.line, Geom.point, Theme(default_color = "red")),
    layer(x = Nvec, y = avg_LABEL, Geom.line, Geom.point, Theme(default_color = "grey")),
    layer(x = Nvec2, y = avg_NODE_LB, Geom.line, Geom.point, Theme(default_color = "tomato"), shape=[Shape.vline]),
    layer(x = Nvec1, y = avg_LABEL_LB, Geom.line, Geom.point, Theme(default_color = "lightgrey"), shape=[Shape.vline]),
    Scale.y_log10,
    Scale.x_continuous(format= :plain),
    Guide.ylabel(" Mean Time to Solve (s)"),
    Guide.xlabel(" Number of Nodes"),
    Guide.title("EUCLIDEAN GRAPH 3D - Straight Line Distance LB"),
    Theme(key_position = :top),
    Guide.manual_color_key("", [" NODE ", " LABEL ", " NODE-EucLB ", " LABEL-EucLB "], ["red","grey", "tomato", "lightgrey"]),
)
set_default_plot_size(12cm, 20cm)
draw(PDF("Figs/LB_comparison_EUC3D.pdf", 12cm, 20cm), plt)  



##plot2: euc2D with LBs
Nvec =  [50:500:2000; 2000:1000:20000]
Nvec1 = [50:500:2000; 2000:1000:4000]
Nvec2 = [50:500:2000; 2000:1000:3000]
get_sol_vec("lattice", "euc_p", algo = "_node", conn = "")
times_NODE, avg_NODE = get_sol_vec(Nvec, "euc_probs2D", algo = "_node", conn = "")
times_NODE_LB, avg_NODE_LB = get_sol_vec(Nvec2, "euc_probs2D", algo = "_node", conn = "", heur = "_eucLB")

times_LABEL, avg_LABEL = get_sol_vec(Nvec, "euc_probs2D", algo = "_label", conn = "")
times_LABEL_LB, avg_LABEL_LB = get_sol_vec(Nvec, "euc_probs2D", algo = "_label", conn = "", heur = "_eucLB")


plt = plot(
    layer(x = Nvec, y = avg_NODE, Geom.line, Geom.point, Theme(default_color = "red")),
    layer(x = Nvec, y = avg_LABEL, Geom.line, Geom.point, Theme(default_color = "grey")),
    layer(x = Nvec2, y = avg_NODE_LB, Geom.line, Geom.point, Theme(default_color = "tomato"), shape=[Shape.vline]),
    layer(x = Nvec, y = avg_LABEL_LB, Geom.line, Geom.point, Theme(default_color = "lightgrey"), shape=[Shape.vline]),
    Scale.y_log10,
    Scale.x_continuous(format= :plain),
    Guide.xlabel(" Number of Nodes"),
    Guide.ylabel(" Mean Time to Solve (s)"),
    Guide.title("EUCLIDEAN GRAPH 2D - Straight Line Distance LB"),
    Theme(key_position = :top),
    Guide.manual_color_key("", [" NODE ", " LABEL ", " NODE-EucLB ", " LABEL-EucLB "], ["red","grey", "tomato", "lightgrey"]),
)

# draw(PDF("Figs/LB_comparison_EUC2D.pdf", 12cm, 20cm), plt) 


#plot3: lattice2D with LBs
Nvec =  5:50

times_NODE, avg_NODE = get_sol_vec(Nvec, "lattice_probs_2D", algo = "_node", conn = "")
times_NODE_LB, avg_NODE_LB = get_sol_vec(Nvec, "lattice_probs_2D", algo = "_node", conn = "", heur = "_eucLB")

times_LABEL, avg_LABEL = get_sol_vec(Nvec, "lattice_probs_2D", algo = "_label", conn = "")
times_LABEL_LB, avg_LABEL_LB = get_sol_vec(Nvec, "lattice_probs_2D", algo = "_label", conn = "", heur = "_eucLB")


plt = plot(
    layer(x = Nvec.^2, y = avg_NODE, Geom.line, Geom.point, Theme(default_color = "red")),
    layer(x = Nvec.^2, y = avg_LABEL, Geom.line, Geom.point, Theme(default_color = "grey")),
    layer(x = Nvec.^2, y = avg_NODE_LB, Geom.line, Geom.point, Theme(default_color = "tomato"), shape=[Shape.vline]),
    layer(x = Nvec.^2, y = avg_LABEL_LB, Geom.line, Geom.point, Theme(default_color = "lightgrey"), shape=[Shape.vline]),
    Scale.y_log10,
    Scale.x_continuous(format= :plain),
    Guide.xlabel(" Number of Nodes"),
    Guide.ylabel(" Mean Time to Solve (s)"),
    Guide.title("LATTICE GRAPH 2D - Straight Line Distance LB"),
    Theme(key_position = :top),
    Guide.manual_color_key("", [" NODE ", " LABEL ", " NODE-EucLB ", " LABEL-EucLB "], ["red","grey", "tomato", "lightgrey"]),
)

draw(PDF("Figs/LB_comparison_LATTICE2D.pdf", 10cm, 8cm), plt) 


#plot4: lattice3D with LBs
Nvec =  5:50
Nvec1 = 5:7
Nvec2 = 5:19
Nvec3 = 5:36
Nvec4 = 5:31
times_NODE, avg_NODE = get_sol_vec(Nvec1, "lattice_probs_disc", algo = "_node", conn = "")
times_NODE_LB, avg_NODE_LB = get_sol_vec(Nvec2, "lattice_probs_disc", algo = "_node", conn = "", heur = "_eucLB")

times_LABEL, avg_LABEL = get_sol_vec(Nvec3, "lattice_probs_disc", algo = "_label", conn = "")
times_LABEL_LB, avg_LABEL_LB = get_sol_vec(Nvec4, "lattice_probs_disc", algo = "_label", conn = "", heur = "_eucLB")


plt = plot(
    layer(x = Nvec1.^2, y = avg_NODE, Geom.line, Geom.point, Theme(default_color = "red")),
    layer(x = Nvec3.^2, y = avg_LABEL, Geom.line, Geom.point, Theme(default_color = "grey")),
    layer(x = Nvec2.^2, y = avg_NODE_LB, Geom.line, Geom.point, Theme(default_color = "tomato"), shape=[Shape.vline]),
    layer(x = Nvec4.^2, y = avg_LABEL_LB, Geom.line, Geom.point, Theme(default_color = "lightgrey"), shape=[Shape.vline]),
    Scale.y_log10,
    Scale.x_continuous(format= :plain),
    Guide.xlabel(" Number of Nodes"),
    Guide.ylabel(" Mean Time to Solve (s)"),
    Guide.title("LATTICE GRAPH3D - Straight Line Distance LB"),
    Theme(key_position = :top),
    Guide.manual_color_key("", [" NODE ", " LABEL ", " NODE-EucLB ", " LABEL-EucLB "], ["red","grey", "tomato", "lightgrey"]),
)

draw(PDF("Figs/LB_comparison_LATTICE3D.pdf", 10cm, 8cm), plt) 



###########################################################
## 2:  box plot for lattice probs and euc probs AND average plot
NvecL = [5:8:13; 18:3:25; 26:2:50]
NvecE = [50:500:2000; 2000:1000:20000]
xvecL = 5*NvecL.^2
times_L, avg_times_L = get_sol_vec(NvecL, "lattice_probs", conn = "")
bpltL = get_boxplot_plt(xvecL, times_L, color = "cyan", xlabel = "Number of Nodes", xmax = 20000)
    
times_e, avg_times_e = get_sol_vec(NvecE, "euc_probs")
bpltE = get_boxplot_plt(NvecE, times_e, color = "red", xlabel = "")


avg_plt = plot(
    layer(x = NvecE, y = avg_times_e, Geom.line, Geom.point, Theme(default_color = "red")),
    layer(x = xvecL, y = avg_times_L, Geom.line, Geom.point, Theme(default_color = "cyan")),
    Scale.y_log10,
    Scale.x_continuous(format= :plain),
    Guide.xlabel(""),
    Guide.ylabel(" Mean Time to Solve (s)"),
    Theme(key_position = :top),
    Guide.manual_color_key("", ["Euclidean", "Lattice"], ["red","cyan"]),
)
set_default_plot_size(15cm, 45cm)
stacked = vstack(avg_plt, bpltE, bpltL)
title(stacked, "")



###########################################################
##  3: Node Selection Results  --EUC--
Nvec =  [50:500:2000; 2000:1000:20000]
times_p, avg_times_p = get_sol_vec(Nvec, "euc_probs_disc", algo = "_node")
bpltNODE = get_boxplot_plt(Nvec, times_p, color = "red", xlabel = "Number of Nodes")
set_default_plot_size(20cm, 12cm)
stack = vstack(bpltNODE)
title(stack, "Time to Solve - Discrete Node Selection")


############################################################








###########################################################
##  4: Pruning results (euc)
Nvec =  [50:500:2000; 2000:1000:20000]
times_p, avg_times_p = get_sol_vec(Nvec, "euc_probs2D_prune")
bplt1 = get_boxplot_plt(Nvec, times_p, color = "grey", xlabel = "")
times, avg_times = get_sol_vec(Nvec, "euc_probs2D")
bplt2 = get_boxplot_plt(Nvec[1:8], times[1:8,:], color = "red", xlabel = "Number of Nodes", xmax = Nvec[end])

plt = plot(layer(x = Nvec, y = avg_times_p, Geom.line, Geom.point, Theme(default_color = "grey")),
           layer(x = Nvec[1:8], y = avg_times[1:8], Geom.line, Geom.point, Theme(default_color = "red")),
           Scale.y_log10,
           Scale.x_continuous(format= :plain),
           Guide.xlabel(""),
           Guide.ylabel(" Mean Time to Solve (s)"),
           Guide.manual_color_key("", ["Pruning", "No Prune"], ["grey", "red"]),
           Theme(key_position = :top)
            )
stacked = vstack(plt, bplt1, bplt2)
set_default_plot_size(15cm, 25cm)
title(stacked, "")

###########################################################
##  5: Plot comparison of cost-to-go for lattice AND euc
NvecE = [50:500:2000; 2000:1000:3000] 
NvecL = 5:22
timesEA, avg_timesEA = get_sol_vec(NvecE, "euc_probs")
timesEE, avg_timesEE = get_sol_vec(NvecE, "euc_probs_euc")

timesLA, avg_timesLA = get_sol_vec(NvecL, "lattice_probs", conn = "")
timesLE, avg_timesLE = get_sol_vec(NvecL, "lattice_probs_euc", conn = "")
set_default_plot_size(15cm, 10cm)
#plot for average time for euc 
pltE = plot(layer(x = NvecE, y = avg_timesEA, Geom.line, Geom.point, Theme(default_color = "red")),
           layer(x = NvecE, y = avg_timesEE, Geom.line, Geom.point, Theme(default_color = "grey")),
           Scale.y_log10,
           Scale.x_continuous(format= :plain),
           Guide.xlabel(""),
           Guide.ylabel(" Mean Time to Solve (s)"),
        #    Theme(key_position = :top),
           Guide.manual_color_key("", ["A* LB","Euclidean LB"], ["red", "grey"])
            )
    
#plot for average time for lattice
pltL = plot(layer(x = 5*NvecL.^2, y = avg_timesLA, Geom.line, Geom.point, Theme(default_color = "cyan")),
           layer(x = 5*NvecL.^2, y = avg_timesLE, Geom.line, Geom.point, Theme(default_color = "grey")),
           Scale.y_log10,
           Scale.x_continuous(format= :plain),
           Guide.xlabel(""),
           Guide.ylabel(" Mean Time to Solve (s)"),
            #   Theme(key_position = :top),
           Guide.manual_color_key("", ["A* LB", "Euclidean LB"], ["cyan", "grey"])
            )

display(pltE)
display(pltL)


###########################################################
## 6: Plot a 2D solution
using Gadfly
N = 550
k  = 5
string1 = "euc_probs2D"
@load "Problems\\$(string1)\\$(N)_$(k)" euc_inst
@load "Solutions\\$(string1)\\$(N)_$(k)_node" tdp cost pathL gen
plt = plot_euc_graph(euc_inst, path = pathL, gen = gen, color_ends = false);
# draw(PDF("Figs/TASE_path_gen_ex.pdf", 10cm, 10cm), plt)  
plt
 

## Plot Generator vs noise restricted zones....
N = 550
k  = 5
string1 = "euc_probs2D"
@load "Problems\\$(string1)\\$(N)_$(k)" euc_inst
@load "Solutions\\$(string1)\\$(N)_$(k)_node" tdp cost pathL gen

plt_gen = gen_plot(pathL, gen, euc_inst)
display(plt_gen)
draw(PDF("Figs/TASE_gen_and_noiseR.pdf", 11cm, 6cm), plt_gen)

###############################################################
## 7: Plot parametrized graph with an edge or few labeling power, time, etc.
N = 50
k  = 2      
string1 = "euc_probs2D"
string2 = "4conn_"
@load "Problems\\$(string1)\\$(N)_$(string2)$(k)" euc_inst
euc_inst.GFlipped .= 0
graph = make_graph(euc_inst)
@load "Solutions\\$(string1)\\$(N)_$(string2)$(k)" tdp cost pathL gen lX
plt = plot_euc_graph_labeled(euc_inst, label_strings = ["P_ij", "t_ij",], label_vals = [10. 15.], label_edge_idxs = [12],  label_units = ["W","s"], color_ends = false);
# display(plt)
draw(PDF("Figs/sol_ex.pdf", 10cm, 10cm), plt)  #if draw, then displaying after is 


###############################################################
## 8: plot 


## testy boy
using GraphPlot
using Graphs
g = smallgraph(:karate)
nodelabel = 1:nv(g)
plt = gplot(g, nodelabel=nodelabel);
display(plt)

## 
using Graphs: smallgraph
using GraphPlot
g = smallgraph(:karate)
gplot(g)



##################################################################
## Misc code.....
##

#load same problem as prior.....
gen_with_0 = [0; gen]
time_pointsMILP = 0:1/(length(gen_with_0)-1):1

noiseR_along_path = Float64.([!euc_inst.GFlipped[pathL[idx-1],pathL[idx]] for idx in 2:length(pathL)])

noiseR_along_path= [0;noiseR_along_path]
gen_repeat = repeat(gen,inner=2)
time_pointsMILP = repeat(time_points, inner=2)

_ = popfirst!(gen_with_0)
_ = popfirst!(time_pointsMILP)
_ = pop!(time_pointsMILP)
_ = popfirst!(gen_with_0)

plt_gen = plot(
            layer(x = time_pointsMILP, y = gen_repeat, Geom.line, Theme(default_color = "black")),
            layer(xmin=time_points[1:end-1], xmax = time_points[2:end].+ 0.001, Geom.vband, color = G[1:end],), #max gen as ribbons....
            Scale.color_discrete_manual("white", "lightcoral"),
            Guide.colorkey(title="", labels = [""]),
            Guide.manual_color_key("",["Noise Resricted   ","Gen Throttle"],["lightcoral","black"]),
            Guide.xlabel("Time (normalized)"),
            Guide.ylabel("Gen Throttle"),
            Theme(grid_line_width = 0pt),
            Coord.cartesian(xmin = 0, xmax = 1),
            Theme(key_position=:top,
            )
)



G = [!euc_inst.GFlipped[pathL[idx-1],pathL[idx]] for idx in 2:length(pathL)] #1 less than time_vec... 

time = Vector(0:length(gen))
timerep = repeat(time, inner = 2)
popfirst!(timerep)
pop!(timerep)
genrep = repeat(gen, inner = 2)

#normalize time vecs
time = time./time[end]
timerep = timerep./timerep[end]

plt_gen = plot(
            layer(x= timerep, y = genrep, Geom.line, Theme(default_color = "black")), 
            layer(xmin=time[1:end-1], xmax = time[2:end].+ 0, Geom.vband, color = G), #max gen as ribbons....
            Scale.color_discrete_manual("white", "lightcoral"),
            Guide.colorkey(title="", labels = [""]),
            Guide.manual_color_key("",["Noise Resricted   ", "Gen Throttle"],["lightcoral", "black"]),
            Guide.xlabel("Time (normalized)"),
            Guide.ylabel("Gen State,  g ∈ {0,1}"),
            Guide.yticks(ticks = [0,1]),
            Coord.cartesian(xmin = 0, xmax = time[end]),
            Theme(key_position=:top)
)

##
set_default_plot_size(15cm, 12cm)
display(plt_gen)
draw(PDF("Figs/TASE_gen_and_noiseR.pdf", 15cm, 12cm), plt_gen)
# compose(context(0,0,1.0, 0.5), render(plt_gen))
#now make random plot for legend to past into plot...
# plt = plot(x = 1:5, y = 1:5, Guide.manual_color_key("",["Gen On", "Gen Off"],["magenta", "cyan"]))
# display(plt)



pathA1 = a_star(graph, euc_inst.S, euc_inst.E, euc_inst.C)
pathA, LB, costA = astar_proc(pathA1, euc_inst.C)
plot_euc_graph(euc_inst)
plain_plt = plot_euc_graph(euc_inst)
plt = plot_euc_graph(euc_inst, path = pathL, gen = gen)
stacked = vstack(plt, plain_plt)



bool = false
(pathA == pathL) && (bool = true)
set_default_plot_size(20cm, 40cm)
title(stacked, "N=$(N) | k=$(k) | Path == SPP: $(bool)")
##
################################
#     Solution Plotting        #
################################
function a_star_plot(size, instance; dim = "")
    @load "Problems/euc_probs$(dim)/"

end

function hybrid_sol_plot(size, instance)

end 

plot_euc_graph(euc_inst, path = path, gen = gen)
plt = plot(x = 1:length(dist), y = dist, Geom.line)
graph = make_graph(euc_inst)
pathA1 = a_star(graph, path[end], euc_inst.E, euc_inst.C)
path_list, LB_vec, cost = astar_proc(astar_out, euc_inst.C)
pathA = [pathA1[i].src for i in 1:length(pathA1)]
append!(pathA, pathA1[end].dst)
plot_euc_graph(euc_inst, path = pathA, gen = zeros(length(pathA)))



## Plot avg time to solve, 2D vs 3D
Nvec2 = 30:40:2000
# Nvec = [50:50:2000; 3000:1000:20000]


tvec4 = zeros(length(Nvec2))
tvecSC = zeros(length(Nvec2))
cnt = 0
for n in Nvec2
    cnt += 1
    @load "Solutions\\euc_probs\\$(n)_4conn" tdp cost path gen lX
    tvec4[cnt] = tdp
    @load "Solutions\\euc_probsSC\\$(n)_4conn" tdp cost path gen lX
    tvecSC[cnt] = tdp
    
end
set_default_plot_size(25cm, 15cm)
Gadfly.push_theme(:default)

plt2 = Gadfly.plot(
    layer(x = Nvec2, y = tvec4, Geom.line, Theme(default_color = "grey")),
    layer(x = Nvec2, y = tvecSC, Geom.line, Theme(default_color = "red")),
    Guide.manual_color_key("", ["SPP", "GPenalty"], ["grey", "red"]),
    Guide.title("Time to Solve"),
    Guide.xlabel("Numer of Nodes"),
    Guide.ylabel("Time (s)"),
    Scale.x_continuous(format= :plain),
    Scale.y_log10()
           )








## Plot Time to Solve for EUC problems....
Nvec2 = 30:40:2000
# Nvec = [50:50:2000; 3000:1000:20000]


tvec4 = zeros(length(Nvec2))
tvecSC = zeros(length(Nvec2))
cnt = 0
for n in Nvec2
    cnt += 1
    @load "Solutions\\euc_probs\\$(n)_4conn" tdp cost path gen lX
    tvec4[cnt] = tdp
    @load "Solutions\\euc_probsSC\\$(n)_4conn" tdp cost path gen lX
    tvecSC[cnt] = tdp
    
end
set_default_plot_size(25cm, 15cm)
Gadfly.push_theme(:default)

plt2 = Gadfly.plot(
    layer(x = Nvec2, y = tvec4, Geom.line, Theme(default_color = "grey")),
    layer(x = Nvec2, y = tvecSC, Geom.line, Theme(default_color = "red")),
    Guide.manual_color_key("", ["SPP", "GPenalty"], ["grey", "red"]),
    Guide.title("Time to Solve"),
    Guide.xlabel("Numer of Nodes"),
    Guide.ylabel("Time (s)"),
    Scale.x_continuous(format= :plain),
    Scale.y_log10()
           )










