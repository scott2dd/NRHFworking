using Revise
import JSON
using Graphs
using HybridUAVPlanning
using MAPFHybrid
using JLD2
using Random
import Statistics.mean
using DataFrames    
using Base.Threads
using CairoMakie
using JLD2
using DataFrames
using GLM
using MAPFHybrid
using ColorSchemes

struct MAPFInstance
    inst_str::String
    starts::Vector{Int64}
    goals::Vector{Int64}
end
function make_MAPF_instance(nDim::Int64, nAgents::Int64; graph_inst::Int64 = 1)
    inst_str = "Problems/lattice_probs_2D/$(nDim)_$(graph_inst)"
    nNodes = nDim^2 #this work only as long as we are doing lattice 2D
    starts = randperm(nNodes)[1:nAgents]
    goals = Vector{Int64}(undef, nAgents)
    for i in 1:nAgents
        possible_goals = setdiff(1:nNodes, [starts; goals[1:i-1]])
        goals[i] = rand(possible_goals)
    end
    
    return MAPFInstance(inst_str, starts, goals)
end


function solve_MAPF(instance::MAPFInstance; time_lim = 1.0)
    inst_str = instance.inst_str
    loaded = load(inst_str) #lattice_inst
    lattice_inst = loaded["lattice_inst"]
    starts = instance.starts
    goals = instance.goals
    initial_states = HybridState[]
    [push!(initial_states, HybridState(nodeIdx=starts[i], time=0, b=0,g=0)) for i in 1:length(starts)]
    
    graph = make_graph(lattice_inst)
    env = HybridEnvironment(nNodes = nv(graph), goals=goals, obstacles=HybridLocation[], state_graph = graph, locs = lattice_inst.locs, orig_graph = graph, euc_inst = lattice_inst)
    
    CBS   = CBSSolver{HybridState,HybridAction,Int64,SumOfCosts,HybridConflict,HybridConstraints,HybridEnvironment}(env=env)
    ECBS = ECBSSolver{HybridState,HybridAction,Int64,SumOfCosts,HybridConflict,HybridConstraints,HybridEnvironment}(env=env, weight=1.3)
    time_CBS = @elapsed sol_nodeCBS, total_nodesCBS, times_subroutineCBS, times_astarCBS = search!(CBS, initial_states, time_lim = time_lim)
    time_ECBS = @elapsed sol_nodeECBS, total_nodesECBS, times_subroutineECBS, times_astarECBS = search!(ECBS, initial_states, time_lim = time_lim)
    return sol_nodeCBS, sol_nodeECBS, total_nodesCBS, total_nodesECBS, time_CBS, time_ECBS, times_subroutineCBS, times_subroutineECBS, times_astarCBS, times_astarECBS
end

function get_paths(solution)
    paths =  Vector{Int64}[]
    nAgents = length(solution)
    for agenti in 1:nAgents
        path = []
        states = solution[agenti].states
        for state_tuple in states
            node = state_tuple[1].nodeIdx
            push!(path, node)
        end
        push!(paths, path)  
    end
    return paths
end

function get_gens(solution)
    gens =  Vector{Int64}[]
    nAgents = length(solution)
    for agenti in 1:nAgents
        gen = solution[agenti].gen
        push!(gens, gen)
    end
    return gens
end

##############################################################################
## Solve TASE instances 
##############################################################################
@load "Problems/MAPF_instances" instances
solve_MAPF(instances[1], time_lim = 2.0); #JIT compile

dflength = length(instances)

df = DataFrame(inst=zeros(Int, dflength), nAgents=zeros(Int, dflength), Dim=zeros(Int, dflength), nNodes=zeros(Int, dflength), nCBSNodes=zeros(Int, dflength), nECBSNodes=zeros(Int, dflength), time_CBS=zeros(dflength), time_ECBS=zeros(dflength), times_subroutineCBS=fill([1.0], dflength), times_subroutineECBS=fill([1.0], dflength), times_astarCBS=fill([1.0], dflength), times_astarECBS=fill([1.0], dflength), cost_CBS=zeros(dflength), cost_ECBS=zeros(dflength))

CBSsol_nodes = Vector{MAPFHybrid.CBSHighLevelNode}(undef, length(instances))
ECBSsol_nodes = Vector{MAPFHybrid.ECBSHighLevelNode}(undef, length(instances))

time_start = time()
solved = zeros(Bool, dflength)
Threads.@threads for inst_idx in 1:dflength
    instance = instances[inst_idx]
    nAgents = length(instance.starts)
    loaded = load(instance.inst_str)
    lattice_inst = loaded["lattice_inst"]
    graph = make_graph(lattice_inst)
    nNodes = nv(graph)
    Dim = Int(sqrt(nNodes))
    time_elapsed = round(time() - time_start)
    sol_nodeCBS, sol_nodeECBS, total_nodesCBS, total_nodesECBS, time_CBS, time_ECBS, times_subroutineCBS, times_subroutineECBS, times_astarCBS, times_astarECBS = 
        solve_MAPF(instance, time_lim = 900.0)
    solved[inst_idx] = true
    printstyled("Solved $(sum(solved)) / $(dflength) instances       \r", color=:green)
    cost_CBS, cost_ECBS = sol_nodeCBS.cost, sol_nodeECBS.cost
    df[inst_idx, :] = (inst_idx, nAgents, Dim, nNodes, total_nodesCBS, total_nodesECBS, time_CBS, time_ECBS, times_subroutineCBS, times_subroutineECBS, times_astarCBS, times_astarECBS, cost_CBS, cost_ECBS)
    CBSsol_nodes[inst_idx] = sol_nodeCBS
    ECBSsol_nodes[inst_idx] = sol_nodeECBS
end


#now save the df and sync with onedrive 
results = df
@save "Solutions/MAPF_results" results CBSsol_nodes ECBSsol_nodes

run(`onedrive --synchronize --single-directory juliaOD/HybridVehicle/Solutions`) 

###########################################################################
## plotting  TASE
###########################################################################

@load "Solutions/MAPF_results" results CBSsol_nodes ECBSsol_nodes

## CBS 2D Plot (lines)
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10,
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    title="CBS Results"
)
#plot lines for time_CBS but only if time to solve is less than 900
group_out = []
for group in groupby(results, :nAgents)
    
    group_out = group
    #take mean of of time_CBS for each nNodes
    lines!(ax, group.nNodes[group.time_CBS .< 900], group.time_CBS[group.time_CBS .< 900], label="Agents = $(group.nAgents[1])")
    scatter!(ax, group.nNodes[group.time_CBS.<900], group.time_CBS[group.time_CBS.<900], markersize = 6)
end
#add a legend
fig[1,2] = Legend(fig, ax, "", framevisible=false)
fig


## 
result = map([1, 2, 3, 4, 5]) do x
    x * 2
end

## ECBS 2D Plot (lines)
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10,
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    title="ECBS Results"
)
#plot lines for time_CBS but only if time to solve is less than 900
for group in groupby(results, :nAgents)
    lines!(ax, group.nNodes[group.time_ECBS .< 900], group.time_ECBS[group.time_ECBS .< 900], label="Agents = $(group.nAgents[1])")
    scatter!(ax, group.nNodes[group.time_ECBS .< 900], group.time_ECBS[group.time_ECBS .< 900], markersize=6)
end

fig[1, 2] = Legend(fig, ax, "", framevisible=false)
fig

## Same plots as above, but subplotted
fig3 = Figure()
ax1 = Axis(fig3[1, 1], yscale=log10,
    ylabel="Time to Solve (s)",
    title="CBS Results",
    xlabel = "",
    limits = (nothing, nothing, 1*10^(-3.5), 10^3.5),
    yminorticksvisible = true,
    yminorticks = [1e3, 1e1, 1e-1, 1e-3]
)
ax2 = Axis(fig3[2, 1], yscale=log10,
    ylabel="Time to Solve (s)",
    title="ECBS Results",
    xlabel = "Number of Nodes",
    limits=(nothing, nothing, 1 * 10^(-3.5), 10^3.5),
    yminorticksvisible = true,
    yminorticks=[1e3, 1e1, 1e-1, 1e-3]
)
#plot lines 
for group in groupby(results, :nAgents)
    lines!(ax1, group.nNodes[group.time_CBS.<900], group.time_CBS[group.time_CBS.<900], label="Agents = $(group.nAgents[1])")
    scatter!(ax1, group.nNodes[group.time_CBS.<900], group.time_CBS[group.time_CBS.<900], markersize=6)

    lines!(ax2, group.nNodes[group.time_ECBS.<900], group.time_ECBS[group.time_ECBS.<900], label="Agents = $(group.nAgents[1])")
    scatter!(ax2, group.nNodes[group.time_ECBS.<900], group.time_ECBS[group.time_ECBS.<900], markersize=6)
end

fig3[1, 2] = Legend(fig3, ax1, "", framevisible=true)
fig3







## CBS 2D Plot (color bars)
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10,
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    title="CBS Results, unsolved problems not shown"
)

colormap = ColorSchemes.amp
colorfunc = x -> colormap[(x-minimum(results.nAgents))/(maximum(results.nAgents)-minimum(results.nAgents))]
for group in groupby(results, :nAgents)
    scatter!(ax, group.nNodes, group.time_CBS, label="Agents = $(group.nAgents[1])", color=colorfunc(group.nAgents[1]))
end

cb = Colorbar(fig[1, 2], limits=(minimum(results.nAgents), maximum(results.nAgents)), label="Number of Agents", colormap=colormap)
fig


## ECBS 2D Plot
fig = Figure()
ax = Axis(fig[1, 1], yscale=log10,
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    title="ECBS Results, unsolved problems not shown"
)

colormap = ColorSchemes.amp
colorfunc = x -> colormap[(x-minimum(results.nAgents))/(maximum(results.nAgents)-minimum(results.nAgents))]
for group in groupby(results, :nAgents)
    scatter!(ax, group.nNodes, group.time_ECBS, label="Agents = $(group.nAgents[1])", color=colorfunc(group.nAgents[1]))
end
cb = Colorbar(fig[1, 2], limits=(minimum(results.nAgents), maximum(results.nAgents)), label="Number of Agents", colormap=colormap)

fig

## Plot A* times for CBS vs subroutine times
import Statistics.mean
fig = Figure(fontsize=10)
ax = Axis(fig[1, 1],
    xlabel="Number of Nodes",
    ylabel="Mean time in subroutine (s)",
    yscale=log10,
    title="A* times vs Subroutine times"
)

colormap_astar = reverse(ColorSchemes.berlin)
colorfunc_astar = x -> colormap_astar[(x-minimum(results.nAgents))/(maximum(results.nAgents)-minimum(results.nAgents))]

colormap_label = ColorSchemes.buda
colorfunc_label = x -> colormap_label[(x-minimum(results.nAgents))/(maximum(results.nAgents)-minimum(results.nAgents))]

for group in groupby(results, :nAgents)
    mean_times_astar = mean.(group.times_astarCBS)
    mean_times_label = mean.(group.times_subroutineCBS)
    scatter!(ax, group.nNodes, mean_times_astar, label="Agents = $(group.nAgents[1])", color=colorfunc_astar(group.nAgents[1]), markersize=6, marker = :circle)
    scatter!(ax, group.nNodes, mean_times_label, label="Agents = $(group.nAgents[1])", color=colorfunc_label(group.nAgents[1]), markersize=12, marker = :hline)
end
cb = Colorbar(fig[1, 2], limits=(minimum(results.nAgents), maximum(results.nAgents)), label="Number of Agents - SPP (A*)", colormap=colormap_astar)
cb = Colorbar(fig[1, 3], limits=(minimum(results.nAgents), maximum(results.nAgents)), label="Number of Agents - NRHFSPP (Label)", colormap=colormap_label)

fig



## now plot ECBS/CBS solution quality.  
# TODO - need to save J* returned from ECBS and CBS!!


###############################################################################
## Making problems... run a bunch and pick out the non-trivial ones.  More than 1 node.  Then save problems.
###############################################################################
using JLD2
#do 5x5 up to 50x50....
instances = MAPFInstance[]
nAgentLB(x) = Int(round(0.16x + 1.5))
nAgentUB(x) = Int(round(0.88x + 6))

nAgentsVec = [5 10 15 30 40 50]
for nDim in 5:50
    # for nAgents in nAgentLB(nDim):nAgentUB(nDim)
    for nAgents in nAgentsVec
        println("nDim = $(nDim), nAgents = $(nAgents)")
        nAgents >= 0.4*nDim^2 && continue
        count = 1
        inst = 1
        while true #loop until we get 10 instances at this size
            count > 10 && (inst += 1)
            instance = make_MAPF_instance(nDim, nAgents, graph_inst=inst)
            sol_nodeCBS, sol_nodeECBS, total_nodesCBS, total_nodesECBS, time_CBS, time_ECBS, times_subroutineCBS, times_subroutineECBS, times_astarCBS, times_astarECBS = solve_MAPF(instance, time_lim=0.1)
            if total_nodesCBS > 1  #if not solved at root and not infeasible at root
                push!(instances, instance)
                count += 1
            end
            count >= 10 && break
        end
    end
end

##
@save "Problems/MAPF_instances" instances












##############################################################################
## ACC CODE!!!
##############################################################################
## Load instances and solve them
##############################################################################
@load "Problems/MAPF_instances_ACC.jld2" instances5 instances10 instances15

sol_node, total_nodes, time_master, times_subroutine, times_astar = solve_MAPF(instances10[1])
#init a DataFrame thatt stores, for each problem we solve, number of agents, number of nodes, number of CBS nodes, total time to solve, times to solve w/ subroutine, times to solve w/ astar
dflength = length(instances5) + length(instances10) + length(instances15)
df = DataFrame(dfidx=1:dflength, inst=zeros(Int, dflength), nAgents=zeros(Int, dflength), Dim=zeros(Int, dflength), nNodes=zeros(Int, dflength), nCBSNodes=zeros(Int, dflength), time_master=zeros(dflength), times_subroutine=fill([1.0], dflength), times_astar=fill([1.0], dflength))

zipped = hcat(1:(length(instances5)+length(instances10)+length(instances15)), [1:length(instances5); 1:length(instances10); 1:length(instances15)], [instances5; instances10; instances15])
#initialize a vector of length 168 of any type
sol_nodes = Vector{Any}(undef, 168)
Threads.@threads for rowk in 1:size(zipped, 1)
    dfidx, inst_idx, inst = zipped[rowk, :]
    nAgents = length(inst.starts)
    inst_str = inst.inst_str
    loaded = load(inst_str) #lattice_inst
    lattice_inst = loaded["lattice_inst"]
    graph = make_graph(lattice_inst)
    nNodes = nv(graph)
    Dim = sqrt(nNodes)
    printstyled("Solving dfidx$(dfidx) | instance $(inst_idx) with $(nAgents) agents, $(Dim) x $(Dim) grid, $(nNodes) nodes\n", color=:lightgreen)
    sol_node, CBSnodes, time_master, times_subroutine, times_astar = solve_MAPF(inst)
    # println("  Total time: $(time_master)")
    # println("  CBS Nodes: $(CBSnodes)")
    df[dfidx, :] = (dfidx, inst_idx, nAgents, Dim, nNodes, CBSnodes, time_master, times_subroutine, times_astar)
    sol_nodes[dfidx] = sol_node
end
printstyled("Done!\n", color=:green)

#now save the df!
results = df
@save "Solutions/MAPF_results_ACC.jld2" results sol_nodes



###########################################################################
# plotting  ACC
###########################################################################
## now plot boxplot with Makie
using CairoMakie
using JLD2
using DataFrames
@load "Solutions/MAPF_results_ACC.jld2" results sol_nodes


#plot boxplot of time to solve vs Dim
f = Figure(fontsize=30, resolution=(1000, 500))
ax1 = Axis(f[1, 1],
    xlabel="Number of Nodes",
    ylabel="Number of CBS Nodes",
    xticks=([25, 100, 225], ["5x5", "10x10", "15x15"]),
    yscale=log10)
#grab nCBSNodes column for Dim = 5
nCBSNodes5 = results[results.Dim.==5, :nCBSNodes]
nCBSNodes10 = results[results.Dim.==10, :nCBSNodes]
nCBSNodes15 = results[results.Dim.==15, :nCBSNodes]
boxplot!(fill(25, length(nCBSNodes5)), nCBSNodes5, width=50, label="4 UAVs")
boxplot!(fill(100, length(nCBSNodes10)), nCBSNodes10, width=50, label="5 UAVs")
boxplot!(fill(225, length(nCBSNodes15)), nCBSNodes15, width=50, label="10 UAVs")
#Now add legen using labels
f[1, 2] = Legend(f, ax1, "", framevisible=false)
#now show plot
CairoMakie.activate!()
save("FigsMAPF/CBSnodes_boxplot.pdf", f, pt_per_unit=4)
display(f)

#get maximum number of CBS nodes 
maxCBS = maximum(results.nCBSNodes)

## now do same with time to solve
f2 = Figure(fontsize=30, resolution=(1000, 500))
ax2 = Axis(f2[1, 1],
    xlabel="Number of Nodes",
    ylabel="Time to Solve (s)",
    xticks=([25, 100, 225], ["5x5", "10x10", "15x15"]),
    yscale=log10)
#grab nCBSNodes column for Dim = 5
time5 = results[results.Dim.==5, :time_master]
time10 = results[results.Dim.==10, :time_master]
time15 = results[results.Dim.==15, :time_master]
boxplot!(fill(25, length(nCBSNodes5)), time5, width=50, label="4 UAVs")
boxplot!(fill(100, length(nCBSNodes10)), time10, width=50, label="5 UAVs")
boxplot!(fill(225, length(nCBSNodes15)), time15, width=50, label="10 UAVs")
#Now add legen using labels
f2[1, 2] = Legend(f2, ax2, "", framevisible=false)
#now show plot
CairoMakie.activate!()
save("FigsMAPF/times_boxplot.pdf", f2, pt_per_unit=2)
display(f2)


## Plot Paths from a soln
using HybridUAVPlanning
using JLD2
using CairoMakie
@load "Solutions/MAPF_results_ACC.jld2" results sol_nodes
@load "Solutions/MAPF_instances_ACC.jld2" instances5 instances10 instances15
instances = vcat(instances5, instances10, instances15)
inst = 13
dim = 10
#get location in data frame for given inst and dim
idx = findfirst((results.inst .== inst) .& (results.Dim .== dim))
resulti = results[idx, :]
P_solved = sol_nodes[idx]
mapf = instances[idx]

P_init, _, _, _, _ = solve_MAPF(mapf, time_lim=0.0)
paths_final = get_paths(P_solved.solution)
paths_init = get_paths(P_init.solution)

#now get lattice inst for plotting....
@load mapf.inst_str lattice_inst

#now plot!
plt_final = plot_MAPF(lattice_inst, paths_final)
plt_init = plot_MAPF(lattice_inst, paths_init)

#now save each plot as pdf
using Gadfly
draw(PDF("FigsMAPF/paths_init.pdf", 10cm, 10cm), plt_init)
draw(PDF("FigsMAPF/paths_final.pdf", 10cm, 10cm), plt_final)
display(plt_init)
display(plt_final)

## plot gen from solution...
#for each agent, get final path and generator pattern

gens_init = gens = get_gens(P_init.solution)
gens = get_gens(P_solved.solution)
plts = []
for k in 1:length(paths_final) #for each, get gen and NR
    noiseRk = noiseR_along_paths[k]
    genk = gens[k]

    #call gen_plot
    legend = false
    # k == length(paths_final) && (legend = true)


    plt = gen_plot(paths_final[k], genk, lattice_inst, legendbool=legend, max_time=maximum(length.(paths_final)), uavi=k)
    push!(plts, plt)
end

#now stack plots into one plot
set_default_plot_size(20cm, 20cm)
plt = vstack(plts...)
#now save plt as pdf
draw(PDF("FigsMAPF/gens.pdf", 20cm, 20cm), plt)


##plot A* times vs labeling algo times, boxplot for each 


f3 = Figure(fontsize=30, resolution=(1000, 500))
ax2 = Axis(f3[1, 1],
    xlabel="Number of Nodes",
    ylabel="Mean Time In Subroutine (s)",
    yscale=log10,
    xticks=([25, 100, 225], ["5x5", "10x10", "15x15"])
)
#grab nCBSNodes column for Dim = 5
astar5 = mean.(results[results.Dim.==5, :times_astar])
astar10 = mean.(results[results.Dim.==10, :times_astar])
astar15 = mean.(results[results.Dim.==15, :times_astar])
label5 = mean.(results[results.Dim.==5, :times_subroutine])
label10 = mean.(results[results.Dim.==10, :times_subroutine])
label15 = mean.(results[results.Dim.==15, :times_subroutine])


boxplot!(fill(25, length(nCBSNodes5)), label5, color=:red, width=50)
boxplot!(fill(100, length(nCBSNodes10)), label10, color=:red, width=50)
boxplot!(fill(225, length(nCBSNodes15)), label15, color=:red, label="NRHF\nSPP", width=50)
boxplot!(fill(25, length(nCBSNodes5)), astar5, color=:cyan, width=50)
boxplot!(fill(100, length(nCBSNodes10)), astar10, color=:cyan, width=50)
boxplot!(fill(225, length(nCBSNodes15)), astar15, color=:cyan, label="SPP", width=50)

#Now add legen using labels
f3[1, 2] = Legend(f3, ax2, "Subroutine", framevisible=false)
#now show plot
display(f3)
CairoMakie.activate!()
save("FigsMAPF/subroutine_boxplot.pdf", f3, pt_per_unit=10)


## junk code from old approach...
using FileIO
#do 5's
nDim = 5
instances5 = MAPFInstance[]
nCBSNodes = Int64[]
nAgents = 4
for probi in 1:100
    inst = rand(1:10)
    instance = make_MAPF_instance(nDim, nAgents, inst)
    
    println("Problem $(probi)")
    println("  Starts: $(instance.starts)")
    println("  Goals: $(instance.goals)")
    #now solve....
    sol_node, total_nodes, time_master, times_subroutine,_ = solve_MAPF(instance)
    push!(nCBSNodes, total_nodes)
    println("  CBS Nodes: $(total_nodes)")
    println("  Time: $(time_master)")
    if sol_node.cost != -1 && total_nodes > 1 #if not timeout and not solved at root
        push!(instances5, instance)
    end
end



nDim = 10
instances10 = MAPFInstance[]
nCBSNodes = Int64[]
nAgents = 5
for probi in 1:100
    inst = rand(1:10)
    instance = make_MAPF_instance(nDim, nAgents, inst)
    
    println("Problem $(probi)")
    println("  Starts: $(instance.starts)")
    println("  Goals: $(instance.goals)")
    #now solve....
    sol_node, total_nodes, time_master, times_subroutine,_ = solve_MAPF(instance)
    push!(nCBSNodes, total_nodes)
    println("  CBS Nodes: $(total_nodes)")
    println("  Time: $(time_master)")
    if sol_node.cost != -1 && total_nodes > 1 #if not timeout and not solved at root
        push!(instances10, instance)
    end
end

##
nDim = 15
instances15 = MAPFInstance[]
nCBSNodes = Int64[]
nAgents = 10
for probi in 1:100
    inst = rand(1:10)
    instance = make_MAPF_instance(nDim, nAgents, inst)
    
    println("Problem $(probi)")
    println("  Starts: $(instance.starts)")
    println("  Goals: $(instance.goals)")
    #now solve....
    sol_node, total_nodes, time_master, times_subroutine,times_astar = solve_MAPF(instance)
    push!(nCBSNodes, total_nodes)
    println("  CBS Nodes: $(total_nodes)")
    println("  Time: $(time_master)")
    println("  total time A*: $(sum(times_astar))")
    if sol_node.cost != -1 && total_nodes > 1 #if not timeout and not solved at root
        push!(instances15, instance)
    end
end

# @save "Problems/MAPF_instances.jld2" instances5 instances10 instances15 


