using Graphs
using HybridUAVPlanning
using MAPFHybrid
using JLD2
using Random
import Statistics.mean
print("\033c")

@load "Problems/lattice_probs_2D/25_3"

nNodes = size(lattice_inst.C,1)
nAgents= 50
starts = randperm(nNodes)[1:nAgents]
perm = randperm(nNodes)
goals  = [x for x in perm if x ∉ starts][1:nAgents]

initial_states = HybridState[]
[push!(initial_states, HybridState(nodeIdx=starts[i], time=0, b=0,g=0)) for i in 1:length(starts)]

graph = make_graph(lattice_inst)
env = HybridEnvironment(nNodes = nv(graph), goals=goals, obstacles=HybridLocation[], state_graph = graph, locs = lattice_inst.locs, orig_graph = graph, euc_inst = lattice_inst)
solver = CBSSolver{HybridState,HybridAction,Int64,SumOfCosts,HybridConflict,HybridConstraints,HybridEnvironment}(env=env)
ecbs_solver = ECBSSolver{HybridState,HybridAction,Int64,SumOfCosts,HybridConflict,HybridConstraints,HybridEnvironment}(env=env, weight = 1.3)


time_master = @elapsed sol_node, total_nodes, times_subroutine = search!(solver, initial_states, time_lim=5*60.);
solution = []
time_ecbs = @elapsed sol_node_ecbs, total_nodes_ecbs, times_subroutine_ecbs = search!(ecbs_solver, initial_states, time_lim=5*60.);
solution = sol_node.solution
solution_ecbs = sol_node_ecbs.solution
if isempty(solution)
    println("Solution Empty from CBS")
else
    cost = 0
    makespan = 0
    for s in solution
        cost += s.cost
        makespan = max(s.cost, makespan)
    end
    print("\033c")
    printstyled("CBS Statistics:\n", bold=true, color=:blue)
    println("Cost: ", cost)
    println("Makespan: ", makespan)
    println("Total Nodes: ", total_nodes)
    println("Time: ", time_master)
    println("Avg subproblem time: ", round(mean(times_subroutine), digits=3))
end

if isempty(solution_ecbs)
    println("Solution Empty from ECBS")
else
    cost = 0
    makespan = 0
    for s in solution_ecbs
        cost += s.cost
        makespan = max(s.cost, makespan)
    end
    printstyled("ECBS Statistics:\n", bold=true, color=:blue)
    println("Cost: ", cost)
    println("Makespan: ", makespan)
    println("Total Nodes: ", total_nodes_ecbs)
    println("Time: ", time_ecbs)
    println("Avg subproblem time: ", round(mean(times_subroutine_ecbs), digits=3))
end


## Let's compare a call label_temporal to a call label_temporal_focal
agent_idx = 1
statei = MAPFHybrid.HybridState(nodeIdx=1, time=0, b=0,g=0)
goali = 2500
focal_state_heuristic(s) = 0
focal_transition_heuristic(s1a, s1b) = 0
constraints = MAPFHybrid.get_empty_constraint(MAPFHybrid.HybridConstraints)

t1 = @elapsed plan1 = MAPFHybrid.label_temporal_low_IQ(ecbs_solver.env, constraints, agent_idx, statei, goali);
println("Plan1 cost: ", plan1.cost)
t2 = @elapsed plan2 = MAPFHybrid.label_temporal(ecbs_solver.env, constraints, agent_idx, statei, goali);
plan2 != nothing && println("Plan2 cost: ", plan2.cost)
 
println("Speed up by: ",round(t1/t2,digits = 1 ), "x !")

t1 = @elapsed plan1 = MAPFHybrid.label_temporal_focal(ecbs_solver.env, constraints, agent_idx, statei, goali, 1.3, focal_state_heuristic, focal_transition_heuristic);
println("Speed up by: ", round(t1 / t2, digits=1), "x !")





##
using DataStructures
a = [10, (-1, 0, 0, 0,)]
b = [10, (1, 0, 0, 0)]
c = [10, (0, 0, 0, 0)]
d = [9, (10, 0,0,0  )]
h = MutableBinaryMinHeap([a])
push!(h, b)
push!(h, c)

topp = pop!(h)
println(topp)
