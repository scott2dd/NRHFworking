#=
each of these lines is a set of problems...... 
1) lattice v euc
2) 2d v 3D
3) Node v Label
4) astar v euc LB
=# 
using HybridUAVPlanning
solve_euc(algo = "node", dims = "2D", heur = "astar", tlim = 2*60)
solve_euc(algo = "label", dims = "2D", heur = "astar", tlim = 2*60)
solve_euc(algo = "node", dims = "3D", heur = "astar", tlim = 2*60)
solve_euc(algo = "label", dims = "3D", heur = "astar", tlim = 2*60)

solve_lattice(algo = "node", dims = "2D", heur = "astar", tlim = 2*60) 
solve_lattice(algo = "label", dim1s = "2D", heur = "astar", tlim = 2*60)


## Lattice 3D
# solve_lattice(algo = "node", dims = "3D", heur = "astar", tlim = 2*60)
solve_lattice(algo = "label", dims = "3D", heur = "astar", tlim = 2*60)
solve_lattice(algo = "node", dims = "3D", heur = "euc", tlim = 2*60) 
solve_lattice(algo = "label", dims = "3D", heur = "euc", tlim = 2*60)

##
solve_lattice(algo = "node", dims = "2D", heur = "euc", tlim = 1*60)
solve_lattice(algo = "label", dims = "2D", heur = "euc", tlim = 1*60)

solve_euc(algo = "node", dims = "2D", heur = "euc", tlim = 1*60) 
solve_euc(algo = "label", dims = "2D", heur = "euc", tlim = 1*60)
solve_euc(algo = "node", dims = "3D", heur = "euc", tlim = 1*60)
solve_euc(algo = "label", dims = "3D", heur = "euc", tlim = 1*60)

## Do our variable connectivity EXPERIMENT
connvec = [4, 5, 8, 10, 12]
for conn in connvec
    solve_euc(algo = "node", dims = "2D", heur = "astar", tlim = 2*60, conn = conn)
    solve_euc(algo = "label", dims = "2D", heur = "astar", tlim = 2*60, conn = conn)
end


## Solve 1, test temporal labeling algo
using JLD2 
using HybridUAVPlanning
using MAPFHybrid
constraints = MAPFHybrid.HybridConstraints(Set{MAPFHybrid.VertexConstraint}(), Set{MAPFHybrid.EdgeConstraint}())
# push!(constraints.vertex_constraints, HybridUAVPlanning.VertexConstraint(1, 31))
push!(constraints.edge_constraints, MAPFHybrid.EdgeConstraint(0, 1, 31))

MAPFHybrid.EdgeConstraint(0, 1, 31) in constraints.edge_constraints
@load "Problems/euc_probs_disc/50_4conn_1"
tdp = @elapsed cost, pathL = hybrid_label_temporal(euc_inst, constraints, 1, 1, 31)


## Test hybrid_label_selection_smart
using HybridUAVPlanning
using JLD2


@load "Problems\\euc_probs_disc\\5000_4conn_1" euc_inst
@time hybrid_label_selection_dumb(euc_inst, heur="astar")
@time hybrid_label_selection(euc_inst, heur="astar")
@time hybrid_node_selection(euc_inst, heur="astar")
@time hybrid_node_selection_dumb(euc_inst, heur="astar")

@load "Problems\\euc_probs2D\\10000_1" euc_inst
@time hybrid_label_selection_dumb(euc_inst, heur="astar")
@time hybrid_label_selection(euc_inst, heur="astar")
@time hybrid_node_selection_dumb(euc_inst, heur="astar")
@time hybrid_node_selection(euc_inst, heur="astar")


@load "Problems\\euc_probs2D\\20000_1" euc_inst
@time hybrid_label_selection_dumb(euc_inst, heur="astar")
@time hybrid_label_selection(euc_inst, heur="astar")
@time hybrid_node_selection_dumb(euc_inst, heur="astar")
@time hybrid_node_selection(euc_inst, heur="astar")


#lattice 2D
@load "Problems\\lattice_probs_2D\\50_1" lattice_inst
@time hybrid_label_selection_dumb(lattice_inst, heur="astar")
@time hybrid_label_selection(lattice_inst, heur="astar")




## Checking heuristics....
using JLD2, SparseArrays
@load "Problems\\lattice_probs_disc\\5_1" lattice_inst
@time hybrid_node_selection(lattice_inst, heur = "euc");

###############################
#0 performance testing.....
# @load "Problems\\euc_probs_disc\\50_4conn_1" euc_inst
# @code_warntype hybrid_label_selection(euc_inst)
#####################################################
## need to rename.... we made an error in euc3D naming DataStructures
Nvec =  [50:500:2000; 2000:1000:20000]
for N in Nvec
    for k in 1:10
        @load "Solutions\\euc_probs_disc\\$(N)_$(k)_4conn_node_eucLB" tdp cost pathL gen
        @save "Solutions\\euc_probs_disc\\$(N)_4conn_$(k)_node_eucLB" tdp cost pathL gen
    end
end


