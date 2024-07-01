#moved functions to module on 5/8/24
include("MakeProblems.jl")
##  Make 1 lattice 
lattice_inst = MakeProblems.make_instance3D([25, 25, 1 ], abs(rand(Int)));
maximum(nonzeros(lattice_inst.C))
graph = make_graph(lattice_inst)
plot_euc_graph(lattice_inst)

## make 1 euc...
euc_inst = make_instance_euc(500, abs(rand(Int)), conn = 4, Dim = 2)
graph = make_graph(euc_inst)
plot_euc_graph(euc_inst)


## Make euc problems....

Nvec = [50:500:2000; 2000:1000:20000]

conn = 4
for N in Nvec
    println(N)  
    for k in 1:10
        println("  $k" )
        euc_inst = make_instance_euc(N, N*1876+k, conn = conn, Dim = 3)
        @save "Problems\\euc_probs_disc\\$(N)_$(conn)conn_$(k)" euc_inst
    end
end

## Make lattice problems... 3D
Nvec = 5:50
for N in Nvec
    println(N)  
    for k in 1:10
        println("  $k" )
        lattice_inst = make_instance3D([N, N, 5], N*1876+k)
        @save "Problems\\lattice_probs_disc\\$(N)_$(k)" lattice_inst
    end
end
    
## Make 2D euc problems....
Nvec = [50:500:2000; 2000:1000:20000]
# Nvec = 5:100

conn = 4
for N in Nvec
    println(N)  
    for k in 1:10
        println("  $k" )
        euc_inst = make_instance_euc(N, N*1876+k, conn = conn, Dim = 2)
        @save "Problems\\euc_probs2D\\$(N)_$(k)" euc_inst #$(conn)conn_
    end
end

## Make 2D lattice problems...
Nvec = 5:100
Nvec = 51:100
for N in Nvec
    println(N)  
    for k in 1:10
        println("  $k" )
        lattice_inst = make_instance3D([N, N, 1], N*1876+k)
        @save "Problems\\lattice_probs_2D\\$(N)_$(k)" lattice_inst
    end
end



##############################################
## EXPERIMENT: make 2D euc problems of varying connectivity
##############################################
Nvec = [50:1000:2000; 2000:2000:20000]

connvec = [3, 4, 5, 8, 10, 12]
for N in Nvec
    # N < 14000 && continue
    printstyled("N =  $N \n",color=:green)
    for conn in connvec
        print("$(conn)conn, " )
        @threads for k in 1:10
            euc_inst = make_instance_euc(N, N*1876+k, conn = conn, Dim = 2)
            @save "Problems\\connectivity_expr\\$(N)_$(k)_$(conn)conn" euc_inst 
        end
    end
    println("")
end

##
import Statistics.mean
connvec = [3, 4, 5, 8, 10, 12]
for conn in connvec
    @load "Problems\\connectivity_expr\\2000_1_$(conn)conn" euc_inst
    euc = euc_inst
    A = euc.Alist
    Aavg = mean([length(A[i]) for i in 1:length(A)])
    println("$(conn)conn, $(Aavg) mean")
end    

## Draw connectivity plots
using Compose
connvec = [4 12]
for conn in connvec
    @load "Problems\\connectivity_expr\\2000_1_$(conn)conn" euc_inst
    #plot graph 
    graph = make_graph(euc_inst, false)
    plt = plot_euc_graph(euc_inst)
    plt
    #now save
    draw(PDF("Figs\\$(conn)conn.pdf", 8inch, 8inch), plt)

    #count total edges
    A = euc_inst.Alist
    Atotal = sum([length(A[i]) for i in 1:length(A)])/2
    println("$(conn)conn, $(Atotal) edges")
end


##############################################################################################
#old, do not need anymore....
##############################################################################################

## make specific instances of grid prolems
xx,yy = 3,3 #5, 10, 25
zz = 1
for k = 1:1
    seedy = k
    prob = make_instance3D(seedy, [xx,yy,zz])
    @save "Problems\\grid_probs3D\\$(xx)x$(yy)x$(zz)_$(seedy)" prob
end

## make and save graph def...
xx,yy,zz = 3,3,1
Alist,A, C, F =  make_ACF3D([xx,yy,zz])
Z = 2 .*C

A, C, F, Z = sparse(A), sparse(C), sparse(F), sparse(Z);
graphdef = GraphDef3D(Alist, A, F, C, Z);

@save "Problems\\grid_probs3D\\GraphDef_$(xx)x$(yy)x$(zz)" graphdef




##
using Graphs
using GraphPlot
using Gadfly
locs = rand(2,10)
g, dists = euclidean_graph(locs, cutoff = 0.5)
gplot(g, locs[1,:], locs[2,:])


## need to rename path to pathL....
Nvec = [50:500:2000; 2000:1000:20000]
for N in Nvec
    for k in 1:10
        try
            @load "Solutions\\euc_probs\\$(N)_4conn_$(k)" tdp cost path gen lX
            pathL = path
            @save "Solutions\\euc_probs\\$(N)_4conn_$(k)" tdp cost pathL gen lX
        catch
            println("No change")
        end
    end
end


## rename RCSPP as no RCSPP in filename.... bonk....

Nvec = [50:500:2000; 2000:1000:20000]

for N in Nvec
    println(N)
    for k in 1:10
        @load "Problems\\euc_probs2D\\$(N)_$(k)_RCSPP" euc_inst
        @save "Problems\\euc_probs2D\\$(N)_$(k)" euc_inst
    end
    sleep(10)
end

# ## make grid problems for 5x5 up to 50x50
# for nn = 5:50
#     seedy = nn + 101
#     def = make_instance(seedy, [nn,nn])
#     @save "Problems\\smooth_probs3D\\$(nn)x$(nn)" def

# end
# struct MapDef
#     Alist::Vector{Vector{Any}}
#     A::Array{Float64}
#     C::Array{Float64}
#     G::Array{Float64}
#     Z::Array{Float64}
#     anchor_list::Vector{Float64}
#     Dim::Vector{Float64} #[n.m]
#     nN::Int64 #number of houses (total N - n*m + nN)
# end

## testy boy....

# xx = 5
# yy = 2
# zz = 3
# top_side = []
# for k in 1:zz
#     a = xx*yy*k - xx + 1
#     b = xx*yy*k
#     println((a,b))
#     append!(top_side, a:b)
# end
# println(top_side)
