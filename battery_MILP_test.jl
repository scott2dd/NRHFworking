using CPLEX
using JuMP
using GLM
include("utils.jl")
include("Label_LB.jl")
Nvec = 5:100
kvec = 1:30
## first define function for using naive battery model

function MILP_RCSPP_naive(def)
    S, E, C = def.S, def.E, def.C
    Alist = def.Alist
    Bstart = 5000*3.6
    N = size(C, 1)
    A = zeros(N,N)
    for i in 1:length(Alist)
        Ai = Alist[i]
        for j in Ai
            A[i,j] = 1
            A[j,i] = 1

        end
    end

    P = rand((150:300), N,N)
    t = rand((20:30), N,N)
    bigM = 99
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0, "CPXPARAM_Emphasis_Memory" => 1, "CPXPARAM_WorkMem" => 9000,  "CPXPARAM_MIP_Pool_Capacity" => 0))
    @variable(m, x[i=1:N,j=1:N], Bin)
    @constraint(m, [i=1:N], x[i,i] == 0)

    @variable(m, Bstart >= b[i=1:N] >= 0)
    @constraint(m, [i = 1:N, j = 1:N; A[i,j] == 0], x[i,j] == 0)
    @constraint(m, sum(x[S,j] for j=1:N) == 1)
    @constraint(m, sum(x[j,E] for j=1:N) == 1)

    @constraint(m, sum(x[j,S] for j=1:N) == 0)
    @constraint(m, sum(x[E,j] for j=1:N) == 0)

    @constraint(m, [i = setdiff(1:N, [S,E])], sum(x[i,j] for j=1:N)-sum(x[j,i] for j=1:N) == 0)


    
    @constraint(m, b[S] == Bstart)
    @constraint(m, [i = 1:N, j = 1:N], b[j] <=  b[i] - P[i,j]*t[i,j]/15 + bigM*(1-x[i,j]))


    @objective(m, Min, sum(x[i,j]*C[i,j] for  i=1:N, j=1:N))

    optimize!(m)
    time = solve_time(m)
    return time
end

#now wrtie function for using linear battery model
function MILP_RCSPP_linear(def)
    S, E, C = def.S, def.E, def.C
    Bstart = 5000*3.6
    Alist = def.Alist

    N = size(C, 1)
    A = zeros(N,N)
    for i in 1:length(Alist)
        Ai = Alist[i]
        for j in Ai
            A[i,j] = 1
            A[j,i] = 1

        end
    end
    
    P = rand((150:300), N,N)
    t = rand((20:30), N,N)
    Req = 0.0012
    OCV, Smin, Smax = get_OCV_func_LiPo()
    onebyVL = get_one_by_VspGLM(OCV, 100:-5:Smin, Req, maximum(P))
    Cc, Aa, Bb = coef(onebyVL)

    bigM = 99
    m = Model(optimizer_with_attributes(CPLEX.Optimizer, "CPX_PARAM_SCRIND" => 0, "CPXPARAM_Emphasis_Memory" => 1, "CPXPARAM_WorkMem" => 9000,  "CPXPARAM_MIP_Pool_Capacity" => 0))
    @variable(m, x[i=1:N,j=1:N], Bin)
    @constraint(m, [i=1:N], x[i,i] == 0)

    @variable(m, Bstart >= b[i=1:N] >= 0)
    @constraint(m, [i = 1:N, j = 1:N; A[i,j] == 0], x[i,j] == 0)
    @constraint(m, sum(x[S,j] for j=1:N) == 1)
    @constraint(m, sum(x[j,E] for j=1:N) == 1)

    @constraint(m, sum(x[j,S] for j=1:N) == 0)
    @constraint(m, sum(x[E,j] for j=1:N) == 0)

    @constraint(m, [i = setdiff(1:N, [S,E])], sum(x[i,j] for j=1:N)-sum(x[j,i] for j=1:N) == 0)


    @constraint(m, b[S] == Bstart)



    
    
    @constraint(m, [i = 1:N, j = 1:N], b[j] <=  b[i] - P[i,j]*(Aa*P[i,j] + Bb*b[i] +  Cc)*t[i,j] + bigM*(1-x[i,j]))


    @objective(m, Min, sum(x[i,j]*C[i,j] for  i=1:N, j=1:N))

    optimize!(m)
    xsol = value.(x)
    time = solve_time(m)
    return time
end

function RCSPP_label(def::EucGraph, model; heur = "astar", params = [1,1,1])
    S, E, C = def.S, def.E, def.C
    Bstart = 5000*3.6
    Qstart = 999
    Alist = def.Alist

    N = size(C, 1)
    locs= def.locs

    
    Power = rand((150:300), N,N)
    t = rand((20:30), N,N)
    Req = 0.0012
    OCV, Smin, Smax = get_OCV_func_LiPo()
    tLR = @elapsed onebyVL = get_one_by_VspGLM(OCV, 100:-5:Smin, Req, maximum(Power))
    Cc, Aa, Bb = coef(onebyVL)
        ####### pass third argument of heuristics as locs
    
    N = length(Alist)
    # println("Solving with $(N) Nodes || B0 = $(Bstart)")
    graph = make_graph(def)
    Fvec = fill(NaN, N)
    heur_astar = get_heur_astar(E, locs, Fvec)
    heur_astar(S)
    if heur == "astar"
        heur_label! = get_heur_label(Fvec, graph, C, E, heur_astar)
    elseif heur == "euc"
        heur_label! = get_heur_label_euc(Fvec, locs, E)
    elseif heur == "manhattan"

    else
        println("invalid heuristic... breaking")
        return tLR
    end
    if heur_label!(S) == Inf     
        println("Infinite Start Heuristic...")
        return tLR
    end
    
    Q = MutableBinaryMinHeap([   (heur_label!(S) + 0, [0.0, Bstart, Qstart, S, 1.0, heur_label!(S) + 0]) ] )
    P = [zeros(0,6) for k = 1:size(C,1)] #set of treated labels | new labels compare to this AND Q
    X = Vector{Int}[ [S] ]  #hold path info
    Y = Vector{Int}[  []   ] #hold gen info
    
    z=0
    look_stack = Int[]
    dist_stack = Float32[]
    f_stack = Float32[]
    prog = ProgressUnknown("Working hard...", spinner=true)
    while true #loop until get to end node, or Q is empty
        isempty(Q) && (printstyled("Q empty, Z  = $z... \n", color=:light_cyan); break)

        #pull minimum cost label....
        next = pop!(Q)
        label_treated = next[2]
        i = Int(label_treated[4])
        K = Int(label_treated[5])
        append!(look_stack, i)
        append!(dist_stack, heur_label!(i))
        append!(f_stack, next[1])
        # println(f_stack[end])
        #now add it to P
        P[i] = vcat(P[i], reshape(label_treated,(1,6)))
        # println(i)
        #if we are at the end, then stop algo and return this label
        if i == E
            opt_cost =  label_treated[1]
            opt_path =  X[K]
            opt_gen =   Y[K]
            return tLR
        end
        for j in Alist[i]
            j==i && continue
            j∈X[K] && continue
            gc = label_treated[1] + C[i,j]
            h =  heur_label!(j) 
            if model == "linear"
                Δ = Power[i,j]*(Aa*Power[i,j] + Bb*label_treated[2] +  Cc)*t[i,j]
            elseif model == "naive"
                Δ = Power[i,j]*t[i,j]/15
            else
                error("Not a valid battery model")
            end
            label = []
            if label_treated[2] - Δ >= 0 
                label =  [gc, label_treated[2]-Δ , label_treated[3],  j, Int(length(X)+1), gc+h]
                    # println("label, $(i)->$(j) GEN")

                if EFF_heap(Q, label) && EFF_P(P, label)
                    push!(Q, (gc+h,label))

                    push!(X, [X[K]; j ])
                    push!(Y, [Y[K]; 1 ])
                end
            end
        end
        z+=1
        z == 200_000 && (printstyled("Z BREAK... @Z=$(z)\n", color=:light_red); break)
        z%500 == 0 && ProgressMeter.next!(prog)

    end
    return tLR
end


############################################################################
@load "Problems/euc_probs2D/5_1_RCSPP" euc_inst
t1 = MILP_RCSPP_naive(euc_inst)
t2 = MILP_RCSPP_linear(euc_inst)
tdp1 = @elapsed RCSPP_label(euc_inst, "linear")
tdp1 = @elapsed RCSPP_label(euc_inst, "naive")


#Now loop...
for N in Nvec
    for k in kvec
        println((N,k))
        @load "Problems\\euc_probs2D\\$(N)_$(k)_RCSPP" euc_inst
        feas_check(euc_inst) == 0 && (tnaive = missing; tlinear = missing)

        try
            a = 1
            # tnaive  = MILP_RCSPP_naive(euc_inst)
            # tlinear = MILP_RCSPP_linear(euc_inst)
        catch e
            println(e)
            tnaive  = missing
            tlinear = missing
        end
        tdp_linear = @elapsed tdpLR_linear = RCSPP_label(euc_inst, "linear")
        tdp_linear += tdpLR_linear
        tdp_naive  = @elapsed tdpLR_naive  = RCSPP_label(euc_inst, "naive")
        tdp_naive += tdpLR_naive

        # @save "Solutions\\battery_compare_MILP\\$(N)_$(k)_naive" tnaive
        # @save "Solutions\\battery_compare_MILP\\$(N)_$(k)_linear" tlinear

        @save "Solutions\\battery_compare_MILP\\$(N)_$(k)_naiveDP" tdp_naive
        @save "Solutions\\battery_compare_MILP\\$(N)_$(k)_linearDP" tdp_linear

    end
end

## Load data into arrays of time to solve... then plot
Mlinear = zeros(Union{Missing, Float64}, length(Nvec), 30)
Mnaive =  zeros(Union{Missing, Float64}, length(Nvec), 30)
MlinearDP = zeros(Union{Missing, Float64}, length(Nvec), 30)
MnaiveDP =  zeros(Union{Missing, Float64}, length(Nvec), 30)

for N in Nvec
    for k in kvec
        @load "Solutions/battery_compare_MILP/$(N)_$(k)_naive" tnaive     
        @load "Solutions/battery_compare_MILP/$(N)_$(k)_linear" tlinear     
        @load "Solutions\\battery_compare_MILP\\$(N)_$(k)_naiveDP" tdp_naive
        @load "Solutions\\battery_compare_MILP\\$(N)_$(k)_linearDP" tdp_linear
        
        Mnaive[N-4,k] = tnaive
        Mlinear[N-4,k] = tlinear
        MnaiveDP[N-4,k] = tdp_naive
        MlinearDP[N-4,k] = tdp_linear

    end
end

mean_linear = [ mean(skipmissing(Mlinear[N-4,:])) for N in Nvec  ]
mean_naive = [ mean(skipmissing(Mnaive[N-4,:])) for N in Nvec  ]

mean_linearDP = [ mean(skipmissing(MlinearDP[N-4,:])) for N in Nvec  ]
mean_naiveDP = [ mean(skipmissing(MnaiveDP[N-4,:])) for N in Nvec  ]


#  plot of time to solve vs problem size

plt = plot(layer(x=Nvec, y=mean_linear, Geom.line, Theme(default_color = "red")),
        layer(x=Nvec, y=mean_naive, Geom.line, Theme(default_color="grey")),
        Guide.xlabel("Number of Nodes"),
        Guide.ylabel("Time to Solve (s)"),
        Guide.manual_color_key("",["Linear Model", "Constant ΔSOC"],["red", "grey"]))

display(plt)
        
pltDP = plot(layer(x=Nvec, y=mean_linearDP, Geom.line, Theme(default_color = "red")),
layer(x=Nvec, y=mean_naiveDP, Geom.line, Theme(default_color="grey")),
Guide.xlabel("Number of Nodes"),
Guide.ylabel("Time to Solve (s)"),
Guide.manual_color_key("",["Linear Model", "Constant ΔSOC"],["red", "grey"]))
display(pltDP)

##
N=10

P = rand((150:300), N,N)
t = rand((20:30), N,N)
Cmax = 5000*3.6
Req = 0.0012
OCV, Smin, Smax = get_OCV_func_LiPo()
onebyVL = get_one_by_VspGLM(OCV, 100:-5:Smin, Req, maximum(P))
Cc, Aa, Bb = coef(onebyVL)



## 
P = 0:5:100
Svec = 0:5:100
Vdata = rand(length(P), length(Svec))
X = repeat(Svec, outer = length(P))
Y = repeat(P, outer=length(Svec))
Z = vec(Vdata)
data = DataFrame`(X=X, Y=Y, Z =Z)
