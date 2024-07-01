module MakeProblems
import Random.seed!
using SparseArrays
using JLD2
using HybridUAVPlanning
import Base.Threads.@threads
import Graphs
#TO DO 
#why are lattice instances taking so much longer than euc instances?
# started this 1/10/23
# Things to check:
#    is the Astar bound accounting for diagonal movemnt?
#    are we doing manhattan/euc accidently with lattice?  
#    were we using the old function loaded in utils.jl?  renamed to OLD 
#        -> check this!!!

get_single_idx3D = HybridUAVPlanning.get_single_idx3D
get_3D_idx = HybridUAVPlanning.get_3D_idx
function make_ACF3D(Dim; init_graph=nothing)
    xx = Dim[1] #width 
    yy = Dim[2] #height?
    zz = Dim[3]
    N = prod(Dim)
    Alist = [[] for k in 1:N]
    A = zeros(N, N)
    C = zeros(N, N)
    F = zeros(N, N)
    for i = 1:N
        if (i - 1) % xx != 0 #if its not left side add neighbor to left
            if init_graph != nothing && !Graphs.has_edge(init_graph, i, i-1)
                #do nothing
            else #if init_graph is nothing OR there is this edge in init_graph, add it!
                push!(Alist[i], i - 1)
                A[i, i-1] = 1
                A[i-1, i] = 1
                C[i, i-1] = 1
                C[i-1, i] = 1
            end
            #left and up
            if i <= (zz - 1) * xx * yy
                move = xx * yy - 1
                push!(Alist[i], i + move)
                A[i, i+move] = 1
                A[i+move, i] = 1
                C[i, i+move] = sqrt(2)
                C[i+move, i] = sqrt(2)
            end
            #left and down
            if i > xx * yy
                move = -xx * yy - 1
                push!(Alist[i], i + move)
                A[i, i+move] = 1
                A[i+move, i] = 1
                C[i, i+move] = sqrt(2)
                C[i+move, i] = sqrt(2)
                F[i, i+move] = 1
            end
        end
        if i % xx != 0 #if its not right side 
            if init_graph != nothing && !Graphs.has_edge(init_graph, i, i + 1)
                #do nothing
            else
                push!(Alist[i], i + 1)
                A[i, i+1] = 1
                A[i+1, i] = 1
                C[i, i+1] = 1
                C[i+1, i] = 1
            end
            #right and up
            if i <= (zz - 1) * xx * yy
                move = xx * yy + 1
                push!(Alist[i], i + move)
                A[i, i+move] = 1
                A[i+move, i] = 1
                C[i, i+move] = sqrt(2)
                C[i+move, i] = sqrt(2)
            end
            #right and down
            if i > xx * yy
                move = -xx * yy + 1
                push!(Alist[i], i + move)
                A[i, i+move] = 1
                A[i+move, i] = 1
                C[i, i+move] = sqrt(2)
                C[i+move, i] = sqrt(2)
                F[i, i+move] = 1
            end
        end



        bottom_side = []
        for k in 1:zz
            a = xx * yy * (k - 1) + 1
            b = xx * yy * (k - 1) + xx
            append!(bottom_side, a:b)
        end
        if i ∉ bottom_side #if its not bottom ("south") side
            if init_graph != nothing && !Graphs.has_edge(init_graph, i, i-xx)
                #do nothing
            else
                push!(Alist[i], i - xx)
                A[i, i-xx] = 1
                A[i-xx, i] = 1
                C[i, i-xx] = 1
                C[i-xx, i] = 1
            end
            #south and up
            if i <= (zz - 1) * xx * yy
                move = xx * yy - xx
                push!(Alist[i], i + move)
                A[i, i+move] = 1
                A[i+move, i] = 1
                C[i, i+move] = sqrt(2)
                C[i+move, i] = sqrt(2)
            end
            #south and down
            if i > xx * yy
                move = -xx * yy - xx
                push!(Alist[i], i + move)
                A[i, i+move] = 1
                A[i+move, i] = 1
                C[i, i+move] = sqrt(2)
                C[i+move, i] = sqrt(2)
                F[i, i+move] = 1
            end

        end


        top_side = []
        for k in 1:zz
            a = xx * yy * k - xx + 1
            b = xx * yy * k
            append!(top_side, a:b)
        end
        if i ∉ top_side #if its not top ("north") side
            if init_graph != nothing && !Graphs.has_edge(init_graph, i, i+xx)
                #do nothing
            else
                push!(Alist[i], i + xx)
                A[i, i+xx] = 1
                A[i+xx, i] = 1
                C[i, i+xx] = 1
                C[i+xx, i] = 1
            end
            #north and up
            if i <= (zz - 1) * xx * yy
                move = xx * yy + xx
                push!(Alist[i], i + move)
                A[i, i+move] = 1
                A[i+move, i] = 1
                C[i, i+move] = sqrt(2)
                C[i+move, i] = sqrt(2)
            end
            #north and down
            if i > xx * yy
                move = -xx * yy + xx
                push!(Alist[i], i + move)
                A[i, i+move] = 1
                A[i+move, i] = 1
                C[i, i+move] = sqrt(2)
                C[i+move, i] = sqrt(2)
                F[i, i+move] = 1
            end
        end

    end
    #first multiply by 100, then round and Int
    #making digonals the same cost as straight lines (sqrt(2) rounded to 1)...
    # cause SLD to OVERestiamate, thus making it an unusable heuristic
    #multiply every entry in C by 100
    C = 100 .* C
    C = round.(C)
    C = Int.(C)
    return Alist, A, C, F
end



function make_G3D(Dim)  #replace this sith the EUC version to keep consistent....
    xx, yy, zz = Dim[1], Dim[2], Dim[3]
    N = prod(Dim)
    G = ones(Bool, N, N)
    for k in 1:(floor(xx / 5)) #place some small blocks...
        width = Int(rand(1:xx/5))
        length = Int(rand(1:yy/5))
        height = Int(rand(1:max(zz / 5, 1)))
        zz == 1 && (height = 0)
        point = [rand(1:(xx-width)), rand(1:(yy-length)), rand(1:(max(zz - height, 1)))]

        list = []
        [push!(list, [i, j, k]) for i in point[1]:point[1]+width, j in point[2]:point[2]+length, k in point[3]:point[3]+height]
        list_i = get_single_idx3D.(list, Ref(Dim))
        G[list_i, :] .= 0
    end

    for k in 1:(ceil(xx / 10)) #place some larger blocks...
        width = Int(rand(1:xx/3))
        length = Int(rand(1:yy/3))
        height = Int(rand(1:max(zz / 5, 1)))
        zz == 1 && (height = 0)


        point = [rand(1:(xx-width)), rand(1:(yy-length)), rand(1:(max(zz - height, 1)))]

        list = []
        [push!(list, [i, j, k]) for i in point[1]:point[1]+width, j in point[2]:point[2]+length, k in point[3]:point[3]+height]
        list_i = get_single_idx3D.(list, Ref(Dim))
        G[list_i, :] .= 0
    end
    return G
end

function get_locs_lattice(Dim)
    if Dim[3] == 0
        zz = 2
    else
        zz = 3
    end
    locs = zeros(prod(Dim), zz)
    for i in 1:size(locs, 1)
        if Dim[3] == 0
            locs[i, :] = get_2D_idx(i, Dim)
        else
            locs[i, :] = get_3D_idx(i, Dim)
        end
    end
    return locs
end

function make_instance3D(Dim, seedy; init_graph=nothing)
    seed!(seedy)
    xx, yy, zz = Dim[1], Dim[2], Dim[3]
    nn = Dim[1]
    N = prod(Dim)
    if Dim[end] == 1
        height = 2
    elseif Dim[end] > 1
        height = 3
    end
    Sijk = [1, 1, 1]                  #start
    Eijk = [xx, yy, zz]      #end
    S = get_single_idx3D(Sijk, Dim)
    E = get_single_idx3D(Eijk, Dim)

    Alist, A, C, F = make_ACF3D(Dim, init_graph=init_graph)
    zz = 2
    Z = 2 .* C ##
    locs = get_locs_lattice(Dim)
    #Make G...
    # G = make_G3D(Dim)
    locs_for_G = 100 * locs ./ maximum(locs)#scale locs to 100x100x100... 
    Gf = get_G_grid(locs_for_G, Alist, Dim=height)
    # GFlipped = .!G
    # GFlipped = sparse(GFlipped)

    # if sum(abs.(Sijk - Eijk)) <= 10
    #     drop = 2
    # elseif sum(abs.(Sijk - Eijk)) <= 25
    #     drop = 4
    # elseif sum(abs.(Sijk - Eijk)) <= 50
    #     drop = 8
    # elseif sum(abs.(Sijk - Eijk)) <= 100
    #     drop = 15
    # else
    #     drop = 25
    # end
    B0 = 5 #sum(abs.(Sijk - Eijk)) - drop #Manhattan distance minus a few...
    nn = 100
    Q0 = nn * (zz) ## how many edges can we recharge along...
    anchor_list = sparse(zeros(N))

    Bmax = B0
    SC = 1
    if minimum(C) == 1 #need SC to be integer, but less than minimum of C...
        C *= 2
        Z *= 2
        Bmax *= 2
        Q0 *= 2
        B0 *= 2
    end

    # prob = GridProb(S, E, GFlipped, B0, G0, Bmax, SC, anchor_list, Dim)
    prob = EucGraphInt(S, E, Alist, F, C, Gf, Z, B0, Q0, Bmax, SC, anchor_list, locs)

    return prob
end


import Statistics.norm
import Random.seed!
function make_instance_euc(N, seedy; conn=4, Dim=3)
    seed!(seedy)
    locs = 100 * rand(N, Dim)
    C = [norm(locs[i, :] - locs[j, :]) for i = 1:size(locs, 1), j = 1:size(locs, 1)]

    Cmax, SE = findmax(C)
    S = SE[1]
    E = SE[2]

    Alist, C, F = get_ACF_euc(locs, conn=conn, Dim=Dim)
    Gf = get_G_grid(locs, Alist, Dim=Dim)

    Z = 2 * C
    anchor_list = sparse(zeros(N))
    Q0 = N
    B0 = Int(round(norm(locs[S, :] - locs[E, :]) / 8))
    Bmax = B0

    SC = 1
    Cnon0 = C[C.>0]
    Cmin = minimum(Cnon0)
    if Cmin == 1
        C *= 2
        Z *= 2
        Bmax *= 2
        Q0 *= 2
        B0 *= 2
    end
    euc_inst = EucGraphInt(S, E, Alist, F, C, Gf, Z, B0, Q0, Bmax, SC, anchor_list, locs)
    return euc_inst
end

function get_G_euc(locs, Alist; Dim=3, Nb=15)
    #take a few "blocks" and make all nodes inside them noise restricted.
    #then all edges off each of those nodes are noise restricted....
    N = size(locs, 1)
    GFlipped = spzeros(Bool, N, N) #0 if can use Gen, 1 otherwise....
    loc_copy = copy(locs)
    for k = 1:Nb
        xyz = 100 .* [rand(), rand(), rand()]
        LWH = 100 .* [0.2rand((0.8:0.01:1)), 0.2rand((0.8:0.01:1)), 0.2rand((0.8:0.01:1))]  #was 100 * .2() .2() .2()....
        if Dim == 2
            LWH[end] = 1
            xyz[3] = 0.0
            if size(locs, 2) == 3
                loc_copy[:, 3] .= 0.5
            end
        end

        nodes_in = in_box(loc_copy, xyz, LWH)

        #got nodes in the box, now restrict edges in/out for each nodes
        if isempty(nodes_in)
            continue
        end
        for i in nodes_in
            Ai = Alist[i]
            for j in Ai
                GFlipped[i, j] = 1
                GFlipped[j, i] = 1

            end
        end
    end
    Gsparse = sparse(GFlipped)
    return Gsparse
end

function get_G_grid(locs, Alist; Dim=3, p_noise=0.15) #THIS NEEDS TO BE LOSER FOR THE MAZES!!! Less NR?
    #take a few "blocks" and make all nodes inside them noise restricted.
    #then all edges off each of those nodes are noise restricted....
    N = size(locs, 1)

    #get grid_size such that #nodes per grid cell is about 3
    grid_size = Int(round((N / 3)^(1 / Dim)))

    GFlipped = zeros(Bool, N, N) #0 if can use Gen, 1 otherwise....
    if size(locs, 2) == 3 #2D lattice has size==3, all z coords are 1
        loc_copy = copy(locs)
    elseif size(locs, 2) == 2 #only happens for 2D euc
        loc_copy = ones(size(locs, 1), 3)
        loc_copy[:, 1:2] = locs
    end
    #make grid (boxes) and determine if they are noise restricted (p_noise chance)

    if Dim == 2
        xyz_array = zeros(grid_size, grid_size, 1, 3)
        stepx, stepy = maximum(locs[:, 1]) / grid_size, maximum(locs[:, 2]) / grid_size
        stepz = 10
        loc_copy[:, 3] .= 1.0
        bool_array = rand(grid_size, grid_size, 1) .< p_noise
    elseif Dim == 3
        xyz_array = zeros(grid_size, grid_size, grid_size, 3)
        stepx, stepy, stepz = maximum(locs[:, 1]) / grid_size, maximum(locs[:, 2]) / grid_size, maximum(locs[:, 3]) / grid_size
        bool_array = rand(grid_size, grid_size, grid_size) .< p_noise
    end

    for i in 1:size(xyz_array, 1)
        for j in 1:size(xyz_array, 2)
            for k in 1:size(xyz_array, 3)
                xyz_array[i, j, k, :] = [(i - 1) * stepx, (j - 1) * stepy, (k - 1) * stepz]
            end
        end
    end



    for i in 1:size(xyz_array, 1), j in 1:size(xyz_array, 2), k in 1:size(xyz_array, 3)
        xyz = xyz_array[i, j, k, :]
        ijk_bool = bool_array[i, j, k]
        if rand() < p_noise #ijk_bool 
            nodes_in = in_box(loc_copy, xyz, [stepx, stepy, stepz])
        else
            nodes_in = []
        end

        #got nodes in the box, now restrict edges in/out for each nodes
        if isempty(nodes_in)
            continue
        end
        for i in nodes_in
            Ai = Alist[i]
            for j in Ai
                GFlipped[i, j] = 1
                GFlipped[j, i] = 1

            end
        end
    end
    Gsparse = sparse(GFlipped)
    return Gsparse
end

function get_ACF_euc(locs; conn=4, Dim=3)
    N = size(locs, 1)
    F = spzeros(Bool, N, N)
    # Cfull = [norm(locs[i,:] - locs[j,:]) for i = 1:N, j = 1:N]
    Csparse = spzeros(N, N)
    Alist = [[] for k in 1:N]
    for i in 1:size(locs, 1)  #for each node, make 4 connections
        Ci = [norm(locs[i, :] - locs[j, :]) for j = 1:N]
        closest = partialsortperm(Ci, 1:conn+1)
        closest = setdiff(closest, i)
        for j in closest
            j ∉ Alist[i] && (push!(Alist[i], j))
            i ∉ Alist[j] && (push!(Alist[j], i))
            Csparse[i, j] = Int(round(norm(locs[i, :] - locs[j, :])))
            Csparse[j, i] = Int(round(norm(locs[i, :] - locs[j, :])))
            (Dim > 2 && locs[i, 3] > locs[j, 3]) && (F[i, j] = 1)
        end
    end
    return Alist, Csparse, F
end

function in_box(locs, xyz, lwh)
    x, y = xyz[1], xyz[2]
    if size(locs, 2) == 2
        z = 0
    else
        z = xyz[3]
    end
    l, w, h = lwh[1], lwh[2], lwh[3]
    L = x
    R = x + l
    F = y
    Bc = y + w
    T = z + h
    Bt = z
    nodes = []
    for i in 1:size(locs, 1)
        xi, yi = locs[i, 1], locs[i, 2]
        if size(locs, 2) == 2
            zi = 0.001
        elseif size(locs, 2) == 3
            zi = locs[i, 3]
        end
        # xi, yi, zi = locs[i,1], locs[i,2], 0.001

        if xi <= R && xi >= L && yi <= Bc && yi >= F && zi <= T && zi >= Bt
            # if xi < R && xi > L && yi < Bc && yi > F 
            push!(nodes, i)
        end
    end
    return nodes
end

end