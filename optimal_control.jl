####################################################
## load and solve 1.... to test...
####################################################
OCP, path, gen = MILP_to_opt_ctrl(10,1, prob = "lattice_probs_2D")
N = 101
@time OCP_solNL = solve_gen_optimal_control(OCP, path, gen, N, linear = false)
OCP_sol = OCP_solNL
@save "Solutions\\OCP\\15_1_to_$(N)_NL" OCP_sol

####################################################
## load one lattice_2D problem, and solve with varying number of time points....
####################################################
using Gadfly, JLD2, HybridUAVPlanning
function solve_one_vartime(;NMILP = 15)
    prob = "lattice_probs_2D"
    OCP, path, gen = MILP_to_opt_ctrl(NMILP,1, prob = prob)
    printstyled("Solving $(prob) $(NMILP)_1 \n", color=:light_green)
    # for i = [20:100; 101:25:500]
    for i = 50:50:10_000
        print("$(i) Time Points:  ")
        
        OCP_sol = solve_gen_optimal_control(OCP, path, gen, i, linear=false)
        @save "Solutions\\OCP\\$(NMILP)_1_to_$(i)_NL"  OCP_sol
        
        print("$(OCP_sol.time_to_solve) (s) \n")
    end
end

solve_one_vartime(NMILP=10)
    


####################################################
## Plot Trajectory with Generator
####################################################
# @load get_tag(550,2, with_path = true) OCP_sol

#need to get normal solution plot, but without a path and gen....
using Interpolations, Gadfly, HybridUAVPlanning, JLD2
@load "Problems\\lattice_probs_2D\\15_1" lattice_inst
@load "Solutions\\OCP\\15_1_to_50_NL" OCP_sol
locs, path, gen = OCP_sol.locs, OCP_sol.path, OCP_sol.gen

gen[gen.<=0] .= 0

time_MILP = LinRange(0,1, length(path))
interpx, interpy = linear_interpolation(time_MILP, locs[path,1]), linear_interpolation(time_MILP, locs[path,2])
x_OCP, y_OCP = interpx.(OCP_sol.timevec), interpy.(OCP_sol.timevec)

plt_graph = plot_euc_graph(lattice_inst)

plt_traj = plot(
                layer(x = x_OCP[1:end-1], y = -y_OCP[1:end-1], xend = x_OCP[2:end], yend = -y_OCP[2:end] , Geom.segment,color = gen, Theme(line_width = 1mm)),
                # layer(x = locs[path[[1,length(path)]],1], y = locs[path[[1,length(path)]],2]), Geom.point, Theme(default_color = "red"), shape=[Shape.xcross],
                Guide.xlabel("X Position"),
                Guide.ylabel("Y Position"),
                Guide.xticks(label=false),
                Guide.yticks(label=false),
                Guide.colorkey(title="u(t)"),
                Theme(grid_line_width = 0pt),
                Scale.color_continuous(colormap = Scale.lab_gradient("grey", "magenta"), maxvalue = 100),
                Guide.manual_color_key("",["Noise Restricted"],["lightcoral"])
)
# draw(PDF("Figs/traj_plain.pdf", 14cm, 12cm), plt_traj)
# draw(PDF("Figs/graph_plain_MECC.pdf", 14cm, 14cm), plt_graph)

####################################################
## Plot generator over time against MILP
####################################################
locs, path, genOCP = OCP_sol.locs, OCP_sol.path, OCP_sol.gen
time_pointsOCP = OCP_sol.timevec
G_OCP = OCP_sol.noise_restrictions
battery_state = OCP_sol.batt_state

@load "Solutions\\lattice_probs_2D\\15_1_label_eucLB" tdp cost pathL gen
noiseR_along_path = Float64.([lattice_inst.GFlipped[pathL[idx-1],pathL[idx]] for idx in 2:length(pathL)])


genMILP = 100 .* gen
pushfirst!(genMILP,0)
time_pointsMILP = 0:1/(length(genMILP)-1):1
gen_repeat = repeat(genMILP,inner=2)
time_pointsMILPrepeat = repeat(time_pointsMILP, inner=2)
_ = popfirst!(gen_repeat)
_ = popfirst!(time_pointsMILPrepeat)
_ = pop!(time_pointsMILPrepeat)
_ = popfirst!(gen_repeat)


plt_gen = plot(layer(x= time_pointsOCP, y = genOCP, Geom.line, Theme(default_color = "black")), 
            # layer(x = time_pointsMILP[1:end-1], y = genMILP[2:end], xend = time_pointsMILP[2:end], yend = genMILP[2:end], Theme(default_color = "blue"), Geom.segment, Geom.line),
            layer(x = x = time_pointsMILPrepeat, y = gen_repeat, Geom.line, Theme(default_color = "blue")),
            layer(xmin=time_pointsOCP[1:end-1], xmax = time_pointsOCP[2:end].+ 0.0, Geom.vband, color = G_OCP[1:end],), #max gen as ribbons....
            Scale.color_discrete_manual("white", "lightcoral"),
            Guide.colorkey(title="", labels = [""]),
            Guide.manual_color_key("",["Noise Resricted   ","MILP gen   ", "OCP gen"],["lightcoral","blue","black"]),
            Guide.xlabel("Time (normalized)"),
            Guide.ylabel("u(t) - Gen State"),
            Theme(grid_line_width = 0pt),
            Coord.cartesian(xmin = 0, xmax = 1),
            Theme(key_position=:top,
            )
)
set_default_plot_size(15cm, 10cm)
display(plt_gen)
# draw(PDF("Figs/gen_and_noiseR.pdf", 15cm, 10cm), plt_gen)

####################################################
## battery and fuel plot stacked
####################################################
time_points = OCP_sol.timevec
battery_state  = OCP_sol.batt_state
fuel_state = OCP_sol.fuel_state              
fuel_perc = fuel_state
# f1, ff = fuel_state[1], fuel_state[end]

# fuel_perc = [(fuel_state[i] - 0.99ff)/(f1-0.99ff)*100 for i = 1:length(fuel_state)]

plt_battery = plot(x = time_points, y = battery_state, Geom.line, Guide.xlabel(""), Guide.ylabel("Battery SOC (%)"), Coord.cartesian(ymin = 0, ymax = 100), Theme(default_color = "green"))

plt_fuel = plot(x = time_points, y = fuel_perc, Geom.line, Guide.xlabel("Time (Normalized)"), Guide.ylabel("Fuel Level (%)"), Coord.cartesian(ymin = 0, ymax = 100), Theme(default_color = "sienna"))

plt_double = plot(
layer(x = time_points[2:end], y = battery_state[2:end], Geom.line, Theme(default_color = "green")),
layer(x = time_points, y = fuel_perc, Geom.line, Theme(default_color = "sienna")),
Guide.manual_color_key("",["Fuel","Battery"],["sienna","green"]),
Guide.xlabel("Time (normalized)"), Guide.ylabel("Energy Level (% of Max)"), Coord.cartesian(ymin = 0, ymax = 100),
Theme(key_position = :top))

                ##
stacked = vstack(plt_battery, plt_fuel)
plt_batt_fuel = title(stacked, "")

set_default_plot_size(15cm, 10cm)
display(plt_double)
draw(PDF("Figs/batt_fuel_opt_ctrl.pdf", 15cm, 10cm), plt_batt_fuel) 
draw(PDF("Figs/double_batt_fuel.pdf", 15cm, 10cm), plt_double) 



####################################################
## J* as numer of time points increases.....
####################################################
Jstar = Float64[]
time_points = Float64[]

# for i in [20:100; 101:25:500]
for i in 20:100
    @load  "Solutions\\OCP\\10_1_to_$(i)_NL"  OCP_sol
    if OCP_sol != 0
        fuel_state = OCP_sol.fuel_state
        f1, ff = fuel_state[1], fuel_state[end]
        fuel_perc = [(fuel_state[i] - 0.9ff)/(f1-0.9ff)*100 for i = 1:length(fuel_state)]
        push!(time_points, length(OCP_sol.timevec))
        push!(Jstar, fuel_perc[1] - fuel_perc[end])
    end
end

plt_Jstar = plot(x = time_points, y = Jstar .+ 39.07215678, Geom.line,
            Guide.xlabel("Number of time points"),
            Guide.ylabel("J* (Final Fuel Level)"),
        Scale.x_continuous(format =   :plain),

            )
draw(PDF("Figs/J_star.pdf", 15cm, 7cm), plt_Jstar) 
set_default_plot_size(15cm, 7cm)

####################################################
## Time-to-solve vs number of time points....
using Gadfly
using LsqFit

Nvec = 50:50:10_000

lengths = Int64[]
timesNL = Float64[]
for i in Nvec 
    for k in 1:10
        @load  "Solutions\\OCP\\15_1_to_$(i)_NL"  OCP_sol
        if OCP_sol != 0
            push!(lengths, length(OCP_sol.timevec))
            push!(timesNL, OCP_sol.time_to_solve)
        end
    end
end
# @. linearfit(x, p) = p[1]*x^2 + p[2]*x + p[3]
@. nonlinearfit(x,p) = p[1]*x
p0 = [0.005]
# lin_line = curve_fit(linearfit, lengths, timesL, p0)
nonlin_line = curve_fit(nonlinearfit, lengths, timesNL, p0)

# trend_linear(x) = linearfit(x, lin_line.param)
trend_nonlinear(x) = nonlinearfit(x, nonlin_line.param)
# trend_nonlinear(x) = 0.005x
plt_time_solve = plot(
    layer(x = lengths, y = trend_nonlinear.(lengths), Geom.line, Theme(default_color = "grey")),
    layer(x = lengths, y = timesNL, Geom.point, Theme(default_color = "red", grid_line_width = 0mm,highlight_width=0.5pt)),
    Guide.manual_color_key("",["Raw Data","y = 0.004x"],["red","grey"]),
    Guide.xlabel("Number of Time Points"),
    Guide.ylabel("Time to Solve OCP (s)"),
    Scale.x_continuous(format = :plain),
    Theme(key_position=:top),
    Scale.y_log10
    )

draw(PDF("Figs/time_to_solve.pdf", 15cm, 10cm), plt_time_solve )
set_default_plot_size(15cm, 10cm)
display(plt_time_solve) 


##############################################################################
## get plot of average path length vs graph size....
Nvec = [50:500:2000; 2000:1000:20000]
Ls, avg_Ls = get_pathlength_vec(Nvec, "euc_probs2D", conn = "", algo = "_node")

plt = plot(layer(x = Nvec, y = avg_Ls, Geom.line, Geom.point),
            layer(x = Nvec, y = 30log10.(Nvec), Geom.line))

###############################################################################
## Scaling up the problem.....
import Statistics.mean
nvec = 10:50:2000
u_star, t_star = gen_optimal_control(nvec[1]);
tvec = Float64[]
for n in nvec
    println("N = $(n)")
    tN = Float64[]
    for k = 1:10
        t = @elapsed u_star, t_star = gen_optimal_control(n);
        push!(tN, t)
    end
    push!(tvec, mean(tN))
end
using Gadfly
plt = Gadfly.plot(x = nvec, y = tvec, Geom.point, Geom.line, Theme(default_color = "red"), Guide.xlabel.("Number of Time Points"), Guide.ylabel("Time to Solve (s)"), Scale.y_log10)
draw(PDF("Figs/opt_control_scaling.pdf", 12cm, 10cm), plt)



## solve all the euc2D label sols aspect_ratio
function solve_OCPs()
    Nvec = [50:500:2000; 2000:1000:20000]
    K = 10
    OCP, path, gen = MILP_to_opt_ctrl(50,1)
    OCP_sol = solve_gen_optimal_control(OCP, path, gen, length(path), linear=false)
    for N in Nvec
        # N < 3000 && continue
        for k in 1:K
            println("($(N),$(k)) Solving...")
            OCP, path, gen = MILP_to_opt_ctrl(N,k)
            if path == [0] #if there was no MILP sol found...
                OCP_sol = 0
                @save "Solutions\\OCP\\$(OCP.tag)_NL" OCP_sol
            else
                OCP_sol = solve_gen_optimal_control(OCP, path, gen, length(path), linear=false)
                @save "Solutions\\OCP\\$(OCP_sol.tag)_NL"  OCP_sol
            end
        end
    end
end
solve_OCPs()


##
using Gadfly
plt = plot(layer(x = 1:10, y = rand(10)), Guide.manual_color_key("",["Nominal Voltage"],["red"]))

#################################################
##  --- Verification Step --- #
# Once the optimization is complete, the job isn't really over. 
# The next step is to code up an ODE solver to actually propagate 
# the equations of motion forward through time form the initial 
# conditions and using the controls found from the optimization. 
# Once propoagated, the optimal solution of the states should 
# match those found from the nonlinear progr am solver. Also, 
# if a transformation was used from one state space to another, 
# the ODE solver can provide the original state space trajectories. 
# In this case this is performed.



# --- propagate the equations of motion --- #
X0 = [x0;y0];   # Initial State
tspan = (0.0,tf)                               # time vector (min,max)
tvec = LinRange(0, tf, 2*n)
prob = DifferentialEquations.ODEProblem(ezDynamics!,X0,tspan,p)      # define the ode problem to be solved
sol = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5())    # solve the ODE
# --- store the solution to the variables --- #
x_ode = sol(tvec)[1,:]
y_ode =  sol(tvec)[2,:]

# --- Verification Plot for the Rda and Rat States --- #
# Plots.plot(sol,vars=(1))
# Plots.plot!(x_ode,y_ode,fmt=:png)
# Plots.plot(value.(x_ode)[:],value.(y_ode)[:])
# Plots.plot!(LinRange(0,n*value.(dt),n),value.(Rda)[:])
# Plots.plot!(LinRange(0,n*value.(dt),n),value.(Rat)[:])
# Plots.plot()
# Plots.plot!(x_star,y_star, seriestype = :scatter, label = "NLP", lw = 1)
# Plots.plot!(x_ode,y_ode, label = "ODE", lw = 2, legend = :outertopleft, aspect_ratio = 1)

# Plots.plot!(x_star,y_star)
# Stuff all of the collocation results into a CollocationResults struct for easy plotting. The plotting of the type CollocationResult is defined in a plot recipe in the OptimalWEZAvoidance module.
τ = LinRange(0, tf, n)
iWEZ = value.(pmax) .- value.(R)
cr = CollocationResult(hcat(x_star, y_star), u_star, tf, τ, iWEZ, X0, [xf,yf], default_params)
# Plots.plot([cr], ["B?"]; LIP=false)
Plots.plot(cr; LIP=false, title="Scenario B, n = $n, tf = $(round(tf; sigdigits=4))")

nframes = 50
anim = Plots.@animate for t in LinRange(sol.t[1], sol.t[end], nframes)
    Plots.plot(sol, vars=(1,2), alpha = .3, legend = false)
    Plots.plot!(sol, vars=(1,2), tspan = (0.0, t))
    Plots.scatter!(sol, vars=(1,2), tspan = (t,t))
end

mygif = gif(anim, "filename.gif", fps = 15)

display(mygif)






##Make dummy solution..  do not need anymore since we are solving MILPS....
xy = hcat(1:10, 1:10)
gen = rand((0:0.01:1), size(xy,1))
time_points = LinRange(0,1, length(gen))
G = rand(Bool, length(gen))
gen[G .== 0] .= 0

battery_state = zeros(length(time_points))
gen_state = zeros(length(time_points))

battery_state[1] = 100
gen_state[1] = 100

#batt if no gen ends at -20
ΔB = 120/length(time_points)
for i in 2:length(battery_state)
    gen_state[i] = gen_state[i-1] - gen[i]*30
    battery_state[i] = battery_state[i-1] - ΔB + 3*ΔB*gen[i]
end

##
using Gadfly, JLD2, HybridUAVPlanning
OCP, path, gen = MILP_to_opt_ctrl(15,1, prob = "lattice_probs_2D")
C,Z,GFlipped = OCP.C, OCP.Z, OCP.GFlipped
N = 1000
Cvec = C_time_discretized(C, path, N)
Zvec = C_time_discretized(Z, path, N)
Pijmax = maximum(Zvec .- Cvec)*1.01
Pnormed(u,Zij,Cij) = 20*(Zij*u - Cij)/Pijmax
gen_split = HybridUAVPlanning.u_discretized(gen, N)


mdot_normed(uu,Zz) = (86uu^2 + 8.3*uu + 0.083)*Zz/94.383 #from http://dx.doi.org/10.1051/matecconf/201925206009
for i in 1:length(Zvec)
    Z_proc = mdot_normed(.80, Zvec[i])
    # println((Z_proc, Zvec[i]))
end

println("N: $(N)  ||  Zsum = $(sum(Zvec))")

##Test nonlinear constraints.... wrt to battery and gen updates....
using Gadfly, JLD2, HybridUAVPlanning

OCP, path, gen = MILP_to_opt_ctrl(15,1, prob = "lattice_probs_2D")
N = 100
Cvec = C_time_discretized(OCP.C, path, N)
Zvec = C_time_discretized(OCP.Z, path, N)

mdot_normed(uu,Zz) = (86uu^2 + 8.3*uu + 0.083)*Zz/94.383 #from http://dx.doi.org/10.1051/matecconf/201925206009
obV = get_one_by_Vsp()
Pijmax = maximum(Zvec .- Cvec)*1.01

Pnormed(u,Zij,Cij) = 20*(Zij*u - Cij)/Pijmax
Λ(b,Pij) = obV(b,abs(Pij))/obV(100,0)
u = 100*rand(N)
b = zeros(N)
g = zeros(N)
g[1], b[1] = 100,50
for timei = 2:N
    println(timei)
    Pij = Pnormed(u[timei]/100, Zvec[timei-1], Cvec[timei-1])
    g[timei] =g[timei-1] - u[timei]/100*Zvec[timei-1]*mdot_normed(u[timei]/100, Zvec[timei-1])
    b[timei] = b[timei-1] + (Zvec[timei-1]*u[timei]/100 - Cvec[timei-1])*Λ(b[timei-1],  Pij)*100/(2*OCP.Bmax)
    println(b[timei]," || ", g[timei])
end

plot(layer(x = 1:N, y=g, Geom.line, Theme(default_color = "brown")), 
    layer(x = 1:N, y=b, Geom.line, Theme(default_color = "blue")),
    layer(x = 1:N, y=u, Geom.line, Theme(default_color = "black"))
    )
2

@load "Solutions\\OCP\\15_1_to_50_NL" OCP_sol
OCP50 = OCP_sol
@load "Solutions\\OCP\\15_1_to_500_NL" OCP_sol
OCP500 = OCP_sol
plt = plot(layer(x = OCP50.timevec, y = OCP50.fuel_state, Theme(default_color = "red"), Geom.line),
            layer(x = OCP500.timevec, y = OCP500.fuel_state, Theme(default_color = "blue"), Geom.line),
            Guide.manual_color_key("", ["50","500"],[colorant"red", colorant"blue"])
)
