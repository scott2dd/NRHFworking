# Test battery....

# Sections:
# 1) Misc and Util functions
# 2) Pulse test function definitions
# 3) plot pulse test for 18650
# 4) plot pulse test for LiPo
# 5) Pull OCV data from CSVs
# 6) Data from pulse test
# 7) Old code for pulse_test1() plotting

using Gadfly
using Compose
using GLM
using Cairo 
include("utils.jl")
Gadfly.push_theme(:default)

#############################################
## 1) Misc and Utils... 
    

OCV, min, max = get_OCV_func()
onebyV = get_one_by_V()

SOC = 6:0.1:99
OCVv = OCV(SOC) 

plt1 = plot(x=SOC,y=OCVv, Geom.line, Coord.cartesian(xmin=0, xmax=100, ymin=0, ymax=4.5), Guide.xlabel("SOC"), Guide.ylabel("OCV (V)"))
plt2 = plot(x=SOC,y=1 ./OCVv, Geom.line, Guide.xlabel("SOC"), Guide.ylabel("1/OCV"))
# display(plt1)
# display(plt2)
#now test
tvec = [0.,2000.]
t0, tf = tvec[1], tvec[2]
N1 = 10
N2 = 50
N3 = 100
Δ1 =  (tf - t0)/N1
Δ2 =  (tf - t0)/N2
Δ3 =  (tf - t0)/N3

S0 = 85.
P = 500.
P = OCV(S0)*2.5 #Discahrge rate of 2.5 Amps
Cmax = 9000.   #1 cell 18650 capacity in couloumbs
Req = 70/1000
R1 = 0.5
R2 = 0.002
C1 = 7200
C2 = 3000 
 
#Propogate forward w/ runge-katte func
# S_vec1, t_vec1 = riemman(OCV, tvec, S0, P, Cmax, N = N1)
# S_vec2, t_vec2 = riemman(OCV, tvec, S0, P, Cmax, N = N2)
S_vec0, t_vec0 = riemman(OCV, tvec, S0, P, Cmax, N = N3)
S_vec2, t_vec2 = riemman2(V_model2, OCV, tvec, S0, P, Cmax, Req, N = N3)
S_vec7, t_vec7 = riemman7(OCV, tvec, S0, P, Cmax, Req, C2, N = N3)

#Now let's do the same for linear....
S_f = S0 - P*onebyV(S0)*(tvec[2]-tvec[1])/(Cmax)*100




pltSOC = plot(
            layer(x= tvec, y= [S0, S_f], Geom.line, Theme(default_color = "grey")),
            layer(x = t_vec0, y = S_vec0, Geom.line, Theme(default_color = "red")),
            layer(x = t_vec2, y = S_vec2, Geom.line, Theme(default_color = "green") ),
            layer(x = t_vec7, y = S_vec7, Geom.line, Theme(default_color = "Blue")),
            Guide.xlabel("Time (s)"),
            Guide.ylabel("SOC (%)"),
            Guide.manual_color_key("", ["Linear", "No Ohmic Drop","Simple Ohmic Drop", "2nd Order RC"], ["grey", "red", "green", "Blue"]),
            Guide.title("-- 18650 Cell @ $(round(Int, P))W"))
pltV = plot(x = t_vec2, y = OCV(S_vec2), Geom.line, Guide.xlabel("Time (s)"), Guide.ylabel("OCV (V)"), Coord.cartesian(xmin=0, xmax=tf, ymin=0, ymax=4.5))
plt_stack = vstack(compose(context(0,0,1.0, 0.5), render(pltSOC)),
                   compose(context(0,0,1.0, 0.5), render(pltV)) )

set_default_plot_size(20cm, 20cm)
display(plt_stack)


#############################################
## 2) Pulse Test definitions 
function pulse_test(Cmax, Req, C2)
    t1 = [0, 200]
    t2 = [200,400]
    t3 = [400,600]
    t4 = [600,800]
    t5 = [800,1000]
    t6 = [1000,1200]
    t7 = [1200,2000]
    tstack = [t1, t2, t3, t4, t5, t6, t7]

    P1 = 10
    P2 = 5
    P3 = 10
    P4 = 5
    P5 = 10
    P6 = 5
    P7 = 10
    N3 = 50
    S0L,S00, S02 ,S07 = 90,90,90,90
    Pstack = [P1,P2,P3,P4,P5,P6,P7]
    S_vecL, t_vecL = Float64[], Float64[] 
    S_vec0, t_vec0 = Float64[], Float64[] 
    S_vec2, t_vec2 = Float64[], Float64[] 
    S_vec7, t_vec7 = Float64[], Float64[] 
    VfL = get_one_by_V()
    OCV, min, max = get_OCV_func()
    for k = 1:length(tstack)
        P = Pstack[k]
        tvec = tstack[k]
        Sv0, t0 = riemman(OCV, tvec, S00, P, Cmax, N = 100)
        Sv2, t2 = riemman2(V_model2, OCV, tvec, S02, P, Cmax, Req, N = 100)
        Sv7, t7 = riemman7(OCV, tvec, S07, P, Cmax, Req, C2, N = 100)

        S_fL = S0L - P*VfL(S0L)*(tvec[2]-tvec[1])/(Cmax)*100
        SvL = [S0L, S_fL]
        tL = tvec[end]
        append!(S_vec0, Sv0); append!(t_vec0, t0)
        append!(S_vec2, Sv2); append!(t_vec2, t2)
        append!(S_vec7, Sv7); append!(t_vec7, t7)
        append!(S_vecL, SvL); append!(t_vecL, tvec)


        S00 = Sv0[end]
        S02 = Sv2[end]
        S07 = Sv7[end]
        S0L = SvL[end]
    end

    return S_vecL, t_vecL, S_vec0, t_vec0, S_vec2, t_vec2, S_vec7, t_vec7, Pstack
end

function pulse_test18650(Cmax, Req, C2)
    t1 = [0, 200]
    t2 = [200,400]
    t3 = [400,600]
    t4 = [600,800]
    t5 = [800,1000]
    t6 = [1000,1200]
    t7 = [1200,1400]
    t8 = [1400, 1600]
    t9 = [1600, 1800]
    tstack = [t1, t2, t3, t4, t5, t6, t7, t8, t9]

    P1 = 10
    P2 = 5
    P3 = 10
    P4 = 5
    P5 = 10
    P6 = 5
    P7 = 10
    P8 = 20
    P9 = 10
    OCV, min, max = get_OCV_func()
    S0L,S00, S02 ,S07 = 90,90,90,90
    Pstack = [P1,P2,P3,P4,P5,P6,P7,P8, P9]
    S_vecL, t_vecL = Float64[], Float64[] 
    S_vec0, t_vec0 = Float64[], Float64[] 
    S_vec2, t_vec2 = Float64[], Float64[] 
    S_vec7, t_vec7 = Float64[], Float64[]
    Cmax = 9000
    onebyVL = get_one_by_Vsp()
    
    for k = 1:length(tstack)
        P = Pstack[k]
        tvec = tstack[k]
        # Sv0, t0 = riemman(OCV, tvec, S00, P, Cmax, N = 100)
        Sv2, t2 = riemman2(V_model2, OCV, tvec, S02, P, Cmax, Req, N = 100)
        Sv7, t7 = riemman7(OCV, tvec, S07, P, Cmax, Req, C2, N = 100)
        
        S_fL = S0L - P*onebyVL(S0L,P)*(tvec[2]-tvec[1])/(Cmax)*100
        SvL = [S0L, S_fL]
        
        S_f0 = S0L - P*onebyVL(S0L,P)*(tvec[2]-tvec[1])/(Cmax*3.5)*100
        Sv0 = [S00, S_f0]

        tL = tvec[end]
        append!(S_vec0, Sv0); append!(t_vec0, tvec)
        append!(S_vec2, Sv2); append!(t_vec2, t2)
        append!(S_vec7, Sv7); append!(t_vec7, t7)
        append!(S_vecL, SvL); append!(t_vecL, tvec)



        S00 = Sv0[end]
        S02 = Sv2[end]
        S07 = Sv7[end]
        S0L = SvL[end]
    end

    return S_vecL, t_vecL, S_vec0, t_vec0, S_vec2, t_vec2, S_vec7, t_vec7, Pstack
end

function pulse_testLiPo()
    R1 = 0.5
    R2 = 0.002
    C1 = 7200
    C2 = 3000
    
    
    Cmax = 5000*3.6
    Req = 4*0.016  #   1 Healthy LiPo is  ~ 15 mΩ 
    t1 = [0, 120]
    t2 = [120,240]
    t3 = [240,360]
    t4 = [360,480]
    t5 = [480, 600]
    t6 = [600,720]
    t7 = [720,840]
    t8 = [840, 960]
    t9 = [960,1080]
    t10 = [1080, 1200]
    P1 = 20 #needs to be amps!
    P2 = 10
    P3 = 20
    P4 = 10
    P5 = 20
    P6 = 10
    P7 = 20
    P8 = 10
    P9 = 20
    P10 = 10

    tstack = [t1, t2, t3, t4, t5, t6,t7]#, t8, t9, t10]
    Pstack = 16 .* [P1,P2,P3,P4,P5,P6,P7]
    
    OCV, Smin, Smax = get_OCV_func_LiPo()
    S0L,S00, S02 ,S07 = 100,100,100,100
    S_vecL, t_vecL = Float64[], Float64[] 
    S_vec0, t_vec0 = Float64[], Float64[] 
    S_vec2, t_vec2 = Float64[], Float64[] 
    S_vec7, t_vec7 = Float64[], Float64[] 
    onebyVL = get_one_by_Vsp(OCV, S0L:-5:Smin, Req, maximum(Pstack))
    for k = 1:length(tstack)
        P = Pstack[k]
        tvec = tstack[k]
        
        # Sv0, t0 = riemman(OCV, tvec, S00, P, Cmax, N = 100)
        Sv2, t2 = riemman2(V_model2, OCV, tvec, S02, P, Cmax, Req, N = 100)
        Sv7, t7 = riemman7(OCV, tvec, S07, P, Cmax, Req, C2, N = 100)
        
        S_fL = S0L - P*onebyVL(S0L,P)*(tvec[2]-tvec[1])/(Cmax)*100
        SvL = [S0L, S_fL]

        S_f0 = S00 - P*(tvec[2]-tvec[1])/(Cmax*15.2)*100
        Sv0 = [S00, S_f0]

        println("S: $(Sv0[end])")
        tL = tvec[end]
        append!(S_vec0, Sv0); append!(t_vec0, tvec)
        append!(S_vec2, Sv2); append!(t_vec2, t2)
        append!(S_vec7, Sv7); append!(t_vec7, t7)
        append!(S_vecL, SvL); append!(t_vecL, tvec)



        S00 = Sv0[end]
        S02 = Sv2[end]
        S07 = Sv7[end]
        S0L = SvL[end]
    end

    return S_vecL, t_vecL, S_vec0, t_vec0, S_vec2, t_vec2, S_vec7, t_vec7, Pstack
end

function pulse_testLiPo_exp(tstack, Pstack)
    R1 = 0.5
    R2 = 0.002
    C1 = 7200
    C2 = 3000
    

    
    Cmax = 5000*3.6
    Req =4*0.016  #   1 Healthy LiPo is  ~ 15 mΩ 


    #need to get power from an average current...
    Pstack = Pmean
    OCV, Smin, Smax = get_OCV_func_LiPo()
    S0L,S00, S02 ,S07 = 100,100,100,100
    S_vecL, t_vecL = Float64[], Float64[] 
    S_vec0, t_vec0 = Float64[], Float64[] 
    S_vec2, t_vec2 = Float64[], Float64[] 
    S_vec7, t_vec7 = Float64[], Float64[] 
    onebyVL = get_one_by_Vsp(OCV, S0L:-5:Smin, Req, maximum(Pstack))
    for k = 1:length(tstack)
        P = Pstack[k]
        tvec = tstack[k]
        
        # Sv0, t0 = riemman(OCV, tvec, S00, P, Cmax, N = 100)

        Sv2, t2 = riemman2(V_model2, OCV, tvec, S02, P, Cmax, Req, N = 100)
        Sv7, t7 = riemman7(OCV, tvec, S07, P, Cmax, Req, C2, N = 100)
        S_fL = S0L - P*onebyVL(S0L,P)*(tvec[2]-tvec[1])/(Cmax)*100
        SvL = [S0L, S_fL]
        
        S_f0 = S00 - P*(tvec[2]-tvec[1])/(Cmax*16)*100
        Sv0 = [S00, S_f0]

        println("S: $(Sv0[end])")
        tL = tvec[end]
        append!(S_vec0, Sv0); append!(t_vec0, tvec) #t0
        append!(S_vec2, Sv2); append!(t_vec2, t2)
        append!(S_vec7, Sv7); append!(t_vec7, t7)
        append!(S_vecL, SvL); append!(t_vecL, tvec)


        S00 = Sv0[end]
        S02 = Sv2[end]
        S07 = Sv7[end]
        S0L = SvL[end]
    end

    return S_vecL, t_vecL, S_vec0, t_vec0, S_vec2, t_vec2, S_vec7, t_vec7, Pstack
end


#############################################
## 3) Pulse test 2 with 18650 
Cmax = 9000
Req = 70/1000
C2 = 3000 
S_vecL, t_vecL, S_vec0, t_vec0, S_vec2, t_vec2, S_vec7, t_vec7, Pstack  = pulse_test18650(Cmax, Req, C2)



set_default_plot_size(15cm, 20cm)
pltPulse = plot(
            layer(x= t_vecL, y= S_vecL, Geom.line, Theme(default_color = "grey")),
            layer(x = t_vec0, y = S_vec0, Geom.line, Theme(default_color = "red")),
            layer(x = t_vec2, y = S_vec2, Geom.line, Theme(default_color = "green")),
            layer(x = t_vec7, y = S_vec7, Geom.line, Theme(default_color = "Blue")),
            Guide.xlabel("Time (s)"),
            Guide.ylabel("SOC (%)"),
            Guide.manual_color_key("", ["Linear Model", "Nominal Voltage Only","Simple Ohmic Drop", "1st Order RC"], ["grey", "red", "green", "Blue"], pos = [0.1w, 0h]),
            Guide.title("18650 Numerical Test"))


Pplot= []
for k in 1:length(Pstack)
    append!(Pplot, [Pstack[k], Pstack[k]])
end
pltP = plot(x = t_vecL, y = Pplot, Geom.line, Guide.xlabel("Time (s)"), Guide.ylabel(" Power (Watts)"), Coord.cartesian(ymin=0))

plt_stack = vstack(compose(context(0,0,1.0, 0.5), render(pltPulse)),
                   compose(context(0,0,1.0, 0.5), render(pltP)) )

display(plt_stack)

#############################################
## 4) Pulse Test for LiPo 4S 
 

S_vecL, t_vecL, S_vec0, t_vec0, S_vec2, t_vec2, S_vec7, t_vec7, Pstack  = pulse_testLiPo()



pltPulse = plot(
            layer(x= t_vecL, y= S_vecL, Geom.line, Theme(default_color = "grey")),
            layer(x = t_vec0, y = S_vec0, Geom.line, Theme(default_color = "red")),
            layer(x = t_vec2, y = S_vec2, Geom.line, Theme(default_color = "green")),
            layer(x = t_vec7, y = S_vec7, Geom.line, Theme(default_color = "Blue")),
            Guide.xlabel("Time (s)"),
            Guide.ylabel("SOC (%)"),
            Guide.manual_color_key("", ["Linear Model", "Nominal Voltage Only","Simple Ohmic Drop", "1st Order RC"], ["grey", "red", "green", "Blue"]),
            Guide.title("Linear Model Performance -- LiPo Numerical Test"))


Pplot= []
for k in 1:length(Pstack)
    append!(Pplot, [Pstack[k], Pstack[k]])
end
pltP = plot(x = t_vecL, y = Pplot, Geom.line, Guide.xlabel("Time (s)"), Guide.ylabel(" Power (Watts)"), Coord.cartesian(ymin=0))



plt_stack = vstack(compose(context(0,0,1.0, 0.5), render(pltPulse)),
                   compose(context(0,0,1.0, 0.5), render(pltP)) )
set_default_plot_size(20cm, 40cm)
display(plt_stack)



#############################################
## 5) collect OCV curve for LiPo

using CSV
using DataFrames

df60_1 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_60_1.csv"))
df60_2 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_60_2.csv"))
df60_3 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_60_3.csv"))
df60_4 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_60_4.csv"))
df60_5 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_60_5.csv"))

df270_1 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_270_1.csv"))
df270_2 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_270_2.csv"))
df270_3 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_270_3.csv"))
df270_4 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_270_4.csv"))

df120_1 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_120_1.csv"))
df120_2 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_120_2.csv"))
df120_3 = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\OCV_120_3.csv"))

df_fin = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\final_rest.csv"))

t60_1 = df60_1."Time (s)"
t60_2 = df60_2."Time (s)"
t60_3 = df60_3."Time (s)"
t60_4 = df60_4."Time (s)"
t60_5 = df60_5."Time (s)"

t270_1 = df270_1."Time (s)"
t270_2 = df270_2."Time (s)"; t270_2[167]="154.5355175";  t270_2 = parse.(Float64, t270_2)
t270_3 = df270_3."Time (s)"; t270_3[126]="118.00917";  t270_3 = parse.(Float64, t270_3)
t270_4 = df270_4."Time (s)"

t120_1 = df120_1."Time (s)"
t120_2 = df120_2."Time (s)"
t120_3 = df120_3."Time (s)"

tf_1 = df60_1."Time (s)"



V60_1 = df60_1."Voltage (V)"
V60_2 = df60_2."Voltage (V)"
V60_3 = df60_3."Voltage (V)"
V60_4 = df60_4."Voltage (V)"
V60_5 = df60_5."Voltage (V)"

V270_1 = df270_1."Voltage (V)"
V270_2 = df270_2."Voltage (V)"
V270_3 = df270_3."Voltage (V)"
V270_4 = df270_4."Voltage (V)"

V120_1 = df120_1."Voltage (V)"
V120_2 = df120_2."Voltage (V)"
V120_3 = df120_3."Voltage (V)"

Vf_1 = df_fin."Voltage (V)"


i60_1 = df60_1."Current (A)"
i60_2 = df60_2."Current (A)"
i60_3 = df60_3."Current (A)"
i60_4 = df60_4."Current (A)"
i60_5 = df60_5."Current (A)"

i270_1 = df270_1."Current (A)"
i270_2 = df270_2."Current (A)"
i270_3 = df270_3."Current (A)"
i270_4 = df270_4."Current (A)"

i120_1 = df120_1."Current (A)"
i120_2 = df120_2."Current (A)"
i120_3 = df120_3."Current (A)"

if_1 = df_fin."Current (A)"


Vr = Float64[]
mah = Float64[]
Vri = mean(V60_1[1:5])
mahi = 5000 
push!(Vr, Vri)
push!(mah, mahi)

i_stack =[]
V_stack = []
t_stack = []
push!(t_stack, t60_1, t60_2, t60_3, t60_4, t60_5)
push!(t_stack, t270_1, t270_2, t270_3, t270_4)
push!(t_stack, t120_1, t120_2, t120_3, tf_1) 

push!(V_stack, V60_1, V60_2, V60_3, V60_4,V60_5)
push!(V_stack, V270_1, V270_2, V270_3, V270_4)
push!(V_stack, V120_1, V120_2, V120_3, Vf_1)

push!(i_stack, i60_1, i60_2, i60_3, i60_4, i60_5)
push!(i_stack, i270_1, i270_2, i270_3, i270_4)
push!(i_stack, i120_1, i120_2, i120_3, if_1)

for i in 1:(length(t_stack)-1)
    Ci = simps(t_stack[i], i_stack[i] .+ 76.3/1000)
    Vri = mean(V_stack[i+1][1:1])
    mahi = mah[end] - Ci*1/3.6
    push!(Vr, Vri )
    push!(mah, mahi)

end

SOC = mah./5000


# Plot OCV curves for modeling paper...
#use Plots.jl because can't do multi axis in Gadfly
using Plots; pyplot()
OCV18650, min, max = get_OCV_func()
OCVLiPo, min, max = get_OCV_func_LiPo()
Plots.plot(OCV18650.itp.knots[1], OCV18650.itp.coefs, ylabel="Open Circuit Voltage (V)", xlabel="SOC (%)", leg=:topleft,  label = "18650 (Left Axis)", right_margin=45Plots.mm)
plt = Plots.plot!(twinx(), OCVLiPo.itp.knots[1], OCVLiPo.itp.coefs,
    c=:red,
    ylabel="",
    label="4S LiPo (Right Axis)",
    leg=:bottomright,
    dpi =500, 
    size = (500,300 ))
Plots.plot!(right_bottom=10mm)
Plots.savefig(plt, "Figs\\OCV_plot.pdf")
display(plt)

#now print the OCV table in LaTeX form 

for i in 1:length(Vr)
    print(" $(round(100*SOC[i],digits=2)) & $(round(Vr[i],digits=4))  \\\\ \\hline \n")
end





######################################################################################
## 6) Pulse test expirmental data
using CSV

df1 =   DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\pulse_1.csv"))
df2 =   DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\pulse_2.csv"))
dfrel = DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\pulse_relaxed.csv"))
df3 =   DataFrame(CSV.File("C:\\Users\\Drew\\OneDrive - University of Cincinnati\\RC_benchmark\\pulse_3.csv"))

#df1 ended early due to vibration error
#df2 is first pulse test (discharged by 20|10|20|10|20|)
#df3 is second pulse test (discharged by 20|10|20|10|20|10)
#dfrel is relaxed voltage after end of second test

s3 = [8 150;
     152 296;
    297  439;
     441 580;
     582 713;
     715 854;
     855 988;
     990 1123]

i3 = df3."Current (A)"
t3 = df3."Time (s)"
V3 = df3."Voltage (V)"
P3 = df3."Electrical Power (W)"
imean= []
Pmean = []
for i in 1:size(s3,1)
    imean_i = mean(i3[s3[i,1]:s3[i,2]] .+ .+ 76.3/1000)
    Pmean_i = mean(P3[s3[i,1]:s3[i,2]])
    push!(imean, imean_i)
    push!(Pmean, Pmean_i)
end

@save "Problems\\pulse_data" imean s3 i3 t3 V3 P3

plot(x=t3, y = P3, Geom.line)
#first get cumulative battery SOC for experimental data....
mah = zeros(length(i3))
mah[1] = 5000
for i in 2:length(i3)
    Δ = simps(t3[1:i], i3[1:i])/3.6
    mah[i] = 5000 - Δ
end

tstack = t3[s3]
tstack = [ tstack[i,:] for i in 1:size(tstack,1) ]
#now run pulse test with experimental data, then compare against experimental data...
S_vecL, t_vecL, S_vec0, t_vec0, S_vec2, t_vec2, S_vec7, t_vec7, Pstack  = pulse_testLiPo_exp(tstack, Pmean)

SOC = 100*mah./5000

#now plot
set_default_plot_size(15cm, 20cm)

pltPulse = plot(
            layer(x= t_vecL, y= S_vecL, Geom.line, Theme(default_color = "grey")),
            layer(x = t_vec0, y = S_vec0, Geom.line, Theme(default_color = "red")),
            layer(x = t_vec2, y = S_vec2, Geom.line, Theme(default_color = "green")),
            layer(x = t_vec7, y = S_vec7, Geom.line, Theme(default_color = "blue")),
            layer(x = t3, y = SOC, Geom.line, Theme(default_color = "cyan")),
            Guide.xlabel("Time (s)"),
            Guide.ylabel("SOC (%)"),
            Guide.manual_color_key("", ["Linear Model", "Nominal Voltage Only","Simple Ohmic Drop", "1st Order RC","Actual SOC"], ["grey", "red", "green", "blue", "cyan"], pos = [0.05w, 0.25h]),
            Coord.cartesian(ymin = 0, ymax = 100),
            Guide.title("4S LiPo Experimental Test"))


Pplot= []
for k in 1:length(Pstack)
    append!(Pplot, [Pstack[k], Pstack[k]])
end
pltP = plot(layer(x = t_vecL, y = Pplot, Geom.line, Theme(default_color = "blue")),
            layer(x = t3, y = P3, Geom.line, Theme(default_color = "red")),
            Guide.manual_color_key("",["Average Pulse Power","Actual Power Draw"],["blue","red"], pos = [0.05w, 0.25h]),
            Guide.xlabel("Time (s)"), Guide.ylabel(" Power (Watts)"), Coord.cartesian(ymin=0)
        )



plt_stack = vstack(compose(context(0,0,1.0, 0.5), render(pltPulse)),
                   compose(context(0,0,1.0, 0.5), render(pltP)) )
display(plt_stack)
draw(PDF("Figs/sol_ex.pdf", 15cm, 20cm), plt_stack)  #if draw, then displaying after is broken?

##
using DataFrames, Gadfly, RDatasets
set_default_plot_size(18cm, 8cm)

labs = ["exp", "sqrt", "log", "winsor", "linear"]
funcs = [x->60*(1.0.-exp.(-0.2*x)), x->sqrt.(x)*10, x->log.(x)*10,
    x->clamp.(x,5,26), x->x*0.6]
x = [1.0:30;]
D = vcat([DataFrame(x=x, y=f(x), linev=l) for (f,l) in zip(funcs, labs)]...)
D[134:136,:y] .= NaN

p1 = plot(D, x=:x, y=:y, linestyle=:linev, Geom.line)
# p2 = plot(dataset("datasets", "CO2"), x=:Conc, y=:Uptake,
    # group=:Plant, linestyle=:Treatment, color=:Type, Geom.line,
    # Scale.linestyle_discrete(order=[2,1]),
    # Theme(key_position=:top, key_title_font_size=-8mm))
display(p1)

#############################################
## 7) Monkey code for pulse_test1 ???

S_vecL, t_vecL, S_vec0, t_vec0, S_vec2, t_vec2, S_vec7, t_vec7, Pstack  = pulse_test(Cmax, Req, C2)



pltPulse = plot(
            layer(x= t_vecL, y= S_vecL, Geom.line, Theme(default_color = "grey")),
            layer(x = t_vec0, y = S_vec0, Geom.line, Theme(default_color = "red")),
            layer(x = t_vec2, y = S_vec2, Geom.line, Theme(default_color = "green")),
            layer(x = t_vec7, y = S_vec7, Geom.line, Theme(default_color = "Blue")),
            Guide.xlabel("Time (s)"),
            Guide.ylabel("SOC (%)"),
            Guide.manual_color_key("", ["My Linear Model", "No Ohmic Drop","Simple Ohmic Drop", "2nd Order RC"], ["grey", "red", "green", "Blue"]),
            Guide.title("Linear OCV Performance -- 18650 Cell @ $(round(Int, P))W"))


Pplot= []
for k in 1:length(Pstack)
    append!(Pplot, [Pstack[k], Pstack[k]])
end
pltP = plot(x = t_vecL, y = Pplot, Geom.line, Guide.xlabel("Time (s)"), Guide.ylabel(" Power (Watts)"), Coord.cartesian(ymin=0))

plt_stack = vstack(compose(context(0,0,1.0, 0.5), render(pltPulse)),
                   compose(context(0,0,1.0, 0.5), render(pltP)) )

display(plt_stack)
