using Plots
using DifferentialEquations
using Sundials
using DataFrames

pyplot()

include("ConservationEquations.jl")
pth_in = pwd() * "/inputs/"
pth_out= pwd() * "/outputs/"

fluidisation_props = YAML.load_file((pth_out * "fluidisation_parameters.yml"))
u_b  = fluidisation_props["bubble velocity"]

ρ_gas = 2.98 

ṁ = 1813.616 * 1000 / 3600 / 2# kg/s
ṅ = ṁ / 30.217 * 1000      # mol/s
V̇ = ṁ / ρ_gas              # m³/s

x_0 = [25.14; 3.550; 21.7; 3.915; 45.906] / 100
n_0 = ṅ .* x_0

# water addition

n_H2O_stoch = n_0[1]
stoch_ratio = 2
n_H2O = n_H2O_stoch * stoch_ratio

n_H2O_add = n_H2O - n_0[4]
n_0[4] = n_0[4] + n_H2O_add

c_0 = n_0 ./ V̇
c_0 = [c_0;
       0.;0.;0.;0.;0.;
       0.;0.;0.;0.;0.;
       (350. + 273); (350. + 273)]

dc_0 = [0.;0.;0.;0.;0.
        0.;0.;0.;0.;0.
        0.;0.;0.;0.;0.;
        0.;0.];

differential_vars = [true, true, true, true, true,
                    false, false, false, false, false,
                    false, false, false, false, false,
                    true, false]

tspan = (0.0,40.0)
p =[1]

println("************************\n DAE SOLVER STARTING \n************************")
prob = DAEProblem(rxn, dc_0, c_0, tspan, p, differential_vars=differential_vars)
sol = solve(prob,IDA(linear_solver=:Dense));
println("\n")
println("************************\n DAE SOLVER COMPLETED \n************************")



co_b = sol[1,:]
h2_b = sol[2,:]
h2o_b = sol[4,:]
co_conv = (c_0[1] .- co_b)./c_0[1]
space = u_b * sol.t

b = plot(space, [(co_b*V̇)/sum(ṅ), (h2_b*V̇)/sum(ṅ), (h2o_b*V̇)/sum(ṅ), co_conv],
     labels= ["CO mole fraction" "H2 mole fraction" "H2O mole fraction" "CO conversion"],
     xaxis  = "length x / m",
     title = "Bubble Phase Results", show = true, reuse = false)


co_c = sol[6,:]
h2_c = sol[7,:]
h2o_c = sol[9,:]
co_conv = (c_0[1] .- co_b)./c_0[1]
space = u_b * sol.t

c = plot(space, [(co_c*V̇)/sum(ṅ), (h2_c*V̇)/sum(ṅ), (h2o_c*V̇)/sum(ṅ), co_conv],
     labels= ["CO mole fraction" "H2 mole fraction" "H2O mole fraction" "CO conversion"],
     xaxis  = "length x / m",
     title = "Cloud Phase Results", show = true, reuse = false)


co_e = sol[11,:]
h2_e = sol[12,:]
h2o_e = sol[14,:]
co_conv = (c_0[1] .- co_b)./c_0[1]
space = u_b * sol.t

e = plot(space, [(co_e*V̇)/sum(ṅ), (h2_e*V̇)/sum(ṅ), (h2o_e*V̇)/sum(ṅ), co_conv],
     labels= ["CO mole fraction" "H2 mole fraction" "H2O mole fraction" "CO conversion"],
     xaxis  = "length x / m",
     title = "Emulsion Phase Results", show = true, reuse = false)     


t = plot(space, [sol[16,:], sol[17,:]],
     labels= ["Bubble Temperature" "Constant Emulsion Temperature at Design Length"],
     xaxis  = "length x / m",
     title = "Temperature Results", show = true, reuse = false)     

