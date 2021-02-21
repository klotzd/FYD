#! usr/bin/env julia

using YAML
using Roots
using JSON

println("************************ \nFLUIDIZATION CALCULATION \n ************************")

cat_props = YAML.load_file("catalystproperties.yml")

d_p = cat_props["diameter"] *10^-6
ρ_cat = cat_props["catalyst density"]
μ = cat_props["gas_viscosity"]
ρ_gas = cat_props["gas_density"]
d_b = cat_props["bubble_diameter"]
D   = cat_props["bed_diameter"]
ṁ   = cat_props["mass_flow"]


trigger = false

# initial flow calculations
V̇ = ṁ / 2/ ρ_gas
u_0 = V̇ / (D^2 / 4 * pi)
g = 9.81

# minimum fluidization
ε_mf = 0.586 * ( μ^2 / (g * (ρ_cat - ρ_gas) * ρ_gas * d_p^3) )^(0.029) * (ρ_gas/ρ_cat)^(0.021)
u_mf = d_p^2/(150 *μ) * g * (ρ_cat - ρ_gas) * (ε_mf^3 / (1-ε_mf))

# transport velocity
Re_mf = u_mf * d_p * ρ_cat / μ

if Re_mf < 1
    # Stoke's Law
    u_t = g * (ρ_cat - ρ_gas) * d_p^2 / (18 * μ)

elseif (Re_mf > 1) && (Re_mf < 1000)
    # Trambouze's Correlation
    f(x) = x-(4*d_p*(ρ_cat-ρ_gas)*g)/(3*ρ_gas*exp(-5.5+69.43/(log(ρ_gas*d_p*x/μ)+7.99)))^(1/2)
    u_t = fzero(f, 3)
else
    println("Warning: Re_mf indicates turbulent \n No drag coefficient method available")
end

# velocity check
println(" \n fluidization check: \n ")

if (u_0 < u_t) && (u_0 > u_mf)
    println("fludisation possible")
else
    println("ERROR: fluidisation not possible")
    trigger = true
end

# bubble velocity
if d_b == "estimate"
    d = 0.00376 * (u_0*100 - u_mf*100)^2 / 100
    d_b = d
elseif d_b == "Geldart A"
    d_ = 0.1
    d_b = d
else
    d_b = d_b
end

u_br = 0.64 * (d_b * g)^(0.5) 
u_b = u_0 - u_mf + u_br      

# phase fractions
if round(u_b / u_mf) > 6
    f_b = (u_0-u_mf)/(u_b)
else
    if d_p > (550 * 10^-6)
        α = 0.4
    else
        α = 0.5
    end
    f_b = (u_0-u_mf)/(u_b-u_mf*(1+α))
end

f_c = (3*u_mf/ε_mf)/(0.711*(g*d_b)^2-u_mf/ε_mf)
f_e = 1 - f_b - f_c
ε = (1-f_b)*(1-ε_mf)

# catalyst distribution
γ_b = f_b * (1 - ε)
γ_c = f_c
γ_e = f_e / (1 - ε_mf)

##### FOR NOW APPROXIMATE
γ_b = .005
γ_c = .3
γ_e = 1.5

# collect results

if trigger
    println("************************ \nFLUIDIZATION ERROR \n ************************")
else
fluidization = Dict("min fluid velocity" => u_mf,
            "min fluid voidage" => ε_mf,
            "transport velocity"=> u_t,
            "superficial velocity" => u_0,
            "bubble velocity"=> u_b,
            "bubble diameter" => d_b,
            "bubble fraction" => f_b,
            "cloud fraction" => f_c,
            "emulsion fraction" => f_e,
            "voidage" => ε,
            "gamma_b" => γ_b,
            "gamma_c" => γ_c,
            "gamma_e" => γ_e)

println("Fluidisation Properties:")
print(json(fluidization,4))

# write results            
YAML.write_file("fluidisation_parameters.yml", fluidization)

println("************************ \nFLUIDIZATION COMPLETED \n ************************")
end