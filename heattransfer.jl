#! usr/bin/env julia
using YAML
using JSON

println("************************ \nHEAT TRANSFER CALCULATION \n************************")
pth_in = pwd() * "/inputs/"
pth_out= pwd() * "/outputs/"

# data import section
species_props = YAML.load_file((pth_in * "species_properties.yml"))
cat_props = YAML.load_file((pth_in * "catalyst_properties.yml"))
feed_props = YAML.load_file((pth_in * "feed_properties.yml"))
fluidisation_props = YAML.load_file((pth_out *"fluidisation_parameters.yml"))

# calculate mean gas props
M_feed = feed_props["CO-frac"] * species_props["CO-M"] +
            feed_props["CO2-frac"] * species_props["CO2-M"] +
            feed_props["H2-frac"] * species_props["H2-M"] +
            feed_props["H2O-frac"] * species_props["H2O-M"] +
            feed_props["N2-frac"] * species_props["N2-M"]

cp_feed = feed_props["CO-frac"] * species_props["CO-cp"] * species_props["CO-M"] +
            feed_props["CO2-frac"] * species_props["CO2-cp"] * species_props["CO2-M"]+
            feed_props["H2-frac"] * species_props["H2-cp"] * species_props["H2-M"] +
            feed_props["H2O-frac"] * species_props["H2O-cp"] * species_props["H2O-M"] +
            feed_props["N2-frac"] * species_props["N2-cp"] * species_props["N2-M"]

k_feed = (feed_props["CO-frac"] * species_props["CO-k"] +
            feed_props["CO2-frac"] * species_props["CO2-k"] +
            feed_props["H2-frac"] * species_props["H2-k"] +
            feed_props["H2O-frac"] * species_props["H2O-k"] +
            feed_props["N2-frac"] * species_props["N2-k"] ) * 1000

T_feed = feed_props["Temperature"]
P_feed = feed_props["Pressure"]
ρ_gas =  cat_props["gas_density"]
D = cat_props["bed_diameter"]
μ = cat_props["gas_viscosity"]
ε = fluidisation_props["voidage"]
u_mf = fluidisation_props["min fluid velocity"]
d_b  = fluidisation_props["bubble diameter"]

g = 9.81

# interfacial coefficient from kunii & levenspiel:
H_if = 1.32 * 4.5 * (u_mf * ρ_gas * cp_feed / d_b) + 5.85 * ((k_feed * cp_feed * ρ_gas)^(1/2) * (g*100)^(1/4) / (d_b*100)^(5/4))

# wall coefficient from walton & levenspiel:
G = (cat_props["mass_flow"] / 2) / (D^2 * pi / 4)
H_w  = cp_feed / 30.32 * G * 0.6 * (cat_props["diameter"] * 10^-6 * (G / ε) / μ)^0.75

# nusselt number sensibility check 
Ga = fluidisation_props["Galileo No"]
Nu_max = (0.09*Ga^0.69)^(0.423) * Ga^0.34 * (0.8)^(1/3)

if (Nu_max < 10) || (Nu_max > 1000)
    println("WARNING: CHECK THE NUSSELT NUMBER \n")
else 
end


heat_transfer = Dict("H_interface" => H_if,
            "H_wall" => H_w,
            "c_p" => cp_feed,
            "Nu_max" => Nu_max)

println("Heat Transfer Properties: \n")
println("\n")
print(json(heat_transfer,4))
println("\n")
# write results            
YAML.write_file((pth_out * "HT_parameters.yml"), heat_transfer)

println("************************ \nHEAT TRANSFER COMPLETED\n************************")