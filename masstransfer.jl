#! usr/bin/env julia
using YAML
using JSON

println("************************ \nDIFFUSION CALCULATION \n************************")
pth_in = pwd() * "/inputs/"
pth_out= pwd() * "/outputs/"
# data import section
species_props = YAML.load_file((pth_in * "species_properties.yml"))
feed_props = YAML.load_file((pth_in * "feed_properties.yml"))
fluidisation_props = YAML.load_file((pth_out * "fluidisation_parameters.yml"))

# calculate mean gas props
M_feed = feed_props["CO-frac"] * species_props["CO-M"] +
            feed_props["CO2-frac"] * species_props["CO2-M"] +
            feed_props["H2-frac"] * species_props["H2-M"] +
            feed_props["H2O-frac"] * species_props["H2O-M"] +
            feed_props["N2-frac"] * species_props["N2-M"]

σ_feed = feed_props["CO-frac"] * species_props["CO-sigma"] +
            feed_props["CO2-frac"] * species_props["CO2-sigma"] +
            feed_props["H2-frac"] * species_props["H2-sigma"] +
            feed_props["H2O-frac"] * species_props["H2O-sigma"] +
            feed_props["N2-frac"] * species_props["N2-sigma"]

T_feed = feed_props["Temperature"]
P_feed = feed_props["Pressure"]

species = ["CO", "CO2", "H2", "H2O", "N2"]

D_bubble = Dict()
for molecule in species
    M_string = molecule * "-M"
    σ_string = molecule * "-sigma"
    D_bubble[molecule] = 1.859 * 10^-3 * T_feed^(3/2) * (1/M_feed + 1/species_props[M_string])^(1/2) / (P_feed / 101.325 * 1/2 * (σ_feed + species_props[σ_string]))
end

u_mf = fluidisation_props["min fluid velocity"]
d_b  = fluidisation_props["bubble diameter"]
ε_mf = fluidisation_props["min fluid voidage"]
u_b  = fluidisation_props["bubble velocity"]
g = 9.81

K = Dict()

for molecule in species
    bubble_string = "K_bubble-" * molecule 
    K[bubble_string] = 4.5 * (u_mf/d_b) + 5.85 * (D_bubble[molecule]^(1/2) * (g*100)^(1/4) / (d_b*100)^(5/4))
end

for molecule in species
    cloud_string = "K_cloud-" * molecule 
    K[cloud_string] = 6.77 * (ε_mf * D_bubble[molecule]/10000 * u_b / (d_b^3))^(1/2)
end

println("Diffusivities: \n")
print(json(D_bubble,4))
println("\n")
println("Interfacial Coefficients: \n")
print(json(K,4))
println("\n")
println("************************ \nDIFFUSION COMPLETED \n************************")
YAML.write_file((pth_out * "MT_parameters.yml"), K)