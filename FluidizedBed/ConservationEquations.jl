#! usr/bin/env julia
using YAML

constants = YAML.load_file("constants.yml")
R = constants["gas constant"]

cat_props = YAML.load_file("catalystproperties.yml")
fluidisation_props = YAML.load_file("fluidisation_parameters.yml")
MT_props = YAML.load_file("MT_parameters.yml")


function rxn(out, dc, c, p, t)
    
    # Temperature
    T = 350 + 273.15 # K
    

    # reaction rate constant
    k = 3.2*10^-6 * exp(-( 112 * 1000) / R * ( 1/T - 1/678) )
    k_cat = cat_props["catalyst density"] * k
    
    # Equilibrium constant
    Z  = (1000/T) - 1;
    K  = exp(Z*(Z*(0.63508-0.29353*Z)+4.1778)+0.31688);
    
    ######
    # variables
    ######
    
    
    # extract concentrations
    # mol/m³
    
    # bubble
    c_b_CO, c_b_H2, c_b_CO2, c_b_H2O, c_b_N2 = c[1:5]
    # cloud
    c_c_CO, c_c_H2, c_c_CO2, c_c_H2O, c_c_N2 = c[6:10]
    # emulsion
    c_e_CO, c_e_H2, c_e_CO2, c_e_H2O, c_e_N2 = c[11:15]
    
    # calculate partial pressures
    # kPa
    
    # bubble
    pi_b_CO  = c_b_CO * R * T / 1000   
    pi_b_H2  = c_b_H2 * R * T / 1000 
    pi_b_CO2 = c_b_CO2 * R * T / 1000 
    pi_b_H2O = c_b_H2O * R * T / 1000
    
    # cloud
    pi_c_CO  = c_c_CO * R * T / 1000   
    pi_c_H2  = c_c_H2 * R * T / 1000 
    pi_c_CO2 = c_c_CO2 * R * T / 1000 
    pi_c_H2O = c_c_H2O * R * T / 1000
    
    # emulsion
    pi_e_CO  = c_e_CO * R * T / 1000   
    pi_e_H2  = c_e_H2 * R * T / 1000 
    pi_e_CO2 = c_e_CO2 * R * T / 1000 
    pi_e_H2O = c_e_H2O * R * T / 1000
    
    # calculate rates in bubble
    # mol/sec m³ bubble
    
    rb = k_cat * exp(-( 112 * 1000) / R * ( 1/T - 1/678) ) * 
        (pi_b_CO * pi_b_H2O - pi_b_CO2 * pi_b_H2 / K)
    
    # calculate rates in cloud
    # mol/sec m³ bubble
    
    rc = k_cat * exp(-( 112 * 1000) / R * ( 1/T - 1/678) ) * 
        (pi_c_CO * pi_c_H2O - pi_c_CO2 * pi_c_H2 / K)
    
    # calculate rates in emulsion
    # mol/sec m³ bubble
    
    re = k_cat * exp(-( 112 * 1000) / R * ( 1/T - 1/678) ) * 
        (pi_e_CO * pi_e_H2O - pi_e_CO2 * pi_e_H2 / K)
    
    
    #####
    # Balances
    #####
    
    # BUBBLE
    
    out[1] = dc[1] - fluidisation_props["gamma_b"] * rb - MT_props["K_bubble-CO"] * (c_b_CO - c_c_CO)        # CO
    out[2] = dc[2] + fluidisation_props["gamma_b"]  * rb - MT_props["K_bubble-H2"] * (c_b_H2 - c_c_H2)       # H2
    out[3] = dc[3] + fluidisation_props["gamma_b"]  * rb - MT_props["K_bubble-CO2"] * (c_b_CO2 - c_c_CO2)    # CO2
    out[4] = dc[4] - fluidisation_props["gamma_b"] * rb - MT_props["K_bubble-H2O"] * (c_b_H2O - c_c_H2O)     # H2O
    out[5] = dc[5] - MT_props["K_bubble-N2"]* (c_b_N2 - c_c_N2)                                              # N2
    
    # CLOUD
    
    out[6]  = -fluidisation_props["gamma_c"] * rc + MT_props["K_cloud-CO"] * (c_c_CO - c_e_CO) - MT_props["K_bubble-CO"] * (c_b_CO - c_c_CO)        # CO
    out[7]  =  fluidisation_props["gamma_c"] * rc + MT_props["K_cloud-H2"] * (c_c_H2 - c_e_H2) - MT_props["K_bubble-H2"] * (c_b_H2 - c_c_H2)        # H2
    out[8]  =  fluidisation_props["gamma_c"]* rc + MT_props["K_cloud-CO2"] * (c_c_CO2 - c_e_CO2) - MT_props["K_bubble-CO2"] * (c_b_CO2 - c_c_CO2)   # CO2
    out[9]  = -fluidisation_props["gamma_c"] * rc + MT_props["K_cloud-H2O"] * (c_c_H2O - c_e_H2O) - MT_props["K_bubble-H2O"] * (c_b_H2O - c_c_H2O)  # H2O
    out[10] =  MT_props["K_cloud-N2"] * (c_c_N2 - c_e_N2) - MT_props["K_bubble-N2"] * (c_b_N2 - c_c_N2)                                             # N2
    
    # EMULSION
    
    out[11] = -fluidisation_props["gamma_c"] * re - MT_props["K_cloud-CO"] * (c_c_CO - c_e_CO)       # CO
    out[12] =  fluidisation_props["gamma_c"] * re - MT_props["K_cloud-H2"] * (c_c_H2 - c_e_H2)       # H2
    out[13] =  fluidisation_props["gamma_c"] * re - MT_props["K_cloud-CO2"] * (c_c_CO2 - c_e_CO2)    # CO2
    out[14] = -fluidisation_props["gamma_c"] * re - MT_props["K_cloud-H2O"] * (c_c_H2O - c_e_H2O)    # H2O
    out[15] =  MT_props["K_cloud-N2"] * (c_c_N2 - c_e_N2)                                            # N2
    
end