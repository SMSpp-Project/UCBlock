# Function used to sample the load and renewable distributions for long period uncertainty using the Point Estimate Method

# scen_s_sample = number of scenarios s
# sigma_load = standard deviation associated to load distribution
# sigma_ren = standard deviation associated to renewable prouction distribution

function pem_extraction(scen_s_sample::Int, sigma_load, sigma_ren)

    if scen_s_sample > 1 # more than a single scenario s to generate
        Distribution_load = truncated(Normal(1.0, sigma_load), 0.0, +Inf) # load distribution
        Distribution_ren = truncated(Normal(1.0, sigma_ren), 0.0, +Inf) # renewable production distribution
        pem_load = pem(Distribution_load, scen_s_sample) # extracted point for load (with their probability)
        pem_ren = pem(Distribution_ren, scen_s_sample) # extracted point for renewable production

        point_load = pem_load.x
        point_ren = pem_ren.x

        scen_probability = pem_load.p # probability associated to each extracted point from load distribution
    else
        point_load = [1.0]
        point_ren = [1.0]
        scen_probability = [1.0]
    end

    return point_load, point_ren, scen_probability
end
