## Test the instances with EnergyCommunity.jl
#
# This file aims to obtain the results of the instances with EnergyCommunity.jl.
# If a stochastic instance is requested, it will be solved with the stochastic version of EnergyCommunity.jl, otherwise it will be solved with the deterministic version of EnergyCommunity.jl.
# A stochastic instance is identified by the presence of "_sto.yml" in the name of the configuration file, otherwise it is considered deterministic.
#
# To run this file you can run in the terminal:
# julia test_instance_with_EC_jl.jl {file_name [optional]}
# the results will be printed in the terminal.
#
# {file_name} is the optional name of the configuration file, which is a YAML file, e.g. `energy_community_model.yaml`.
#
# Each instance is defined by:
# - configuration file: configuration file of the instance, which contains the parameters of the model. 
#   The configuration file is a YAML file, e.g. `energy_community_model.yaml`
# - data file(s): data file of the instance, which contains the data of the model. 
#   The data file(s) are CSV files, e.g. `input_resources.csv` as defined in the 
#   configuration file under general->optional_datasets
#
# For more details on the configuration files, please see the documentation of EnergyCommunity.jl:
# https://spsunipi.github.io/EnergyCommunity.jl/dev/configuration/configuration/
#
# This file is structured as follows:
# 1. Imports the necessary packages
# 2. Defines the path to the configuration file
# 3. Create the model and solve it
# 4. Print the results

## 1. Imports the necessary packages

# Package manager to setup the environment
import Pkg
Pkg.activate(".")
# `Pkg.instantiate()` is intentionally NOT called here: the script later
# does `Pkg.add(url=..., rev=...)` to pin EnergyCommunity to the right branch
# (`main` vs `stochastic`), which would conflict with a pre-existing manifest
# that does not list EnergyCommunity yet on a fresh checkout.

## 2. Defines the path to the configuration file

# Load the data
fconfig = "energy_community_model_new.yml"  # default value
no_thermal = false  # --no-thermal: fix every thermal generator install to 0
no_asset   = false  # --no-asset: fix every installable asset to 0 (NA case)
for a in ARGS
    global fconfig, no_thermal, no_asset
    if a == "--no-thermal"
        no_thermal = true
    elseif a == "--no-asset"
        no_asset = true
    elseif !startswith(a, "--")
        fconfig = string(a)
    end
end
println("Using configuration file: ", fconfig)
no_thermal && println("Fixing thermal generator installs to 0 (--no-thermal).")
no_asset   && println("Fixing all installable assets to 0 (--no-asset).")

# define if network is stochastic: if the configuration file contains "_sto.yml"
# it is considered stochastic, otherwise it is deterministic

is_stochastic = occursin("_sto.yml", fconfig)

# Ensure environment is set up with the correct version of EnergyCommunity.jl;
# this MUST happen before the `using EnergyCommunity` below, otherwise Julia
# would fail to find the package.
if is_stochastic
    Pkg.add(url="https://github.com/SPSUnipi/EnergyCommunity.jl", rev="stochastic")
else
    Pkg.add(url="https://github.com/SPSUnipi/EnergyCommunity.jl", rev="main")
end

using EnergyCommunity
using Gurobi  # Commercial solver; if you don't have a license you can use HiGHS
using HiGHS
using JuMP
using Random
# `Distributions` and `StochasticPrograms` are needed by the sampler files
# we `include` below (they declare `@sampler` / `@define_scenario` macros
# and the `MvNormal`/`Normal`/`truncated`/`pem` calls).
using Distributions
using StochasticPrograms
using PointEstimateMethod

# Use the same RNG seed as `csv2nc4.jl` so that the (s, eps) scenario draws
# coincide between the SMS++ TSSB pipeline and this script — necessary for
# the FO comparison to be meaningful.
Random.seed!(123)

# default solver
optimizer = Gurobi.Optimizer  # if you have a Gurobi license, otherwise use HiGHS
# optimizer = HiGHS.Optimizer  # Uncomment this line to use HiGHS instead of Gurobi

## 3. Create the model and solve it

obj_value = nothing
optimal_design = nothing

if is_stochastic
    println("The model is stochastic.")

    # Mirrors `examples/RunStochModel(TBD).jl` of
    # https://github.com/SPSUnipi/EnergyCommunity.jl branch `stochastic`,
    # so the YAML must follow that schema (keys `general.n_s`, `general.n_eps`,
    # `general.uncertain_var`, per-component `profile.std`, single
    # `market.profile`).

    # Re-read the data already loaded above so we can pull the stochastic
    # parameters out of the YAML; `read_input` is provided by EnergyCommunity.
    data = read_input(fconfig)
    (gen_data, users_data, market_data) = explode_data(data)

    scen_s_sample = field(gen_data, "n_s")
    scen_eps_sample = field(gen_data, "n_eps")
    unc_var = field(gen_data, "uncertain_var")

    # Long-period uncertainty parameters. Defaults match the upstream example;
    # override via `general.sigma_load`, `general.sigma_pv`, `general.sigma_wind`,
    # `general.mean_pv`, `general.mean_wind` if present in the YAML.
    sigma_load = get(gen_data, "sigma_load", 0.3)
    mean_pv    = get(gen_data, "mean_pv",    1.0)
    sigma_pv   = get(gen_data, "sigma_pv",   0.1)
    mean_wind  = get(gen_data, "mean_wind",  0.95)
    sigma_wind = get(gen_data, "sigma_wind", 0.15)

    # Globals used by the repo-side sampler (peak_set / time_set / user_set
    # are referenced as free variables inside `scen_eps_sampler.jl`).
    init_step = field(gen_data, "init_step")
    final_step = field(gen_data, "final_step")
    time_set = init_step:final_step
    user_set = user_names(gen_data, users_data)
    peak_set = unique(profile(market_data, "peak_categories")[time_set])

    # Repo-side helpers used by `scen_eps_sampler.jl` that are NOT exported
    # by EC.jl@stochastic. Defined here at top-level so the included file
    # can call them directly.
    profile_d_(d, name, default=nothing) = let p = get(d, "profile", nothing)
        isnothing(p) ? default : get(p, name, default)
    end
    profile_component_d_(d, c_name, name, default=nothing) =
        haskey(d, c_name) ? profile_d_(d[c_name], name, default) : default
    # The included sampler refers to these as `profile_d` / `profile_component_d`.
    profile_d           = profile_d_
    profile_component_d = profile_component_d_

    # Override the EC.jl-shipped sampling (`pem_extraction` and
    # `scenarios_generator`) with the repo versions so that `csv2nc4.jl`
    # — which writes the SMS++ TSSB — and this script draw bit-identical
    # scenarios when both seed Random with the same value. The
    # `Scenario_Load_Renewable` struct itself is NOT re-included from the
    # repo: `StochasticEC(...)` dispatches on `EnergyCommunity.Scenario_*`,
    # so we keep the EC.jl-side definition (its fields match the repo one).
    include(joinpath(@__DIR__, "pem_extraction.jl"))
    include(joinpath(@__DIR__, "scen_eps_sampler.jl"))

    # Reset the RNG immediately before sampling so the sequence is the same
    # one consumed by `csv2nc4.jl`.
    Random.seed!(123)

    # Sample the long-period (`s`) points and their probabilities.
    (point_s_load, point_s_pv, point_s_wind, scen_probability) =
        pem_extraction(scen_s_sample, sigma_load,
                       mean_pv, sigma_pv,
                       mean_wind, sigma_wind,
                       unc_var)

    # Build the full (s, eps) scenario set. The repo `scenarios_generator`
    # has signature
    #   scenarios_generator(data, point_s_load, point_s_pv, point_s_wind,
    #                       n_s, n_eps, unc_var; point_probability,
    #                       control_risimulation)
    # which matches the EC.jl@stochastic spec used by the rest of this script.
    sampled_scenarios = scenarios_generator(
        data,
        point_s_load, point_s_pv, point_s_wind,
        scen_s_sample, scen_eps_sample, unc_var;
        point_probability=scen_probability,
    )

    # Cooperative (CO) stochastic model.
    model = StochasticEC(fconfig, EnergyCommunity.GroupCO(), optimizer,
                         sampled_scenarios, scen_s_sample, scen_eps_sample)
    build_specific_model!(EnergyCommunity.GroupCO(), model, optimizer)

    # Emulate the "no thermal" and "no-asset" (NA) variants in the stochastic
    # flow, where the YAML schema has no explicit knob for them. The
    # first-stage design variable `x_us[u, a]` is forced to 0 for every asset
    # we want to remove. EC.jl@stochastic builds the deterministic equivalent,
    # so reaching the proxy first-stage variables goes through
    # `proxy(stochasticprogram, 1)`.
    if no_thermal || no_asset
        proxy_model = StochasticPrograms.proxy(model.model, 1)
        x_us = proxy_model[:x_us]
        for u in user_set
            for a in device_names(users_data[u])
                if no_asset ||
                   ( no_thermal &&
                     EnergyCommunity.asset_type( users_data[u] , a ) ==
                       EnergyCommunity.THER )
                    JuMP.fix( x_us[u, a] , 0.0 )
                end
            end
        end
    end

    # Solver tuning. `set_parameters_ECmodel!` upstream only knows about
    # CPLEX / Gurobi (it dispatches via `occursin(name, opt)`); for
    # other solvers (HiGHS in particular) calling it raises and resets the
    # scheduler, so we skip it and let the solver use its defaults.
    if occursin("CPLEX", string(optimizer)) ||
       occursin("Gurobi", string(optimizer))
        set_parameters_ECmodel!(model, 1e-6, 60 * 60, Threads.nthreads(), 1)
    end

    # Solves the deterministic equivalent and stores the solution into
    # `model.results`.
    optimize_deterministic_ECmodel(model)

    obj_value = objective_value(model.model)
    optimal_design = model.results[:x_us].data
else
    println("The model is deterministic.")

    # Create the model
    model = ModelEC(fconfig, EnergyCommunity.GroupCO(), optimizer)

    # build the model
    build_model!(model)

    # Solve the model
    optimize!(model)

    obj_value = objective_value(model)
    optimal_design = value.(model.results[:x_us])
end


## 4. Print the results

println("Optimal value: ", obj_value)
println("Optimal installed capacity by user: ", optimal_design)