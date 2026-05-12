## Test the instances with EnergyCommunity.jl
#
# Solve an Energy Community instance with EnergyCommunity.jl and print its
# objective value. The instance is selected by the YAML configuration file
# passed on the command line: a `_sto.yml` filename triggers the stochastic
# branch (using EnergyCommunity.jl@stochastic), any other filename triggers
# the deterministic branch (using EnergyCommunity.jl@main).
#
# Usage:
#
#     julia test_instance_with_EC_jl.jl [yml] [--no-thermal | --no-asset]
#
# - `[yml]` — YAML configuration file (default `energy_community_model_new.yml`).
# - `--no-thermal` — fix every thermal generator install to 0 (matches the
#   SMS++ pipeline run without `--with-thermal-blocks`).
# - `--no-asset`   — fix every installable asset to 0 (NA case, replaces the
#   removed `*_NA*.yml` files).
#
# For more details on the YAML schema see
# https://spsunipi.github.io/EnergyCommunity.jl/dev/configuration/configuration/.

## 1. Imports the necessary packages

import Pkg
Pkg.activate(".")
# `Pkg.instantiate()` is intentionally NOT called here: this script later
# does `Pkg.add(url=..., rev=...)` to pin EnergyCommunity to the right branch
# (`main` vs `stochastic`), which would conflict with a pre-existing manifest
# that does not list EnergyCommunity yet on a fresh checkout.

## 2. Parse command-line arguments

fconfig    = "energy_community_model_new.yml"  # default YAML
no_thermal = false                              # --no-thermal flag
no_asset   = false                              # --no-asset flag
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

# Stochastic vs deterministic flow selection by filename convention.
is_stochastic = occursin("_sto.yml", fconfig)

# Pin the right EnergyCommunity.jl revision before `using EnergyCommunity`.
# Set `SKIP_PKG_ADD=1` in the environment to skip this network-dependent
# step on local re-runs (the package must already be in the manifest).
if get(ENV, "SKIP_PKG_ADD", "0") != "1"
    if is_stochastic
        Pkg.add(url="https://github.com/SPSUnipi/EnergyCommunity.jl", rev="stochastic")
    else
        Pkg.add(url="https://github.com/SPSUnipi/EnergyCommunity.jl", rev="main")
    end
end

using EnergyCommunity
using Gurobi  # commercial solver; if not licensed, switch `optimizer` to HiGHS
using HiGHS
using JuMP
using Random
# `Distributions` and `StochasticPrograms` are required by `pem_extraction.jl`
# / `scen_eps_sampler.jl` (the `@sampler` / `@define_scenario` macros, the
# `MvNormal` / `Normal` / `truncated` / `pem` calls).
using Distributions
using StochasticPrograms
using PointEstimateMethod

# Same RNG seed used by `csv2nc4.jl` so that the `(s, eps)` scenario draws
# coincide between the SMS++ TSSB pipeline and this script.
Random.seed!(123)

optimizer = Gurobi.Optimizer
# optimizer = HiGHS.Optimizer  # uncomment to use HiGHS instead

## 3. Build the model and solve it

obj_value      = nothing
optimal_design = nothing

if is_stochastic
    println("The model is stochastic.")

    # YAML schema mirrors `examples/RunStochModel(TBD).jl` of
    # https://github.com/SPSUnipi/EnergyCommunity.jl branch `stochastic`:
    # `general.n_s`, `general.n_eps`, `general.uncertain_var`, per-component
    # `profile.std`, single `market.profile`.

    data = read_input(fconfig)
    (gen_data, users_data, market_data) = explode_data(data)

    scen_s_sample   = field(gen_data, "n_s")
    scen_eps_sample = field(gen_data, "n_eps")
    unc_var         = field(gen_data, "uncertain_var")

    # Long-period uncertainty parameters; defaults match the upstream example.
    # Override via `general.{sigma_load,sigma_pv,sigma_wind,mean_pv,mean_wind}`
    # in the YAML.
    sigma_load = get(gen_data, "sigma_load", 0.3)
    mean_pv    = get(gen_data, "mean_pv",    1.0)
    sigma_pv   = get(gen_data, "sigma_pv",   0.1)
    mean_wind  = get(gen_data, "mean_wind",  0.95)
    sigma_wind = get(gen_data, "sigma_wind", 0.15)

    # Globals expected by the repo-side sampler (referenced as free variables
    # inside `scen_eps_sampler.jl`).
    init_step  = field(gen_data, "init_step")
    final_step = field(gen_data, "final_step")
    time_set   = init_step:final_step
    user_set   = user_names(gen_data, users_data)
    peak_set   = unique(profile(market_data, "peak_categories")[time_set])

    # Repo-side helpers used by `scen_eps_sampler.jl` that are not exported
    # by EnergyCommunity.jl@stochastic. Defined at top level so the included
    # file can call them directly.
    profile_d_(d, name, default=nothing) = let p = get(d, "profile", nothing)
        isnothing(p) ? default : get(p, name, default)
    end
    profile_component_d_(d, c_name, name, default=nothing) =
        haskey(d, c_name) ? profile_d_(d[c_name], name, default) : default
    profile_d           = profile_d_
    profile_component_d = profile_component_d_

    # Use the repo-side sampling (`pem_extraction` + `scenarios_generator`)
    # so that this script and `csv2nc4.jl` draw bit-identical scenarios when
    # both seed `Random` with the same value. The `Scenario_Load_Renewable`
    # struct itself is left to EnergyCommunity.jl@stochastic — `StochasticEC`
    # dispatches on `EnergyCommunity.Scenario_*`.
    include(joinpath(@__DIR__, "pem_extraction.jl"))
    include(joinpath(@__DIR__, "scen_eps_sampler.jl"))

    # Reset the RNG immediately before sampling so the sequence is the same
    # one consumed by `csv2nc4.jl`.
    Random.seed!(123)

    (point_s_load, point_s_pv, point_s_wind, scen_probability) =
        pem_extraction(scen_s_sample, sigma_load,
                       mean_pv, sigma_pv,
                       mean_wind, sigma_wind,
                       unc_var)

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

    # `--no-thermal` / `--no-asset` are not knobs of the YAML schema, so we
    # emulate them by fixing the corresponding first-stage `x_us[u, a]` to 0
    # on the deterministic-equivalent's first-stage proxy.
    if no_thermal || no_asset
        proxy_model = StochasticPrograms.proxy(model.model, 1)
        x_us = proxy_model[:x_us]
        for u in user_set
            for a in device_names(users_data[u])
                if no_asset ||
                   ( no_thermal &&
                     EnergyCommunity.asset_type(users_data[u], a) ==
                       EnergyCommunity.THER )
                    JuMP.fix(x_us[u, a], 0.0)
                end
            end
        end
    end

    # Solver tuning. `set_parameters_ECmodel!` upstream only knows about
    # CPLEX / Gurobi (it dispatches via `occursin(name, opt)`); for any other
    # solver (HiGHS in particular) calling it raises and resets the
    # scheduler, so we skip it and let the solver use its defaults.
    if occursin("CPLEX", string(optimizer)) ||
       occursin("Gurobi", string(optimizer))
        set_parameters_ECmodel!(model, 1e-6, 60 * 60, Threads.nthreads(), 1)
    end

    optimize_deterministic_ECmodel(model)

    # The SMS++ TSSB test compares against the LP-relaxation value (the
    # LagrangianDualSolver in BSPar-2S.txt converges to the LP bound, and
    # the MILPSolver runs with `intRelaxIntVars=1`), so the EC.jl
    # reference must also be the LP optimum. The bulk
    # `StochasticPrograms.relax_integrality` entry point skips Decision
    # variables, and `JuMP.relax_integrality` on the materialized
    # `deterministic_model` doesn't propagate either. The pure-API way
    # is to iterate every Decision (stage 1 first-stage + per-scenario
    # stage > 1) plus every regular JuMP variable of the DEP and unset
    # the integer / binary attribute one-by-one.
    let sp = model.model, det = model.deterministic_model
        # First-stage Decisions.
        for dvar in StochasticPrograms.all_decision_variables(sp, 1)
            JuMP.is_integer(dvar) && JuMP.unset_integer(dvar)
            JuMP.is_binary(dvar)  && JuMP.unset_binary(dvar)
        end
        # Stage > 1 Decisions (scenario-dependent).
        n_scen = StochasticPrograms.num_scenarios(sp)
        for stage in 2:StochasticPrograms.num_stages(sp)
            for dvar in StochasticPrograms.all_decision_variables(sp, stage)
                for s in 1:n_scen
                    JuMP.is_integer(dvar, s) && JuMP.unset_integer(dvar, s)
                    JuMP.is_binary(dvar, s)  && JuMP.unset_binary(dvar, s)
                end
            end
        end
        # Plain JuMP variables on the DEP that are not Decisions.
        for var in JuMP.all_variables(det)
            JuMP.is_integer(var) && JuMP.unset_integer(var)
            JuMP.is_binary(var)  && JuMP.unset_binary(var)
        end
        JuMP.optimize!(det)
    end

    obj_value      = JuMP.objective_value(model.deterministic_model)
    optimal_design = model.results[:x_us].data
else
    println("The model is deterministic.")

    model = ModelEC(fconfig, EnergyCommunity.GroupCO(), optimizer)
    build_model!(model)

    # Deterministic equivalent of `--no-thermal` / `--no-asset`. In the
    # deterministic EC.jl model, `x_us` is an `@expression` defined as
    # `n_us[u, a] * nom_capacity` (see EnergyCommunity/base_model.jl); the
    # actual decision is `n_us`, so fixing `n_us` to 0 is equivalent to
    # fixing `x_us = 0`. `force=true` is required because `n_us` already
    # has a default lower bound of 0.
    if no_thermal || no_asset
        n_us = model.model[:n_us]
        for u in model.user_set
            for a in device_names(model.users_data[u])
                if no_asset ||
                   ( no_thermal &&
                     EnergyCommunity.asset_type(model.users_data[u], a) ==
                       EnergyCommunity.THER )
                    JuMP.fix(n_us[u, a], 0.0; force=true)
                end
            end
        end
    end

    # Solver tuning. The deterministic branch can't use
    # set_parameters_ECmodel! (that helper expects the stochastic
    # deterministic_model field), so we set the attributes directly on
    # the underlying JuMP Model.
    if occursin("Gurobi", string(optimizer))
        JuMP.set_optimizer_attribute(model.model, "TimeLimit", 60 * 60)
        JuMP.set_optimizer_attribute(model.model, "Threads", Threads.nthreads())
        JuMP.set_optimizer_attribute(model.model, "OutputFlag", 1)
    elseif occursin("CPLEX", string(optimizer))
        JuMP.set_optimizer_attribute(model.model, "CPX_PARAM_TILIM", 60 * 60)
        JuMP.set_optimizer_attribute(model.model, "CPX_PARAM_THREADS", Threads.nthreads())
        JuMP.set_optimizer_attribute(model.model, "CPX_PARAM_SCRIND", 1)
    end

    # The SMS++ LDS_UC pipeline solves the LP relaxation (BSPar-2S.txt
    # MILPSolver has intRelaxIntVars=1, and LagrangianDualSolver
    # converges to the LP bound), so the EC.jl reference must also be
    # the LP optimum, not the MILP one. Today the LP and MILP coincide
    # on CO/NC because the optimum installs zero thermal, but relaxing
    # the integer/binary variables here keeps the refs valid even if a
    # future YAML makes the install > 0.
    JuMP.relax_integrality(model.model)

    optimize!(model)

    obj_value      = objective_value(model)
    optimal_design = value.(model.results[:x_us])
end

## 4. Print the results
#
# EC.jl maximizes Social Welfare (`@objective(model, Max, SW)`), while SMS++
# minimizes cost on the same instance; for a community whose welfare ends
# up negative (typical when grid imports + investment dominate exports),
# `objective_value` is the negative of the SMS++ cost. Print `-obj_value`
# so the number is directly comparable to (and can be copied into) the
# SMS++ batch-ec reference values.

println("Optimal value (SMS++ cost convention, = -welfare): ", -obj_value)
println("Optimal installed capacity by user: ", optimal_design)
