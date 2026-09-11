#!/usr/bin/env julia
#=
  json2nc4.jl  --  Convert single-bus thermal UC instances to SMS++ netCDF.

  Usage:
    julia json2nc4.jl  <input.json>  [output.nc4]

  If output is omitted, replaces .json with .nc4 in the same directory.

  This converter targets the SINGLE-BUS, thermal-only UC instances produced by
  `generate_single_bus_uc` in DatasetGenerator_UCOTS.jl (schema
  "v4.0-singlebus-reserve"). The network / OTS part of the framework is
  intentionally ignored: instances have a single node and a single nodal
  balance, so no DCNetworkBlock / OTSNetworkBlock and no lines are written.

  The produced file conforms to the SMS++ UCBlock format with:
    - one ThermalUnitBlock per generator (all on node 0);
    - ActivePowerDemand on the single node;
    - optional primary (FCR) and secondary (aFRR) reserve: a single reserve
      zone, PrimaryDemand / SecondaryDemand at UCBlock level and PrimaryRho /
      SecondaryRho at ThermalUnitBlock level.
    - optional reactive power: ReactivePowerDemand at UCBlock level and the
      commitment-gated capability MaxReactivePowerOn / MinReactivePowerOn at
      ThermalUnitBlock level (the off-state box defaults to {0}).

  All quantities are assumed to be ALREADY discretised to the instance
  time-step: minimum up/down times in PERIODS, ramps in MW/period, costs per
  period. The converter copies them verbatim.

  Author: Fabrizio Lacalandra / Antonio Frangioni
  Date:   2026
=#

using JSON
using NCDatasets

# =========================================================================
#  Main conversion function
# =========================================================================

"""
    convert_json_to_nc4(json_path::String, nc_path::String)

Read a single-bus thermal UC JSON instance and write the SMS++ netCDF file.
"""
function convert_json_to_nc4(json_path::String, nc_path::String;
                             nuclear::Bool = false, mod_time::Int = 8,
                             mod_frac::Float64 = 0.25)

    data = JSON.parsefile(json_path)

    meta = data["metadata"]
    th   = data["generators"]["thermal"]

    n_gen     = Int(meta["n_gen"])
    n_periods = Int(meta["n_periods"])

    # --- single-bus sanity check -----------------------------------------
    n_bus = Int(get(meta, "n_bus", 1))
    n_bus == 1 || error("json2nc4: this converter only handles single-bus " *
                        "instances (n_bus = $n_bus). Network/OTS is out of scope.")

    # --- thermal generator data ------------------------------------------
    p_min         = Float64.(th["p_min"])
    p_max         = Float64.(th["p_max"])
    c_lin         = Float64.(th["c_lin"])
    c_quad        = Float64.(th["c_quad"])
    c_fixed       = Float64.(th["c_fixed"])
    startup_cost  = Float64.(th["startup_cost"])
    min_up_time   = Int.(th["min_up_time"])
    min_down_time = Int.(th["min_down_time"])
    ramp_up       = Float64.(th["ramp_up"])
    ramp_down     = Float64.(th["ramp_down"])
    ramp_up_str   = Float64.(th["ramp_up_str"])
    ramp_down_str = Float64.(th["ramp_down_str"])
    pt0           = Float64.(th["pt0"])
    storia0       = Int.(th["storia0"])

    # reserve participation factors (optional)
    primary_rho   = haskey(th, "primary_rho")   ? Float64.(th["primary_rho"])   : zeros(n_gen)
    secondary_rho = haskey(th, "secondary_rho") ? Float64.(th["secondary_rho"]) : zeros(n_gen)

    # --- demand: single load = total system demand -----------------------
    # loads.profile is a 1-element list whose entry is the length-n_periods
    # total demand series.
    demand = Float64.(data["loads"]["profile"][1])
    length(demand) == n_periods ||
        error("json2nc4: demand length $(length(demand)) != n_periods $n_periods")

    # --- reserve demand (optional) ---------------------------------------
    has_reserve     = haskey(data, "reserve")
    has_primary     = has_reserve && haskey(data["reserve"], "primary_demand")
    has_secondary   = has_reserve && haskey(data["reserve"], "secondary_demand")
    primary_demand   = has_primary   ? Float64.(data["reserve"]["primary_demand"])   : Float64[]
    secondary_demand = has_secondary ? Float64.(data["reserve"]["secondary_demand"]) : Float64[]

    # --- reactive power (optional) ---------------------------------------
    # A "reactive" block carries the system reactive demand; each thermal unit
    # then offers reactive power within a [min, max] capability band.
    has_reactive    = haskey(data, "reactive")
    reactive_demand = has_reactive ? Float64.(data["reactive"]["reactive_demand"]) : Float64[]
    max_reactive    = has_reactive ? Float64.(th["max_reactive_power"]) : Float64[]
    min_reactive    = has_reactive ? Float64.(th["min_reactive_power"]) : Float64[]
    if has_reactive
        length(reactive_demand) == n_periods ||
            error("json2nc4: reactive demand length $(length(reactive_demand)) " *
                  "!= n_periods $n_periods")
    end

    # =====================================================================
    #  Write netCDF file
    # =====================================================================
    isfile(nc_path) && rm(nc_path)

    NCDataset(nc_path, "c") do ds

        # --- Root: SMS++ Block file header -------------------------------
        ds.attrib["SMS++_file_type"] = Int64(1)   # Block file

        # --- Block_0: the UCBlock ----------------------------------------
        blk = defGroup(ds, "Block_0")
        blk.attrib["type"] = "UCBlock"

        # --- Dimensions --------------------------------------------------
        defDim(blk, "TimeHorizon", n_periods)
        defDim(blk, "NumberUnits", n_gen)
        defDim(blk, "NumberNodes", 1)
        defDim(blk, "NumberElectricalGenerators", n_gen)

        # reserve zones: a single zone each; PrimaryZones/SecondaryZones are
        # NOT written, so by default every node belongs to zone 0.
        if has_primary
            defDim(blk, "NumberPrimaryZones", 1)
        end
        if has_secondary
            defDim(blk, "NumberSecondaryZones", 1)
        end

        # --- GeneratorNode: every generator on node 0 --------------------
        v_gn = defVar(blk, "GeneratorNode", UInt32, ("NumberElectricalGenerators",))
        v_gn[:] = fill(UInt32(0), n_gen)

        # --- ActivePowerDemand(NumberNodes, TimeHorizon) in C order ------
        # NCDatasets is column-major, so list dims reversed and store a
        # (TimeHorizon, NumberNodes) = (n_periods, 1) array.
        v_apd = defVar(blk, "ActivePowerDemand", Float64,
                       ("TimeHorizon", "NumberNodes"))
        v_apd[:, :] = reshape(demand, n_periods, 1)

        # --- PrimaryDemand(NumberPrimaryZones, TimeHorizon) in C order ---
        if has_primary
            v_pd = defVar(blk, "PrimaryDemand", Float64,
                          ("TimeHorizon", "NumberPrimaryZones"))
            v_pd[:, :] = reshape(primary_demand, n_periods, 1)
        end
        if has_secondary
            v_sd = defVar(blk, "SecondaryDemand", Float64,
                          ("TimeHorizon", "NumberSecondaryZones"))
            v_sd[:, :] = reshape(secondary_demand, n_periods, 1)
        end

        # --- ReactivePowerDemand(NumberNodes, TimeHorizon) in C order ----
        if has_reactive
            v_rpd = defVar(blk, "ReactivePowerDemand", Float64,
                           ("TimeHorizon", "NumberNodes"))
            v_rpd[:, :] = reshape(reactive_demand, n_periods, 1)
        end

        # =================================================================
        #  UnitBlock sub-groups (one ThermalUnitBlock per generator)
        # =================================================================
        for g in 1:n_gen
            ug = defGroup(blk, "UnitBlock_$(g-1)")
            ug.attrib["type"] = nuclear ? "NuclearUnitBlock" : "ThermalUnitBlock"

            defVar(ug, "MinPower",      Float64, ())[:] = p_min[g]
            defVar(ug, "MaxPower",      Float64, ())[:] = p_max[g]
            defVar(ug, "LinearTerm",    Float64, ())[:] = c_lin[g]
            defVar(ug, "QuadTerm",      Float64, ())[:] = c_quad[g]
            defVar(ug, "ConstTerm",     Float64, ())[:] = c_fixed[g]
            defVar(ug, "StartUpCost",   Float64, ())[:] = startup_cost[g]
            defVar(ug, "DeltaRampUp",   Float64, ())[:] = ramp_up[g]
            defVar(ug, "DeltaRampDown", Float64, ())[:] = ramp_down[g]
            defVar(ug, "StartUpLimit",  Float64, ())[:] = min(ramp_up_str[g],   p_max[g])
            defVar(ug, "ShutDownLimit", Float64, ())[:] = min(ramp_down_str[g], p_max[g])
            defVar(ug, "MinUpTime",     UInt32,  ())[:] = UInt32(min_up_time[g])
            defVar(ug, "MinDownTime",   UInt32,  ())[:] = UInt32(min_down_time[g])
            defVar(ug, "InitialPower",  Float64, ())[:] = pt0[g]
            defVar(ug, "InitUpDownTime", Int32,  ())[:] = Int32(storia0[g])

            # reserve participation factors (scalar, constant over time).
            # Written only when the corresponding reserve zone exists AND the
            # unit actually offers that reserve (rho > 0).
            if has_primary && primary_rho[g] > 0.0
                defVar(ug, "PrimaryRho", Float64, ())[:] = primary_rho[g]
            end
            if has_secondary && secondary_rho[g] > 0.0
                defVar(ug, "SecondaryRho", Float64, ())[:] = secondary_rho[g]
            end

            # reactive-power capability, gated by the commitment: a stopped
            # unit produces essentially no reactive power, so the json capability
            # band is the on-state coefficient (^on) and the off-state box is
            # {0}. This realises the state-dependent bound Qoff + Qon*u_t of
            # ThermalUnitBlock; the ^off terms default to 0 (unit off => q = 0).
            if has_reactive
                defVar(ug, "MaxReactivePowerOn", Float64, ())[:] = max_reactive[g]
                defVar(ug, "MinReactivePowerOn", Float64, ())[:] = min_reactive[g]
            end

            # nuclear modulation data (only when emitting NuclearUnitBlocks).
            # Constant modulation ramps = mod_frac * the thermal ramps; the unit
            # is free to modulate from the start (InitModulation = ModulationTime).
            # The time-varying incentive to ramp comes from the UCBlock demand.
            if nuclear
                defVar(ug, "ModulationTime",  UInt32, ())[:] = UInt32(mod_time)
                defVar(ug, "InitModulation",  UInt32, ())[:] = UInt32(mod_time)
                defVar(ug, "ModulationDeltaRampUp",   Float64, ())[:] =
                    mod_frac * ramp_up[g]
                defVar(ug, "ModulationDeltaRampDown", Float64, ())[:] =
                    mod_frac * ramp_down[g]
            end
        end

        # --- Informational metadata on Block_0 ---------------------------
        blk.attrib["source_json"] = basename(json_path)
        blk.attrib["n_gen"]       = n_gen
        blk.attrib["n_periods"]   = n_periods
        blk.attrib["dt_hours"]    = Float64(get(meta, "dt_hours", 1.0))
        haskey(meta, "horizon") && (blk.attrib["horizon"] = meta["horizon"])
        haskey(meta, "seed")    && (blk.attrib["seed"]    = meta["seed"])
    end  # NCDataset

    println("Written: $nc_path")
    println("  Generators: $n_gen, Periods: $n_periods, " *
            "Primary reserve: $(has_primary ? "yes" : "no"), " *
            "Secondary reserve: $(has_secondary ? "yes" : "no"), " *
            "Reactive power: $(has_reactive ? "yes" : "no")")
end


# =========================================================================
#  Standalone single-unit ThermalUnitBlock instances (for testing TUDPS)
# =========================================================================

"""
    emit_thermal_single_tubs(json_path, outdir; n_units)

For the first `n_units` thermal generators of the JSON instance, write a
*standalone* ThermalUnitBlock netCDF file (one per unit) directly testable by
the ThermalUnitBlock_Solver tester.

Each file carries the unit's thermal data and a *time-varying* `LinearTerm`
synthesised from the system demand profile so that the unit is incentivised
to follow the load and therefore to ramp and cycle: a standalone unit with
constant costs would just sit at a constant power (or off) and the solvers
under comparison would only meet trivial schedules.

The unconstrained per-period optimum of `a p^2 + b_t p` is `p* = -b_t/(2a)`;
we set `b_t = -2 a * p_target_t` so that `p*` tracks `p_target_t`, a
load-following target swinging between `MinPower` and `MaxPower` with the
(normalised) demand.

Returns the list of written file paths.
"""
function emit_thermal_single_tubs(json_path::String, outdir::String;
                                  n_units::Int = 5)
    data = JSON.parsefile(json_path)
    meta = data["metadata"]
    th   = data["generators"]["thermal"]
    n_gen     = Int(meta["n_gen"])
    n_periods = Int(meta["n_periods"])

    p_min         = Float64.(th["p_min"])
    p_max         = Float64.(th["p_max"])
    c_quad        = Float64.(th["c_quad"])
    startup_cost  = Float64.(th["startup_cost"])
    min_up_time   = Int.(th["min_up_time"])
    min_down_time = Int.(th["min_down_time"])
    ramp_up       = Float64.(th["ramp_up"])
    ramp_down     = Float64.(th["ramp_down"])
    ramp_up_str   = Float64.(th["ramp_up_str"])
    ramp_down_str = Float64.(th["ramp_down_str"])
    pt0           = Float64.(th["pt0"])
    storia0       = Int.(th["storia0"])

    demand = Float64.(data["loads"]["profile"][1])
    dmin, dmax = extrema(demand)
    dn = dmax > dmin ? (demand .- dmin) ./ (dmax - dmin) : fill(0.5, n_periods)

    mkpath(outdir)
    base = splitext(basename(json_path))[1]
    K = min(n_units, n_gen)
    written = String[]

    for g in 1:K
        # ensure a strictly positive quadratic so the vertex p* is well defined
        a = c_quad[g] > 1e-9 ? c_quad[g] : 0.01
        p_target = p_min[g] .+ dn .* (p_max[g] - p_min[g])
        lin = -2.0 .* a .* p_target           # places p* on the load-following target

        nc_path = joinpath(outdir, "$(base)_TUB$(g-1).nc4")
        isfile(nc_path) && rm(nc_path)
        NCDataset(nc_path, "c") do ds
            ds.attrib["SMS++_file_type"] = Int64(1)   # Block file
            blk = defGroup(ds, "Block_0")
            blk.attrib["type"] = "ThermalUnitBlock"

            defDim(blk, "TimeHorizon", n_periods)
            defDim(blk, "NumberIntervals", n_periods)

            defVar(blk, "MinPower",       Float64, ())[:] = p_min[g]
            defVar(blk, "MaxPower",       Float64, ())[:] = p_max[g]
            defVar(blk, "DeltaRampUp",    Float64, ())[:] = ramp_up[g]
            defVar(blk, "DeltaRampDown",  Float64, ())[:] = ramp_down[g]
            defVar(blk, "QuadTerm",       Float64, ())[:] = a
            defVar(blk, "LinearTerm",     Float64, ("NumberIntervals",))[:] = lin
            # ConstTerm = 0 so that being on is *profitable* (the per-period
            # cost at the load-following target is -a*p_target^2 < 0): the
            # unit then wants to stay on and track the demand, so ramps and
            # start-up/shut-down limits actually bind
            defVar(blk, "ConstTerm",      Float64, ())[:] = 0.0
            defVar(blk, "StartUpCost",    Float64, ())[:] = startup_cost[g]
            defVar(blk, "StartUpLimit",   Float64, ())[:] = min(ramp_up_str[g],   p_max[g])
            defVar(blk, "ShutDownLimit",  Float64, ())[:] = min(ramp_down_str[g], p_max[g])
            defVar(blk, "MinUpTime",      UInt32,  ())[:] = UInt32(min_up_time[g])
            defVar(blk, "MinDownTime",    UInt32,  ())[:] = UInt32(min_down_time[g])
            defVar(blk, "InitialPower",   Float64, ())[:] = pt0[g]
            defVar(blk, "InitUpDownTime", Int32,   ())[:] = Int32(storia0[g])

            blk.attrib["source_json"] = basename(json_path)
            blk.attrib["unit"]        = g - 1
            blk.attrib["n_periods"]   = n_periods
        end
        push!(written, nc_path)
        println("Written: $nc_path  (ThermalUnitBlock, $n_periods periods)")
    end
    return written
end

# =========================================================================
#  Standalone single-unit NuclearUnitBlock instances (for testing NUDPS)
# =========================================================================

"""
    emit_nuclear_single_tubs(json_path, outdir; mod_time, mod_frac, n_units)

For the first `n_units` thermal generators of the JSON instance, write a
*standalone* NuclearUnitBlock netCDF file (one per unit) directly testable by
the ThermalUnitBlock_Solver tester (a NuclearUnitBlock is-a ThermalUnitBlock).

Each file carries the unit's thermal data plus the modulation data, and a
*time-varying* `LinearTerm` synthesised from the system demand profile so that
the unit is incentivised to follow the load and therefore to ramp: with the
small modulation ramp `mod_frac * Delta`, big moves require a modulation, which
the window constraint limits. This is what exercises the new DP logic; a
standalone unit with constant costs would just sit at a constant power and never
modulate.

The unconstrained per-period optimum of `a p^2 + b_t p` is `p* = -b_t/(2a)`; we
set `b_t = -2 a * p_target_t` so that `p*` tracks `p_target_t`, a load-following
target swinging between `MinPower` and `MaxPower` with the (normalised) demand.

Returns the list of written file paths.
"""
function emit_nuclear_single_tubs(json_path::String, outdir::String;
                                  mod_time::Int = 8, mod_frac::Float64 = 0.25,
                                  n_units::Int = 5)
    data = JSON.parsefile(json_path)
    meta = data["metadata"]
    th   = data["generators"]["thermal"]
    n_gen     = Int(meta["n_gen"])
    n_periods = Int(meta["n_periods"])

    p_min         = Float64.(th["p_min"])
    p_max         = Float64.(th["p_max"])
    c_lin         = Float64.(th["c_lin"])
    c_quad        = Float64.(th["c_quad"])
    c_fixed       = Float64.(th["c_fixed"])
    startup_cost  = Float64.(th["startup_cost"])
    min_up_time   = Int.(th["min_up_time"])
    min_down_time = Int.(th["min_down_time"])
    ramp_up       = Float64.(th["ramp_up"])
    ramp_down     = Float64.(th["ramp_down"])
    ramp_up_str   = Float64.(th["ramp_up_str"])
    ramp_down_str = Float64.(th["ramp_down_str"])
    pt0           = Float64.(th["pt0"])
    storia0       = Int.(th["storia0"])

    demand = Float64.(data["loads"]["profile"][1])
    dmin, dmax = extrema(demand)
    dn = dmax > dmin ? (demand .- dmin) ./ (dmax - dmin) : fill(0.5, n_periods)

    mkpath(outdir)
    base = splitext(basename(json_path))[1]
    K = min(n_units, n_gen)
    written = String[]

    for g in 1:K
        # ensure a strictly positive quadratic so the vertex p* is well defined
        a = c_quad[g] > 1e-9 ? c_quad[g] : 0.01
        p_target = p_min[g] .+ dn .* (p_max[g] - p_min[g])
        lin = -2.0 .* a .* p_target           # places p* on the load-following target

        nc_path = joinpath(outdir, "$(base)_NUB$(g-1).nc4")
        isfile(nc_path) && rm(nc_path)
        NCDataset(nc_path, "c") do ds
            ds.attrib["SMS++_file_type"] = Int64(1)   # Block file
            blk = defGroup(ds, "Block_0")
            blk.attrib["type"] = "NuclearUnitBlock"

            defDim(blk, "TimeHorizon", n_periods)
            defDim(blk, "NumberIntervals", n_periods)

            defVar(blk, "MinPower",       Float64, ())[:] = p_min[g]
            defVar(blk, "MaxPower",       Float64, ())[:] = p_max[g]
            defVar(blk, "DeltaRampUp",    Float64, ())[:] = ramp_up[g]
            defVar(blk, "DeltaRampDown",  Float64, ())[:] = ramp_down[g]
            defVar(blk, "QuadTerm",       Float64, ())[:] = a
            defVar(blk, "LinearTerm",     Float64, ("NumberIntervals",))[:] = lin
            # ConstTerm = 0 so that being on is *profitable* (the per-period cost
            # at the load-following target is -a*p_target^2 < 0): the unit then
            # wants to stay on and track the demand, which makes the modulation
            # constraints bind. With the unit's real (positive) fixed cost it
            # would just stay off and the test would not exercise modulation.
            defVar(blk, "ConstTerm",      Float64, ())[:] = 0.0
            defVar(blk, "StartUpCost",    Float64, ())[:] = startup_cost[g]
            defVar(blk, "StartUpLimit",   Float64, ())[:] = min(ramp_up_str[g],   p_max[g])
            defVar(blk, "ShutDownLimit",  Float64, ())[:] = min(ramp_down_str[g], p_max[g])
            defVar(blk, "MinUpTime",      UInt32,  ())[:] = UInt32(min_up_time[g])
            defVar(blk, "MinDownTime",    UInt32,  ())[:] = UInt32(min_down_time[g])
            defVar(blk, "InitialPower",   Float64, ())[:] = pt0[g]
            defVar(blk, "InitUpDownTime", Int32,   ())[:] = Int32(storia0[g])

            # modulation data
            defVar(blk, "ModulationTime", UInt32, ())[:] = UInt32(mod_time)
            defVar(blk, "InitModulation", UInt32, ())[:] = UInt32(mod_time)
            defVar(blk, "ModulationDeltaRampUp",   Float64, ())[:] = mod_frac * ramp_up[g]
            defVar(blk, "ModulationDeltaRampDown", Float64, ())[:] = mod_frac * ramp_down[g]

            blk.attrib["source_json"] = basename(json_path)
            blk.attrib["unit"]        = g - 1
            blk.attrib["n_periods"]   = n_periods
            blk.attrib["modulation_time"] = mod_time
            blk.attrib["modulation_frac"] = mod_frac
        end
        push!(written, nc_path)
        println("Written: $nc_path  (NuclearUnitBlock, $n_periods periods, " *
                "ModulationTime=$mod_time, mod ramp=$(round(mod_frac, digits=3))*Delta)")
    end
    return written
end

# =========================================================================
#  CLI entry point
# =========================================================================

const USAGE = """
Usage:
  julia json2nc4.jl <input.json> [output.nc4]            # thermal UCBlock (default)
  julia json2nc4.jl <input.json> [output.nc4] --nuclear  # UCBlock of NuclearUnitBlocks
  julia json2nc4.jl <input.json> --nuclear-single <dir>  # standalone NuclearUnitBlock
                                                         # files, one per unit
  julia json2nc4.jl <input.json> --thermal-single <dir>  # standalone ThermalUnitBlock
                                                         # files, one per unit
Options (nuclear modes):
  --mod-time N    modulation interval tau^M (>=2), default 8
  --mod-frac F    modulation ramp = F * thermal ramp (0..1), default 0.25
  --n-units K     (single modes only) number of units to emit, default 5
"""

function main()
    if length(ARGS) < 1
        print(stderr, USAGE)
        exit(1)
    end

    # split flags from positional arguments
    positional = String[]
    nuclear = false
    single_dir = ""
    thermal_single_dir = ""
    mod_time = 8
    mod_frac = 0.25
    n_units = 5
    i = 1
    while i <= length(ARGS)
        a = ARGS[i]
        if a == "--nuclear"
            nuclear = true
        elseif a == "--nuclear-single"
            single_dir = ARGS[i += 1]
        elseif a == "--thermal-single"
            thermal_single_dir = ARGS[i += 1]
        elseif a == "--mod-time"
            mod_time = parse(Int, ARGS[i += 1])
        elseif a == "--mod-frac"
            mod_frac = parse(Float64, ARGS[i += 1])
        elseif a == "--n-units"
            n_units = parse(Int, ARGS[i += 1])
        elseif startswith(a, "--")
            println(stderr, "Unknown option: $a"); print(stderr, USAGE); exit(1)
        else
            push!(positional, a)
        end
        i += 1
    end

    if isempty(positional)
        print(stderr, USAGE); exit(1)
    end
    json_path = positional[1]
    if !isfile(json_path)
        println(stderr, "Error: file not found: $json_path")
        exit(1)
    end

    # mode 1a: standalone single-unit ThermalUnitBlock files
    if !isempty(thermal_single_dir)
        emit_thermal_single_tubs(json_path, thermal_single_dir;
                                 n_units = n_units)
        return
    end

    # mode 1: standalone single-unit NuclearUnitBlock files
    if !isempty(single_dir)
        emit_nuclear_single_tubs(json_path, single_dir;
                                 mod_time = mod_time, mod_frac = mod_frac,
                                 n_units = n_units)
        return
    end

    # default / --nuclear: one UCBlock netCDF
    nc_path = length(positional) >= 2 ? positional[2] :
              (replace(json_path, r"\.json$" => ".nc4") == json_path ?
               json_path * ".nc4" : replace(json_path, r"\.json$" => ".nc4"))

    convert_json_to_nc4(json_path, nc_path;
                        nuclear = nuclear, mod_time = mod_time,
                        mod_frac = mod_frac)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
