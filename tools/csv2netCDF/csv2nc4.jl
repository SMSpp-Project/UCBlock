using Pkg
Pkg.activate(".")  # Activate environment from Project.toml
Pkg.instantiate()

using YAML
# the official repo, i.e., https://github.com/JuliaGeo/NetCDF.jl,
# does not support (yet) the concept of group :(
using NCDatasets
using DataStructures

using Parameters
using DataFrames
using XLSX
using JLD2
using CSV

using Distributions
using PointEstimateMethod

using StochasticPrograms

using Random

# include additional useful functions, i.e., main type definitions and read data
include("utils.jl")

# Define the scenario
include("scenario_definition.jl")

# Include the samplers for long period uncertainty
include("pem_extraction.jl")

# Include the sampler for distributions associated to short period uncertainty and a function to generate scenarios
include("scen_eps_sampler.jl")

# setting the seed
Random.seed!(123)


# YAML layout assumed by this driver (the only supported format):
#   * EC-wide profiles (`time_res`, `energy_weight`, `reward_price`,
#     `peak_categories`) under `general.profile`.
#   * Pricing fields (`buy_price`, `sell_price`, `consumption_price`,
#     `peak_tariff`, `peak_weight`) under per-tariff blocks selected by
#     `users.<u>.tariff_name`.

@inline ec_profile(name) = profile(gen_data, name)
@inline ec_profile_d(name, default) = profile_d(gen_data, name, default)
@inline user_market_data(u) = field(market_data, field(users_data[u], "tariff_name"))


function csvEC2nc4(
    deterministic::Bool=false,
    sampled_scenarios::Union{Nothing,Vector{Scenario_Load_Renewable}}=nothing,
)

    middle = "_"
    if occursin("_CO", file_name)
        middle = string(middle, "CO_")
    elseif occursin("_NA", file_name)
        middle = string(middle, "NA_")
    elseif occursin("_NC", file_name)
        middle = string(middle, "NC_")
    end

    last = ""
    if "--with-thermal-blocks" in OPTION_ARGS && !occursin("_NA", file_name)
        last = string(last, "_TUB")
    end
    if "--with-network-blocks" in OPTION_ARGS
        last = string(last, "_NB")
    end

    # The mode "c" stands for creating a new file (clobber)
    ds = NCDataset(string("../../data/nc4/EC_Data/EC", middle, "Test", last, ".nc4"), "c", attrib=OrderedDict("SMS++_file_type" => 1))
    block = defGroup(ds, "Block_0", attrib=OrderedDict("id" => "0", "type" => "UCBlock"))

    # Store the number of nodes
    n_users = length(user_set)
    defDim(block, "NumberNodes", n_users)

    # Store the number of time steps/horizons
    defDim(block, "TimeHorizon", n_steps)

    # Store the number of `ECNetworkBlock`(s), i.e., the number of peak periods/categories
    peak_categories = ec_profile("peak_categories")[time_set]
    peak_set = unique(peak_categories)
    n_peaks = length(peak_set)
    defDim(block, "NumberNetworks", n_peaks)

    # Create buy, sell, reward, and consumption, i.e., the constant term, price data arrays.
    # UCBlock applies a single scalar BuyPrice/SellPrice/PeakTariff to the community-wide
    # imports/exports/peaks, so the cost is `price · Σ_u flow_u`. To stay consistent with
    # EnergyCommunity.jl — which computes Σ_u price_u · flow_u with the same per-user data —
    # we read the price ONCE from a representative tariff (all `tariff_name` blocks point to
    # the same CSV column in the current setup); summing across users would scale every
    # cost term by the number of users.
    project_lifetime = field(gen_data, "project_lifetime")
    year_set = 1:project_lifetime
    discount_factor = sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)
    energy_weight_profile = ec_profile("energy_weight")
    time_res_profile = ec_profile("time_res")
    ref_market = user_market_data(first(user_set))

    # `BuyPrice`, i.e., the tariff that user pay to buy electricity at each time horizon
    buy_price_data = [profile(ref_market, "buy_price")[t] *
                      energy_weight_profile[t] *
                      time_res_profile[t]
                      for t in time_set] * discount_factor

    # `SellPrice`, i.e., the tariff that user gain to sell electricity at each time horizon
    sell_price_data = [profile(ref_market, "sell_price")[t] *
                       energy_weight_profile[t] *
                       time_res_profile[t]
                       for t in time_set] * discount_factor

    # `RewardPrice`, i.e., the reward awarded to the community
    reward_price_data = [ec_profile_d("reward_price", fill(0.0, n_steps))[t] *
                         energy_weight_profile[t] *
                         time_res_profile[t]
                         for t in time_set] * discount_factor

    # `PeakTariff`, i.e., the peak tariff cost
    peak_tariff_data = [profile(ref_market, "peak_tariff")[w] *
                        profile(ref_market, "peak_weight")[w]
                        for w in peak_set] * discount_factor

    # `ConstantTerm`, i.e., the consumption price applied to total user load
    const_term_data = [profile(ref_market, "consumption_price")[t] *
                       sum(profile_component(users_data[u], l, "load")[t]
                           for u in user_set for l in asset_names(users_data[u], LOAD)) *
                       energy_weight_profile[t] *
                       time_res_profile[t]
                       for t in time_set] * discount_factor

    if (!("--with-network-blocks" in OPTION_ARGS) &&
        allequal(sell_price_data) &&
        allequal(buy_price_data) &&
        allequal(peak_tariff_data) &&
        allequal(reward_price_data))

        n_intervals = [count(x -> x == w, peak_categories) for w in peak_set]
        @assert length(unique(n_intervals)) == 1 "The values of n_intervals are not all equal"
        defDim(block, "NumberIntervals", n_intervals[1])

        # Store the specific classname of the NetworkBlock, i.e., `ECNetworkBlock` and `ECNetworkData`, to
        # inform UCBlock about the specific type of network (since it deals with both transmission and
        # community networks)
        network_block_classname = defVar(block, "NetworkBlockClassname", String, ())
        network_block_classname[1] = "ECNetworkBlock"
        network_data_classname = defVar(block, "NetworkDataClassname", String, ())
        network_data_classname[1] = "ECNetworkData"

        # `ActivePowerDemand`, i.e., the electricity demand of each node/user at each time horizon
        ## A T T E N T I O N: The data is stored in the NetCDF file in the same order as they are
        ## stored in memory. As Julia uses the column-major ordering for arrays, the order of dimensions
        ## will appear reversed when the data is loaded in languages or programs using row-major
        ## ordering such as C/C++, Python/NumPy or the tools ncdump/ncgen.
        ## To store the demand in the correct shape, i.e., NumberNodes x TimeHorizon, we need to store
        ## it transposed, i.e., TimeHorizon x NumberNodes.
        power_demand = defVar(block, "ActivePowerDemand", Float64, ("TimeHorizon", "NumberNodes")) # ("NumberNodes", "TimeHorizon"))
        power_demand[:, :] = [profile_component(users_data[u], "load", "load")[t]
                              for t in time_set, u in user_set] # for u in user_set, t in time_set]

        # `SellPrice`, i.e., the tariff that user gain to sell electricity at each time horizon
        sell_price = defVar(block, "SellPrice", Float64, ())
        sell_price[:] = sell_price_data[1]

        # `BuyPrice`, i.e., the tariff that user pay to buy electricity at each time horizon
        buy_price = defVar(block, "BuyPrice", Float64, ())
        buy_price[:] = buy_price_data[1]

        # `RewardPrice`, i.e., the reward awarded to the community, if any
        if reward_price_data[1] > 0
            reward_price = defVar(block, "RewardPrice", Float64, ())
            reward_price[:] = reward_price_data[1]
        end

        # `PeakTariff`, i.e., the peak tariff cost
        peak_tariff = defVar(block, "PeakTariff", Float64, ())
        peak_tariff[:] = peak_tariff_data[1]

        # `NetworkConstantTerms`, i.e., the constant term of each ECNetworkBlock
        const_term = defVar(block, "NetworkConstantTerms", Float64, ("NumberNetworks",))
        last_t = 1
        for (i_w, w) in enumerate(peak_set)
            last_i = findlast(x -> x == w, peak_categories)
            const_term[i_w] = sum(const_term_data[last_t:last_i])
            n_intervals = count(x -> x == w, peak_categories)
            last_t += n_intervals
        end

    else

        # Store the specific classname of the NetworkData to inform UCBlock
        # about the specific type of network (since it deals with both
        # transmission and community networks)
        network_data_classname = defVar(block, "NetworkDataClassname", String, ())
        network_data_classname[1] = "NetworkData"

        # `NetworkConstantTerms` MUST be written at the top-level UCBlock even when
        # individual `NetworkBlock_n` groups carry their own scalar `ConstantTerm`:
        # `UCBlock::deserialize` (UCBlock.cpp ~408) overwrites every NetworkBlock's
        # constant_term with `v_network_constant_terms[n]`. Without the top-level
        # vector that array is zero-filled and the per-block ConstantTerm is lost.
        nct = defVar(block, "NetworkConstantTerms", Float64, ("NumberNetworks",))
        last_t = 1
        for (i_w, w) in enumerate(peak_set)
            last_i = findlast(x -> x == w, peak_categories)
            nct[i_w] = sum(const_term_data[last_t:last_i])
            n_intervals = count(x -> x == w, peak_categories)
            last_t += n_intervals
        end

        # Create w `ECNetworkBlock`(s) for each peak period/category, each of them span w_t time steps/horizons
        last_t = 1
        for (i_w, w) in enumerate(peak_set)

            ecnb = defGroup(block, "NetworkBlock_$(i_w - 1)", attrib=OrderedDict("type" => "ECNetworkBlock"))

            # `NumberIntervals`, i.e., the number of sub time horizons spanned by each peak period, i.e., an `ECNetworkBlock`
            n_intervals = count(x -> x == w, peak_categories)
            defDim(ecnb, "NumberIntervals", n_intervals)

            # Store the number of nodes in each NetworkBlock
            n_users = length(user_set)
            defDim(ecnb, "NumberNodes", n_users)

            last_i = findlast(x -> x == w, peak_categories)

            # `ActiveDemand`, i.e., the electricity demand of each node/user at each intervals
            ## A T T E N T I O N: The data is stored in the NetCDF file in the same order as they are
            ## stored in memory. As Julia uses the column-major ordering for arrays, the order of dimensions
            ## will appear reversed when the data is loaded in languages or programs using row-major
            ## ordering such as C/C++, Python/NumPy or the tools ncdump/ncgen.
            ## To store the demand in the correct shape, i.e., NumberIntervals x NumberNodes, we need to store
            ## it transposed, i.e., NumberNodes x NumberIntervals.
            power_demand = defVar(ecnb, "ActiveDemand", Float64, ("NumberNodes", "NumberIntervals")) # ("NumberIntervals", "NumberNodes"))
            power_demand[:, :] = [profile_component(users_data[u], "load", "load")[t]
                                  for u in user_set, t in last_t:last_i] # for t in last_t:last_i, u in user_set]

            # `BuyPrice`, i.e., the tariff that user pay to buy electricity at each time horizon
            if (allequal(buy_price_data[last_t:last_i]))
                buy_price = defVar(ecnb, "BuyPrice", Float64, ())
                buy_price[:] = buy_price_data[last_t]
            else
                buy_price = defVar(ecnb, "BuyPrice", Float64, ("NumberIntervals",))
                buy_price[:] = buy_price_data[last_t:last_i]
            end

            # `SellPrice`, i.e., the tariff that user gain to sell electricity at each time horizon
            if (allequal(sell_price_data[last_t:last_i]))
                sell_price = defVar(ecnb, "SellPrice", Float64, ())
                sell_price[:] = sell_price_data[last_t]
            else
                sell_price = defVar(ecnb, "SellPrice", Float64, ("NumberIntervals",))
                sell_price[:] = sell_price_data[last_t:last_i]
            end

            # `RewardPrice`, i.e., the reward awarded to the community...
            if (allequal(reward_price_data[last_t:last_i]))
                if reward_price_data[last_t] > 0 # ... if any
                    reward_price = defVar(ecnb, "RewardPrice", Float64, ())
                    reward_price[:] = reward_price_data[last_t]
                end
            else
                reward_price = defVar(ecnb, "RewardPrice", Float64, ("NumberIntervals",))
                reward_price[:] = reward_price_data[last_t:last_i]
            end

            # `PeakTariff`, i.e., the peak tariff cost
            peak_tariff = defVar(ecnb, "PeakTariff", Float64, ())
            peak_tariff[:] = peak_tariff_data[i_w]

            # `ConstantTerm`, i.e., the consumption price
            const_term = defVar(ecnb, "ConstantTerm", Float64, ())
            const_term[:] = sum(const_term_data[last_t:last_i])

            last_t += n_intervals
        end
    end

    # --------------------------------------------------------------------------------------- #

    # Create g `UnitBlock`(s) for each electrical generator/device

    n_devices = reduce(+, [d != "generator" ? 1 :
                          div(field_component(users_data[u], d, "max_capacity"), field_component(users_data[u], d, "nom_capacity"))
                          for u in user_set
                          for d in asset_names(users_data[u], SMSPP_DEVICES)], init=0)
    # number of UnitBlock
    defDim(block, "NumberUnits", n_devices)

    # AbstractPath
    if !deterministic # stochastic model
        path_dim = 0
        path_group_idx_data = String[]
        path_element_idx_data = Int[]
        path_range_idx_data = Int[]
        # IntermittentUnitBlock(s) whose MaxPower changes per scenario.
        # For each entry we record:
        #   ub_idx       = UnitBlock_<ub_idx> position inside the UCBlock
        #   user, asset  = keys to look up scen.Ren[user][asset]
        #   max_capacity = scaling factor applied to the per-unit profile
        intermittent_units = Tuple{Int,String,String,Float64}[]
    end

    if n_devices > 0

        # each UnitBlock has just one electrical generator
        defDim(block, "NumberElectricalGenerators", n_devices)

        # `GeneratorNode` is a 1D variable that represent the node/user owner
        # of each electrical generator/device
        generator_node = defVar(block, "GeneratorNode", UInt32, ("NumberElectricalGenerators",))

        last_g = 0
        for (i_u, u) in enumerate(user_set)

            for g in asset_names(users_data[u], SMSPP_DEVICES)

                if g in ("PV", "wind")

                    ub = defGroup(block, "UnitBlock_$(last_g)", attrib=OrderedDict("type" => "IntermittentUnitBlock"))

                    # store the maximum installable capacity of the pv/wind asset
                    max_capacity = defVar(ub, "MaxCapacity", Float64, ())
                    max_capacity[:] = field_component(users_data[u], g, "max_capacity")

                    # store the maximum power of the pv/wind asset
                    max_power_data = [field_component(users_data[u], g, "max_capacity") *
                                      profile_component(users_data[u], g, "ren_pu")[t]
                                      for t in time_set]
                    if (allequal(max_power_data))
                        max_power = defVar(ub, "MaxPower", Float64, ())
                        max_power[:] = max_power_data[1]
                    else
                        max_power = defVar(ub, "MaxPower", Float64, ("TimeHorizon",))
                        max_power[:] = max_power_data[:]
                    end

                    # store the Net Present Value of the pv/wind asset
                    investment_cost = defVar(ub, "InvestmentCost", Float64, ())
                    investment_cost[:] = sum(y == 0 ? field_component(users_data[u], g, "CAPEX_lin") : # investment cost of the component
                                             ((field_component(users_data[u], g, "OEM_lin") + # operation and maintenance cost of the component
                                               ((mod(y, field_component(users_data[u], g, "lifetime_y")) == 0 && y != project_lifetime) ?
                                                field_component(users_data[u], g, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                               ((mod(y, field_component(users_data[u], g, "lifetime_y")) != 0 && y == project_lifetime) ?
                                                field_component(users_data[u], g, "CAPEX_lin") *
                                                (1.0 - mod(y, field_component(users_data[u], g, "lifetime_y")) /
                                                       field_component(users_data[u], g, "lifetime_y")) : 0.0)) * # residual value of the component
                                              (1 / (1 + field(gen_data, "d_rate"))^y)) for y in append!([0], year_set)) *
                                         field_component(users_data[u], g, "max_capacity")

                    if !deterministic # stochastic model
                        path_dim += 1
                        append!(path_group_idx_data, [string(last_g), "x_intermittent"]) # i.e., last_g wrt B, V x_intermittent
                        append!(path_element_idx_data, [typemax(UInt32), 0]) # i.e., _ wrt B, 0 wrt V x_intermittent
                        append!(path_range_idx_data, [typemax(UInt32), 1]) # i.e., _ wrt B, 1 or _ wrt V x_intermittent
                        push!(intermittent_units, (last_g, u, g, field_component(users_data[u], g, "max_capacity")))
                    end

                    last_g += 1
                    generator_node[last_g] = i_u - 1 # assign the ownership of the current pv/wind asset to the respective user

                elseif g == "batt"

                    ub = defGroup(block, "UnitBlock_$(last_g)", attrib=OrderedDict("type" => "BatteryUnitBlock"))

                    # ----------- Battery -----------

                    # store the maximum installable capacity of the battery
                    batt_max_capacity = defVar(ub, "BatteryMaxCapacity", Float64, ())
                    batt_max_capacity[:] = field_component(users_data[u], g, "max_capacity")

                    # store the maximum power of the battery
                    batt_max_power = defVar(ub, "MaxPower", Float64, ())
                    batt_max_power[:] = field_component(users_data[u], g, "max_capacity")

                    # store the maximum C-rate of the battery in charge
                    max_C_ch = field_component(users_data[u], g, "max_C_ch")
                    if max_C_ch > 1
                        batt_max_C_ch = defVar(ub, "MaxCRateCharge", Float64, ())
                        batt_max_C_ch[:] = max_C_ch
                    end

                    # store the maximum C-rate of the battery in discharge
                    max_C_dch = field_component(users_data[u], g, "max_C_dch")
                    if max_C_dch > 1
                        batt_max_C_dch = defVar(ub, "MaxCRateDischarge", Float64, ())
                        batt_max_C_dch[:] = max_C_dch
                    end

                    # set a negative initial power negative to use the cyclical notation
                    initial_storage = defVar(ub, "InitialStorage", Float64, ())
                    initial_storage[:] = -1

                    # store the minimum storage of the battery
                    min_storage_data = [field_component(users_data[u], g, "min_SOC") /
                                        time_res_profile[t] # energy (kWh), i.e., power * time, to power (kW), i.e., energy / time
                                        for t in time_set] * field_component(users_data[u], g, "max_capacity")
                    if (allequal(min_storage_data))
                        min_storage = defVar(ub, "MinStorage", Float64, ())
                        min_storage[:] = min_storage_data[1]
                    else
                        min_storage = defVar(ub, "MinStorage", Float64, ("TimeHorizon",))
                        min_storage[:] = min_storage_data[:]
                    end

                    # store the maximum storage of the battery
                    max_storage_data = [field_component(users_data[u], g, "max_SOC") /
                                        time_res_profile[t] # energy (kWh), i.e., power * time, to power (kW), i.e., energy / time
                                        for t in time_set] * field_component(users_data[u], g, "max_capacity")
                    if (allequal(max_storage_data))
                        max_storage = defVar(ub, "MaxStorage", Float64, ())
                        max_storage[:] = max_storage_data[1]
                    else
                        max_storage = defVar(ub, "MaxStorage", Float64, ("TimeHorizon",))
                        max_storage[:] = max_storage_data[:]
                    end

                    # store the Net Present Value of the battery
                    batt_investment_cost = defVar(ub, "BatteryInvestmentCost", Float64, ())
                    batt_investment_cost[:] = (sum(y == 0 ? field_component(users_data[u], g, "CAPEX_lin") : # investment cost of the component
                                                   ((field_component(users_data[u], g, "OEM_lin") + # operation and maintenance cost of the component
                                                     ((mod(y, field_component(users_data[u], g, "lifetime_y")) == 0 && y != project_lifetime) ?
                                                      field_component(users_data[u], g, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                                     ((mod(y, field_component(users_data[u], g, "lifetime_y")) != 0 && y == project_lifetime) ?
                                                      field_component(users_data[u], g, "CAPEX_lin") *
                                                      (1.0 - mod(y, field_component(users_data[u], g, "lifetime_y")) /
                                                             field_component(users_data[u], g, "lifetime_y")) : 0.0)) * # residual value of the component
                                                    (1 / (1 + field(gen_data, "d_rate"))^y)) for y in append!([0], year_set)) *
                                               field_component(users_data[u], g, "max_capacity"))

                    # ---------- Converter ----------

                    g_conv = field_component(users_data[u], g, "corr_asset") # corresponding converter, i.e., "conv"

                    # store the maximum installable capacity of the converter
                    conv_max_capacity = defVar(ub, "ConverterMaxCapacity", Float64, ())
                    conv_max_capacity[:] = field_component(users_data[u], g_conv, "max_capacity")

                    # store the maximum power of the converter
                    conv_max_power = defVar(ub, "ConverterMaxPower", Float64, ())
                    conv_max_power[:] = field_component(users_data[u], g_conv, "max_capacity")

                    # store the intake roundtrip efficiency of the battery
                    intake_coeff = defVar(ub, "ExtractingBatteryRho", Float64, ())
                    intake_coeff[:] = 1 / (sqrt(field_component(users_data[u], g, "eta")) *
                                           field_component(users_data[u], g_conv, "eta")) # corresponding converter, i.e., "conv"

                    # store the outtake roundtrip efficiency of the battery
                    outtake_coeff = defVar(ub, "StoringBatteryRho", Float64, ())
                    outtake_coeff[:] = sqrt(field_component(users_data[u], g, "eta")) *
                                       field_component(users_data[u], g_conv, "eta") # corresponding converter, i.e., "conv"

                    # store the Net Present Value of the converter
                    conv_investment_cost = defVar(ub, "ConverterInvestmentCost", Float64, ())
                    conv_investment_cost[:] = (sum(y == 0 ? field_component(users_data[u], g_conv, "CAPEX_lin") : # investment cost of the component
                                                   ((field_component(users_data[u], g_conv, "OEM_lin") + # operation and maintenance cost of the component
                                                     ((mod(y, field_component(users_data[u], g_conv, "lifetime_y")) == 0 && y != project_lifetime) ?
                                                      field_component(users_data[u], g_conv, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                                     ((mod(y, field_component(users_data[u], g_conv, "lifetime_y")) != 0 && y == project_lifetime) ?
                                                      field_component(users_data[u], g_conv, "CAPEX_lin") *
                                                      (1.0 - mod(y, field_component(users_data[u], g_conv, "lifetime_y")) /
                                                             field_component(users_data[u], g_conv, "lifetime_y")) : 0.0)) * # residual value of the component
                                                    (1 / (1 + field(gen_data, "d_rate"))^y)) for y in append!([0], year_set)) *
                                               field_component(users_data[u], g_conv, "max_capacity"))

                    if !deterministic # stochastic model
                        path_dim += 2
                        append!(path_group_idx_data, [string(last_g), "x_battery", string(last_g), "x_converter"]) # i.e., last_g wrt B, V x_battery, x_converter
                        append!(path_element_idx_data, [typemax(UInt32), 0, typemax(UInt32), 0]) # i.e., _ wrt B, 0 wrt V x_battery, x_converter
                        append!(path_range_idx_data, [typemax(UInt32), 1, typemax(UInt32), 1]) # i.e., _ wrt B, 1 or _ wrt V x_battery, x_converter
                    end

                    last_g += 1
                    generator_node[last_g] = i_u - 1 # assign the ownership of the current battery to the respective user

                elseif g == "generator"

                    for _ in 1:div(field_component(users_data[u], g, "max_capacity"), field_component(users_data[u], g, "nom_capacity"))

                        ub = defGroup(block, "UnitBlock_$(last_g)", attrib=OrderedDict("type" => "ThermalUnitBlock"))

                        # store the installable capacity of the thermal
                        thermal_capacity = defVar(ub, "Capacity", Float64, ())
                        thermal_capacity[:] = field_component(users_data[u], g, "nom_capacity")

                        # store the minimum power of the thermal
                        thermal_min_power = defVar(ub, "MinPower", Float64, ())
                        thermal_min_power[:] = (field_component(users_data[u], g, "min_technical") *
                                                field_component(users_data[u], g, "nom_capacity"))

                        # store the maximum power of the thermal
                        thermal_max_power = defVar(ub, "MaxPower", Float64, ())
                        thermal_max_power_val =
                            field_component(users_data[u], g, "max_technical") *
                            field_component(users_data[u], g, "nom_capacity")
                        thermal_max_power[:] = thermal_max_power_val

                        # store the start-up limit
                        thermal_start_up_limit = defVar(ub, "StartUpLimit", Float64, ())
                        thermal_start_up_limit[:] = thermal_max_power_val

                        # store the shut-down limit
                        thermal_shut_up_limit = defVar(ub, "ShutDownLimit", Float64, ())
                        thermal_shut_up_limit[:] = thermal_max_power_val

                        # store the Net Present Value of the thermal
                        investment_cost = defVar(ub, "InvestmentCost", Float64, ())
                        investment_cost[:] = sum(y == 0 ? field_component(users_data[u], g, "CAPEX_lin") : # investment cost of the component
                                                 ((((mod(y, field_component(users_data[u], g, "lifetime_y")) == 0 && y != project_lifetime) ?
                                                    field_component(users_data[u], g, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                                   ((mod(y, field_component(users_data[u], g, "lifetime_y")) != 0 && y == project_lifetime) ?
                                                    field_component(users_data[u], g, "CAPEX_lin") *
                                                    (1.0 - mod(y, field_component(users_data[u], g, "lifetime_y")) /
                                                           field_component(users_data[u], g, "lifetime_y")) : 0.0)) * # residual value of the component
                                                  (1 / (1 + field(gen_data, "d_rate"))^y)) for y in append!([0], year_set)) *
                                             field_component(users_data[u], g, "nom_capacity")

                        # store the linear term of the thermal
                        linear_term_data = sum([(field_component(users_data[u], g, "fuel_price") * # fuel consumption wrt the slope of the piece-wise linear cost function
                                                 field_component(users_data[u], g, "slope_map")) *
                                                energy_weight_profile[t] *
                                                time_res_profile[t]
                                                for t in time_set] *
                                               (1 / (1 + field(gen_data, "d_rate"))^y) for y in year_set)
                        if (allequal(linear_term_data))
                            linear_term = defVar(ub, "LinearTerm", Float64, ())
                            linear_term[:] = linear_term_data[1]
                        else
                            linear_term = defVar(ub, "LinearTerm", Float64, ("TimeHorizon",))
                            linear_term[:] = linear_term_data[:]
                        end

                        # store the constant term of the thermal
                        const_term_data = sum([(field_component(users_data[u], g, "OEM_lin") + # operation and maintenance cost of the component
                                                (field_component(users_data[u], g, "fuel_price") * # fuel consumption wrt the intercept of the piece-wise linear cost function
                                                 field_component(users_data[u], g, "inter_map"))) *
                                               energy_weight_profile[t] *
                                               time_res_profile[t]
                                               for t in time_set] *
                                              (1 / (1 + field(gen_data, "d_rate"))^y) for y in year_set) *
                                          field_component(users_data[u], g, "nom_capacity")
                        if (allequal(const_term_data))
                            const_term = defVar(ub, "ConstTerm", Float64, ())
                            const_term[:] = const_term_data[1]
                        else
                            const_term = defVar(ub, "ConstTerm", Float64, ("TimeHorizon",))
                            const_term[:] = const_term_data[:]
                        end

                        if !deterministic # stochastic model
                            path_dim += 1
                            append!(path_group_idx_data, [string(last_g), "x_thermal"]) # i.e., last_g wrt B, V x_thermal
                            append!(path_element_idx_data, [typemax(UInt32), 0]) # i.e., _ wrt B, 0 wrt V x_thermal
                            append!(path_range_idx_data, [typemax(UInt32), 1]) # i.e., _ wrt B, 1 or _ wrt V x_thermal
                        end

                        last_g += 1
                        generator_node[last_g] = i_u - 1 # assign the ownership of the current therms generator to the respective user
                    end
                end
            end
        end
    end

    close(ds)

    if !deterministic # stochastic model

        # The mode "c" stands for creating a new file (clobber).
        # The TSSB instance only references the deterministic UCBlock written above
        # via the `filename` attribute on its inner Block — the deterministic block
        # itself is NOT embedded inline.
        tssb_ds = NCDataset(string("../../data/nc4/EC_Data/TSSB_EC", middle, "Test", last, ".nc4"), "c", attrib=OrderedDict("SMS++_file_type" => 1))
        tssb = defGroup(tssb_ds, "Block_0", attrib=OrderedDict("id" => "0", "type" => "TwoStageStochasticBlock"))

        ## Number of scenarios in the TwoStageStochasticBlock.
        ## We use the number of sampled_scenarios, which already encodes
        ## the (s, eps) combinations returned by scenarios_generator.
        n_scen = length(sampled_scenarios)
        defDim(tssb, "NumberScenarios", n_scen)

        # ----------------------------------------------------------------
        # Stochastic-price detection.
        # ----------------------------------------------------------------
        # `scen_eps_sampler.jl` perturbs market-level prices whenever the
        # YAML market profile defines `std_<name>` (`std_buy_price`,
        # `std_sell_price`, `std_consumption_price`, `std_penalty_price`,
        # `std_peak_tariff`); without those entries the corresponding price
        # is left deterministic. We scan the sampled scenarios to find
        # which price fields actually vary, then emit one
        # `SimpleDataMapping` per (varying field, peak period) targeting the
        # matching ECNetworkBlock setter (registered in
        # `ECNetworkBlock::static_initialization`):
        #   buy_price          -> ECNetworkBlock::set_buy_price
        #   sell_price         -> ECNetworkBlock::set_sell_price
        #   peak_tariff        -> ECNetworkBlock::set_peak_tariff   (scalar)
        #   consumption_price  -> ECNetworkBlock::set_const_term    (scalar)
        #   penalty_price      -> ECNetworkBlock::set_penalty_price
        function _scenario_field_varies(field_extractor)
            isempty(sampled_scenarios) && return false
            ref = field_extractor(sampled_scenarios[1])
            return any(s -> field_extractor(s) != ref, sampled_scenarios)
        end
        varying_price_fields = String[]
        _scenario_field_varies(s -> s.buy_price)         && push!(varying_price_fields, "buy_price")
        _scenario_field_varies(s -> s.sell_price)        && push!(varying_price_fields, "sell_price")
        _scenario_field_varies(s -> s.consumption_price) && push!(varying_price_fields, "consumption_price")
        _scenario_field_varies(s -> s.penalty_price)     && push!(varying_price_fields, "penalty_price")
        _scenario_field_varies(s -> s.peak_tariff)       && push!(varying_price_fields, "peak_tariff")

        # peak_tariff and consumption_price feed scalar setters (slice length 1);
        # the other three feed per-time vectors of length equal to the peak's
        # number of intervals.
        ec_setter_for = Dict(
            "buy_price"         => "ECNetworkBlock::set_buy_price",
            "sell_price"        => "ECNetworkBlock::set_sell_price",
            "peak_tariff"       => "ECNetworkBlock::set_peak_tariff",
            "consumption_price" => "ECNetworkBlock::set_const_term",
            "penalty_price"     => "ECNetworkBlock::set_penalty_price",
        )
        scalar_price_field(name) = name == "peak_tariff" || name == "consumption_price"

        n_intervals_per_peak = [count(x -> x == w, peak_categories) for w in peak_set]

        # price_mappings[k] = (field_name, peak_index_1based, slice_length)
        price_mappings = Tuple{String,Int,Int}[]
        for field_name in varying_price_fields
            for i_w in 1:n_peaks
                len = scalar_price_field(field_name) ? 1 : n_intervals_per_peak[i_w]
                push!(price_mappings, (field_name, i_w, len))
            end
        end
        N_price_tail = isempty(price_mappings) ? 0 : sum(m[3] for m in price_mappings)
        # ----------------------------------------------------------------

        # DiscreteScenarioSet
        #
        # This group will be read by DiscreteScenarioSet::deserialize().
        # Expected layout in C++ (row-major):
        #
        #   dim NumberScenarios
        #   dim ScenarioSize
        #   var Scenarios(NumberScenarios, ScenarioSize)
        #   var PoolWeights(NumberScenarios)
        #
        # As Julia stores arrays in column-major order, the data is written
        # as (ScenarioSize, NumberScenarios) so that C++ will see it as
        # [NumberScenarios][ScenarioSize].
        #
        dss = defGroup(
            tssb,
            "DiscreteScenarioSet",
            attrib=OrderedDict("type" => "DiscreteScenarioSet"),
        )

        # ScenarioSize = number of entries in each scenario vector.
        # Layout (in scenario-vector order):
        #   1. ActivePowerDemand: n_steps * n_users entries, flattened (t, u)
        #      consistent with UCBlock::ActivePowerDemand[t, u].
        #   2. For each IntermittentUnitBlock (PV / wind), in the order they
        #      were emitted into the UCBlock: n_steps entries with
        #      max_capacity * scen.Ren[user][asset][t]. These feed
        #      IntermittentUnitBlock::set_maximum_power on the corresponding
        #      UnitBlock_<ub_idx>.
        #   3. Price tail (only when `varying_price_fields` is non-empty):
        #      one slice per (varying field, peak period) ordered as in
        #      `price_mappings`, of length 1 for the scalar setters
        #      (peak_tariff, consumption_price) and `n_intervals_per_peak[i_w]`
        #      for the vector setters (buy_price, sell_price, penalty_price).
        n_users = length(user_set)
        n_intermittent = length(intermittent_units)
        N_dem = n_steps * n_users
        N_mp  = n_steps
        scenario_size = N_dem + n_intermittent * N_mp + N_price_tail

        defDim(dss, "NumberScenarios", n_scen)
        defDim(dss, "ScenarioSize", scenario_size)

        ## A T T E N T I O N: The data is stored in the NetCDF file in the
        ## same order as they are stored in memory. As Julia uses the
        ## column-major ordering for arrays, the order of dimensions
        ## will appear reversed when the data is loaded in languages or
        ## programs using row-major ordering such as C/C++, Python/NumPy
        ## or the tools ncdump/ncgen.
        ## To store the scenario set in the correct shape, i.e.,
        ## NumberScenarios x ScenarioSize in C++, we store it here as
        ## ScenarioSize x NumberScenarios in Julia.
        scen_mat = Array{Float64}(undef, scenario_size, n_scen)
        weights = Array{Float64}(undef, n_scen)

        for (k, scen) in enumerate(sampled_scenarios)
            vec = Array{Float64}(undef, scenario_size)
            idx = 1

            # Section 1: ActivePowerDemand in (time, user) order, consistent
            # with how ActivePowerDemand is written as [t, u].
            for t in time_set
                for u in user_set
                    vec[idx] = scen.Load[u][t]
                    idx += 1
                end
            end

            # Section 2: per-IntermittentUnitBlock max_power time series.
            for (_, u, asset, max_cap) in intermittent_units
                ren_profile = scen.Ren[u][asset]
                for t in time_set
                    vec[idx] = max_cap * ren_profile[t]
                    idx += 1
                end
            end

            # Section 3: per-(price_field, peak) tail. For each varying price
            # field we lay one slice per peak period, in the same order as
            # `price_mappings`. The aggregation mirrors the deterministic
            # scaling (energy_weight, time_res, peak_weight, discount_factor)
            # so the C++ obj coefficient matches the deterministic case when
            # the perturbation amplitude is zero.
            if !isempty(price_mappings)
                last_t_per_peak = let lt = 1, out = Int[]
                    for nint in n_intervals_per_peak
                        push!(out, lt)
                        lt += nint
                    end
                    out
                end
                for (field_name, i_w, _len) in price_mappings
                    w = peak_set[i_w]
                    last_t = last_t_per_peak[i_w]
                    last_i = last_t + n_intervals_per_peak[i_w] - 1
                    if field_name == "buy_price"
                        for tt in last_t:last_i
                            t = time_set[tt]
                            vec[idx] = scen.buy_price[t] *
                                       energy_weight_profile[t] *
                                       time_res_profile[t] * discount_factor
                            idx += 1
                        end
                    elseif field_name == "sell_price"
                        for tt in last_t:last_i
                            t = time_set[tt]
                            vec[idx] = scen.sell_price[t] *
                                       energy_weight_profile[t] *
                                       time_res_profile[t] * discount_factor
                            idx += 1
                        end
                    elseif field_name == "penalty_price"
                        for tt in last_t:last_i
                            t = time_set[tt]
                            vec[idx] = scen.penalty_price[t] *
                                       energy_weight_profile[t] *
                                       time_res_profile[t] * discount_factor
                            idx += 1
                        end
                    elseif field_name == "consumption_price"
                        # ConstantTerm aggregates Σ_t consumption_price[t] · Σ_u Load_u[t]
                        # over the peak's time window (mirrors deterministic const_term_data).
                        acc = 0.0
                        for tt in last_t:last_i
                            t = time_set[tt]
                            acc += scen.consumption_price[t] *
                                   sum(scen.Load[u][t] for u in user_set) *
                                   energy_weight_profile[t] *
                                   time_res_profile[t]
                        end
                        vec[idx] = acc * discount_factor
                        idx += 1
                    elseif field_name == "peak_tariff"
                        vec[idx] = scen.peak_tariff[w] *
                                   profile(ref_market, "peak_weight")[w] *
                                   discount_factor
                        idx += 1
                    else
                        error("Unhandled stochastic price field: $field_name")
                    end
                end
            end

            scen_mat[:, k] = vec
            weights[k] = probability(scen)
        end

        # Scenarios: stored as (ScenarioSize, NumberScenarios) in Julia
        # so that C++ sees [NumberScenarios][ScenarioSize].
        scenarios_var = defVar(
            dss,
            "Scenarios",
            Float64,
            ("ScenarioSize", "NumberScenarios"),
        )
        scenarios_var[:, :] = scen_mat

        # Scenario weights (probabilities)
        pool_weights_var = defVar(
            dss,
            "PoolWeights",
            Float64,
            ("NumberScenarios",),
        )
        pool_weights_var[:] = weights

        # AbstractPath
        ap = defGroup(tssb, "StaticAbstractPath")

        defDim(ap, "PathDim", path_dim)

        path_length = 2 # 1 B (UnitBlock_*) + 1 V (x_design) for each path
        total_length = path_length * path_dim
        defDim(ap, "TotalLength", total_length)

        path_start = defVar(ap, "PathStart", UInt32, ("PathDim",))
        path_start[:] = collect(0:path_length:total_length-1)[:] # range from 0 to total_length each path_length

        path_node_types = defVar(ap, "PathNodeTypes", Char, ("TotalLength",))
        path_node_types[:] = collect("BV"^path_dim)[:] # repeat BV path_dim times

        path_group_idx = defVar(ap, "PathGroupIndices", String, ("TotalLength",))
        path_group_idx[:] = path_group_idx_data[:]

        path_element_idx = defVar(ap, "PathElementIndices", UInt32, ("TotalLength",))
        path_element_idx[:] = path_element_idx_data[:]

        path_range_idx = defVar(ap, "PathRangeIndices", UInt32, ("TotalLength",))
        path_range_idx[:] = path_range_idx_data[:]

        # StochasticBlock
        sb = defGroup(tssb, "StochasticBlock", attrib=OrderedDict("type" => "StochasticBlock"))

        # SimpleDataMapping
        #
        # One mapping per stochastic quantity. The scenario vector is split as
        # described above: Section 1 (demand), then one slice per intermittent
        # unit (Section 2), then one slice per (varying price field, peak
        # period) (Section 3, only when price perturbations are detected).
        # Each mapping declares:
        #
        #   * which C++ setter consumes the slice
        #   * SetSize=(0,0) i.e. Range/Range mode
        #   * SetElements=[fromStart, fromEnd, toStart, toEnd] picking the
        #     scenario slice and pushing it into the target's full-length
        #     argument [0, N_target).
        #
        # Mapping 0 targets the inner UCBlock itself (empty AbstractPath);
        # mappings 1..n_intermittent target UnitBlock_<ub_idx> via a single
        # "B" hop carrying the UnitBlock index; the price mappings target
        # NetworkBlock_<i_w-1> via a single "B" hop carrying the NetworkBlock
        # index (= n_devices + (i_w-1) in UCBlock's sub-Block ordering).

        number_mappings = 1 + n_intermittent + length(price_mappings)

        defDim(sb, "NumberDataMappings", number_mappings)
        defDim(sb, "SetSize_dim", 2 * number_mappings)
        defDim(sb, "SetElements_dim", 4 * number_mappings)

        v_FunctionName = defVar(sb, "FunctionName", String, ("NumberDataMappings",))
        v_DataType     = defVar(sb, "DataType",     Char,   ("NumberDataMappings",))
        v_Caller       = defVar(sb, "Caller",       Char,   ("NumberDataMappings",))
        v_SetSize      = defVar(sb, "SetSize",      UInt32, ("SetSize_dim",))
        v_SetElements  = defVar(sb, "SetElements",  UInt32, ("SetElements_dim",))

        function_names = String["UCBlock::set_active_power_demand";
                                fill("IntermittentUnitBlock::set_maximum_power", n_intermittent)]
        for (field_name, _i_w, _len) in price_mappings
            push!(function_names, ec_setter_for[field_name])
        end
        v_FunctionName[:] = function_names
        v_DataType[:]     = fill('D', number_mappings)
        v_Caller[:]       = fill('B', number_mappings)

        # Range/Range for every mapping
        v_SetSize[:] = fill(UInt32(0), 2 * number_mappings)

        set_elements = UInt32[]
        # Mapping 0: scenario [0, N_dem) -> demand argument [0, N_dem)
        append!(set_elements, UInt32[0, N_dem, 0, N_dem])
        # Mapping i: scenario slice for the i-th intermittent unit -> [0, T)
        for i in 1:n_intermittent
            offset = N_dem + (i - 1) * N_mp
            append!(set_elements, UInt32[offset, offset + N_mp, 0, N_mp])
        end
        # Mappings for the per-(price_field, peak) tail.
        let tail_offset = N_dem + n_intermittent * N_mp
            for (_field_name, _i_w, len) in price_mappings
                append!(set_elements, UInt32[tail_offset, tail_offset + len, 0, len])
                tail_offset += len
            end
        end
        v_SetElements[:] = set_elements

        # AbstractPath nested in StochasticBlock: tells the deserializer how
        # to navigate from the inner Block (the loaded UCBlock) to each setter
        # target. Empty path = the UCBlock itself; "B" + UInt32 index = enter
        # the sub-Block at that group position. UCBlock orders its sub-blocks
        # as [UnitBlock_0..n_devices-1, NetworkBlock_0..n_peaks-1], so the
        # NetworkBlock for peak (i_w-1) lives at index n_devices + (i_w-1).
        ap = defGroup(sb, "AbstractPath", attrib=OrderedDict("type" => "AbstractPath"))

        # Concatenated steps across all paths: one 'B' per max_power mapping,
        # plus one 'B' per (price_field, peak) mapping.
        total_length_inner = n_intermittent + length(price_mappings)

        defDim(ap, "PathDim", number_mappings)
        defDim(ap, "TotalLength", total_length_inner)

        v_PathStart        = defVar(ap, "PathStart",         UInt32, ("PathDim",))
        v_PathNodeTypes    = defVar(ap, "PathNodeTypes",     Char,   ("TotalLength",))
        v_PathGroupIndices = defVar(ap, "PathGroupIndices",  UInt32, ("TotalLength",))
        v_PathElementIdx   = defVar(ap, "PathElementIndices", UInt32, ("TotalLength",))
        v_PathRangeIdx     = defVar(ap, "PathRangeIndices",   UInt32, ("TotalLength",))

        # PathStart[k] is the start position in the concatenated TotalLength
        # array for the k-th mapping; the k-th path covers indices
        # [PathStart[k], PathStart[k+1]) (with the last path running to the
        # end). Mapping 0 (demand) has length 0 (empty path); each subsequent
        # mapping (intermittent + price) is a single 'B' hop.
        path_starts_inner = UInt32[0]
        # n_intermittent + length(price_mappings) single-step paths follow
        for i in 1:(n_intermittent + length(price_mappings))
            push!(path_starts_inner, UInt32(i - 1))
        end
        v_PathStart[:] = path_starts_inner

        if total_length_inner > 0
            node_types = Char[]
            group_indices = UInt32[]
            for unit in intermittent_units
                push!(node_types, 'B')
                push!(group_indices, UInt32(unit[1]))
            end
            for (_field_name, i_w, _len) in price_mappings
                push!(node_types, 'B')
                push!(group_indices, UInt32(n_devices + (i_w - 1)))
            end
            v_PathNodeTypes[:]    = node_types
            v_PathGroupIndices[:] = group_indices
            v_PathElementIdx[:]   = fill(typemax(UInt32), total_length_inner)
            v_PathRangeIdx[:]     = fill(typemax(UInt32), total_length_inner)
        end

        # The deterministic UCBlock is referenced via filename — not embedded inline.
        defGroup(sb, "Block", attrib=OrderedDict("id" => "0", "filename" => string("EC", middle, "Test", last, ".nc4[0]")))

        close(tssb_ds)
    end
end

## Parameters

@assert 0 <= length(ARGS) <= 3

NO_OPTION_ARGS = filter(arg -> !startswith(arg, "--"), ARGS)
@assert 0 <= length(NO_OPTION_ARGS) <= 1

OPTION_ARGS = setdiff(ARGS, NO_OPTION_ARGS)
@assert 0 <= length(OPTION_ARGS) <= 2
@assert issubset(OPTION_ARGS, ["--with-thermal-blocks", "--with-network-blocks"])

file_name = !isempty(NO_OPTION_ARGS) ?
            string(NO_OPTION_ARGS[1], endswith(NO_OPTION_ARGS[1], ".yml") ? "" : ".yml") :
            "energy_community_model_CO_sto.yml"

## Initialization

data = read_input(file_name)

gen_data, users_data, market_data = data["general"], data["users"], data["market"]

user_set = user_names(gen_data, users_data)

init_step = field(gen_data, "init_step")
final_step = field(gen_data, "final_step")
time_set = init_step:final_step
n_steps = length(time_set)

# Peak set is required as a global by Scen_eps_sampler / scenario_definition.
peak_set = unique(ec_profile("peak_categories")[time_set])

# number of scenarios to be extracted
scen_s_sample = field(gen_data, "scen_s_sample")
scen_eps_sample = field(gen_data, "scen_eps_sample")

is_det = false
if scen_s_sample == 1 && scen_eps_sample == 1
    is_det = true
else
    # standard deviation associated with load and renewable production in long period uncertainty
    sigma_load = field(gen_data, "sigma_load")
    sigma_ren = field(gen_data, "sigma_ren")
end

# converters, i.e., CONV, are modeled with the corresponding BatteryUnitBlock in SMS++
SMSPP_DEVICES = setdiff(DEVICES, "--with-thermal-blocks" in OPTION_ARGS ? [CONV] : [CONV, THER])  # devices codes in SMS++


# Preprocessing to create the data structure (sampled_scenarios) for stochastic applications.
# If the model is deterministic, sampled_scenarios is nothing.
sampled_scenarios = nothing
if !is_det

    scen_s_set = 1:scen_s_sample
    scen_eps_set = 1:scen_eps_sample

    # Extraction of the point used to sample the distributions associated to the long period uncertainty
    (point_s_load,
        point_s_ren,
        scen_probability) = pem_extraction(scen_s_sample, sigma_load, sigma_ren)

    # To define an empty stochastic model we have to declare previously the scenarios.
    # sampled_scenarios is a list of Scenario_Load_Renewable defined in scenario_definition.jl;
    # see definition for more information.
    # Notable quantities are:
    #   sampled_scenarios[i].scen_s : scenario s
    #   sampled_scenarios[i].scen_eps : scenario epsilon
    #   probability(sampled_scenarios[1]) : denotes the probability of the scenario
    #   sampled_scenarios[i].Load["user1"][1] : load of user1 at time 1
    #   sampled_scenarios[i].Ren["user1"]["PV"][1] : PV production of user1 at time 1
    sampled_scenarios = scenarios_generator(data, point_s_load, point_s_ren, scen_probability, scen_s_sample, scen_eps_sample)
end

## Data aggregation and netCDF files generation
csvEC2nc4(is_det, sampled_scenarios)
