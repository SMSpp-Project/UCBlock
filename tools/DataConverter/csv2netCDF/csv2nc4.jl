using YAML
# the official repo, i.e., https://github.com/JuliaGeo/NetCDF.jl, 
# does not support (yet) the concept of group :(
using NCDatasets
using DataStructures

using Parameters
using DataFrames
using XLSX
using JLD2
using YAML
using CSV

using Distributions
using PointEstimateMethod

using StochasticPrograms

using Random

# include additional useful functions, i.e., main type definitions and read data
include("utils.jl")

# setting the seed
Random.seed!(123)

function csvEC2nc4(deterministic::Bool=false)

    middle = "_"
    if occursin("_CO", file_name)
        middle = "CO_"
    elseif occursin("_NA", file_name)
        middle = "NA_"
    elseif occursin("_NC", file_name)
        middle = "NC_"
    end

    last = ""
    if "--with-thermal-blocks" in OPTION_ARGS && !occursin("_NA", file_name)
        last = string(last, "_TUB")
    end
    if "--with-network-blocks" in OPTION_ARGS
        last = string(last, "_NB")
    end

    # The mode "c" stands for creating a new file (clobber)
    ds = NCDataset(string("../../../netCDF_files/EC_Data/EC", middle, "Test", last, ".nc4"), "c", attrib=OrderedDict("SMS++_file_type" => 1))
    block = defGroup(ds, "Block_0", attrib=OrderedDict("id" => "0", "type" => "UCBlock"))

    # Store the number of nodes
    n_users = length(user_set)
    defDim(block, "NumberNodes", n_users)

    # Store the number of time steps/horizons
    defDim(block, "TimeHorizon", n_steps)

    # Store the number of `ECNetworkBlock`(s), i.e., the number of peak periods/categories
    peak_categories = profile(market_data, "peak_categories")[time_set]
    peak_set = unique(peak_categories)
    n_peaks = length(peak_set)
    defDim(block, "NumberNetworks", n_peaks)

    # Create buy, sell, reward, and consumption, i.e., the constant term, price data arrays
    project_lifetime = field(gen_data, "project_lifetime")
    year_set = 1:project_lifetime

    # `BuyPrice`, i.e., the tariff that user pay to buy electricity at each time horizon
    buy_price_data = [profile(market_data, "buy_price")[t] *
                      profile(market_data, "energy_weight")[t] *
                      profile(market_data, "time_res")[t]
                      for t in time_set] *
                     sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    # `SellPrice`, i.e., the tariff that user gain to sell electricity at each time horizon
    sell_price_data = [profile(market_data, "sell_price")[t] *
                       profile(market_data, "energy_weight")[t] *
                       profile(market_data, "time_res")[t]
                       for t in time_set] *
                      sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    # `RewardPrice`, i.e., the reward awarded to the community
    reward_price_data = [profile(market_data, "reward_price")[t] *
                         profile(market_data, "energy_weight")[t] *
                         profile(market_data, "time_res")[t]
                         for t in time_set] *
                        sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    # `PenaltyPrice`, i.e., the penalty price for energy squilibrium
    #= penalty_price_data = [profile(market_data, "penalty_price")[t] *
                            profile(market_data, "energy_weight")[t] *
                            profile(market_data, "time_res")[t]
                            for t in time_set] *
                            sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set) =#

    # `PeakTariff`, i.e., the peak tariff cost
    peak_tariff_data = [profile(market_data, "peak_tariff")[w] *
                        profile(market_data, "peak_weight")[w]
                        for w in peak_set] *
                       sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    # `ConstantTerm`, i.e., the consumption price
    const_term_data = [sum(profile(market_data, "consumption_price")[t] *
                           profile_component(users_data[u], l, "load")[t]
                           for u in user_set for l in asset_names(users_data[u], LOAD)) *
                       profile(market_data, "energy_weight")[t] *
                       profile(market_data, "time_res")[t]
                       for t in time_set] *
                      sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

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
        # path_group_idx_data = Int[]
        path_group_idx_data = String[]
        path_element_idx_data = Int[]
        path_range_idx_data = Int[]
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
                        # append!(path_group_idx_data, [last_g, 0]) # i.e., last_g wrt B, 0 wrt V x_intermittent
                        append!(path_group_idx_data, [string(last_g), "x_intermittent"]) # i.e., last_g wrt B, V x_intermittent
                        append!(path_element_idx_data, [typemax(UInt32), 0]) # i.e., _ wrt B, 0 wrt V x_intermittent
                        append!(path_range_idx_data, [typemax(UInt32), 1]) # i.e., _ wrt B, 1 or _ wrt V x_intermittent
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
                                        profile(market_data, "time_res")[t] # energy (kWh), i.e., power * time, to power (kW), i.e., energy / time
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
                                        profile(market_data, "time_res")[t] # energy (kWh), i.e., power * time, to power (kW), i.e., energy / time
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
                        # append!(path_group_idx_data, [last_g, 0, last_g, 1]) # i.e., last_g wrt B, 0 wrt V x_battery, 1 wrt V x_converter
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
                        thermal_max_power[:] = (field_component(users_data[u], g, "max_technical") *
                                                field_component(users_data[u], g, "nom_capacity"))

                        # store the start-up limit
                        thermal_start_up_limit = defVar(ub, "StartUpLimit", Float64, ())
                        thermal_start_up_limit[:] = thermal_max_power

                        # store the shut-down limit
                        thermal_shut_up_limit = defVar(ub, "ShutDownLimit", Float64, ())
                        thermal_shut_up_limit[:] = thermal_max_power

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
                                                profile(market_data, "energy_weight")[t] *
                                                profile(market_data, "time_res")[t]
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
                                               profile(market_data, "energy_weight")[t] *
                                               profile(market_data, "time_res")[t]
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
                            # append!(path_group_idx_data, [last_g, 0]) # i.e., last_g wrt B, 0 wrt V x_thermal
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

        # The mode "c" stands for creating a new file (clobber)
        tssb_ds = NCDataset(string("../../../netCDF_files/EC_Data/TSSB_EC", middle, "Test", last, ".nc4"), "c", attrib=OrderedDict("SMS++_file_type" => 1))
        tssb = defGroup(tssb_ds, "Block_0", attrib=OrderedDict("id" => "0", "type" => "TwoStageStochasticBlock"))

        defDim(tssb, "NumberScenarios", scen_s_sample)

        # ScenarioGenerator
        # dss = defGroup(tssb, "ScenarioGenerator", attrib=OrderedDict("type" => "DiscreteScenarioSet"))

        # defDim(dss, "NumberScenarios", scen_s_sample)
        # defDim(dss, "ScenarioSize", )

        ## A T T E N T I O N: The data is stored in the NetCDF file in the same order as they are
        ## stored in memory. As Julia uses the column-major ordering for arrays, the order of dimensions
        ## will appear reversed when the data is loaded in languages or programs using row-major
        ## ordering such as C/C++, Python/NumPy or the tools ncdump/ncgen.
        ## To store the scenario set in the correct shape, i.e., NumberScenarios x ScenarioSize, we need to store
        ## it transposed, i.e., ScenarioSize x NumberScenarios.
        # scenario_set = defVar(block, "ScenarioSet", Float64, ("ScenarioSize", "NumberScenarios")) # ("NumberScenarios", "ScenarioSize"))
        # scenario_set[:, :] = [ ]

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

        # path_group_idx = defVar(ap, "PathGroupIndices", UInt32, ("TotalLength",))
        path_group_idx = defVar(ap, "PathGroupIndices", String, ("TotalLength",))
        path_group_idx[:] = path_group_idx_data[:]

        path_element_idx = defVar(ap, "PathElementIndices", UInt32, ("TotalLength",))
        path_element_idx[:] = path_element_idx_data[:]

        path_range_idx = defVar(ap, "PathRangeIndices", UInt32, ("TotalLength",))
        path_range_idx[:] = path_range_idx_data[:]

        # StochasticBlock
        sb = defGroup(tssb, "StochasticBlock", attrib=OrderedDict("type" => "StochasticBlock"))

        # UCBlock nc4 file
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
            "energy_community_model_CO.yml"

## Initialization

data = read_input(file_name)

gen_data, users_data, market_data = data["general"], data["users"], data["market"]

user_set = user_names(gen_data, users_data)

init_step = field(gen_data, "init_step")
final_step = field(gen_data, "final_step")
time_set = init_step:final_step
n_steps = length(time_set)

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

## Data aggregation and netCDF files generation
csvEC2nc4(is_det)