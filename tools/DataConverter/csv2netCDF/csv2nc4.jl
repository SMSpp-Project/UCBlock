using YAML
# the official repo, i.e., https://github.com/JuliaGeo/NetCDF.jl, 
# does not support (jet) the concept of group :(
using NCDatasets
using DataStructures

include("utils.jl")

function csvEC2nc4()

    n_users = length(user_set)

    # The mode "c" stands for creating a new file (clobber)
    ds = NCDataset(string("../../../netCDF_files/EC_Test.nc4"), "c", attrib=OrderedDict("SMS++_file_type" => 1))

    block = defGroup(ds, "Block_0", attrib=OrderedDict("id" => "0", "type" => "UCBlock"))

    defDim(block, "NumberNodes", n_users)

    n_timesteps = length(time_set)
    defDim(block, "TimeHorizon", n_timesteps)

    # Store the specific classname of the NetworkBlock, i.e., `ECNetworkBlock` and `ECNetworkData`, to
    # inform UCBlock about the specific type of network (since it deals with both transmission and 
    # community networks)
    network_block_classname = defVar(block, "NetworkBlockClassname", String, ())
    network_block_classname[1] = "ECNetworkBlock"
    network_data_classname = defVar(block, "NetworkDataClassname", String, ())
    network_data_classname[1] = "ECNetworkData"

    # Store the number of `ECNetworkBlock`(s), i.e., the number of peak period/category
    peak_categories = profile(market_data, "peak_categories")[time_set]
    peak_set = unique(peak_categories)
    n_peaks = length(peak_set)
    defDim(block, "NumberNetworks", n_peaks)

    if ("-store-demand-in-father" in ARGS)
        # `NumberIntervals`, i.e., the number of sub time horizon spanned by each peak period, i.e., an `ECNetworkBlock`
        n_intervals = count(x -> x == peak_set[1], peak_categories)
        defDim(block, "NumberIntervals", n_intervals)

        # Store the first index (-1 since in C++ the array's indexing starts from
        # zero) of each peak period/category, i.e., of each `(EC)NetworkBlock`
        peak_start_idx = defVar(block, "StartNetworkIntervals", UInt32, ("NumberNetworks",))
        peak_start_idx[:] = [findfirst(x -> x == w, peak_categories) - 1
                             for w in peak_set]

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
    end

    # Create buy, sell, reward, and consumption price data arrays
    project_lifetime = field(gen_data, "project_lifetime")
    year_set = 1:project_lifetime

    # `BuyPrice`, i.e., the tariff that user pay to buy electricity at each time horizon
    buy_price_data = [profile(market_data, "energy_weight")[t] *
                      profile(market_data, "time_res")[t] *
                      profile(market_data, "buy_price")[t]
                      for t in time_set] *
                     sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    # `SellPrice`, i.e., the tariff that user gain to sell electricity at each time horizon
    sell_price_data = [profile(market_data, "energy_weight")[t] *
                       profile(market_data, "time_res")[t] *
                       profile(market_data, "sell_price")[t]
                       for t in time_set] *
                      sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    # `RewardPrice`, i.e., the reward awarded to the community
    reward_price_data = [profile(market_data, "energy_weight")[t] *
                         profile(market_data, "time_res")[t] *
                         profile(market_data, "reward_price")[t]
                         for t in time_set] *
                        sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    # `PeakTariff`, i.e., the peak tariff cost
    peak_tariff_data = [(profile(market_data, "peak_weight")[w] *
                         profile(market_data, "peak_tariff")[w])
                        for w in peak_set] *
                       sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    constant_term = [sum(profile(market_data, "energy_weight")[t] *
                         profile(market_data, "time_res")[t] *
                         (profile(market_data, "consumption_price")[t] *
                          sum(Float64[
                             profile_component(users_data[u], l, "load")[t]
                             for l in asset_names(users_data[u], LOAD)]))
                         for u in user_set)
                     for t in time_set] *
                    sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    # Create w `ECNetworkBlock`(s) for each peak period/category, each of them span w_t time step/horizon
    last_t = 1
    for (i_w, w) in enumerate(peak_set)

        ecnb = defGroup(block, "NetworkBlock_$(i_w-1)", attrib=OrderedDict("type" => "ECNetworkBlock"))

        # `NumberIntervals`, i.e., the number of sub time horizon spanned by each peak period, i.e., an `ECNetworkBlock`
        n_intervals = count(x -> x == w, peak_categories)
        defDim(ecnb, "NumberIntervals", n_intervals)

        last_i = findlast(x -> x == w, peak_categories)

        if !("-store-demand-in-father" in ARGS)
            # Store the number of nodes in each NetworkBlock
            n_users = length(user_set)
            defDim(ecnb, "NumberNodes", n_users)

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
        end

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

        # `RewardPrice`, i.e., the reward awarded to the community
        if (allequal(reward_price_data[last_t:last_i]))
            reward_price = defVar(ecnb, "RewardPrice", Float64, ())
            reward_price[:] = reward_price_data[last_t]
        else
            reward_price = defVar(ecnb, "RewardPrice", Float64, ("NumberIntervals",))
            reward_price[:] = reward_price_data[last_t:last_i]
        end

        # `PeakTariff`, i.e., the peak tariff cost
        peak_tariff = defVar(ecnb, "PeakTariff", Float64, ())
        peak_tariff[:] = peak_tariff_data[i_w]

        # `ConstTerm`, i.e., the consumption price
        const_term = defVar(ecnb, "ConstTerm", Float64, ())
        const_term[:] = sum(constant_term[last_t:last_i])

        # # `MaxNodeInjection` to bound the node injection
        # max_injection = defVar(ecnb, "MaxNodeInjection", Float64, ("NumberNodes", "NumberIntervals")) # ("NumberIntervals", "NumberNodes"))
        # # `reduce(+, itr; init)`, i.e., sum() over (possible) empty collection
        # max_injection[:, :] = [reduce(+, [field_component(users_data[u], r, "max_capacity") *
        #                                   profile_component(users_data[u], r, "ren_pu")[t]
        #                                   for r in asset_names(users_data[u], REN)], init=0.0) +
        #                        reduce(+, [field_component(users_data[u], b, "max_capacity")
        #                                   for b in asset_names(users_data[u], BATT)], init=0.0)
        #                        for u in user_set, t in last_t:last_i]

        last_t += n_intervals
    end

    # --------------------------------------------------------------------------------------- #

    # Create g `(Battery/Intermittent)UnitBlock`(s) for each electrical generator/device

    devices = [d for u in user_set
               for d in asset_names(users_data[u], SMSPP_DEVICES)]
    n_devices = length(devices)
    # number of UnitBlock
    defDim(block, "NumberUnits", n_devices)
    # each UnitBlock has just one electrical generator
    defDim(block, "NumberElectricalGenerators", n_devices)

    # `GeneratorNode` is a 1D variable that represent the node/user owner
    # of each electrical generator/device
    generator_node = defVar(block, "GeneratorNode", UInt32, ("NumberElectricalGenerators",))

    last_g = 1
    for (i_u, u) in enumerate(user_set)

        for g in intersect(device_names(users_data[u]), devices)

            if g in ("PV", "wind")

                ub = defGroup(block, "UnitBlock_$(last_g - 1)", attrib=OrderedDict("type" => "IntermittentUnitBlock"))

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

                # Net Present Value of the component
                investment_cost = defVar(ub, "InvestmentCost", Float64, ())
                investment_cost[:] = sum(y == 0 ? field_component(users_data[u], g, "CAPEX_lin") : # investment cost of the component
                                         (field_component(users_data[u], g, "OEM_lin") + # operation and maintenance cost of the component
                                          ((mod(y, field_component(users_data[u], g, "lifetime_y")) == 0 && y != project_lifetime) ?
                                           field_component(users_data[u], g, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                          ((mod(y, field_component(users_data[u], g, "lifetime_y")) != 0 && y == project_lifetime) ?
                                           field_component(users_data[u], g, "CAPEX_lin") *
                                           (1.0 - mod(y, field_component(users_data[u], g, "lifetime_y")) /
                                                  field_component(users_data[u], g, "lifetime_y")) : 0.0)) * # residual value of the component
                                         (1 / ((1 + field(gen_data, "d_rate"))^y))
                                         for y in append!([0], year_set)) * field_component(users_data[u], g, "max_capacity")

            elseif g == "batt"

                ub = defGroup(block, "UnitBlock_$(last_g - 1)", attrib=OrderedDict("type" => "BatteryUnitBlock"))

                # ----------- Battery -----------

                # store the maximum installable capacity of the battery
                batt_max_capacity = defVar(ub, "BatteryMaxCapacity", Float64, ())
                batt_max_capacity[:] = field_component(users_data[u], g, "max_capacity")

                # store the maximum power of the battery
                batt_max_power = defVar(ub, "MaxPower", Float64, ())
                batt_max_power[:] = field_component(users_data[u], g, "max_capacity")

                # store the minimum storage of the battery
                min_storage_data = [(field_component(users_data[u], g, "min_SOC") *
                                     field_component(users_data[u], g, "max_capacity")) /
                                    profile(market_data, "time_res")[t] # energy (kWh), i.e., power * time, to power (kW), i.e., energy / time
                                    for t in time_set]
                if (allequal(min_storage_data))
                    min_storage = defVar(ub, "MinStorage", Float64, ())
                    min_storage[:] = min_storage_data[1]
                else
                    min_storage = defVar(ub, "MinStorage", Float64, ("TimeHorizon",))
                    min_storage[:] = min_storage_data[:]
                end

                # store the maximum storage of the battery
                max_storage_data = [(field_component(users_data[u], g, "max_SOC") *
                                     field_component(users_data[u], g, "max_capacity")) /
                                    profile(market_data, "time_res")[t] # energy (kWh), i.e., power * time, to power (kW), i.e., energy / time
                                    for t in time_set]
                if (allequal(max_storage_data))
                    max_storage = defVar(ub, "MaxStorage", Float64, ())
                    max_storage[:] = max_storage_data[1]
                else
                    max_storage = defVar(ub, "MaxStorage", Float64, ("TimeHorizon",))
                    max_storage[:] = max_storage_data[:]
                end

                # Net Present Value of the battery
                batt_investment_cost = defVar(ub, "BatteryInvestmentCost", Float64, ())
                batt_investment_cost[:] = (sum(y == 0 ? field_component(users_data[u], g, "CAPEX_lin") : # investment cost of the component
                                               (field_component(users_data[u], g, "OEM_lin") + # operation and maintenance cost of the component
                                                ((mod(y, field_component(users_data[u], g, "lifetime_y")) == 0 && y != project_lifetime) ?
                                                 field_component(users_data[u], g, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                                ((mod(y, field_component(users_data[u], g, "lifetime_y")) != 0 && y == project_lifetime) ?
                                                 field_component(users_data[u], g, "CAPEX_lin") *
                                                 (1.0 - mod(y, field_component(users_data[u], g, "lifetime_y")) /
                                                        field_component(users_data[u], g, "lifetime_y")) : 0.0)) * # residual value of the component
                                               (1 / ((1 + field(gen_data, "d_rate"))^y))
                                               for y in append!([0], year_set)) * field_component(users_data[u], g, "max_capacity"))

                # ---------- Converter ----------

                g_conv = field_component(users_data[u], g, "corr_asset") # corresponding converter, i.e., "conv"

                # store the maximum installable capacity of the converter
                conv_max_capacity = defVar(ub, "ConverterMaxCapacity", Float64, ())
                conv_max_capacity[:] = field_component(users_data[u], g_conv, "max_capacity")

                # store the maximum power of the converter
                conv_max_power = defVar(ub, "ConverterMaxPower", Float64, ())
                conv_max_power[:] = field_component(users_data[u], g_conv, "max_capacity")

                # store the intake roundtrip efficency of the battery
                intake_coeff_data = [1 / (sqrt(field_component(users_data[u], g, "eta")) *
                                          # corresponding converter, i.e., "conv"
                                          field_component(users_data[u], g_conv, "eta"))
                                     for t in time_set]
                if (allequal(intake_coeff_data))
                    intake_coeff = defVar(ub, "ExtractingBatteryRho", Float64, ())
                    intake_coeff[:] = intake_coeff_data[1]
                else
                    intake_coeff = defVar(ub, "ExtractingBatteryRho", Float64, ("TimeHorizon",))
                    intake_coeff[:] = intake_coeff_data[:]
                end

                # store the outtake roundtrip efficency of the battery
                outtake_coeff_data = [sqrt(field_component(users_data[u], g, "eta")) *
                                      # corresponding converter, i.e., "conv"
                                      field_component(users_data[u], g_conv, "eta")
                                      for t in time_set]
                if (allequal(outtake_coeff_data))
                    outtake_coeff = defVar(ub, "StoringBatteryRho", Float64, ())
                    outtake_coeff[:] = outtake_coeff_data[1]
                else
                    outtake_coeff = defVar(ub, "StoringBatteryRho", Float64, ("TimeHorizon",))
                    outtake_coeff[:] = outtake_coeff_data[:]
                end

                # Net Present Value of the converter
                conv_investment_cost = defVar(ub, "ConverterInvestmentCost", Float64, ())
                conv_investment_cost[:] = (sum(y == 0 ? field_component(users_data[u], g_conv, "CAPEX_lin") : # investment cost of the component
                                               (field_component(users_data[u], g_conv, "OEM_lin") + # operation and maintenance cost of the component
                                                ((mod(y, field_component(users_data[u], g_conv, "lifetime_y")) == 0 && y != project_lifetime) ?
                                                 field_component(users_data[u], g_conv, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                                ((mod(y, field_component(users_data[u], g_conv, "lifetime_y")) != 0 && y == project_lifetime) ?
                                                 field_component(users_data[u], g_conv, "CAPEX_lin") *
                                                 (1.0 - mod(y, field_component(users_data[u], g_conv, "lifetime_y")) /
                                                        field_component(users_data[u], g_conv, "lifetime_y")) : 0.0)) * # residual value of the component
                                               (1 / ((1 + field(gen_data, "d_rate"))^y))
                                               for y in append!([0], year_set)) * field_component(users_data[u], g_conv, "max_capacity"))
            end

            generator_node[last_g] = i_u - 1 # assign the ownership of the current electrical generator to the respective user
            last_g += 1
        end
    end

    close(ds)
end

## Parameters

file_name = "energy_community_model.yml"

## Initialization

data = read_input(file_name)

gen_data, users_data, market_data = data["general"], data["users"], data["market"]

user_set = user_names(gen_data, users_data)

init_step = field(gen_data, "init_step")
final_step = field(gen_data, "final_step")
time_set = init_step:final_step

# converters, i.e., CONV, are modeled with the corresponding BatteryUnitBlock in SMS++
SMSPP_DEVICES = setdiff(DEVICES, [CONV])  # devices codes in SMS++

## Data aggregation and netCDF files generation
csvEC2nc4()