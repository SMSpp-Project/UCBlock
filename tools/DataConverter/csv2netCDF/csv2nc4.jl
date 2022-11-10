using YAML
# the official repo, i.e., https://github.com/JuliaGeo/NetCDF.jl, 
# does not support (jet) the concept of group :(
using NCDatasets
using DataStructures

include("utils.jl")

function csvEC2nc4()

    # The mode "c" stands for creating a new file (clobber)
    ds = NCDataset("-with-network-blocks" in ARGS ? "../../../netCDF_files/EC_Test_NB.nc4" :
                   "../../../netCDF_files/EC_Test.nc4", "c", attrib=OrderedDict("SMS++_file_type" => 1))

    block = defGroup(ds, "Block_0", attrib=OrderedDict("id" => "0", "type" => "UCBlock"))

    n_users = length(user_set)
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

    # --------------------------------------------------------------------------------------- #

    # Create w `(EC)NetworkBlock`(s) for each peak period/category, each of them span t time step/horizon

    # Store the number of `(EC)NetworkBlock`(s), i.e., the number of peak period/category
    peak_categories = profile(market_data, "peak_categories")[time_set]
    peak_set = unique(peak_categories)
    n_peaks = length(peak_set)
    defDim(block, "NumberNetworks", n_peaks)

    # Create (sell/buy/consumption) price data arrays
    project_lifetime = field(gen_data, "project_lifetime")
    year_set = 1:project_lifetime

    sell_price_data = [profile(market_data, "energy_weight")[t] *
                       profile(market_data, "time_res")[t] *
                       profile(market_data, "sell_price")[t]
                       for t in time_set] *
                      sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    buy_price_data = [profile(market_data, "energy_weight")[t] *
                      profile(market_data, "time_res")[t] *
                      profile(market_data, "buy_price")[t]
                      for t in time_set] *
                     sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    constant_term = [sum(profile(market_data, "energy_weight")[t] *
                         profile(market_data, "time_res")[t] *
                         (
                             profile(market_data, "consumption_price")[t] *
                             sum(Float64[
                                 profile_component(users_data[u], l, "load")[t]
                                 for l in asset_names(users_data[u], LOAD)])
                         )
                         for u in user_set, t in time_set) *
                     sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)]

    peak_tariff = [profile(market_data, "peak_weight")[w] *
                   profile(market_data, "peak_tariff")[w] for w in peak_set] *
                  sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    reward_price_data = [profile(market_data, "energy_weight")[t] *
                         profile(market_data, "time_res")[t] *
                         profile(market_data, "reward_price")[t]
                         for t in time_set] *
                        sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    if ("-with-network-blocks" in ARGS) &&
       allequal([count(x -> x == w, peak_categories)
                 for w in peak_set]) &&
       allequal(sell_price_data) && allequal(buy_price_data) && allequal(peak_tariff) && allequal(reward_price_data)

        # `NumberIntervals`, i.e., the number of sub time horizon spanned by each peak period, i.e., an `ECNetworkBlock`
        n_intervals = count(x -> x == peak_set[1], peak_categories)
        defDim(block, "NumberIntervals", n_intervals)

        # Store the first index (-1 since in C++ the array's indexing starts from
        # zero) of each peak period/category, i.e., of each `(EC)NetworkBlock`
        # peak_start_idx = defVar(block, "StartNetworkIntervals", UInt32, ("NumberNetworks",))
        # peak_start_idx[:] = [findfirst(x -> x == w, peak_categories) - 1
        #                      for w in peak_set]

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

        # `BuyPrice`, i.e., the tariff that user pay to buy electricity at each time horizon
        buy_price = defVar(block, "BuyPrice", Float64, ())
        buy_price[:] = buy_price_data[1]

        # `SellPrice`, i.e., the tariff that user gain to sell electricity at each time horizon
        sell_price = defVar(block, "SellPrice", Float64, ())
        sell_price[:] = sell_price_data[1]

        # `RewardPrice`, i.e., the reward awarded to the community
        reward_price = defVar(block, "RewardPrice", Float64, ())
        reward_price[:] = reward_price_data[1]

        # `RenewableProduction` to bound the node injection
        max_injection = defVar(block, "MaxInjection", Float64, ("NumberNodes", "TimeHorizon")) # ("TimeHorizon", "NumberNodes"))
        # `reduce(+, itr; init)`, i.e., sum() over (possible) empty collection
        max_injection[:, :] = [reduce(+, [field_component(users_data[u], r, "max_capacity") *
                                          (profile_component(users_data[u], r, "ren_pu")[t] != 0 ?
                                           profile_component(users_data[u], r, "ren_pu")[t] : 1.0)
                                          for r in asset_names(users_data[u], REN)], init=0.0) +
                               reduce(+, [field_component(users_data[u], b, "max_capacity")
                                          for b in asset_names(users_data[u], BATT)], init=0.0)
                               for u in user_set, t in time_set]

        # `MaxTariff`, i.e., the peak tariff cost
        max_tariff = defVar(block, "MaxTariff", Float64, ())
        max_tariff[:] = peak_tariff[1]

        # `ConstantTerm`
        const_term = defVar(block, "ConstTerm", Float64, ())
        const_term[:] = sum(constant_term)

    else

        # create one `(EC)NetworkBlock` for each peak period/category
        last_t = 1
        for (i_w, w) in enumerate(peak_set)

            ecnb = defGroup(block, "NetworkBlock_$(i_w-1)", attrib=OrderedDict("type" => "ECNetworkBlock"))

            # Store the number of nodes in each NetworkBlock
            n_users = length(user_set)
            defDim(ecnb, "NumberNodes", n_users)

            # `NumberIntervals`, i.e., the number of sub time horizon spanned by each peak period, i.e., an `ECNetworkBlock`
            n_intervals = count(x -> x == w, peak_categories)
            defDim(ecnb, "NumberIntervals", n_intervals)

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
            buy_price = defVar(ecnb, "BuyPrice", Float64, ("NumberIntervals",))
            buy_price[:] = buy_price_data[last_t:last_i]

            # `SellPrice`, i.e., the tariff that user gain to sell electricity at each time horizon
            sell_price = defVar(ecnb, "SellPrice", Float64, ("NumberIntervals",))
            sell_price[:] = sell_price_data[last_t:last_i]

            # `RewardPrice`, i.e., the reward awarded to the community
            reward_price = defVar(ecnb, "RewardPrice", Float64, ("NumberIntervals",))
            reward_price[:] = reward_price_data[last_t:last_i]

            # `RenewableProduction` to bound the node injection
            max_injection = defVar(ecnb, "MaxInjection", Float64, ("NumberNodes", "NumberIntervals")) # ("NumberIntervals", "NumberNodes"))
            # `reduce(+, itr; init)`, i.e., sum() over (possible) empty collection
            max_injection[:, :] = [reduce(+, [field_component(users_data[u], r, "max_capacity") *
                                              (profile_component(users_data[u], r, "ren_pu")[t] != 0 ?
                                               profile_component(users_data[u], r, "ren_pu")[t] : 1.0)
                                              for r in asset_names(users_data[u], REN)], init=0.0) +
                                   reduce(+, [field_component(users_data[u], b, "max_capacity")
                                              for b in asset_names(users_data[u], BATT)], init=0.0)
                                   for u in user_set, t in last_t:last_i]

            last_t += n_intervals

            # `MaxTariff`, i.e., the peak tariff cost
            max_tariff = defVar(ecnb, "MaxTariff", Float64, ())
            max_tariff[:] = peak_tariff[i_w]

            # `ConstantTerm`
            const_term = defVar(ecnb, "ConstTerm", Float64, ())
            const_term[:] = sum(constant_term)
        end
    end

    # --------------------------------------------------------------------------------------- #

    # Let's create g `(Battery/Intermittent)UnitBlock`(s) for each electrical generator/device

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

                # store the maximum power, i.e., the maximum capacity, of the pv/wind device
                max_power = defVar(ub, "MaxPower", Float64, ("TimeHorizon",))
                max_power[:] = [field_component(users_data[u], g, "max_capacity") *
                                (profile_component(users_data[u], g, "ren_pu")[t] != 0 ?
                                 profile_component(users_data[u], g, "ren_pu")[t] : 1.0)
                                for t in time_set]

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

                # store the maximum power, i.e., the maximum capacity, of the converter related to the battery
                max_power = defVar(ub, "MaxPower", Float64, ())
                max_power[:] = (field_component(users_data[u], g, "max_C_dch") *
                                field_component(users_data[u], g, "max_capacity"))

                # store the minimum and maximum storage of the battery
                min_storage = defVar(ub, "MinStorage", Float64, ())
                min_storage[:] = (field_component(users_data[u], g, "min_SOC") *
                                  field_component(users_data[u], g, "max_capacity"))

                max_storage = defVar(ub, "MaxStorage", Float64, ())
                max_storage[:] = (field_component(users_data[u], g, "max_SOC") *
                                  field_component(users_data[u], g, "max_capacity"))

                # ASSUMPTION: if there is a battery there is ALWAYS also a converter!

                # store the intake roundtrip efficency of the battery
                intake_coeff = defVar(ub, "ExtractingBatteryRho", Float64, ("TimeHorizon",))
                intake_coeff[:] = [1 / # profile(market_data, "time_res")[t] /
                                   (sqrt(field_component(users_data[u], g, "eta")) *
                                    # corresponding converter, i.e., "conv"
                                    field_component(users_data[u], field_component(users_data[u], g, "corr_asset"), "eta"))
                                   for t in time_set]

                # store the outtake roundtrip efficency of the battery
                outtake_coeff = defVar(ub, "StoringBatteryRho", Float64, ("TimeHorizon",))
                outtake_coeff[:] = [# profile(market_data, "time_res")[t] *
                    (sqrt(field_component(users_data[u], g, "eta")) *
                     # corresponding converter, i.e., "conv"
                     field_component(users_data[u], field_component(users_data[u], g, "corr_asset"), "eta"))
                    for t in time_set]

                g_conv = field_component(users_data[u], g, "corr_asset") # corresponding converter, i.e., "conv"

                # Net Present Value of the component (both for battery and converter)
                battery_investment_cost = defVar(ub, "InvestmentCost", Float64, ())
                battery_investment_cost[:] = (sum(y == 0 ? field_component(users_data[u], g, "CAPEX_lin") : # investment cost of the component
                                                  (field_component(users_data[u], g, "OEM_lin") + # operation and maintenance cost of the component
                                                   ((mod(y, field_component(users_data[u], g, "lifetime_y")) == 0 && y != project_lifetime) ?
                                                    field_component(users_data[u], g, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                                   ((mod(y, field_component(users_data[u], g, "lifetime_y")) != 0 && y == project_lifetime) ?
                                                    field_component(users_data[u], g, "CAPEX_lin") *
                                                    (1.0 - mod(y, field_component(users_data[u], g, "lifetime_y")) /
                                                           field_component(users_data[u], g, "lifetime_y")) : 0.0)) * # residual value of the component
                                                  (1 / ((1 + field(gen_data, "d_rate"))^y))
                                                  for y in append!([0], year_set)) +
                                              sum(y == 0 ? field_component(users_data[u], g_conv, "CAPEX_lin") : # investment cost of the component
                                                  (field_component(users_data[u], g_conv, "OEM_lin") + # operation and maintenance cost of the component
                                                   ((mod(y, field_component(users_data[u], g_conv, "lifetime_y")) == 0 && y != project_lifetime) ?
                                                    field_component(users_data[u], g_conv, "CAPEX_lin") : 0.0) - # replacement cost of the component
                                                   ((mod(y, field_component(users_data[u], g_conv, "lifetime_y")) != 0 && y == project_lifetime) ?
                                                    field_component(users_data[u], g_conv, "CAPEX_lin") *
                                                    (1.0 - mod(y, field_component(users_data[u], g_conv, "lifetime_y")) /
                                                           field_component(users_data[u], g_conv, "lifetime_y")) : 0.0)) * # residual value of the component
                                                  (1 / ((1 + field(gen_data, "d_rate"))^y))
                                                  for y in append!([0], year_set))) * field_component(users_data[u], g_conv, "max_capacity")

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