using YAML
# the official repo, i.e., https://github.com/JuliaGeo/NetCDF.jl, 
# does not support (jet) the concept of group :(
using NCDatasets

include("utils.jl")

function csvEC2nc4()

    # The mode "c" stands for creating a new file (clobber)
    ds = NCDataset("../../../netCDF_files/EC_Data/EC_Test.nc4", "c")

    block = defGroup(ds, "Block_0") # UCBlock
    n_users = length(user_set)
    defDim(block, "NumberNodes", n_users)
    n_timesteps = length(time_set)
    defDim(block, "TimeHorizon", n_timesteps)
    devices = [d for u in user_set
               for d in asset_names(users_data[u], SMSPP_DEVICES)]
    n_devices = length(devices)
    defDim(block, "NumberElectricalGenerators", n_devices)
    # w `ECNetworkBlock` for each time peak period/category; each of then span t time step/horizon + 
    # g `(Battery/Intermittent)UnitBlock` for each electrical generator/device
    peak_categories = profile(market_data, "peak_categories")
    peak_set = unique(peak_categories)
    n_peaks = length(peak_set)
    defDim(block, "NumberUnits", n_peaks + n_devices)

    # `ActivePowerDemand` is a 2D variable that represent the electricity
    # demand for each node/user wrt each time step/horizon
    power_demand = defVar(block, "ActivePowerDemand", Float64, ("NumberNodes", "TimeHorizon"))
    power_demand[:, :] = [profile_component(users_data[u], "load", "load")[t]
                          for u in user_set, t in time_set]

    # Let's create `EnergyCommunityNetworkBlock`(s)

    project_lifetime = field(gen_data, "project_lifetime")
    year_set = 1:project_lifetime

    sell_price_data = [profile(market_data, "energy_weight")[t] *
                       profile(market_data, "time_res")[t] *
                       profile(market_data, "sell_price")[t]
                       for t in time_set] * sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    buy_price_data = [profile(market_data, "energy_weight")[t] *
                      profile(market_data, "time_res")[t] *
                      profile(market_data, "buy_price")[t]
                      for t in time_set] * sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)

    _consumption_price_data = [profile(market_data, "energy_weight")[t] *
                               profile(market_data, "time_res")[t] *
                               (
                                   profile(market_data, "consumption_price")[t] *
                                   sum(Float64[
                                       profile_component(users_data[u], l, "load")[t]
                                       for l in asset_names(users_data[u], LOAD)])
                               )
                               for u in user_set, t in time_set]
    consumption_price_data = [sum(_consumption_price_data[u, t]
                                  for (u, _) in enumerate(user_set), t in time_set) *
                              sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)]

    last_t = 1
    for (i_w, w) in enumerate(peak_set)

        ecnb = defGroup(block, "NetworkBlock_$(i_w-1)")

        # `NumberIntervals`, i.e., the number of sub time horizon spanned by each peak period, i.e., an `ECNetworkBlock`
        # n_intervals = findlast(x -> x == w, peak_categories)
        n_intervals = count(x -> x == w, peak_categories)
        defDim(ecnb, "NumberIntervals", n_intervals)

        # Vector variables

        # `BuyPrice`, i.e., the tariff that user pay to buy electricity at each time horizon
        buy_price = defVar(ecnb, "BuyPrice", Float64, ("NumberIntervals",))
        buy_price = buy_price_data[last_t:n_intervals]

        # `SellPrice`, i.e., the tariff that user gain to sell electricity at each time horizon
        sell_price = defVar(ecnb, "SellPrice", Float64, ("NumberIntervals",))
        sell_price = sell_price_data[last_t:n_intervals]

        last_t += n_intervals

        # Scalar variables

        # `MaxTariff`, i.e., the peak tariff cost
        peak_tariff = defVar(ecnb, "MaxTariff", Float64, ())
        peak_tariff = profile(market_data, "peak_weight")[w] * profile(market_data, "peak_tariff")[w]

        # `ConstantTerm`
        constant_term = defVar(ecnb, "ConstantTerm", Float64, ())
        constant_term = sum(consumption_price_data)
    end

    # Let's create `(Battery/Intermittent)UnitBlock`(s)

    # `GeneratorNode` is a 1D variable that represent the node/user owner
    # of each electrical generator/device
    generator_node = defVar(block, "GeneratorNode", UInt32, ("NumberElectricalGenerators",))

    last_g = 1
    for (i_u, u) in enumerate(user_set)

        for (_, g) in enumerate(
            intersect(device_names(users_data[u]),
                devices))

            ub = defGroup(block, "UnitBlock_$(last_g-1)")

            generator_node[last_g] = i_u # assign the ownership of the current electrical generator to the respective user

            last_g += 1



        end
    end

    close(ds)
end

function csvNC2nc4()

end

## Parameters

file_name = "energy_community_model.yml"
model = "EC"

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

if model == "EC"
    csvEC2nc4()
elseif model == "NC"
    csvNC2nc4()
else
    throw(TaskFailedException("model unknown"))
end