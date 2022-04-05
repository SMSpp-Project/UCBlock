using Parameters, JuMP, DataFrames
import XLSX, JLD2, YAML, CSV


@enum ASSET_TYPE LOAD = 0 REN = 1 BATT = 2 CONV = 3
ANY = collect(instances(ASSET_TYPE))  # all assets code
DEVICES = setdiff(ANY, [LOAD])  # devices codes
GENS = [REN]  # generator codes

type_codes = Base.Dict(
    "renewable" => REN,
    "battery" => BATT,
    "converter" => CONV,
    "load" => LOAD,
)

@with_kw mutable struct USER_TYPE
    user_name::String = ""  # Name of the user
    d_rate::Float64 = 0.0  # Discount rate
    max_opt_out::Float64 = 0.0  # Maximum opt-out value
    eta_line::Float64 = 1  # Efficiency of the line between the user and the junction box

    time_steps::Array{Int,1} = Array{Int,1}(undef, 0) # time step id
    asset_names::Array{String,1} = Array{String,1}(undef, 0)  # Name representing the assets [Name of the asset]

    asset_type::Dict{String,ASSET_TYPE} = Dict{String,ASSET_TYPE}()  # Type of the asset [See enumerator ASSET_TYPE]
    CAPEX_cost::Dict{String,Float64} = Dict{String,Float64}() # Investment cost of the assets [€/base unit]
    OEM_cost::Dict{String,Float64} = Dict{String,Float64}()  # Operating and maintenance cost per unit asset [€/base unit]
    lifetime::Dict{String,Float64} = Dict{String,Float64}()  # Lifetime of each asset [y]
    eta::Dict{String,Float64} = Dict{String,Float64}()  # Efficiency of the asset (roundtrip for batteries)
    max_dispatch::Dict{String,Float64} = Dict{String,Float64}()  # Maximum dispatch of each asset [pu]
    min_dispatch::Dict{String,Float64} = Dict{String,Float64}()  # Minimum dispatch of each asset [pu]
    corr_asset::Dict{String,String} = Dict{String,String}()  # Corresponding element, when applicable: i.e. battery correspond to a converter and viceversa
    max_capacity::Dict{String,Float64} = Dict{String,Float64}()  # Estimated max capacity [base unit]


    load::Dict{Int,Float64} = Dict{Int,Float64}()  # 1D array representing the load of the user
    ren_pu::Dict{String,Dict{Int,Float64}} = Dict{String,Dict{Int,Float64}}()  # 2D array representing the available renewable production for every renewable asset
end

@with_kw mutable struct AGGREGATOR_TYPE
    sigma_min::Float64 = 0  # Minimum value of the user surplus kept by the aggregator
    sigma_max::Float64 = 0  # Maximum value of the user surplus kept by the aggregator
    sigma_steps::Int64 = 0  # Numebr of steps to discretize the sigma value
    project_lifetime::Int64 = 20  # Number of years of the project
    #res_val_coeff::Int64 = 0.8  # Coefficient to model the remaining value of an assets

    discount_names::Vector{String} = Vector{String}()  # Names of the discounts
    discount_asset::Dict{String,String} = Dict{String,String}()  # Assets corresponding to the current discount type
    discount_max::Dict{String,Float64} = Dict{String,Float64}()  # max discount
    discount_min::Dict{String,Float64} = Dict{String,Float64}()  # min discount
    discount_steps::Dict{String,Int64} = Dict{String,Int64}()  # number of steps to test each discount
end

@with_kw mutable struct MARKET_TYPE
    time_steps::Array{Int,1} = Array{Int,1}(undef, 0)  # time step id

    peak_period_names::Vector{String} = Vector{String}()  # name of each peak period
    peak_weight::Dict{String,Float64} = Dict{String,Float64}()  # weight of the peak power period
    peak_tariff::Dict{String,Float64} = Dict{String,Float64}()  # peak tariff for every peak period
    init_peak_period::Dict{String,Int} = Dict{String,Int}()  # Initial time step of each peak period
    final_peak_period::Dict{String,Int} = Dict{String,Int}()  # Final time step of each peak period

    sell_price::Dict{Int,Float64} = Dict{Int,Float64}()  # 1D array representing the variable seller price to the public market [€/MWh] (variable component)
    buy_price::Dict{Int,Float64} = Dict{Int,Float64}()  # 1D array representing the buyer price to the public market [€/MWh]
    consumption_price::Dict{Int,Float64} = Dict{Int,Float64}()  # 1D array representing the price on consumption to the public market [€/MWh] (fixed component)
    reward_price::Dict{Int,Float64} = Dict{Int,Float64}()  # 1D array representing the reward awarded to the community [€/MWh]
    time_weight::Dict{Int,Float64} = Dict{Int,Float64}()  # 1D array representing a weight for each time step, except for time resolution
    # it means how many time steps it is representing: example: a year represented by 12 days -> time_weight = 365/12
end

# # Return whether time_step belong to the peak period peak_period
# @inline inPeakPeriod(market_data::MARKET_TYPE, time_step::Int, peak_period::String) =
#     (time_step >= market_data.init_peak_period[peak_period]) && (time_step <= market_data.final_peak_period[peak_period])
# # Return the peak period that contains the requested time step
# @inline getPeakPeriod(market_data::MARKET_TYPE, time_step::Int) = [peak_period for peak_period in market_data.peak_period_names if inPeakPeriod(market_data, time_step, peak_period)]
# Return the time periods related to the current peak period
@inline getTimePeriod(market_data, peak_period::String) = [
    time_step for time_step in market_data.time_steps if
    inPeakPeriod(market_data, time_step, peak_period)
]

@with_kw struct DATA_TYPE
    n_users::Int = 0  # Number of users to read
    init_step::Int = 0  # Initial time step
    final_step::Int = 0  # Final time step
    n_steps::Int = 0 # Simulation intervals
    time_res::Float64 = 0.0 # Time resolution [h]
    aggregator_data::AGGREGATOR_TYPE = AGGREGATOR_TYPE()  # Data of the aggregator
    users_data::Array{USER_TYPE,1} = Array{USER_TYPE,1}(undef, 0)  # Vector of the data of users
    market_data::MARKET_TYPE = MARKET_TYPE()  # Market data
end

# Get the previous time step, with circular time step
@inline pre(time_step::Int, gen_data::Dict) =
    if (time_step > field(gen_data, "init_step"))
        time_step - 1
    else
        field(gen_data, "final_step")
    end
@inline pre(time_step::V, time_set::UnitRange{V}) where {V<:Int} =
    if (time_step > time_set[1])
        time_step - 1
    else
        time_set[end]
    end


"Function to safely get a field of a dictionary with default value"
@inline field_d(d::AbstractDict, field, default = nothing) =
    (field in keys(d) ? d[field] : default)
@inline field_i(d, field) = field_d(d, field, 0)
@inline field_f(d, field) = field_d(d, field, 0.0)
"Function get field that throws an error if the field is not found"
@inline function field(d, field, desc = nothing)
    if d isa AbstractDict && field in keys(d)
        return d[field]
    else
        msg =
            isnothing(desc) ?
            "Field $field not found in dictionary $(keys(d))" : desc
        throw(KeyError(msg))
    end
end

"Function to get the general parameters"
general(d) = field(d, "general")
"Function to get the users configuration"
users(d) = field(d, "users")
"Function to get the market configuration"
market(d) = field(d, "market")
"Function to get the profile dictionary"
profiles(d) = field_d(d, "profile")
"Function to get the components list of a dictionary"
components(d) = d
"Function to get the components value of a dictionary"
component(d, c_name) = field(components(d), c_name)
"Function to get the components value of a dictionary"
field_component(d, c_name, f_name) = field(component(d, c_name), f_name)
"Function to get a specific profile"
function profile(d, profile_name)
    profile_block = profiles(d)
    return field(profile_block, profile_name)
end
"Function to get a specific profile"
function profile_component(d, c_name, profile_name)
    profile_block = profiles(component(d, c_name))
    return field(profile_block, profile_name)
end

"Function to get the asset type of a component"
asset_type(d, comp_name) = type_codes[field(component(d, comp_name), "type")]

"Function to get the list of the assets for a user"
function asset_names(d, a_type::ASSET_TYPE)
    comps = components(d)
    return [at for at in keys(comps) if asset_type(comps, at) == a_type]
end

"Function to get the list of the assets for a user"
asset_names(d) = collect(keys(components(d)))

"Function to get the list of the assets for a user in a list of elements"
function asset_names(d, a_types::Vector{ASSET_TYPE})
    comps = components(d)
    return [at for at in keys(comps) if asset_type(comps, at) ∈ a_types]
end

"Function to get the list of the assets for a user in a list of elements except a list of given types"
function asset_names_ex(d, ex::Vector{ASSET_TYPE})
    comps = components(d)
    accepted_types = set_diff(ANY, ex)
    return asset_names(d, accepted_types)
end

"Function to get the list of devices for a user"
device_names(d) = asset_names(d, DEVICES)

"Get the list of users"
function user_names(gen_data, users_data)
    # get the list of users if set
    user_list = field_d(gen_data, "user_list")
    if isnothing(user_list)
        @info "List of users not specified: all users selected"
        user_list = collect(keys(users_data))
    elseif isempty(user_list)
        throw(ErrorException("Input user list is empty"))
    elseif !(user_list isa AbstractVector)
        throw(ErrorException("Input user list is not a vector"))
    end
    return sort!(user_list)
end


"""
Function to parse a string value of a profile to load the corresponding dataframe
"""
function parse_dataprofile(
    gen_config,
    data,
    profile_name,
    profile_value::AbstractString,
)

    if profile_value in names(data)
        return data[!, profile_value]
    else
        throw(
            KeyError(
                "Profile name $profile_value not found in available dataframes",
            ),
        )
    end
end

"""
Function to parse a string value of a profile to load the corresponding dataframe
"""
function parse_dataprofile(
    gen_config,
    data,
    profile_name,
    profile_value::AbstractVector{T},
) where {T<:Integer}

    if length(profile_value) <
       gen_config["final_step"] - gen_config["init_step"] + 1
        throw(
            ErrorException(
                "Not enough profile data available for the current profile list",
            ),
        )
    end

    return profile_value
end

"""
Function to parse a personalized processing to generate the data
When profile_value is a dictionary, then the user is asking a custom processing of data by a function
"""
function parse_dataprofile(gen_config, data, profile_name, profile_value::Dict)

    func_name = field(profile_value, "function")
    inputs = field(profile_value, "inputs")

    # load input data for the function
    input_data = []
    for i_data in inputs
        push!(
            input_data,
            parse_dataprofile(
                gen_config,
                data,
                profile_name * " inputs",
                i_data,
            ),
        )
    end

    # prepare the execution of the function
    cmd_expr = Expr(
        :call,
        Symbol(func_name),
        gen_config,
        data,
        profile_name,
        input_data...,
    )

    # execute the function
    ret_value = eval(cmd_expr)

    return ret_value
end


"""
Function to parse a string value of a profile to load the corresponding dataframe
"""
function parse_dataprofile(
    gen_config,
    data,
    profile_name,
    profile_value::T,
) where {T<:Real}

    n_steps = gen_config["final_step"] - gen_config["init_step"] + 1

    return fill(convert(Float64, profile_value), n_steps)
end

"""
Function to throw error for unformatted data
"""
function parse_dataprofile(gen_config, data, profile_name, profile_value::Any)
    @error "Data descriptor for $profile_name not accepted"
end

"""
Function to read the input of the optimization model described as a yaml file
"""
function read_input(file_name::AbstractString)

    data = YAML.load_file(file_name)

    gen_data = general(data)

    # optional data by csv files
    opt_data = DataFrame()
    opt_files = field(gen_data, "optional_datasets")
    if !isnothing(opt_files)  # datasets are available
        for f_name in opt_files
            # read dataset and join to the original dataset
            d = CSV.read(f_name, DataFrame)
            if isempty(opt_data)
                opt_data = d
            else
                opt_data = innerjoin(opt_data, d, on = "time")
            end
        end
    end

    # process main fields of the assets to populate the dictionary with all filled data
    market_data = market(data)
    users_data = users(data)

    function change_profile!(d_dict, data_profiles)
        profile_dict = profiles(d_dict)
        if !isnothing(profile_dict) && length(profile_dict) > 0
            for p_name in keys(profile_dict)
                profile_dict[p_name] = parse_dataprofile(
                    gen_data,
                    opt_data,
                    p_name,
                    profile_dict[p_name],
                )
            end
        end
    end

    change_profile!(market_data, opt_data)

    for u_name in keys(users_data)
        comp_dict = components(users_data[u_name])
        if !isnothing(comp_dict)
            for c_name in keys(comp_dict)
                change_profile!(comp_dict[c_name], opt_data)
            end
        end
    end

    return data
end

"Return main data elements of the dataset: general parameters, users data and market data"
function explode_data(data)
    return general(data), users(data), market(data)
end

parse_to_float(x::AbstractString) = parse(Float64, x)
parse_to_float(x::Any) = Float64(x)


"Function to parse the peak power categories and tariff"
function parse_peak_quantity_by_time_vectors(
    gen_config,
    data,
    profile_name,
    peak_categories,
    peak_tariffs,
)
    # initialization of output dictionary
    peak_tariffs_by_category = Dict{String,Float64}()

    for (p_cat, t_value) in zip(peak_categories, peak_tariffs)
        if p_cat ∈ keys(peak_tariffs_by_category) &&
           peak_tariffs_by_category[p_cat] != t_value
            throw(
                ErrorException(
                    "Peak tariff category $p_cat corresponds multiple prices (e.g. $t_value and $(peak_tariffs_by_category[p_cat])",
                ),
            )
        elseif p_cat ∉ keys(peak_tariffs_by_category)
            peak_tariffs_by_category[p_cat] = parse_to_float(t_value)
        end
    end

    return peak_tariffs_by_category
end


"""
    jump_to_dict

Function to turn a JuMP model to a dictionary
"""
function jump_to_dict(model::Model)
    results = Dict{Symbol,Any}()

    for key_model in keys(model.obj_dict)
        push!(results, key_model => value.(model[key_model]))
    end

    results
end



"""
    calculate_energy_ratios(users_data, _P_ren_us, user_agg_set, agg_id, time_set, _P_tot_us, _x_us)

Calculate energy ratios
'''
# Outputs
- PV_frac
- PV_frac_tot
- wind_frac
- wind_frac_tot
'''
"""
function calculate_energy_ratios(
    users_data,
    user_agg_set,
    agg_id,
    time_set,
    _P_ren_us,
    _P_tot_us,
    _x_us,
)

    user_set = setdiff(user_agg_set, agg_id)

    # PV fraction of the aggregate case noagg
    PV_frac_tot =
        sum(
            (
                length(users_data[u].asset_type) == 0 ||
                !any(
                    a_type == REN for (name, a_type) in users_data[u].asset_type
                )
            ) ? 0.0 :
            (_P_ren_us[u, t] <= 0) ? 0.0 :
            _P_ren_us[u, t] * sum(
                Float64[
                    users_data[u].ren_pu[pv][t] * _x_us[u, pv] for
                    pv in users_data[u].asset_names if occursin("pv", pv)
                ],
            ) / sum(
                users_data[u].ren_pu[r][t] * _x_us[u, r] for
                r in users_data[u].asset_names if
                users_data[u].asset_type[r] == REN
            ) for u in user_set, t in time_set
        ) / sum(users_data[u].load[t] for u in user_set, t in time_set)

    # fraction of PV production with respect to demand by user (agg) noagg case
    PV_frac = JuMP.Containers.DenseAxisArray(
        vcat(
            PV_frac_tot,
            [
                (
                    length(users_data[u].asset_type) == 0 ||
                    !any(
                        a_type == REN for
                        (name, a_type) in users_data[u].asset_type
                    )
                ) ? 0.0 :
                sum(
                    (_P_ren_us[u, t] <= 0) ? 0.0 :
                    _P_ren_us[u, t] * sum(
                        Float64[
                            users_data[u].ren_pu[pv][t] * _x_us[u, pv] for
                            pv in users_data[u].asset_names if
                            occursin("pv", pv)
                        ],
                    ) / sum(
                        users_data[u].ren_pu[r][t] * _x_us[u, r] for
                        r in users_data[u].asset_names if
                        users_data[u].asset_type[r] == REN
                    ) for t in time_set
                ) / sum(users_data[u].load[t] for t in time_set) for
                u in user_set
            ],
        ),
        user_agg_set,
    )

    # wind fraction of the aggregate case noagg
    wind_frac_tot =
        sum(
            (
                length(users_data[u].asset_type) == 0 ||
                !any(
                    a_type == REN for (name, a_type) in users_data[u].asset_type
                )
            ) ? 0.0 :
            (_P_ren_us[u, t] <= 0) ? 0.0 :
            _P_ren_us[u, t] * sum(
                Float64[
                    users_data[u].ren_pu[w][t] * _x_us[u, w] for
                    w in users_data[u].asset_names if occursin("wind", w)
                ],
            ) / sum(
                users_data[u].ren_pu[r][t] * _x_us[u, r] for
                r in users_data[u].asset_names if
                users_data[u].asset_type[r] == REN
            ) for u in user_set, t in time_set
        ) / sum(users_data[u].load[t] for u in user_set, t in time_set)

    # fraction of wind production with respect to demand by user (agg)noagg case
    wind_frac = JuMP.Containers.DenseAxisArray(
        vcat(
            wind_frac_tot,
            [
                (
                    length(users_data[u].asset_type) == 0 ||
                    !any(
                        a_type == REN for
                        (name, a_type) in users_data[u].asset_type
                    )
                ) ? 0.0 :
                sum(
                    (_P_ren_us[u, t] <= 0) ? 0.0 :
                    _P_ren_us[u, t] * sum(
                        Float64[
                            users_data[u].ren_pu[pv][t] * _x_us[u, pv] for
                            pv in users_data[u].asset_names if
                            occursin("wind", pv)
                        ],
                    ) / sum(
                        users_data[u].ren_pu[r][t] * _x_us[u, r] for
                        r in users_data[u].asset_names if
                        users_data[u].asset_type[r] == REN
                    ) for t in time_set
                ) / sum(users_data[u].load[t] for t in time_set) for
                u in user_set
            ],
        ),
        user_agg_set,
    )

    return PV_frac, wind_frac
end


"""
    calculate_grid_ratios_noagg(users_data, user_agg_set, agg_id, time_set, _P_tot_us_noagg)

Calculate energy ratios
'''
# Outputs
- grid_frac_noagg
- grid_frac_tot_noagg
'''
"""
function calculate_grid_ratios_noagg(
    users_data,
    user_agg_set,
    agg_id,
    time_set,
    _P_tot_us_noagg,
)

    user_set = setdiff(user_agg_set, agg_id)

    # fraction of grid resiliance of the aggregate case noagg
    grid_frac_tot_noagg =
        sum(max(-_P_tot_us_noagg[u, t], 0) for u in user_set, t in time_set) /
        sum(users_data[u].load[t] for u in user_set, t in time_set)

    # fraction of grid reliance with respect to demand by user noagg case
    grid_frac_noagg = JuMP.Containers.DenseAxisArray(
        vcat(
            grid_frac_tot_noagg,
            [
                sum(max(-_P_tot_us_noagg[u, t], 0) for t in time_set) /
                sum(users_data[u].load[t] for t in time_set) for
                u in user_set
            ],
        ),
        user_agg_set,
    )

    return grid_frac_noagg
end

"""
    calculate_grid_ratios_noagg(users_data, user_agg_set, agg_id, time_set, _P_tot_us_agg)

Calculate energy ratios
'''
# Outputs
- grid_frac_agg
- grid_frac_tot_agg
'''
"""
function calculate_grid_ratios_agg(
    users_data,
    user_agg_set,
    agg_id,
    time_set,
    _P_tot_us_agg,
)

    user_set = setdiff(user_agg_set, agg_id)

    # fraction of grid reliance with respect to demand of the aggregate agg case
    grid_frac_tot_agg =
        sum(
            max(-sum(_P_tot_us_agg[u, t] for u in user_set), 0) for
            t in time_set
        ) / sum(users_data[u].load[t] for u in user_set, t in time_set)

    # fraction of grid reliance with respect to demand by user agg case
    grid_frac_agg = JuMP.Containers.DenseAxisArray(
        vcat(
            grid_frac_tot_agg,
            [
                sum(max(-_P_tot_us_agg[u, t], 0) for t in time_set) /
                sum(users_data[u].load[t] for t in time_set) for
                u in user_set
            ],
        ),
        user_agg_set,
    )

    return grid_frac_agg
end


"""
    calculate_shared_energy(users_data, user_agg_set, agg_id, time_set,
        _P_ren_us, _P_tot_us, shared_en_frac, shared_cons_frac)

Calculate the shared produced energy (en) and the shared consumption (cons) ratios

'''
# Outputs
- shared_en_frac_us_agg
- shared_en_tot_frac_agg
- shared_cons_frac_us_agg
- shared_cons_tot_frac_agg
'''
"""
function calculate_shared_energy_agg(
    users_data,
    user_agg_set,
    agg_id,
    time_set,
    _P_tot_us_agg,
    _P_ren_us_agg,
)

    # total sum of power sold by users in agg case
    _P_tot_P_sum_us_agg = JuMP.Containers.DenseAxisArray(
        [sum(max(_P_tot_us_agg[u, t], 0) for u in user_set) for t in time_set],
        time_set,
    )
    # total sum of power bought by users in agg case
    _P_tot_N_sum_us_agg = JuMP.Containers.DenseAxisArray(
        [sum(max(-_P_tot_us_agg[u, t], 0) for u in user_set) for t in time_set],
        time_set,
    )

    # fraction of shared energy with respect to total sold energy by users in the agg case
    shared_en_frac_agg = JuMP.Containers.DenseAxisArray(
        [
            (_P_tot_P_sum_us_agg[t] > 0) ?
            (
                _P_tot_P_sum_us_agg[t] -
                max(sum(_P_tot_us_agg[u, t] for u in user_set), 0)
            ) / _P_tot_P_sum_us_agg[t] : 0.0 for t in time_set
        ],
        time_set,
    )

    # fraction of shared consumption with respect to total bought energy by users in the agg case
    shared_cons_frac_agg = JuMP.Containers.DenseAxisArray(
        [
            (_P_tot_N_sum_us_agg[t] > 0) ?
            (
                _P_tot_P_sum_us_agg[t] -
                max(sum(_P_tot_us_agg[u, t] for u in user_set), 0)
            ) / _P_tot_N_sum_us_agg[t] : 0.0 for t in time_set
        ],
        time_set,
    )

    # fraction of energy share reliance with respect to demand by user aggnoagg case
    shared_en_tot_frac_agg =
        sum(
            (
                isempty(asset_names(users_data[u])) ||
                isempty(asset_names(users_data[u], GENS))
            ) ? 0.0 : shared_en_frac_agg[t] * max(_P_tot_us_agg[u, t], 0) for
            u in user_set, t in time_set
        ) / sum(_P_ren_us_agg[u, t] for u in user_set, t in time_set)

    # fraction of energy share reliance with respect to demand by user aggnoagg case
    shared_en_frac_agg_tot = JuMP.Containers.DenseAxisArray(
        vcat(
            shared_en_tot_frac_agg,
            [
                (
                    isempty(asset_names(users_data[u])) ||
                    isempty(asset_names(users_data[u], GENS))
                ) ? 0.0 :
                sum(
                    shared_en_frac_agg[t] * max(_P_tot_us_agg[u, t], 0) for
                    t in time_set
                ) / sum(_P_ren_us_agg[u, t] for t in time_set) for
                u in user_set
            ],
        ),
        user_agg_set,
    )

    # fraction of shared consumption reliance with respect to demand by user aggnoagg case
    shared_cons_tot_frac_agg =
        sum(
            shared_cons_frac_agg[t] * max(-_P_tot_us_agg[u, t], 0) for
            u in user_set, t in time_set
        ) / sum(users_data[u].load[t] for u in user_set, t in time_set)

    # fraction of shared consumption reliance with respect to demand by user aggnoagg case
    shared_cons_frac_tot_agg = JuMP.Containers.DenseAxisArray(
        vcat(
            shared_cons_tot_frac_agg,
            [
                sum(
                    shared_cons_frac_agg[t] * max(-_P_tot_us_agg[u, t], 0)
                    for t in time_set
                ) / sum(users_data[u].load[t] for t in time_set) for
                u in user_set
            ],
        ),
        user_agg_set,
    )

    return shared_en_frac_agg_tot, shared_cons_frac_tot_agg
end

"""
    calculate_shared_energy_abs_agg(users_data, user_set, time_set,
        _P_ren_us, _P_tot_us, shared_en_frac, shared_cons_frac)

Calculate the absolute shared produced energy (en), the shared consumption (cons) and self consumption (self cons)

'''
# Outputs
- shared_en_frac_us_agg
- shared_en_tot_frac_agg
- shared_cons_frac_us_agg
- shared_cons_tot_frac_agg
- self_cons_frac_us_agg
- self_cons_tot_frac_agg
'''
"""
function calculate_shared_energy_abs_agg(
    users_data,
    user_set,
    time_set,
    _P_tot_us_agg,
    _P_ren_us_agg,
)

    _P_tot_P_sum_us_agg = JuMP.Containers.DenseAxisArray(
        [sum(max(_P_tot_us_agg[u, t], 0) for u in user_set) for t in time_set],
        time_set,
    )
    # total sum of power bought by users in agg case
    _P_tot_N_sum_us_agg = JuMP.Containers.DenseAxisArray(
        [sum(max(-_P_tot_us_agg[u, t], 0) for u in user_set) for t in time_set],
        time_set,
    )

    # fraction of shared energy with respect to total sold energy by users in the agg case
    shared_en_frac_agg = JuMP.Containers.DenseAxisArray(
        [
            (_P_tot_P_sum_us_agg[t] > 0) ?
            (
                _P_tot_P_sum_us_agg[t] -
                max(sum(_P_tot_us_agg[u, t] for u in user_set), 0)
            ) / _P_tot_P_sum_us_agg[t] : 0.0 for t in time_set
        ],
        time_set,
    )

    # fraction of shared consumption with respect to total bought energy by users in the agg case
    shared_cons_frac_agg = JuMP.Containers.DenseAxisArray(
        [
            (_P_tot_N_sum_us_agg[t] > 0) ?
            (
                _P_tot_P_sum_us_agg[t] -
                max(sum(_P_tot_us_agg[u, t] for u in user_set), 0)
            ) / _P_tot_N_sum_us_agg[t] : 0.0 for t in time_set
        ],
        time_set,
    )

    # energy share reliance with respect to demand by user aggnoagg case
    shared_en_us_agg = JuMP.Containers.DenseAxisArray(
        [
            (
                isempty(asset_names(users_data[u])) ||
                isempty(asset_names(users_data[u], GENS))
            ) ? 0.0 :
            sum(
                shared_en_frac_agg[t] * max(_P_tot_us_agg[u, t], 0) for
                t in time_set
            ) for u in user_set
        ],
        user_set,
    )

    # energy share reliance with respect to demand by user aggnoagg case
    shared_en_tot_agg = sum(
        (
            isempty(asset_names(users_data[u])) ||
            isempty(asset_names(users_data[u], GENS))
        ) ? 0.0 : shared_en_frac_agg[t] * max(_P_tot_us_agg[u, t], 0) for
        u in user_set, t in time_set
    )

    # shared consumption reliance with respect to demand by user aggnoagg case
    shared_cons_us_agg = JuMP.Containers.DenseAxisArray(
        [
            sum(
                shared_cons_frac_agg[t] * max(-_P_tot_us_agg[u, t], 0) for
                t in time_set
            ) for u in user_set
        ],
        user_set,
    )

    # shared consumption reliance with respect to demand by user aggnoagg case
    shared_cons_tot_agg = sum(
        shared_cons_frac_agg[t] * max(-_P_tot_us_agg[u, t], 0) for
        u in user_set, t in time_set
    )

    return shared_en_us_agg,
    shared_en_tot_agg,
    shared_cons_us_agg,
    shared_cons_tot_agg,
    shared_en_frac_agg,
    shared_cons_frac_agg
end