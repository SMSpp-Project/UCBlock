# Scenario definition

@define_scenario Scenario_Load_Renewable = begin
    scen_s::Int
    scen_eps::Int
    peak_tariff::Dict{String,Float64}
    buy_price::Dict{Int,Float64}
    consumption_price::Dict{Int,Float64}
    sell_price::Dict{Int,Float64}
    reward_price::Dict{Int,Float64}
    penalty_price::Dict{Int,Float64}
    Load::Dict{String,Dict{Int,Float64}} # load_demand of each user per each scenario (s,eps)
    Ren::Dict{String,Dict{String,Dict{Int,Float64}}} # renewable production of each users plant per each scenario (s,eps)
    @zero begin
        n = 1
        empdict = Dict{Int,Float64}()
        a = Dict{String,Dict{Int,Float64}}()
        b = Dict{String,Dict{String,Dict{Int,Float64}}}()
        for u in user_set
            a[u] = Dict{Int,Float64}()
            b[u] = Dict{String,Dict{Int,Float64}}()
        end
        return Scenario_Load_Renewable(n, n, Dict{String,Float64}(), empdict, empdict, empdict, empdict, empdict, a, b)
    end

    @expectation begin
        p_t = Dict{String,Float64}()
        b_V = Dict{Int,Float64}()
        b_F = Dict{Int,Float64}()
        s_p = Dict{Int,Float64}()
        r_p = Dict{Int,Float64}()
        p_p = Dict{Int,Float64}()
        p = Dict{String,Dict{Int,Float64}}()
        q = Dict{String,Dict{String,Dict{Int,Float64}}}
        for u in user_set
            p[u] = Dict{Int,Float64}()
            q[u] = Dict{String,Dict{Int,Float64}}()
        end
        control = false
        for t in time_set
            get!(b_V, t, sum([probability(sc) * sc.buy_price[t] for sc in scenarios]))
            get!(b_F, t, sum([probability(sc) * sc.consumption_price[t] for sc in scenarios]))
            get!(s_p, t, sum([probability(sc) * sc.sell_price[t] for sc in scenarios]))
            get!(r_p, t, sum([probability(sc) * sc.reward_price[t] for sc in scenarios]))
            get!(p_p, t, sum([probability(sc) * sc.penalty_price[t] for sc in scenarios]))
        end

        for peak in peak_set
            get!(p_t, peak, sum([probability(sc) * sc.peak_tariff[peak] for sc in scenarios]))
        end

        length_scenarios = length(scenarios)

        utils = Array{Dict{Int,Float64}}(undef, length_scenarios + 1) # variable used to store values
        for scen = 1:length_scenarios+1
            utils[scen] = Dict{Int,Float64}()
        end
        for n in time_set
            get!(utils[1], n, 0)
        end
        for u in user_set
            count = 1
            for sc in scenarios
                count = count + 1
                prob_scen = probability(sc)
                load_scen = sc.Load[u]
                utils2 = Dict{Int,Float64}() # variable used to store values
                for n in time_set
                    utils2[n] = prob_scen * load_scen[n]
                end
                utils[count] = merge(+, utils[count-1], utils2)
            end
            p[u] = utils[length_scenarios+1]
            for name = asset_names(users_data[u], REN)
                temp1 = Array{Dict{Int,Float64}}(undef, length_scenarios + 1)
                for scen = 1:length_scenarios+1
                    temp1[scen] = Dict{Int,Float64}()
                end
                for n in time_set
                    get!(temp1[1], n, 0)
                end
                count2 = 1
                for sc in scenarios
                    count2 = count2 + 1
                    temp2 = Dict{Int,Float64}() # temporary dictionary to store the expected scenario for renewable production of each asset "name"
                    ren_scen_ass = sc.Ren[u][name]
                    prob_scen = probability(sc)
                    for n in time_set
                        temp2[n] = prob_scen * ren_scen_ass[n]
                    end
                    temp1[count2] = merge(+, temp1[count2-1], temp2)
                end
                get!(q[u], name,
                    temp1[length_scenarios+1])
            end
        end
        return Scenario_Load_Renewable(1, 1, p_t, b_V, b_F, s_p, r_p, p_p, p, q)
    end
end

""" Sampler for distribution associated to short period uncertainty
	
Sampler function

INPUT: users data: dictionary containing all the user data information

OUTPUT: point_load_demand: sampled point from the normalized distribution associated with load demand
        point_ren_production: sampled point from the normalized distribution associated with renewable production

"""
function scenario_eps_point_sampler(data_user; deterministic::Bool=false)

    point_load_demand = Dict{String,Array{Float64}}() # extracted point for each user
    point_ren_production = Dict{String,Dict{String,Array{Float64}}}() # extracted point for each user and asset

    n_step = length(profile_component(data_user["user1"], "load", "load"))

    for u in user_set

        if deterministic == true
            array_n_load = ones(n_step)
            point_load_demand[u] = array_n_load

            point_ren_production[u] = Dict{String,Array{Float64}}()
            for name = asset_names(data_user[u], REN)
                point_ren = ones(n_step)
                get!(point_ren_production[u], name, point_ren)
            end
        else

            # calculate the normalized std for load
            std_n_load = (profile_component(data_user[u], "load", "std") .* 2) ./ profile_component(data_user[u], "load", "load")

            # define load distribution for short period uncertainty
            load_distribution = MvNormal(ones(n_steps), std_n_load)

            # load extraction
            array_n_load = broadcast(abs, rand(load_distribution))

            # control when the extracted point are < 0
            for t in time_set
                if array_n_load[t] < 0
                    array_n_load[t] = 0
                end
            end

            # new load demand of user u in the specific scenario s considered
            point_load_demand[u] = array_n_load

            point_ren_production[u] = Dict{String,Array{Float64}}()
            for name = asset_names(data_user[u], REN)
                # calculate the normalized std for ren production
                std_n_ren = ((name == "PV") ? profile_component(data_user[u], name, "std") .* 2 : profile_component(data_user[u], name, "std"))

                # define load distribution for short period uncertainty
                ren_distribution = MvNormal(ones(n_steps), std_n_ren)

                # renewable extraction
                point_ren = broadcast(abs, rand(ren_distribution))

                # control to set 0 the extracted production when initial renewable production was 0 or when the extracted point are < 0
                for t in time_set
                    if profile_component(data_user[u], name, "ren_pu")[t] == 0 || point_ren[t] < 0
                        point_ren[t] = 0
                    end
                end
                get!(point_ren_production[u], name, point_ren)
            end
        end
    end
    return point_load_demand, point_ren_production
end

""" Scenarios generator function

INPUT: data: dictionary containing all the data information
       point_s_load: extracted point for long period uncertainty distribution in load demand
       point_s_ren: extracted point for long period uncertainty distribution in renewable production
       n_scen_s: number of scenarios s to generate
       n_scen_eps: number of scenarios eps to generate
       first_stage: boolean value to know if we are in the first stage
       second_stage: boolean value to know if we are in the second stage
       point_eps_load: previously extracted point in short period uncertainty distributions, used only in the second phase
       point_eps_ren: previously extracted point in short period uncertainty distributions, used only in the second phase

OUTPUT: sampled_scenarios: array with the generaated scenarios
        point_load_sampled (only FIRST STAGE): array containing the point sampled for short period uncertainty, to be reused in the second phase
        point_ren_sampled (only FIRST STAGE): array containing the point sampled for short period uncertainty, to be reused in the second phase
"""

function scenarios_generator(
    data::Dict{Any,Any},
    point_s_load::Vector{Float64}, # extracted point for long period uncertainty distribution in load demand
    point_s_ren::Vector{Float64}, # extracted point for long period uncertainty distribution in renewable production
    n_scen_s::Int, # number of scenarios s to generate
    n_scen_eps::Int; # number of scenarios eps to generate
    point_probability::Vector{Float64}=Float64[], # probability of each point
    deterministic::Bool=false
)

    data_user = users(data)
    data_market = market(data)

    n_scen = n_scen_s * n_scen_eps # total number of scenarios

    # array containing each scenario
    sampled_scenarios = Array{Scenario_Load_Renewable}(undef, n_scen)
    for scen = 1:n_scen
        sampled_scenarios[scen] = zero(Scenario_Load_Renewable) # initialize an empty scenario 
    end

    point_load_sampled = Array{Any}(undef, n_scen_eps) # array used to store all the extracted point for load demand
    point_ren_sampled = Array{Any}(undef, n_scen_eps) # array used to store all the extracted point for renewable production

    # extract epsilon point and store them in the output variables
    for eps = 1:n_scen_eps
        (point_load_sampled[eps], point_ren_sampled[eps]) = scenario_eps_point_sampler(users(data), deterministic=deterministic)
    end

    for s = 1:n_scen_s
        for eps = 1:n_scen_eps
            scen = (s - 1) * n_scen_eps + eps

            # we have now to evaluate the real value sampled for load and renewable production
            load_demand = Dict{String,Dict{Int,Float64}}()
            ren_production = Dict{String,Dict{String,Dict{Int,Float64}}}()

            for u in user_set
                load_demand[u] = array2dict(point_load_sampled[eps][u] .* (point_s_load[s] * profile_component(data_user[u], "load", "load")))

                ren_production[u] = Dict{String,Dict{Int,Float64}}()
                for name = asset_names(data_user[u], REN)
                    temp = array2dict(point_ren_sampled[eps][u][name] .* (point_s_ren[s] * profile_component(data_user[u], name, "ren_pu")))

                    get!(ren_production[u], name, temp)
                end
            end

            sampled_scenarios[scen] = Scenario_Load_Renewable(
                s,
                eps,
                profile(data_market, "peak_tariff"),
                array2dict(profile(data_market, "buy_price")),
                array2dict(profile(data_market, "consumption_price")),
                array2dict(profile(data_market, "sell_price")),
                array2dict(profile(data_market, "reward_price")),
                array2dict(profile(data_market, "penalty_price")),
                load_demand,
                ren_production,
                probability=point_probability[s] / n_scen_eps)
        end
    end

    return sampled_scenarios

end
