using Formatting, JuMP, Plots, DataStructures
import Gurobi, XLSX, DataFrames

include("utils.jl")

function build_model(data, user_set::Vector=Vector(); output_divider=1000)

    # get main parameters
    gen_data, users_data, market_data = explode_data(data)
    n_users = length(users_data)
    init_step = field(gen_data, "init_step")
    final_step = field(gen_data, "final_step")
    n_steps = final_step - init_step + 1
    project_lifetime = field(gen_data, "project_lifetime")
    peak_categories = profile(market_data, "peak_categories")

    # Set definitions

    year_set = 1:project_lifetime
    time_set = init_step:final_step
    peak_set = unique(peak_categories)

    # Set definition when optional value is not included
    if isempty(user_set)
        user_set = user_names(gen_data, users_data)
    end

    ## Model definition

    # Definition of JuMP model
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))

    ## Constant expressions

    @variable(
        model,
        0 <=
        E_batt_us[
            u=user_set,
            b=asset_names(users_data[u], BATT),
            t=time_set,
        ] <=
        field_component(users_data[u], b, "max_capacity")
    )  # Energy stored in the battery

    @variable(
        model,
        0 <=
        P_conv_P_us[
            u=user_set,
            c=asset_names(users_data[u], CONV),
            time_set,
        ] <=
        field_component(users_data[u], c, "max_capacity")
    )  # Converter dispatch positive when supplying to AC

    @variable(
        model,
        0 <=
        P_conv_N_us[
            u=user_set,
            c=asset_names(users_data[u], CONV),
            time_set,
        ] <=
        field_component(users_data[u], c, "max_capacity")
    )  # Converter dispatch positive when absorbing from AC

    @variable(
        model,
        0 <=
        P_ren_us[u=user_set, time_set] <=
        sum(
            Float64[
                field_component(users_data[u], r, "max_capacity") for
                r in asset_names(users_data[u], REN)
            ],
        )
    )  # Dispath of renewable assets

    @variable(
        model,
        0 <=
        x_us[u=user_set, a=device_names(users_data[u])] <=
        field_component(users_data[u], a, "max_capacity")
    )  # Design of assets of the user


    # NPV of the aggregator
    @variable(model, NPV_agg >= 0)

    ## Expressions

    # CAPEX by user and asset
    @expression(
        model,
        CAPEX_us[u in user_set, a in device_names(users_data[u])],
        x_us[u, a] * field_component(users_data[u], a, "CAPEX_lin")
    )
    # CAPEX tot by user
    @expression(
        model,
        CAPEX_tot_us[u in user_set],
        sum(CAPEX_us[u, a] for a in device_names(users_data[u]))
    )

    # Maintenance cost by user and asset
    @expression(
        model,
        C_OEM_us[u in user_set, a in device_names(users_data[u])],
        x_us[u, a] * field_component(users_data[u], a, "OEM_lin")
    )
    # Maintenance cost by user
    @expression(
        model,
        C_OEM_tot_us[u in user_set],
        sum(C_OEM_us[u, a] for a in device_names(users_data[u]))
    )

    # Replacement cost by year, user and asset
    @expression(
        model,
        C_REP_us[
            y in year_set,
            u in user_set,
            a in device_names(users_data[u]),
        ],
        (
            mod(y, field_component(users_data[u], a, "lifetime_y")) == 0 &&
            y != project_lifetime
        ) ? CAPEX_us[u, a] : 0.0
    )
    # Replacement cost by year and user
    @expression(
        model,
        C_REP_tot_us[y in year_set, u in user_set],
        (y != project_lifetime) ?
        sum(
            GenericAffExpr{Float64,VariableRef}[
                C_REP_us[y, u, a] for a in device_names(users_data[u])
            ],
        ) : 0.0
    )
    # Replacement cost by year, user and asset
    @expression(
        model,
        C_RV_us[y in year_set, u in user_set, a in device_names(users_data[u])],
        (
            y == project_lifetime &&
            mod(y, field_component(users_data[u], a, "lifetime_y")) != 0
        ) ?
        CAPEX_us[u, a] * (
            1.0 -
            mod(y, field_component(users_data[u], a, "lifetime_y")) /
            field_component(users_data[u], a, "lifetime_y")
        ) : 0.0
    )
    # Replacement cost by year and user
    @expression(
        model,
        R_RV_tot_us[y in year_set, u in user_set],
        sum(
            GenericAffExpr{Float64,VariableRef}[
                C_RV_us[y, u, a] for a in device_names(users_data[u])
            ],
        )
    )

    # Economic balance of each user with respect to the public market
    @expression(
        model,
        R_Energy_us[u in user_set, t in time_set],
        profile(market_data, "energy_weight")[t] *
        profile(market_data, "time_res")[t] *
        (
            profile(market_data, "sell_price")[t] *
            (P_public_P_us[u, t] + P_micro_P_us[u, t]) -
            profile(market_data, "buy_price")[t] *
            (P_public_N_us[u, t] + P_micro_N_us[u, t]) -
            profile(market_data, "consumption_price")[t] * sum(
                Float64[
                    profile_component(users_data[u], l, "load")[t] for
                    l in asset_names(users_data[u], LOAD)
                ],
            )
        )  # economic flow with the market
    )

    # Total reward awarded to the community at each time step
    @expression(
        model,
        R_Reward_agg[t in time_set],
        profile(market_data, "energy_weight")[t] *
        profile(market_data, "time_res")[t] *
        profile(market_data, "reward_price")[t] *
        sum(P_micro_N_us[u, t] for u in user_set)
    )

    # Total reward awarded to the community by year
    @expression(model, R_Reward_agg_tot, sum(R_Reward_agg))

    # Total reward awarded to the community in NPV terms
    @expression(
        model,
        R_Reward_agg_NPV,
        R_Reward_agg_tot *
        sum(1 / ((1 + field(gen_data, "d_rate"))^y) for y in year_set)
    )

    # Economic balance of the group of each user with respect to the public market
    @expression(
        model,
        R_Energy_tot_us[u in user_set],
        sum(R_Energy_us[u, t] for t in time_set)
    )

    # Economic balance of the aggregation with respect to the public market in every hour
    @expression(
        model,
        R_Energy_agg_time[t in time_set],
        sum(R_Energy_us[u, t] for u in user_set)
    )

    # Yearly revenue of the user
    @expression(
        model,
        yearly_rev[u=user_set],
        R_Energy_tot_us[u] - C_OEM_tot_us[u]
    )

    # Cash flow
    @expression(
        model,
        Cash_flow_us[u in user_set, y in append!([0], year_set)],
        (y == 0) ? 0.0 - CAPEX_tot_us[u] : # the investment costs CAPEX
        (
            # in `ECNetworkBlock` (done!)
            R_Energy_tot_us[u] - # net economic balance wrt the public market
            C_Peak_tot_us[u] - # the costs due to the peak power
            # in `UnitBlock`
            C_OEM_tot_us[u] - # the fuel costs
            C_REP_tot_us[y, u] + # the replacement costs of the assets
            R_RV_tot_us[y, u] # the residual value of the assets
        )
    )

    @expression(
        model,
        NPV_us[u in user_set],
        sum(
            Cash_flow_us[u, y] / ((1 + field(gen_data, "d_rate"))^y) for
            y in append!([0], year_set)
        )
    )

    @expression(
        model,
        Cash_flow_agg[y in append!([0], year_set)],
        (y == 0) ? 0.0 : R_Reward_agg_tot
    )

    @expression(
        model,
        Cash_flow_tot[y in append!([0], year_set)],
        sum(Cash_flow_us[u, y] for u in user_set) + Cash_flow_agg[y]
    )

    # Social welfare of the entire aggregation
    @expression(model, SW, sum(NPV_us) + R_Reward_agg_NPV)

    # Social welfare of the users
    @expression(model, SW_us, SW - NPV_agg)

    # Other expression

    # Total converters dispatch when supplying to the grid
    @expression(
        model,
        P_conv_P_tot_us[u=user_set, t=time_set],
        sum(P_conv_P_us[u, c, t] for c in asset_names(users_data[u], CONV))
    )

    # Total converters dispatch when absorbing from the grid
    @expression(
        model,
        P_conv_N_tot_us[u=user_set, t=time_set],
        sum(P_conv_N_us[u, c, t] for c in asset_names(users_data[u], CONV))
    )

    # Total converters dispatch by user
    @expression(
        model,
        P_conv_tot_us[u=user_set, t=time_set],
        P_conv_P_tot_us[u, t] - P_conv_N_tot_us[u, t]
    )

    # Total converters dispatch by user and type of component
    @expression(
        model,
        P_conv_us[
            u=user_set,
            c=asset_names(users_data[u], CONV),
            t=time_set,
        ],
        P_conv_P_us[u, c, t] - P_conv_N_us[u, c, t]
    )

    # Total converters dispatch
    @expression(
        model,
        P_conv_tot_agg[t=time_set],
        sum(P_conv_tot_us[u, t] for u in user_set)
    )

    # Total energy available in the batteries
    @expression(
        model,
        E_batt_tot_us[u=user_set, t=time_set],
        sum(E_batt_us[u, b, t] for b in asset_names(users_data[u], BATT))
    )

    ## Constraints

    ## Inequality constraints

    # Set the renewable energy dispatch to be no greater than the actual available energy
    @constraint(
        model,
        con_us_ren_dispatch[u in user_set, t in time_set],
        P_ren_us[u, t] <= sum(
            GenericAffExpr{Float64,VariableRef}[
                profile_component(users_data[u], r, "ren_pu")[t] * x_us[u, r]
                for r in asset_names(users_data[u], REN)
            ],
        )
    )

    # Set the maximum hourly dispatch of converters not to exceed their capacity
    @constraint(
        model,
        con_us_converter_capacity[
            u in user_set,
            c in asset_names(users_data[u], CONV),
            t in time_set,
        ],
        P_conv_P_us[u, c, t] + P_conv_N_us[u, c, t] <= x_us[u, c]
    )

    # Set the maximum hourly dispatch of converters not to exceed the C-rate of the battery in discharge
    @constraint(
        model,
        con_us_converter_capacity_crate_dch[
            u in user_set,
            c in asset_names(users_data[u], CONV),
            t in time_set,
        ],
        P_conv_P_us[u, c, t] <=
        x_us[u, field_component(users_data[u], c, "corr_asset")] *
        field_component(
            users_data[u],
            field_component(users_data[u], c, "corr_asset"),
            "max_C_dch",
        )
    )

    # Set the maximum hourly dispatch of converters not to exceed the C-rate of the battery in charge
    @constraint(
        model,
        con_us_converter_capacity_crate_ch[
            u in user_set,
            c in asset_names(users_data[u], CONV),
            t in time_set,
        ],
        P_conv_N_us[u, c, t] <=
        x_us[u, field_component(users_data[u], c, "corr_asset")] *
        field_component(
            users_data[u],
            field_component(users_data[u], c, "corr_asset"),
            "max_C_ch",
        )
    )

    # Set the minimum level of the energy stored in the battery to be proportional to the capacity
    @constraint(
        model,
        con_us_min_E_batt[
            u in user_set,
            b in asset_names(users_data[u], BATT),
            t in time_set,
        ],
        E_batt_us[u, b, t] >=
        x_us[u, b] * field_component(users_data[u], b, "min_SOC")
    )

    # Set the maximum level of the energy stored in the battery to be proportional to the capacity
    @constraint(
        model,
        con_us_max_E_batt[
            u in user_set,
            b in asset_names(users_data[u], BATT),
            t in time_set,
        ],
        E_batt_us[u, b, t] <=
        x_us[u, b] * field_component(users_data[u], b, "max_SOC")
    )

    ## Equality constraints

    # Set the balance at each battery system
    @constraint(
        model,
        con_us_bat_balance[
            u in user_set,
            b in asset_names(users_data[u], BATT),
            t in time_set,
        ],
        E_batt_us[u, b, t] - E_batt_us[u, b, pre(t, time_set)] +  # Difference between the energy level in the battery. Note that in the case of the first time step, the last id is used
        + profile(market_data, "time_res")[t] *
        P_conv_P_us[u, field_component(users_data[u], b, "corr_asset"), t] / (
            sqrt(field_component(users_data[u], b, "eta")) * field_component(
                users_data[u],
                field_component(users_data[u], b, "corr_asset"),
                "eta",
            )
        ) -
        profile(market_data, "time_res")[t] *
        P_conv_N_us[u, field_component(users_data[u], b, "corr_asset"), t] *
        (
            sqrt(field_component(users_data[u], b, "eta")) * field_component(
                users_data[u],
                field_component(users_data[u], b, "corr_asset"),
                "eta",
            )
        ) == 0
    )

    # Set the objective of maximizing the profit of the aggregator
    # @objective(model, Max, model[:NPV_agg] / output_divider)
    @objective(model, Max, model[:SW] / output_divider)

    return model
end