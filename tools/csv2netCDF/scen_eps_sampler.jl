# Sampler for distribution associated to short period uncertainty.
# Schema aligned with EC.jl@stochastic.

# ---------------------------------------------------------------------------
# Optional stochastic perturbation of market-level prices.
# A price profile becomes scenario-dependent if the YAML market profile
# contains the matching `std_<name>` entry (or `std_peak_tariff` for the
# Dict-typed peak tariff). Without it the price stays deterministic.
# ---------------------------------------------------------------------------

"""
    perturb_price_vector(data_market, name, time_set)

Return a per-time-step price vector for `name` (e.g. "buy_price"). When
`std_<name>` is present in the market profile the deterministic mean is
perturbed via `MvNormal(mean, std)` then folded to non-negative values;
otherwise the deterministic mean is returned untouched.
"""
function perturb_price_vector(data_market, name::AbstractString, time_set)
	mean_vec = profile_d(data_market, name, fill(0.0, length(time_set)))
	std_vec = profile_d(data_market, "std_" * name, nothing)
	isnothing(std_vec) && return mean_vec
	return broadcast(abs, rand(MvNormal(mean_vec, std_vec)))
end

"""
    perturb_peak_tariff(data_market)

Return the peak-tariff Dict (category → tariff). When the market profile
contains `std_peak_tariff` (Dict{String,Float64}) each category is sampled
independently from `Normal(mean, std)`; otherwise the deterministic Dict is
returned untouched.
"""
function perturb_peak_tariff(data_market)
	mean_dict = profile_d(data_market, "peak_tariff", Dict{String, Float64}())
	std_dict = profile_d(data_market, "std_peak_tariff", nothing)
	(isnothing(std_dict) || isempty(std_dict)) && return mean_dict
	out = Dict{String, Float64}()
	for (cat, mu) in mean_dict
		sigma = get(std_dict, cat, 0.0)
		out[cat] = sigma > 0 ? abs(rand(Normal(mu, sigma))) : mu
	end
	return out
end

### Sampler definition

@sampler Scenario_eps_Sampler = begin
	data_market::Dict{Any, Any}      # data of the market
	data_user::Dict{Any, Any}        # data of users
	point_s_load::Array{Any}         # long-period (s) point for load demand
	point_s_pv::Array{Any}           # long-period (s) point for PV
	point_s_wind::Array{Any}         # long-period (s) point for wind
	scen_s::Int                      # current `s` index

	Scenario_eps_Sampler(data_market, data_user, point_s_load, point_s_pv,
	                     point_s_wind, scen_s) =
		new(data_market, data_user, point_s_load, point_s_pv, point_s_wind,
		    scen_s)

	@sample Scenario_Load_Renewable begin
		scen_s        = sampler.scen_s
		data_user     = sampler.data_user
		data_market   = sampler.data_market
		point_s_load  = sampler.point_s_load
		point_s_pv    = sampler.point_s_pv
		point_s_wind  = sampler.point_s_wind

		load_demand = Dict{String,Dict{Int,Float64}}()
		ren_production = Dict{String,Dict{String,Dict{Int,Float64}}}()

		for u in user_set

			# Load — short-period perturbation around `point_s_load[scen_s] * mean`.
			load_mean = profile_component(data_user[u], "load", "load")
			load_scenario_s = point_s_load[scen_s] * load_mean
			# `profile.std` is now mandatory in the EC.jl@stochastic schema; fall
			# back to `|mean|*sigma_load` (a `general.sigma_load` global, default
			# 0.3) when the user-specific column is absent.
			load_std = profile_component_d(data_user[u], "load", "std",
			                               abs.(load_mean) .* sigma_load)
			load_distribution = MvNormal(load_scenario_s, load_std)
			load_demand[u] = array2dict(broadcast(abs, rand(load_distribution)))

			ren_production[u] = Dict{String,Dict{Int,Float64}}()
			for name = asset_names(data_user[u], REN)
				# Pick the long-period centre point for this asset's family.
				asset_type = lowercase(name)  # "pv" / "wind" / ...
				point_centre = if occursin("pv", asset_type)
					point_s_pv[scen_s]
				elseif occursin("wind", asset_type)
					point_s_wind[scen_s]
				else
					1.0  # unknown renewable type → deterministic
				end

				ren_mean = profile_component(data_user[u], name, "ren_pu")
				ren_scenario_s = point_centre * ren_mean
				ren_std = profile_component_d(data_user[u], name, "std",
				                              abs.(ren_mean) .* sigma_pv)
				ren_distribution = MvNormal(ren_scenario_s, ren_std)

				array_ren = broadcast(abs, rand(ren_distribution))
				temp = array2dict(array_ren)

				# Force zero production at time-steps where the deterministic
				# profile is zero (e.g. solar at night).
				for t in time_set
					if profile_component(data_user[u], name, "ren_pu")[t] == 0
						temp[t] = 0
					end
				end
				get!(ren_production[u], name, temp)
			end
		end

		# Market-level price fields: each is perturbed with truncated-Normal
		# noise when the YAML market profile defines the matching `std_<field>`
		# entry (`std_buy_price`, `std_sell_price`, `std_consumption_price`,
		# `std_penalty_price`, `std_peak_tariff`); without `std_*` the price is
		# left deterministic.
		buy_price_arr         = perturb_price_vector(data_market, "buy_price",         time_set)
		sell_price_arr        = perturb_price_vector(data_market, "sell_price",        time_set)
		consumption_price_arr = perturb_price_vector(data_market, "consumption_price", time_set)
		penalty_price_arr     = perturb_price_vector(data_market, "penalty_price",     time_set)
		peak_tariff_dict      = perturb_peak_tariff(data_market)

		return Scenario_Load_Renewable(1,
		                               1,
		                               peak_tariff_dict,
		                               array2dict(buy_price_arr),
		                               array2dict(consumption_price_arr),
		                               array2dict(sell_price_arr),
		                               array2dict(penalty_price_arr),
		                               load_demand,
		                               ren_production)
	end
end

function scenarios_generator(
	data::Dict{Any, Any},
	point_s_load::Vector{Float64},
	point_s_pv::Vector{Float64},
	point_s_wind::Vector{Float64},
	n_scen_s::Int,
	n_scen_eps::Int,
	unc_var::AbstractString;
	point_probability::Vector{Float64} = ones(n_scen_s),
	control_risimulation::Bool = false,
)

	if !control_risimulation
		n_scen = n_scen_s * n_scen_eps
		sampled_scenarios = Array{Scenario_Load_Renewable}(undef, n_scen)
		for scen = 1:n_scen
			sampled_scenarios[scen] = zero(Scenario_Load_Renewable)
		end

		for s = 1:n_scen_s
			# Force a deterministic per-(s, eps) seed so csv2nc4 and the
			# test_instance harness consume bit-identical rand sequences
			# during scenario sampling, regardless of the prior RNG state
			# (which may diverge between the two callers due to subtle
			# differences in earlier macro-generated code, package init,
			# or Dict iteration order). With one seed per scenario the
			# scenarios are still pseudo-random but reproducible across
			# pipelines.
			Random.seed!(123 + s)
			sampler_eps = Scenario_eps_Sampler(market(data), users(data),
			                                   point_s_load, point_s_pv,
			                                   point_s_wind, s)
			for eps = 1:n_scen_eps
				scen = (s - 1) * n_scen_eps + eps
				scenario_sampled = sampler_eps()
				sampled_scenarios[scen] = Scenario_Load_Renewable(
					s, eps,
					scenario_sampled.peak_tariff,
					scenario_sampled.buy_price,
					scenario_sampled.consumption_price,
					scenario_sampled.sell_price,
					scenario_sampled.penalty_price,
					scenario_sampled.Load,
					scenario_sampled.Ren,
					probability = point_probability[s] / n_scen_eps,
				)
			end
		end
		return sampled_scenarios
	else
		# Re-simulation: `n_scen_s` is the index of the single `s` to consider.
		sampled_scenarios = Array{Scenario_Load_Renewable}(undef, n_scen_eps)
		for scen = 1:n_scen_eps
			sampled_scenarios[scen] = zero(Scenario_Load_Renewable)
		end

		s = n_scen_s
		sampler_eps = Scenario_eps_Sampler(market(data), users(data),
		                                   point_s_load, point_s_pv,
		                                   point_s_wind, s)
		for eps = 1:n_scen_eps
			scenario_sampled = sampler_eps()
			sampled_scenarios[eps] = Scenario_Load_Renewable(
				1, eps,
				scenario_sampled.peak_tariff,
				scenario_sampled.buy_price,
				scenario_sampled.consumption_price,
				scenario_sampled.sell_price,
				scenario_sampled.penalty_price,
				scenario_sampled.Load,
				scenario_sampled.Ren,
				probability = 1 / n_scen_eps,
			)
		end
		return sampled_scenarios
	end
end
