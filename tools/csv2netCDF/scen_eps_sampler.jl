# Short-period (eps) sampler used by both `csv2nc4.jl` (when writing the SMS++
# TSSB netCDF) and `test_instance_with_EC_jl.jl` (when computing reference
# objective values via EnergyCommunity.jl@stochastic). It produces a
# `Scenario_Load_Renewable` for each (s, eps) draw, perturbing the relevant
# components around the long-period centre points provided by `pem_extraction`.
#
# Schema aligned with EnergyCommunity.jl@stochastic; see
# `examples/RunStochModel(TBD).jl` of
# https://github.com/SPSUnipi/EnergyCommunity.jl branch `stochastic`.

"""
    perturb_price_vector(data_market, name, time_set)

Per-time-step sample of the market-level price `name` (e.g. `"buy_price"`).
When the market profile defines `std_<name>` the deterministic mean is
perturbed via `MvNormal(mean, std)` and folded to non-negative values;
without `std_<name>` the deterministic mean is returned untouched.
"""
function perturb_price_vector(data_market, name::AbstractString, time_set)
	mean_vec = profile_d(data_market, name, fill(0.0, length(time_set)))
	std_vec = profile_d(data_market, "std_" * name, nothing)
	isnothing(std_vec) && return mean_vec
	return broadcast(abs, rand(MvNormal(mean_vec, std_vec)))
end

"""
    perturb_peak_tariff(data_market)

Sample the peak-tariff Dict (peak category → tariff). When the market profile
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

# `unc_var` is a string of letters in {"L","P","W"} selecting which
# short-period sources of uncertainty to perturb (Load / PV / Wind). For a
# letter not in `unc_var` the sample equals the long-period centre point —
# so that the per-scenario data fed to the LP coincides with the baseline
# profile that the SMS++ pipeline uses to set static bounds (e.g.
# `IntermittentUnitBlock.MaxPower` and the derived `max_node_injection`).

@sampler Scenario_eps_Sampler = begin
	data_market::Dict{Any, Any}      # market data
	data_user::Dict{Any, Any}        # users data
	point_s_load::Array{Any}         # long-period (s) centre points for load
	point_s_pv::Array{Any}           # long-period (s) centre points for PV
	point_s_wind::Array{Any}         # long-period (s) centre points for wind
	scen_s::Int                      # current `s` index
	unc_var::String                  # uncertain sources, subset of "LPW"

	Scenario_eps_Sampler(data_market, data_user, point_s_load, point_s_pv,
	                     point_s_wind, scen_s, unc_var) =
		new(data_market, data_user, point_s_load, point_s_pv, point_s_wind,
		    scen_s, unc_var)

	@sample Scenario_Load_Renewable begin
		scen_s        = sampler.scen_s
		data_user     = sampler.data_user
		data_market   = sampler.data_market
		point_s_load  = sampler.point_s_load
		point_s_pv    = sampler.point_s_pv
		point_s_wind  = sampler.point_s_wind
		unc_var       = sampler.unc_var

		sample_L = occursin("L", unc_var)
		sample_P = occursin("P", unc_var)
		sample_W = occursin("W", unc_var)

		load_demand    = Dict{String, Dict{Int, Float64}}()
		ren_production = Dict{String, Dict{String, Dict{Int, Float64}}}()

		for u in user_set

			# Load: scaled by the long-period centre point, then optionally
			# perturbed via `MvNormal(mean, std)` (folded to non-negative).
			# The standard deviation comes from the user-specific `std`
			# profile when present; otherwise it falls back to
			# `|mean| * sigma_load` (a `general.sigma_load` global, default 0.3).
			load_mean = profile_component(data_user[u], "load", "load")
			load_scenario_s = point_s_load[scen_s] * load_mean
			if sample_L
				load_std = profile_component_d(data_user[u], "load", "std",
				                               abs.(load_mean) .* sigma_load)
				load_distribution = MvNormal(load_scenario_s, load_std)
				load_demand[u] = array2dict(broadcast(abs, rand(load_distribution)))
			else
				load_demand[u] = array2dict(load_scenario_s)
			end

			# Renewables (PV / wind): same structure as load, with the
			# long-period centre and the std/sigma falling back per family.
			ren_production[u] = Dict{String, Dict{Int, Float64}}()
			for name = asset_names(data_user[u], REN)
				asset_type = lowercase(name)  # "pv" / "wind" / ...
				is_pv      = occursin("pv",   asset_type)
				is_wind    = occursin("wind", asset_type)

				point_centre = is_pv   ? point_s_pv[scen_s] :
				               is_wind ? point_s_wind[scen_s] :
				                         1.0  # unknown family → deterministic

				ren_mean = profile_component(data_user[u], name, "ren_pu")
				ren_scenario_s = point_centre * ren_mean

				perturb = (is_pv && sample_P) || (is_wind && sample_W)
				if perturb
					ren_std = profile_component_d(data_user[u], name, "std",
					                              abs.(ren_mean) .* sigma_pv)
					ren_distribution = MvNormal(ren_scenario_s, ren_std)
					array_ren = broadcast(abs, rand(ren_distribution))
				else
					array_ren = ren_scenario_s
				end
				temp = array2dict(array_ren)

				# Force production to zero where the deterministic profile is
				# zero (e.g. solar at night).
				for t in time_set
					if profile_component(data_user[u], name, "ren_pu")[t] == 0
						temp[t] = 0
					end
				end
				get!(ren_production[u], name, temp)
			end
		end

		# Market-level prices: each is perturbed via truncated-Normal noise
		# when the YAML market profile defines the matching `std_<field>`
		# entry, otherwise it is left deterministic.
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

"""
    scenarios_generator(data, point_s_load, point_s_pv, point_s_wind,
                        n_scen_s, n_scen_eps, unc_var;
                        point_probability=ones(n_scen_s),
                        control_risimulation=false)

Build the full `(s, eps)` scenario set. For each long-period `s` a
`Scenario_eps_Sampler` is instantiated and queried `n_scen_eps` times to
produce the short-period draws; the returned scenarios carry probability
`point_probability[s] / n_scen_eps`.

A deterministic per-`s` seed (`Random.seed!(123 + s)`) is set immediately
before instantiating the sampler so that `csv2nc4.jl` and
`test_instance_with_EC_jl.jl` consume bit-identical random sequences when
both seed `Random` with the same value, regardless of any prior RNG
divergence between the two callers.

When `control_risimulation` is `true`, `n_scen_s` is interpreted as the
single `s` index to re-simulate; the generator returns `n_scen_eps`
scenarios all carrying that index, with uniform probability `1/n_scen_eps`.
"""
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
			Random.seed!(123 + s)
			sampler_eps = Scenario_eps_Sampler(market(data), users(data),
			                                   point_s_load, point_s_pv,
			                                   point_s_wind, s, unc_var)
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
	end

	# Re-simulation path: a single `s = n_scen_s` index, `n_scen_eps` draws,
	# uniform probability.
	sampled_scenarios = Array{Scenario_Load_Renewable}(undef, n_scen_eps)
	for scen = 1:n_scen_eps
		sampled_scenarios[scen] = zero(Scenario_Load_Renewable)
	end

	s = n_scen_s
	sampler_eps = Scenario_eps_Sampler(market(data), users(data),
	                                   point_s_load, point_s_pv,
	                                   point_s_wind, s, unc_var)
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
