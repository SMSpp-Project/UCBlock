# Sampler for distribution associated to short period uncertainty
	
### Sampler definition

@sampler Scenario_eps_Sampler = begin
	data_market::Dict{Any, Any} # data of the market
	data_user::Dict{Any, Any} # data of users
	point_s_load::Array{Any} # extracted point for long period uncertainty distribution in load demand
	point_s_ren::Array{Any} # extracted point for long period uncertainty distribution in renewable production
	scen_s::Int # Scenario s considered
	
	Scenario_eps_Sampler(data_market, data_user, point_s_load, point_s_ren, scen_s) = new(data_market, data_user, point_s_load, point_s_ren, scen_s)
	
	@sample Scenario_Load_Renewable begin
		scen_s = sampler.scen_s
		data_user = sampler.data_user
		data_market = sampler.data_market
		point_s_load = sampler.point_s_load
		point_s_ren = sampler.point_s_ren

		load_demand = Dict{String,Dict{Int,Float64}}()
		ren_production = Dict{String,Dict{String,Dict{Int,Float64}}}()

		for u in user_set

			# Define load distribution for short period uncertainty
			load_mean = profile_component(data_user[u], "load", "load")
			load_scenario_s = point_s_load[scen_s] * load_mean # mean load in the scenario s considered
			# Fall back to a multiplicative noise based on the long-period sigma when the
			# profile does not define an explicit std (new-format YAMLs do not).
			load_std = profile_component_d(data_user[u], "load", "std", abs.(load_mean) .* sigma_load)
			load_distribution = MvNormal(load_scenario_s, load_std)

			# Load extraction
			array_load = broadcast(abs, rand(load_distribution))

			# New load demand of user u
			load_demand[u] = array2dict(array_load)

            ren_production[u] = Dict{String,Dict{Int,Float64}}()
			for name = asset_names(data_user[u], REN)
				# Define renewable distribution for short period uncertainty
				ren_mean = profile_component(data_user[u], name, "ren_pu")
				ren_scenario_s = point_s_ren[scen_s] * ren_mean
				ren_std = profile_component_d(data_user[u], name, "std", abs.(ren_mean) .* sigma_ren)
				ren_distribution = MvNormal(ren_scenario_s, ren_std)

				# Renewable extraction
				array_ren = broadcast(abs, rand(ren_distribution))
				temp = array2dict(array_ren)

				# Control to set 0 the extracted production when initial renewable production was 0 
				for t in time_set
					if profile_component(data_user[u], name, "ren_pu")[t] == 0
						temp[t] = 0
					end
				end
				get!(ren_production[u],name,temp)
			end
		end
		# Market-level price fields are stored in the scenario object but are not propagated
		# to the netCDF DiscreteScenarioSet, so any field absent from the YAML (typical of the
		# new format, which lacks `penalty_price`) is filled with empty defaults.
		empty_dict = Dict{Int, Float64}()
		return Scenario_Load_Renewable(1,
						1,
						profile_d(data_market, "peak_tariff", Dict{String, Float64}()),
						array2dict(profile_d(data_market, "buy_price",          fill(0.0, length(time_set)))),
						array2dict(profile_d(data_market, "consumption_price",  fill(0.0, length(time_set)))),
						array2dict(profile_d(data_market, "sell_price",         fill(0.0, length(time_set)))),
						array2dict(profile_d(data_market, "penalty_price",      fill(0.0, length(time_set)))),
						load_demand,
						ren_production)
	end
end

function scenarios_generator(
	data::Dict{Any, Any},
	point_s_load::Vector{Float64}, # extracted point for long period uncertainty distribution in load demand
	point_s_ren::Vector{Float64}, # extracted point for long period uncertainty distribution in renewable production
	point_probability::Vector{Float64}, # probability of each point
	n_scen_s::Int, # number of scenarios s to generate
	n_scen_eps::Int; # number of scenarios eps to generate
	control_risimulation::Bool = false
	)
	
	if control_risimulation == false
		n_scen = n_scen_s*n_scen_eps # total number of scenarios

		# Array containing each scenario
		sampled_scenarios = Array{Scenario_Load_Renewable}(undef,n_scen)
		for scen = 1:n_scen
			sampled_scenarios[scen] = zero(Scenario_Load_Renewable) # initialize an empty scenario 
		end
		
		for s = 1:n_scen_s
			sampler_eps = Scenario_eps_Sampler(market(data),users(data),point_s_load,point_s_ren,s)
			for eps = 1:n_scen_eps
				scen = (s-1)*n_scen_eps+eps
				scenario_sampled = sampler_eps()
				sampled_scenarios[scen] = Scenario_Load_Renewable(
					s,
					eps,
					scenario_sampled.peak_tariff,
					scenario_sampled.buy_price,
					scenario_sampled.consumption_price,
					scenario_sampled.sell_price,
					scenario_sampled.penalty_price,
					scenario_sampled.Load,
					scenario_sampled.Ren,
					probability = point_probability[s]/n_scen_eps)
			end
		end
		return sampled_scenarios
	else # building scenarios for the risimulation, now scen_s doesn't represent the number of scenarios s to be generated
		
		# Array containing each scenario (only one scenario s)
		sampled_scenarios = Array{Scenario_Load_Renewable}(undef,n_scen_eps)
		for scen = 1:n_scen_eps
			sampled_scenarios[scen] = zero(Scenario_Load_Renewable) # initialize an empty scenario 
		end

		s = n_scen_s # scenario s actually considered
		sampler_eps = Scenario_eps_Sampler(market(data),users(data),point_s_load,point_s_ren,s)
			for eps = 1:n_scen_eps
				scen = eps
				scenario_sampled = sampler_eps()
				sampled_scenarios[scen] = Scenario_Load_Renewable(
					1,
					eps,
					scenario_sampled.peak_tariff,
					scenario_sampled.buy_price,
					scenario_sampled.consumption_price,
					scenario_sampled.sell_price,
					scenario_sampled.penalty_price,
					scenario_sampled.Load,
					scenario_sampled.Ren,
					probability = 1/n_scen_eps)
			end
		end
		return sampled_scenarios
end
