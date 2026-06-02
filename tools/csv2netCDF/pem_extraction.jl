# Sample the load / PV / wind distributions for long-period uncertainty using
# the Point Estimate Method. Schema aligned with EC.jl@stochastic.
#
# scen_s_sample: number of long-period (s) scenarios.
# sigma_load:    std deviation of the load distribution (mean = 1.0).
# mean_pv, sigma_pv:   parameters of the PV distribution.
# mean_wind, sigma_wind: parameters of the wind distribution.
# unc_var:       string of letters in {"L","P","W"} selecting which sources
#                of long-period uncertainty are sampled. For each unselected
#                letter the corresponding `point_*` vector is filled with the
#                deterministic mean.
#
# Returns a 4-tuple (point_load, point_pv, point_wind, scen_probability).

using QuadGK

# Central moment of order k of d, by adaptive quadrature. PEM needs central
# moments up to order 2N; estimating those by Monte Carlo (the package default
# for distributions without a direct `moment` method) makes the order >= ~2N-1
# terms pure sampling noise, so the over-determined probability step turns
# infeasible for some N/seed/environment combinations. Quadrature is accurate to
# ~1e-10 for every order, making the PEM construction deterministic and robust.
# Having this method on the distribution type sends `pem` down its direct branch.
function central_moment_quad(d::UnivariateDistribution, k::Int)
    μ = mean(d)
    lo, hi = extrema(d)                 # truncated normal -> (0.0, Inf)
    val, _ = quadgk(x -> (x - μ)^k * pdf(d, x), lo, hi; rtol = 1e-10)
    return val
end

# Run PEM the package's default way (Monte Carlo moments, kept so the draws stay
# bit-identical with test_instance_with_EC_jl.jl) and fall back to deterministic
# quadrature moments only when the MC path makes the probability step infeasible.
# The cases that need the fallback (typically larger N) have no EC.jl reference,
# so the divergence is harmless there.
function pem_robust(d::UnivariateDistribution, N::Int)
    try
        return pem(d, N)
    catch
        @warn "pem_extraction: Monte Carlo PEM failed for N=$N on $d; " *
              "retrying with quadrature moments (no longer bit-identical with EC.jl)."
        return pem(d, N; central_moment_fun = central_moment_quad)
    end
end

function pem_extraction(scen_s_sample::Int,
                        sigma_load,
                        mean_pv, sigma_pv,
                        mean_wind, sigma_wind,
                        unc_var::AbstractString)

    sample_L = occursin("L", unc_var)
    sample_P = occursin("P", unc_var)
    sample_W = occursin("W", unc_var)

    if scen_s_sample > 1 && (sample_L || sample_P || sample_W)
        # Build the truncated normal for each sampled letter and run PEM.
        # Distributions that are not sampled get their deterministic mean
        # replicated `scen_s_sample` times.
        local pem_load = nothing
        local pem_pv = nothing
        local pem_wind = nothing

        if sample_L
            d_load = truncated(Normal(1.0, sigma_load), 0.0, +Inf)
            pem_load = pem_robust(d_load, scen_s_sample)
        end
        if sample_P
            d_pv = truncated(Normal(mean_pv, sigma_pv), 0.0, +Inf)
            pem_pv = pem_robust(d_pv, scen_s_sample)
        end
        if sample_W
            d_wind = truncated(Normal(mean_wind, sigma_wind), 0.0, +Inf)
            pem_wind = pem_robust(d_wind, scen_s_sample)
        end

        # The scenario probability is taken from the first sampled letter
        # (in priority order L > P > W). This matches the EC.jl@stochastic
        # convention: probabilities of the unsampled distributions are
        # implicitly 1, so they don't change the joint mass.
        if !isnothing(pem_load)
            scen_probability = pem_load.p
        elseif !isnothing(pem_pv)
            scen_probability = pem_pv.p
        else
            scen_probability = pem_wind.p
        end

        point_load = isnothing(pem_load) ? fill(1.0, scen_s_sample) : pem_load.x
        point_pv   = isnothing(pem_pv)   ? fill(mean_pv,   scen_s_sample) : pem_pv.x
        point_wind = isnothing(pem_wind) ? fill(mean_wind, scen_s_sample) : pem_wind.x
    else
        point_load = [1.0]
        point_pv   = [mean_pv]
        point_wind = [mean_wind]
        scen_probability = [1.0]
    end

    return (point_load, point_pv, point_wind, scen_probability)
end
