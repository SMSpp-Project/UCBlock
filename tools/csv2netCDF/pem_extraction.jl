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
using LinearAlgebra

# Central moment of order k of d, by adaptive quadrature, accurate to ~1e-10 at
# every order (the PointEstimateMethod package falls back to Monte Carlo for
# distributions without a direct `moment` method, which is pure noise at the
# high orders PEM needs).
function central_moment_quad(d::UnivariateDistribution, k::Int)
    μ = mean(d)
    lo, hi = extrema(d)                 # truncated normal -> (0.0, Inf)
    val, _ = quadgk(x -> (x - μ)^k * pdf(d, x), lo, hi; rtol = 1e-10)
    return val
end

# N-point PEM (Gaussian quadrature) for d, computed solver-free via Golub-Welsch
# on the exact central moments. The package builds the points by solving an
# over-determined feasibility LP with HiGHS; on some machines that LP is reported
# infeasible (-> 0 solutions, a thrown error) or returns a degenerate all-zero
# point (-> weights summing to 0), both for larger N. Golub-Welsch is the
# numerically stable construction: nodes are the eigenvalues of the Jacobi matrix
# obtained from the Cholesky factor of the Hankel moment matrix, weights are the
# squared first eigenvector components. Returns the same (x, p) shape as `pem`.
function gauss_pem(d::UnivariateDistribution, N::Int)
    μ = mean(d)
    m = [k == 0 ? 1.0 : central_moment_quad(d, k) for k in 0:2N]
    M = [m[i + j + 1] for i in 0:N, j in 0:N]          # Hankel, (N+1)x(N+1)
    R = cholesky(Symmetric(M)).U                       # M = R'R, upper triangular
    α = [R[k, k+1] / R[k, k] - (k > 1 ? R[k-1, k] / R[k-1, k-1] : 0.0) for k in 1:N]
    β = [R[k+1, k+1] / R[k, k] for k in 1:N-1]
    E = eigen(SymTridiagonal(α, β))
    ord = sortperm(E.values)
    return (x = E.values[ord] .+ μ, p = (vec(E.vectors[1, :]) .^ 2)[ord])
end

# Probabilities must be (numerically) a distribution: non-negative and summing
# to 1. The package returns garbage on a failed inner solve without throwing.
_valid_pem(r) = all(>=(-1e-9), r.p) && abs(sum(r.p) - 1) < 1e-6

# Run PEM the package's default way (Monte Carlo moments, kept so the draws stay
# bit-identical with test_instance_with_EC_jl.jl on the cases that work) and fall
# back to the deterministic Golub-Welsch construction whenever the package throws
# or returns invalid weights. The cases needing the fallback have no EC.jl
# reference, so the divergence is harmless there.
function pem_robust(d::UnivariateDistribution, N::Int)
    res = try
        r = pem(d, N)
        _valid_pem(r) ? r : nothing
    catch
        nothing
    end
    isnothing(res) || return res
    @warn "pem_extraction: package PEM unusable for N=$N on $d; " *
          "using Golub-Welsch quadrature (no longer bit-identical with EC.jl)."
    return gauss_pem(d, N)
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
