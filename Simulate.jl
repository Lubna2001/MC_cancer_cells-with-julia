using Random

include("MC_cancer_cell_model.jl")

############################################################
# 1. Simulation
############################################################

T_MUT = 90   # optional mutation time


function simulate(seed=1; snapshot_times=SNAPSHOTS)
    rng = MersenneTwister(seed)
    lat = init_lattice(N)
    oxy = init_oxygen(N)

    Nc_hist = Int[]
    Nn_hist = Int[]
    snaps = Dict{Int,Array{Int,2}}()

    mut_pos = nothing  # <-- local to this simulation run

    for t in 1:TMAX
        step!(lat, oxy, rng)

        if t == T_MUT && mut_pos === nothing
            ok, mi, mj = introduce_mutation!(lat, rng)
            mut_pos = ok ? (mi, mj) : nothing
        end

        Nc, Nn = count_states(lat)
        push!(Nc_hist, Nc)
        push!(Nn_hist, Nn)

        if t in snapshot_times
            snaps[t] = copy(lat)
        end
    end

    return Nc_hist, Nn_hist, snaps, mut_pos
end


############################################################
# 2. Monte Carlo ensemble
############################################################

function run_ensemble(seeds)
    Nc_all = Matrix{Int}(undef, length(seeds), TMAX)
    survival = falses(length(seeds))

    for (k, s) in enumerate(seeds)
        Nc, _, _, _ = simulate(s)
        Nc_all[k, :] = Nc
        survival[k] = (Nc[end] > 0)
    end

    Nc_mean = vec(mean(Nc_all, dims=1))
    Nc_std  = vec(std(Nc_all, dims=1))
    P_surv  = mean(survival)

    return Nc_mean, Nc_std, P_surv
end

############################################################
# 3. Derived quantities
############################################################

req_from_Nc(Nc) = sqrt.(Nc ./ pi)

############################################################
# 4. Survival statistics
############################################################

function survival_stats(seeds)
    surv = falses(length(seeds))
    for (k, s) in enumerate(seeds)
        Nc, _, _, _ = simulate(s)
        surv[k] = (Nc[end] > 0)
    end
    k = count(surv)
    n = length(surv)
    p = k / n
    return k, n, p
end

############################################################
# 5. Wilson confidence interval
############################################################

function wilson_ci(k, n; z=1.96)
    p̂ = k / n
    denom = 1 + z^2/n
    center = (p̂ + z^2/(2n)) / denom
    half = (z * sqrt(p̂*(1-p̂)/n + z^2/(4n^2))) / denom
    return center - half, center + half
end



using Plots

# Map lattice state to a small integer that can be colored
# EMPTY=0, NORMAL=1, DEAD=2, CANCER=3 in your model
# We will plot these as categories 0..3.