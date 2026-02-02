using Random
using Plots

include("MC_cancer_cell_model.jl")

############################################################
# 1. Simulation
############################################################

T_MUT = 120   # optional mutation time

function simulate(seed=1; snapshot_times=SNAPSHOTS)
    rng = MersenneTwister(seed)
    lat = init_lattice(N)
    oxy = init_oxygen(N)

    Nc_hist = Int[]
    Nn_hist = Int[]
    snaps = Dict{Int,Array{Int,2}}()

    for t in 1:TMAX
        step!(lat, oxy, rng)

        if t == T_MUT
            introduce_mutation!(lat, rng)
        end

        Nc, Nn = count_states(lat)
        push!(Nc_hist, Nc)
        push!(Nn_hist, Nn)

        if t in snapshot_times
            snaps[t] = copy(lat)
        end
    end

    return Nc_hist, Nn_hist, snaps
end

############################################################
# 2. Monte Carlo ensemble
############################################################

function run_ensemble(seeds)
    Nc_all = Matrix{Int}(undef, length(seeds), TMAX)
    survival = falses(length(seeds))

    for (k, s) in enumerate(seeds)
        Nc, _, _ = simulate(s)
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
        Nc, _, _ = simulate(s)
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

function plot_snapshot(lat; title_txt="", outpath=nothing)
    # heatmap expects a numeric matrix, lat is already Int
    # transpose for conventional x horizontal, y vertical look
    p = heatmap(
        permutedims(lat),
        aspect_ratio = :equal,
        axis = false,
        framestyle = :none,
        title = title_txt,
        color = cgrad([:white, :deepskyblue, :gray70, :red], 4, categorical=true),
        clims = (0, 3),
        colorbar = false
    )

    if outpath !== nothing
        savefig(p, outpath)
    end
    return p
end

function plot_two_snapshots(snaps::Dict{Int,Array{Int,2}}, t1::Int, t2::Int; outpath=nothing)
    @assert haskey(snaps, t1) "Snapshot time t1=$t1 not found in snaps."
    @assert haskey(snaps, t2) "Snapshot time t2=$t2 not found in snaps."

    p1 = plot_snapshot(snaps[t1]; title_txt="t = $t1")
    p2 = plot_snapshot(snaps[t2]; title_txt="t = $t2")

    p = plot(p1, p2, layout=(1,2), size=(900,450))

    if outpath !== nothing
        savefig(p, outpath)
    end
    return p
end
