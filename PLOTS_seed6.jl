using Plots
using Printf

include("simulate.jl")   # model + simulate(...)
include("plots.jl")      # snapshot plotting functions (must include plot_two_snapshots)

# -------------------------------
# Settings
# -------------------------------
seed = 6
Δt = 80
tmut_list = [220,]
cases = [(t, t + Δt) for t in tmut_list]

# Oxygen field is static in time in this model
oxy = init_oxygen(N)

# -------------------------------
# Helper: count empty neighbors (optional, for printing only)
# -------------------------------
function empty_neighbors(lat, i, j)
    n = 0
    for (dx, dy) in NEIGH
        x, y = i + dx, j + dy
        if 1 <= x <= size(lat, 1) && 1 <= y <= size(lat, 2) && lat[x, y] == EMPTY
            n += 1
        end
    end
    return n
end

# -------------------------------
# Run all cases and save snapshot figures
# -------------------------------
for (t_mut, t_late) in cases
    @assert t_late <= TMAX "t_late=$t_late exceeds TMAX=$TMAX"

    global T_MUT = t_mut

    Nc, Nn, snaps, mut_pos = simulate(seed; snapshot_times=[t_mut, t_late])

    mi, mj = mut_pos
    Omut = oxy[mi, mj]

    # optional diagnostic printout
    Enb = empty_neighbors(snaps[t_mut], mi, mj)
    Nc_gain = Nc[t_late] - Nc[t_mut]
    Nn_gain = Nn[t_late] - Nn[t_mut]

    println(@sprintf(
        "t_mut=%3d  t_late=%3d  mut_pos=(%3d,%3d)  O_mut=%.3f  empty_neighbors=%d  Nc_gain=%d  Nn_gain=%d",
        t_mut, t_late, mi, mj, Omut, Enb, Nc_gain, Nn_gain
    ))

    outname = @sprintf("seed%d_tmut%d_t%d.pdf", seed, t_mut, t_late)

    # This should create a 1x2 panel plot:
    # left: snapshot at t_mut with mutation marked and O_mut shown
    # right: snapshot at t_late with same mutation marked
    plot_two_snapshots(
        snaps, t_mut, t_late;
        mut_pos = mut_pos,
        Omut = Omut,
        t_mut = t_mut,
        outpath = outname
    )

    println(@sprintf("Saved %s", outname))
end
