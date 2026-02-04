using Plots
using Printf

include("simulate.jl")  # also includes the model
include("plots.jl")   # contains plot_snapshot and plot_two_snapshots


seed = 6
Δt = 80
tmut_list = [90, 120, 140]
cases = [(t, t + Δt) for t in tmut_list]

oxy = init_oxygen(N)

for (t_mut, t_late) in cases
    @assert t_late <= TMAX "t_late=$t_late exceeds TMAX=$TMAX"

    global T_MUT = t_mut

    Nc, Nn, snaps, mut_pos = simulate(seed; snapshot_times=[t_mut, t_late])

    mi, mj = mut_pos
    Omut = oxy[mi, mj]

    outname = @sprintf("seed%d_tmut%d_t%d.pdf", seed, t_mut, t_late)

    plot_two_snapshots(
        snaps, t_mut, t_late;
        mut_pos = mut_pos,
        Omut    = Omut,
        t_mut   = t_mut,
        outpath = outname
    )

    println(@sprintf(
        "Saved %s   mut_pos=(%d,%d)   O_mut=%.3f",
        outname, mi, mj, Omut
    ))
end
