using Plots
using Printf
gr()

include("simulate.jl")  # this includes MC_cancer_cell_model.jl already

# ---------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------

function plot_snapshot(lat; title_txt="", outpath=nothing, mut_pos=nothing, annotation_txt=nothing)
    p = heatmap(
        permutedims(lat),
        aspect_ratio = :equal,
        axis = false,
        framestyle = :none,
        title = title_txt,
        color = cgrad([:white, :deepskyblue, :gray70, :red], 4, categorical=true),
        clims = (0, 3),
        colorbar = false,              # keep as in your older plots
        legend = :topright,
        legendfontsize = 6
    )

    # Overlay mutation marker if provided
    if mut_pos !== nothing
        mi, mj = mut_pos
        scatter!(p, [mi], [mj],
            markersize = 6,
            markercolor = :yellow,
            markerstrokecolor = :black,
            label = "mutation"
        )
    end

    # Legend entries (dummy points)
    scatter!(p, [NaN], [NaN], markercolor=:deepskyblue, markersize=5, label="normal")
    scatter!(p, [NaN], [NaN], markercolor=:red,        markersize=5, label="cancer")
    scatter!(p, [NaN], [NaN], markercolor=:gray70,     markersize=5, label="dead")
    scatter!(p, [NaN], [NaN], markercolor=:white,      markersize=5, label="empty")

    # Oxygen annotation (small text in the corner)
    if annotation_txt !== nothing
        annotate!(p, 5, size(lat,1) - 20,  text(annotation_txt, 9, :black, :left))
    end

    if outpath !== nothing
        savefig(p, outpath)
    end
    return p
end


function plot_two_snapshots(snaps::Dict{Int,Array{Int,2}}, t1::Int, t2::Int;
                            outpath=nothing, mut_pos=nothing, Omut=nothing, t_mut=nothing)

    @assert haskey(snaps, t1) "Snapshot time t1=$t1 not found in snaps."
    @assert haskey(snaps, t2) "Snapshot time t2=$t2 not found in snaps."

    ann_left = nothing
    if Omut !== nothing && t_mut !== nothing && t1 == t_mut
        ann_left = "O_mut = $(round(Omut, digits=3))"
    end

    p1 = plot_snapshot(snaps[t1]; title_txt="t = $t1", mut_pos=mut_pos, annotation_txt=ann_left)
    p2 = plot_snapshot(snaps[t2]; title_txt="t = $t2", mut_pos=mut_pos)

    p = plot(p1, p2, layout=(1,2), size=(900,450))

    if outpath !== nothing
        savefig(p, outpath)   # THIS saves the combined figure
    end
    return p
end

# ---------------------------------------------------------
# Final figure generation
# ---------------------------------------------------------

# Mutation time and evaluation time
T_MUT = 120
t_mut  = T_MUT
t_late = 200

# Make sure we save both snapshots
snap_times = [t_mut, t_late]

# Oxygen field is static, compute once
oxy = init_oxygen(N)

seeds = [1, 4, 6]

for s in seeds
    Nc, Nn, snaps, mut_pos = simulate(s; snapshot_times=snap_times)

    mi, mj = mut_pos
    Omut = oxy[mi, mj]

    outname = @sprintf("seed%d_tmut%d_t%d.pdf", s, t_mut, t_late)

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

