using Plots
include("test.jl")

# Data from your test (explicit, no recomputation)
Omut = [0.237, 0.226, 0.226, 0.311]
fn   = [0.0352, 0.0390, 0.0335, 0.0268]
fc   = [0.0153, 0.0140, 0.0163, 0.0037]
labels = ["t_mut=5", "t_mut=10", "t_mut=40", "t_mut=90"]

p = plot(
    xlabel = "oxygen concentration at mutation site",
    ylabel = "mean per-step death fraction",
    legend = :topright,
    title  = "Differential hypoxia tolerance"
)

plot!(p, Omut, fn,
      seriestype = :scatter,
      marker = (:circle, 7),
      label = "normal cells")

plot!(p, Omut, fc,
      seriestype = :scatter,
      marker = (:diamond, 7),
      label = "cancer cells")

# Optional: connect points to guide the eye
plot!(p, Omut, fn, label = nothing, lw = 1, ls = :dash)
plot!(p, Omut, fc, label = nothing, lw = 1)

savefig(p, "hypoxia_tolerance_death_fraction.png")
display(p)


# old try_division function
function try_division!(lat, x, y, rng)
    cell = lat[x,y]
    (cell == NORMAL || cell == CANCER) || return false

    # 1) Try direct placement into an empty neighbor (8-neighborhood)
    neigh = shuffle(rng, copy(NEIGH))
    for (dx,dy) in neigh
        nx, ny = x + dx, y + dy
        if inbounds(lat, nx, ny) && lat[nx,ny] == EMPTY
            lat[nx,ny] = cell
            return true
        end
    end
    # No empty neighbor: stop here if pushing is disabled
     # ENABLE_PUSHING || return false

    # 2) No empty neighbor: only push with probability p_push[cell]
    if rand(rng) >= get(p_push, cell, 0.0)
        return false
    end

    # 3) Pushing attempt: try directions until an EMPTY is found along the ray
    dirs = shuffle(rng, copy(NEIGH))
    for (dx,dy) in dirs
        x1, y1 = x + dx, y + dy
        inbounds(lat, x1, y1) || continue

        chain = Tuple{Int,Int}[]
        cx, cy = x1, y1

        while inbounds(lat, cx, cy)
            if lat[cx,cy] == EMPTY
                # shift chain forward toward empty site
                for (px,py) in reverse(chain)
                    lat[px+dx, py+dy] = lat[px,py]
                end
                # place daughter cell next to parent
                lat[x1,y1] = cell
                return true
            end
            push!(chain, (cx,cy))
            cx += dx
            cy += dy
        end
    end

    return false
end

#old step function

function step!(lat, oxy, rng)
    cells = [(i,j) for i in 2:N-1, j in 2:N-1 if lat[i,j] in (NORMAL,CANCER)]
    shuffle!(rng, cells)

    for (i,j) in cells
        cell = lat[i,j]
        p_div, p_death, _ = fate_probs(cell, oxy[i,j])
        r = rand(rng)

        if r < p_death
            lat[i,j] = DEAD
        elseif r < p_death + p_div
            try_division!(lat, i, j, rng)
        end
    end
  clear_dead!(lat, rng)  # important: do it every step

end