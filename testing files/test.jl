const NEIGH = [(dx,dy) for dx in -1:1, dy in -1:1 if !(dx==0 && dy==0)]


function count_empty_neighbors(lat, i, j)
           n = 0
           for (dx,dy) in NEIGH
               ni, nj = i + dx, j + dy
               if 1 ≤ ni ≤ size(lat,1) && 1 ≤ nj ≤ size(lat,2)
                   n += (lat[ni,nj] == EMPTY)
               end
           end
           return n
       end


seed = 4  # random seed for RNG
ts = T_MUT:(T_MUT+10)

Nc, Nn, snaps, mut_pos = simulate(seed; snapshot_times=collect(ts))

@show mut_pos
mut_pos === nothing && error("Mutation failed (mut_pos = nothing).")

mi, mj = mut_pos


for t in ts
    Nc_t, Nn_t = count_states(snaps[t])
    println("t=$t  Nc=$Nc_t  E1=$(count_empty_neighbors(snaps[t], mi, mj))")
end

Nc, Nn, snaps, mut_pos = simulate(seed; snapshot_times=collect(ts))

@show mut_pos

mut_pos === nothing && error("Mutation failed (mut_pos = nothing).")

mi, mj = mut_pos

for t in ts
    println("t=$t  E1=$(count_empty_neighbors(snaps[t], mi, mj))")
end
