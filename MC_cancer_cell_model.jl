############################################################
# Agent-Based Monte Carlo Tumor Growth Model (2D, Oxygen)
# Single-file educational implementation
############################################################

using Random
using Statistics

############################################################
# 1. Parameters
############################################################

const EMPTY   = 0 # empty site
const NORMAL  = 1 # normal cell
const DEAD    = 2 # dead cell
const CANCER  = 3 # cancer cell

N      = 301          # lattice size, you can adjust
TMAX   = 300        # time steps
L_MAX  = 6            # max pushing distance
SNAPSHOTS = [1, 50, 100, 200, 300]
p_clear_dead = 0.02   # per step probability that a DEAD site becomes EMPTY

# const FOUNDER = 4   # cancer founder marker (for plotting only) to locate the muation site
# const ENABLE_PUSHING = false   # set to false for hypoxia-only tests


###################################################################
##################################################################
# =========================================
# Quick push test counters (temporary)
# =========================================
global D_total = 0 # total divisions (space + push)
global D_space = 0 # divisions that succeeded by finding empty space
global D_push  = 0 # divisions that succeeded by pushing
global P_att   = 0 # total push attempts
global P_succ  = 0 # successful pushes

function reset_push_counters!()
    global D_total = 0 
    global D_space = 0 
    global D_push  = 0 
    global P_att   = 0 
    global P_succ  = 0  
    return nothing
end

function push_summary()
    I_push = D_total == 0 ? 0.0 : D_push / D_total   # fraction of divisions that relied on pushing
    E_push = P_att   == 0 ? 0.0 : P_succ / P_att  # pushing efficiency
    return (; D_total, D_space, D_push, P_att, P_succ, I_push, E_push)
end
######################################################################
#######################################################################


# Monte Carlo parameters (effective, qualitative)
params = Dict(
    :normal => Dict(:div => 0.30, :death => 0.05),
    :cancer => Dict(:div => 0.40, :death => 0.02)
)

# Pushing probabilities
p_push = Dict(
    NORMAL => 0.01,   # normal cells almost never push
    CANCER => 0.25    # cancer can push sometimes
)


############################################################
# 2. Initialization
############################################################

function init_lattice(N)
    lat = fill(EMPTY, N, N)
    c = div(N, 2) + 1
    lat[c, c] = NORMAL
    return lat
end

function init_oxygen(N; L=30.0)
    O = zeros(Float64, N, N)
    cx, cy = div(N, 2) + 1, div(N, 2) + 1
    # low in center, increases toward boundary
    Omin = 0.35 # baseline oxygen everywhere (tune 0.2–0.5)

    for i in 1:N, j in 1:N
        r = sqrt((i-cx)^2 + (j-cy)^2)
        O[i,j] = Omin + (1 - Omin) * (1.0 - exp(-r / L))
    end
    return O
end

############################################################
# 3. Oxygen-dependent fate probabilities
############################################################

function fate_probs(cell, oxy)
    if cell == NORMAL
        p_div   = params[:normal][:div] * oxy
        p_death = params[:normal][:death] * (1 - oxy)
    elseif cell == CANCER
        p_div   = params[:cancer][:div] * (0.5 + 0.5*oxy)
        p_death = params[:cancer][:death] * (1 - oxy)
    else
        return 0.0, 0.0, 1.0   # dead/empty should never divide
    end
    p_div = clamp(p_div, 0, 1)
    p_death = clamp(p_death, 0, 1)
    p_quies = max(0.0, 1 - p_div - p_death)
    return p_div, p_death, p_quies
end

############################################################
# 4. Neighborhood and pushing
############################################################

const NEIGH = [(dx,dy) for dx in -1:1, dy in -1:1 if !(dx==0 && dy==0)]

@inline function inbounds(lat, x, y)
    1 <= x <= size(lat,1) && 1 <= y <= size(lat,2)
end

# ----------------------------------------------------------
# Instrumented overload (does NOT change the division logic)
# It only records metrics for cancer cells.
# ----------------------------------------------------------
function try_division!(lat, x, y, rng)
    cell = lat[x,y]
    (cell == NORMAL || cell == CANCER) || return false

    # 1) Try direct placement into an empty neighbor (8-neighborhood)
    neigh = shuffle(rng, copy(NEIGH))
    for (dx,dy) in neigh
        nx, ny = x + dx, y + dy
        if inbounds(lat, nx, ny) && lat[nx,ny] == EMPTY
            lat[nx,ny] = cell

            # counter: space division
            if cell == CANCER
                global D_total += 1
                global D_space += 1
            end

            return true
        end
    end

    # 2) No empty neighbor: only push with probability p_push[cell]
    if rand(rng) >= get(p_push, cell, 0.0)
        return false
    end

    # counter: push attempt
    if cell == CANCER
        global P_att += 1
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

                # counters: push success + push-dependent division
                if cell == CANCER
                    global P_succ += 1
                    global D_total += 1
                    global D_push  += 1
                end

                return true
            end
            push!(chain, (cx,cy))
            cx += dx
            cy += dy
        end
    end

    return false
end

############################################################
# 5. One Monte Carlo step
############################################################

function clear_dead!(lat, rng; p_clear=p_clear_dead)
    for i in 2:size(lat,1)-1, j in 2:size(lat,2)-1
        if lat[i,j] == DEAD && rand(rng) < p_clear
            lat[i,j] = EMPTY
        end
    end
end

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


############################################################
# 6. Metrics
############################################################

function count_states(lat)
    Nc = count(==(CANCER), lat)
    Nn = count(==(NORMAL), lat)
    return Nc, Nn
end

############################################################
# 7. Mutation introduction
############################################################

function introduce_mutation!(lat, rng)
    candidates = Tuple{Int,Int}[]
    for i in 2:N-1, j in 2:N-1
        if lat[i,j] == NORMAL
            # boundary = has at least one empty neighbor
            for (dx,dy) in NEIGH
                if lat[i+dx, j+dy] == EMPTY
                    push!(candidates, (i,j))
                    break
                end
            end
        end
    end
    isempty(candidates) && return (false, 0, 0)
    (i,j) = rand(rng, candidates)
    lat[i,j] = CANCER
    return (true, i, j)
end