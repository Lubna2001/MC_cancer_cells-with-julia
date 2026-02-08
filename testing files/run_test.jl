# =========================================
# run_tests.jl
# load everything (in the correct order) and run seed tests
# =========================================


# 1) metrics first (defines PushStats type)
include("push_metrics.jl")
using .PushMetrics

# 2) model next (defines step!, try_division!, etc.)
include("MC_cancer_cell_model.jl")

# 3) simulation wrapper last (defines simulate_instrumented)
include("simulate_instrumented.jl")

# -----------------------------------------
# run a small test
# -----------------------------------------
seeds = [1, 4, 6, 7, 9, 10]

for s in seeds
    Nc, Nn, snaps, mut_pos, st = simulate_instrumented(s; snapshot_times=[T_MUT, 200])
    println("seed = $s   mut_pos = $mut_pos")
    println(metrics_summary(st))
    println()
end
