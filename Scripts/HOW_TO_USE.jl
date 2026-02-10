############################################################
# HOW TO USE (interactive REPL workflow)
#
# Folder layout:
#   You must be in the repository root folder named:
#     MC_CANCER_MODEL
#
# Steps:
# 1) Open a terminal in the folder MC_CANCER_MODEL
# 2) Start the Julia REPL (Ctrl + shift + p) then type:
#       julia
# you should see a prompt like: julia>
# 3) First time only, install required packages:
#       import Pkg
#       Pkg.add("Plots")
#
# After that, you can 
# copy and paste sections into the Julia REPL interactively
############################################################

# Load simulation code
include("src/simulate.jl")

# Load plotting utilities
include("src/Plots.jl")

# to Load the plotting of seed 6 
include("src/PLOTS_seed6.jl") # you could test a different seed by modifying this file


using Plots

############################################################
# Simulation parameters
############################################################

seed = 1                     # random seed (any integer)
snapshot_times = [150, 200]  # time steps at which snapshots are saved

############################################################
# Run simulation
############################################################

Nc, Nn, snaps = simulate(seed; snapshot_times=snapshot_times)

println("Simulation completed.")
println("Seed = $seed")
println("Snapshots at t = $snapshot_times")

############################################################
# Plot results
############################################################

# Plot individual snapshots
plot(snaps[1], title="t = $(snapshot_times[1])")
plot(snaps[2], title="t = $(snapshot_times[2])")

# Plot two snapshots side by side and save
outname = "snapshots_$(snapshot_times[1])_$(snapshot_times[2])_seed$(seed).png"
plot_two_snapshots(
    snaps,
    snapshot_times[1],
    snapshot_times[2];
    outpath=outname
)

println("Saved figure: $outname")
