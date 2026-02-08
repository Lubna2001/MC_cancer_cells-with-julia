# first you need to start the REPL (ctrl + shift + p)

# then start julia in it by writing   julia 

# make sure to download the package Plots by writing  
 import Pkg; Pkg.add("Plots") 
# in the same folder where you have the julia files

# then load the file Simulation by writing  
 include("Simulate.jl") 
 #  you don't need to load MC_cancer_cell_model, since it is already loaded by loading Simulate.jl, 

# because include("MC_cancer_cell_model.jl") is written in Simulate.jl. 

# write 
using Plots

# write 
Nc, Nn, snaps = simulate(1; snapshot_times=[150, 200])
# where 1 is the number of the seed, you can choose any number for exmple 6 

#  or You can change the seed number  to any integer you want

# the 150, 200 are the times you want to have the snapshots of the lattice, you can check different times, 

# in the simulte.jl file the mutation time was introduced at 60, you can change it in the Simulate.jl file if you want to see how it affects the results
# then you can plot the results by writing 
plot(snaps[1])  # plot first snapshot
plot(snaps[2])  # plot second snapshot

# or plot two snaps next to each other by writing
plot_two_snapshots(snaps, 150, 200; outpath="snapshots_150_200_seed4.png") 

# and saving them with the name snapshots_150_200_seed4.png
# you can change the name of the output file as you want



# notes from latex 
Removing oxygen eliminates spatial heterogeneity in division and death probabilities, making growth depend primarily on local
space availability and stochastic effects. In the current parameter regime, oxygen acts as a weak modulator rather than
a dominant constraint.
