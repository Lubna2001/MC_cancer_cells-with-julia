using Plots
include("Simulate.jl")


## plot individual  Monte C arlo runs
t = 1:TMAX 
p = plot(xlabel="time", ylabel="cancer cells",
         title="Individual Monte Carlo runs")

for s in seeds
    Nc, _, _ = simulate(s)
    plot!(p, t, Nc, label=false, alpha=0.7)
end

display(p) 

function req_from_Nc(Nc)
    sqrt.(Nc ./ pi) * 2
end  # approximate tumor radius from Nc 

#############################################################
# Plotting function for ensemble results
#############################################################
 function plot_req_ensemble(t, Req_mean, Req_std;
                           title_str="Equivalent radius (ensemble)",
                           outfile=nothing)
    p = plot(t, Req_mean,
             ribbon=Req_std,
             xlabel="time",
             ylabel="R_eq",
             label="mean ± std",
             title=title_str)
    if outfile !== nothing
        savefig(p, outfile)
    end
    return p
end





