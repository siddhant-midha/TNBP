using ITensors, ITensorMPS, Plots, LaTeXStrings
using ProgressMeter,Graphs
push!(LOAD_PATH, "functions/")
using BP, Ising2D

L_arr = [4,8,16]  # Different lattice sizes
β_arr = collect(0.1:0.02:1.)
h = 0.0  # External magnetic field (assuming zero if not specified)

# Create the plot
plt = plot(xlabel=L"\beta", 
          ylabel="Absolute Error",
          title="BP Approximation Error vs Temperature for Different Lattice Sizes",
          grid=true)

# Colors and markers for different L values
colors = [:blue, :red, :green, :orange]
markers = [:circle, :square, :diamond, :triangle]

# Loop over different lattice sizes
for (i, L) in enumerate(L_arr)
    println("Computing for L = $L")
    
    # Initialize arrays for this L
    free_energy_arr = Float64[]
    exact_free_energy = Float64[]
    β_valid = Float64[]  # Store only valid β values
    
    # Progress meter to track computation
    @showprogress "Computing free energy for L=$L..." for β in β_arr
        tensors = Ising2D.get_ising_tn(L, β; h=h)
        BPadj_mat, BPedges, BPlinks = BP.get_adj_mat(tensors)
        messages = BP.get_messages(tensors, BPedges, BPlinks;random_part=.1) 
        messages = BP.message_passing(tensors, messages, BPedges, BPadj_mat; α=1.0, max_iters=500, diagnose=false, normalise=true)
        Z_l = BP.get_fixed_point_list(tensors, messages, BPadj_mat)
        
        # Use log-sum trick for numerical stability: log(prod(Z_l)) = sum(log.(Z_l))
        log_Z_sum = sum(log.(real.(Z_l)))    
        free_energy = - log_Z_sum / (2 * L^2)
        push!(free_energy_arr, free_energy)
        push!(exact_free_energy, Ising2D.free_energy(β))
        push!(β_valid, β)
    end
    
    # Calculate absolute error for this L
    absolute_error = abs.(free_energy_arr .- exact_free_energy)
    
    # Add to plot
    plot!(plt, β_valid, absolute_error,
          linewidth=2,
          marker=markers[i],
          markersize=4,
          label="L = $L",
          color=colors[i])
end

# Display the plot and keep it open
display(plt)

# Optionally save the plot
# savefig(plt, "ising_free_energy.png")

# Keep the plot window open (for interactive environments)
readline()
