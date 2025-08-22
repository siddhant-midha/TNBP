include("randompeps.jl")

w_list = [4,6,7,8,9,10]

L= 5 
η_arr = 0.5:0.5:2.
nsamples = 100
wmax = 10 
ti = true 
orthog = false 
annealing = 0.75 
maxiter = 1000 

# Use only the highest weight data which contains all loops up to that weight
cluster_data = load_latest_cluster_file(L,wmax)
loop_objects = cluster_data.all_loops  # This contains loops of all weights up to 10
all_loops = [loop_object.edges for loop_object in loop_objects]
# Storage for results across different η values
results_by_eta = Dict()

N = L  # PEPS dimensions
T = L

for η in η_arr
    println("Processing η = $η")
    
    # Storage for this η across all samples
    eta_contributions_by_weight = Dict{Int, Vector{Float64}}()
    for w in w_list
        eta_contributions_by_weight[w] = Float64[]
    end
    
    # Average over multiple samples for each η
    for sample in 1:nsamples
        if sample % 20 == 0
            println("  Sample $sample/$nsamples")
        end
        
        # Get random PEPS tensor network for this η
        tensors, _ = peps_controllable(N, T; η=η, ti=ti, orthog=orthog)
        adj_mat, edges, links = BP.get_adj_mat(tensors)
        messages = BP.get_messages(tensors, edges, links) 
        messages = BP.message_passing(tensors, messages, edges, adj_mat; α=annealing, max_iters=maxiter, diagnose=false, normalise=true)
        Z_list = BP.get_fixed_point_list(tensors, messages, adj_mat)
        tensors = BP.normalize_tensors(tensors, Z_list)
        
        # Calculate loop contributions and organize by actual loop weight (length)
        loop_contributions_by_weight = Dict{Int, Vector{Float64}}()

        for loop in all_loops
            if all([e in edges for e in loop])
                loop_weight = length(loop)  # Actual weight is the loop length
                loop_contri = scalar(BP.loop_contribution(loop, messages, tensors, edges, links, adj_mat))
                
                if !haskey(loop_contributions_by_weight, loop_weight)
                    loop_contributions_by_weight[loop_weight] = Float64[]
                end
                
                push!(loop_contributions_by_weight[loop_weight], abs(loop_contri))  # Take absolute value
            end
        end
        
        # Add sample averages to eta storage
        for w in w_list
            if haskey(loop_contributions_by_weight, w) && !isempty(loop_contributions_by_weight[w])
                push!(eta_contributions_by_weight[w], mean(loop_contributions_by_weight[w]))
            end
        end
    end
    
    # Store averaged results for this η
    results_by_eta[η] = Dict(w => mean(eta_contributions_by_weight[w]) for w in w_list if !isempty(eta_contributions_by_weight[w]))
end

# Create plot with all η values
even_weights = w_list
colors = palette(:viridis, length(η_arr))

p = plot(xlabel="Loop Weight", 
         ylabel="Average |Loop Contribution|",
         title="Loop Contribution Decay vs Weight for Different η\n(L=$L, Random PEPS, $nsamples samples)",
         grid=true,
         legend=:topright,
         size=(800, 600))

for (i, η) in enumerate(η_arr)
    avg_contribs = [results_by_eta[η][w] + 1e-50 for w in even_weights if haskey(results_by_eta[η], w)]
    valid_weights = [w for w in even_weights if haskey(results_by_eta[η], w)]
    
    if !isempty(avg_contribs)
        plot!(p, valid_weights, avg_contribs,
              linewidth=2,
              markershape=:circle,
              yscale=:ln,
              markersize=4,
              label="η = $η",
              color=colors[i])
    end
end

display(p)
readline()