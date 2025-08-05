"""
test_bp_simple.jl

Simple test to isolate the BP message passing issue
"""

include("../dependencies.jl")
include("../functions/BP.jl")
include("../functions/Ising2D.jl")

using ITensors

function test_bp_minimal()
    """Test BP on minimal system."""
    
    println("🧪 Testing BP on minimal system...")
    
    # Very small system
    L = 3
    β = 0.4
    h = 0.0
    
    println("Creating $(L)x$(L) Ising system...")
    
    try
        # Create tensors
        T = Ising2D.get_ising_tn(L, β; h=h)
        println("✅ Created tensors")
        
        # Get adjacency
        adj_mat, edges, links = BP.get_adj_mat(T)
        println("✅ Got adjacency: $(length(edges)) edges")
        
        # Initialize messages
        messages = BP.get_messages(T, edges, links)
        println("✅ Initialized $(length(messages)) messages")
        
        # Check message structure
        println("Message structure:")
        for i in 1:min(5, length(messages))
            msg = messages[i]
            println("  Message $i: $(typeof(msg)), dims: $(dims(msg))")
        end
        
        # Try message passing with very few iterations
        println("Attempting message passing...")
        try
            messages_new = BP.message_passing(T, messages, edges, links, adj_mat; max_iters=1)
            println("✅ First BP iteration successful")
            
            # Try more iterations
            messages_new = BP.message_passing(T, messages, edges, links, adj_mat; max_iters=10)
            println("✅ 10 BP iterations successful")
            
        catch e
            println("❌ Error in message passing: $e")
            println("  Error type: $(typeof(e))")
            
            # Try to get more info about the error
            if isa(e, UndefRefError)
                println("  This is an UndefRefError - some reference is not initialized")
            end
            
            return false
        end
        
    catch e
        println("❌ Error in setup: $e")
        return false
    end
    
    return true
end

function test_bp_without_cluster_expansion()
    """Test just the BP part without cluster expansion."""
    
    println("\n🧪 Testing complete BP workflow without clusters...")
    
    L = 4
    β = 0.4
    h = 0.0
    
    try
        # Create system
        T = Ising2D.get_ising_tn(L, β; h=h)
        adj_mat, edges, links = BP.get_adj_mat(T)
        
        # BP messages
        messages = BP.get_messages(T, edges, links)
        messages = BP.message_passing(T, messages, edges, links, adj_mat; max_iters=100)
        
        # Partition function
        N = L^2
        Z_bp = BP.mean_free_partition_fn(1:N, T, messages, edges, links, adj_mat)
        println("✅ BP partition function: $Z_bp")
        
        if isnan(Z_bp) || isinf(Z_bp) || Z_bp == 0
            println("❌ Invalid partition function: $Z_bp")
            return false
        end
        
        # Free energy
        free_energy = -log(Z_bp) / (2 * N)
        println("✅ BP free energy: $free_energy")
        
        if isnan(free_energy) || isinf(free_energy)
            println("❌ Invalid free energy: $free_energy")
            return false
        end
        
        # Compare with Onsager
        exact_f = Ising2D.free_energy(β)
        println("✅ Exact free energy: $exact_f")
        println("✅ Error: $(abs(free_energy - exact_f))")
        
        return true
        
    catch e
        println("❌ Error in BP workflow: $e")
        return false
    end
end

# Run tests
if abspath(PROGRAM_FILE) == @__FILE__
    println("🚀 Testing BP in isolation...")
    
    success1 = test_bp_minimal()
    success2 = test_bp_without_cluster_expansion()
    
    if success1 && success2
        println("\n✅ BP tests passed - issue is likely in cluster expansion")
    else
        println("\n❌ BP itself has issues")
    end
end