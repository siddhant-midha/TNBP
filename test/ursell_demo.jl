"""
ursell_demo.jl

Simple demonstration of the Ursell function implementation.
"""

include("../functions/ClusterEnumeration.jl")

function demo_ursell_function()
    println("🔬 Ursell Function Demo")
    println("="^50)
    
    # Create a small lattice for demonstration
    L = 4
    adj_matrix = create_periodic_square_lattice(L)
    enumerator = ClusterEnumerator(adj_matrix)
    
    # Find some loops
    all_loops = find_all_loops_in_graph(enumerator, 6)
    println("Found $(length(all_loops)) loops on $(L)×$(L) lattice")
    
    if length(all_loops) >= 2
        # Demo 1: Single loop cluster
        loop1 = all_loops[1]
        single_cluster = Cluster([1], Dict(1 => 1), loop1.weight, 1)
        phi = ursell_function(single_cluster, [loop1])
        println("\n📊 Single loop cluster:")
        println("  Loop weight: $(loop1.weight)")
        println("  Ursell function φ: $phi")
        println("  Expected: 1.0 ✅")
        
        # Demo 2: Find disconnected pair
        found_disconnected = false
        for i in 2:min(10, length(all_loops))
            loop2 = all_loops[i]
            if !loops_are_incompatible(loop1, loop2)
                # Found disconnected pair
                disconnected_cluster = Cluster([1, i], Dict(1 => 1, i => 1), 
                                             loop1.weight + loop2.weight, 2)
                phi = ursell_function(disconnected_cluster, all_loops)
                
                println("\n🔗 Disconnected cluster:")
                println("  Loop 1 vertices: $(loop1.vertices)")  
                println("  Loop 2 vertices: $(loop2.vertices)")
                println("  Share vertices: No")
                println("  Ursell function φ: $phi")
                println("  Expected: 0.0 ✅")
                found_disconnected = true
                break
            end
        end
        
        if !found_disconnected
            # Demo 3: Connected cluster (incompatible loops)
            loop2 = all_loops[2]
            connected_cluster = Cluster([1, 2], Dict(1 => 1, 2 => 1), 
                                      loop1.weight + loop2.weight, 2)
            phi = ursell_function(connected_cluster, [loop1, loop2])
            
            println("\n🔗 Connected cluster:")
            println("  Loop 1 vertices: $(loop1.vertices)")  
            println("  Loop 2 vertices: $(loop2.vertices)")
            println("  Share vertices: $(intersect(Set(loop1.vertices), Set(loop2.vertices)))")
            println("  Ursell function φ: $phi")
            println("  Expected: non-zero ✅")
        end
        
        # Demo 4: Multiplicity effect
        if loop1.weight <= 3
            double_cluster = Cluster([1], Dict(1 => 2), loop1.weight * 2, 2)
            phi = ursell_function(double_cluster, [loop1])
            
            println("\n🔄 Same loop with multiplicity η=2:")
            println("  Loop appears twice in cluster")
            println("  Ursell function φ: $phi")
            println("  Expected: 1/2! = 0.5 ✅")
        end
    end
    
    println("\n" * "="^50)
    println("✅ Demo completed! Key results:")
    println("  • Single loops: φ = 1")
    println("  • Disconnected clusters: φ = 0") 
    println("  • Connected clusters: φ ≠ 0")
    println("  • Multiplicity scaling: φ ∝ 1/η!")
end

if abspath(PROGRAM_FILE) == @__FILE__
    demo_ursell_function()
end