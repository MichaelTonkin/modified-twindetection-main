using findpoint

# Test function: Sphere function
f(x) = sum(x.^2)

# Problem bounds
bounds = [-5.12 5.12; -5.12 5.12; -5.12 5.12]

options = (testpoint,
allcirclecentre,
allcirclenormal,
allcircleradius,
badreflloc_normed,
0.1)

# Run PSO
best_position = particle_swarm_optimization(findpoint(options), bounds)
println("Best position found: ", best_position)
println("Best value found: ", f(best_position))

function particle_swarm_optimization(f, bounds, num_particles=30, max_iter=100)
    dim = length(bounds[:, 1])
    pos = [rand(dim) .* (bounds[:, 2] - bounds[:, 1]) .+ bounds[:, 1] for _ in 1:num_particles]
    vel = [rand(dim) for _ in 1:num_particles]
    pbest_pos = copy(pos)
    pbest_val = [f(pos[i]) for i in 1:num_particles]
    gbest_pos = pos[argmin(pbest_val)]
    
    w = 0.7
    c1 = 1.4
    c2 = 1.4
    
    for _ in 1:max_iter
        for i in 1:num_particles
            vel[i] = w .* vel[i] .+ c1 .* rand() .* (pbest_pos[i] - pos[i]) .+ c2 .* rand() .* (gbest_pos - pos[i])
            pos[i] = pos[i] .+ vel[i]

            # enforce bounds
            pos[i] = min.(max.(pos[i], bounds[:, 1]), bounds[:, 2])

            if f(pos[i]) < pbest_val[i]
                pbest_val[i] = f(pos[i])
                pbest_pos[i] = pos[i]
            end
        end
        
        gbest_pos = pbest_pos[argmin(pbest_val)]
    end

    return gbest_pos
end


#I think all we need to do is define the bounds as a unit sphere and then map the rotations to their given coordinates. 