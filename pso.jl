using LinearAlgebra, StaticArrays, Rotations, BenchmarkTools, Profile, PProf, Plots, CSV, DataFrames

include("findpointpairs.jl")

# Basis vectors for reciprocal crystal lattice

astar = [1.0,0,0]
bstar = [0,1.0,0]
cstar = [0,0,1.0]


astar = [3.37989705e-02, 0.00000000e+00,0.00000000e+00]
bstar = [-2.06959005e-18, 9.69161288e-02,0.00000000e+00]
cstar = [1.00372045e-03,-1.55853022e-17,1.14499668e-01]


recipbasis = [astar bstar cstar] 
# norm(23astar+2bstar+5cstar)

# Maximum distance from origin to be considered
# rmax = 3.0
rmax = 1.2

# List of hkl values for bad reflections to be investigated
badreflectionhkl = Vector{SVector{3,Int64}}()
badreflectionhkl = SVector{3,Int64}.([[  4,   1,   1],
[  6,   1,   0],
[ -8,   4,   1],
[  2,   2,   0],
[  9,   2,   5],
[ 11,   1,   5],
[  8,   2,   1],
[ 14,   2,   1],
[  3,   3,   0],
[ 13,   2,   5],
[  2,   0,   0],
[ 13,   5,   5],
[  6,   2,   0],
[  3,   0,   0],
[  7,   2,   0],
[ -1,   6,   6],
[  1,   3,   0],
[  7,   5,   6],
[ -7,   1,   1],
[ 22,   1,   5],
[ -2,   5,   6],
[  8,   1,   0],
[ -6,   1,   1],
[ -2,   3,   5],
[ -3,   1,   6],
[ 23,   2,   5],
[  6,   3,   0],
[ 17,   4,   5],
[ 18,   4,   0],
[ -4,   5,   2],
[ 10,   0,   0],
[ 10,   1,   0],
[  7,   8,   6],
[  2,   5,   6],
[ 12,   3,   5],
[ -5,   5,   5],
[  8,   9,   4],
[  5,   4,   0],
[ -1,   2,   6],
[ 10,   1,   5],
[-18,   0,   6],
[-11,   5,   1],
[  8,   5,   6],
[-17,   1,   6],
[ 11,   2,   1],
[ 11,   0,   0],
[ 17,   4,   0],
[ -1,   8,   4],
[ 10,   3,   6],
[  1,   6,   6]])


# Distance (in reciprocal Angstroms) between real reflection and twinned reflection that could lead to confusion (a bad reflection)
#= NOTE: Various things will probably need to be reworked if deltamax can depend on the destination 
(bad) reflection, or worse, on both source and destination =#
deltamax = 0.002 # Final refinement quite slow even with symmetries excluded. Some nontrivial solutions.
# deltamax = 0.000001 # shows only symmetries
# deltamax = 0.001 # makes the final refinement very slow if symmetries aren't excluded.

# Run the makeallhkl function 
# THIS IS REALLY INEFFICIENT WHEN GIVEN BADREFLECTIONHKL

(allloc,allr,allhkl,badreflectionindices) = make_allhkl(recipbasis,rmax,deltamax,badreflectionhkl)

allloc
badreflectionindices

# Create lists of indices for destinations and sources
(destindices,sourceindices) = findpairindices(badreflectionindices,allr,deltamax)

destindices
[sourceindices destindices]

# Normalised locations (useful for computing and testing rotations)
allloc_normed = normalizelocs(allloc) 
badreflloc_normed = allloc_normed[badreflectionindices]

# List of differences in distance from the origin for paired points
rdiff = reflpair_rdiff(destindices,sourceindices,allr)

# Compute all of the vectors and parameters needed for constructing the circles
(allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius) = reflpair_circledetails(destindices,sourceindices,allloc_normed)

allu
allloc_normed

# ROTATIONS MATERIAL BELOW HERE

include("findrotations.jl")
include("findpoint.jl")

# Resolution factor at each round
prevresolution = 2.0
resolution = 0.2

# Number of overlaps required to keep an mrp point as being relevant
keepcriterion = 20 
#keepcriterion = length(badreflectionindices)

allr_dest = allr[destindices]
allr_badrefl = allr[badreflectionindices]
centre = SVector{3,Float64}(0.0,0.0,0.0)
particle_trajectory = []


function particle_swarm_optimization(f, bounds, allcirclecentre,
    allcirclenormal,
    allcircleradius,
    badreflloc_normed,
    allr_dest, 
    allr_badrefl, num_particles=100, max_iter=100, threshold=20)

    options = Dict(
        "allcirclecentre" => allcirclecentre,
        "allcirclenormal" => allcirclenormal,
        "allcircleradius" => allcircleradius,
        "badreflloc_normed" => badreflloc_normed,
        "allr_dest" => allr_dest,
        "allr_badrefl" => allr_badrefl
    )

    dim = length(bounds[:, 1])
    center = zeros(dim)
    radius = 1
    pos = [rand(dim) .* (bounds[:, 2] - bounds[:, 1]) .+ bounds[:, 1] for _ in 1:num_particles]

    tabu_list = []
    proximity_threshold = 0.1 

    vel = [rand(dim) for _ in 1:num_particles]
    pbest_pos = copy(pos)
    pbest_val = [f(pos[i], options) for i in 1:num_particles]
    values_above_threshold = []
    w = 0.7 #inertia weight
    c1 = 1.4 #cognitive constant
    c2 = 0 #social constant (redundant)

    for j in 1:max_iter
        for i in 1:num_particles
            vel[i] = w .* vel[i] .+ c1 .* rand()# .* (pbest_pos[i] - pos[i]) .+ c2 .* rand() .* (gbest_pos - pos[i])
            pos[i] = pos[i] .+ vel[i]

            # enforce bounds
            distance_to_center = norm(pos[i] - center)
            if distance_to_center > radius
                pos[i] = center + (pos[i] - center) / distance_to_center * radius
            end

            # check if the function value exceeds the threshold
            current_val = f(pos[i], options)

            # check if the new position is close to any position in the tabu list
            for tabu_pos in tabu_list
                if norm(pos[i] - tabu_pos) <= proximity_threshold
                    current_val = 0
                    break
                end
            end

            if current_val > threshold
                push!(values_above_threshold, (pos[i], current_val))
                push!(tabu_list, pos[i])
            end

            if (current_val < pbest_val[i]) && !(in(current_val, values_above_threshold))#and it is not in the values above threshold list
                pbest_val[i] = current_val
                pbest_pos[i] = pos[i]
            end 

            push!(particle_trajectory, (copy(pos[i]), i))
            
        end
    end
    
    return values_above_threshold
end

function plot_pso(trajectory, num_particles)

    plotly()  # use plotly backend
    
    p = plot(legend = false)  # create a plot with no legend
    
    for i in 1:num_particles
        particle_trajectory = [pos for (pos, id) in trajectory if id == i]
        # Calculate the data range
        xs = [pos[1] for pos in particle_trajectory]
        ys = [pos[2] for pos in particle_trajectory]
        zs = [pos[3] for pos in particle_trajectory]

        plot!(p, xs, ys, zs, color = :blue, label = false)  # plot trajectories
        marker_size = 3  # size based on range of plot
        scatter!(p, xs, ys, zs, m = (marker_size, 0.5, :blue), label = false)
    end

    # create a grid of points representing the sphere surface
    theta = range(0, stop=2pi, length=100)
    phi = range(0, stop=pi, length=100)
    x = [cos(t)*sin(p) for t in theta, p in phi]
    y = [sin(t)*sin(p) for t in theta, p in phi]
    z = [cos(p) for t in theta, p in phi]

    # add the sphere to the plot
    surface!(p, x, y, z, alpha=0.3, color=:red, light=true)

    display(p)

end

function run_pso(num_particles, max_iter)

    bounds = [-1 1; -1 1; -1 1]

    result = particle_swarm_optimization(findpoint, bounds, allcirclecentre,
    allcirclenormal,
    allcircleradius,
    badreflloc_normed,
    allr_dest, 
    allr_badrefl,
    num_particles, max_iter)

    println("Num of values above threshold: ", length(result))
    return length(result)
end


function test_pso_particles()
    # Assuming run_pso is already defined in your code
    # function run_pso(x, y)
    #     # Your function implementation here
    # end

    # Initialize an empty DataFrame
    df = DataFrame(
        Particles = Int[], 
        Iterations = Int[],
        Run = Int[], 
        Result = Int[], 
        Time_Taken = Float64[], 
        Bytes_Allocated = Float64[]
    )

    # Initialize an empty DataFrame for the averages
    df_average = DataFrame(
        Particles = Int[], 
        Iterations = Int[],
        Result_Average = Float64[], 
        Time_Taken_Average = Float64[], 
        Bytes_Allocated_Average = Float64[]
    )

    pso_parameters = [15]  #, 75, 105, 135]
    iterations = [10, 20, 30, 40, 50]
    num_runs = 10

    # Run your function for each parameter value
    for param in pso_parameters
        for iter in iterations
            for i in 1:num_runs
                result, time_taken, bytes_allocated, gc_flag = @timed run_pso(param, iter)
                push!(df, Dict(
                    :Particles => param, 
                    :Iterations => iter,
                    :Run => i, 
                    :Result => result, 
                    :Time_Taken => time_taken, 
                    :Bytes_Allocated => bytes_allocated
                ))
            end

            # Compute the averages
            average_row = Dict(
                :Particles => param, 
                :Iterations => iter,
                :Result_Average => ceil(mean(filter(row -> row.Particles == param && row.Iterations == iter, df).Result)),
                :Time_Taken_Average => mean(filter(row -> row.Particles == param && row.Iterations == iter, df).Time_Taken),
                :Bytes_Allocated_Average => mean(filter(row -> row.Particles == param && row.Iterations == iter, df).Bytes_Allocated)
            )

            # Append the new row to df_average
            push!(df_average, average_row)
        end
    end

    # Save the detailed DataFrame to a CSV file
    CSV.write("results_detailed.csv", df)

    # Save the averages DataFrame to a CSV file
    CSV.write("results_average.csv", df_average)
end





    #println("Testing run_pso with 30, 1000")
    #result, time_taken, bytes_allocated, gc_flag = @timed run_pso(30, 1000)
    #println("Time taken: ", time_taken, " seconds")


test_pso_particles()
#println("Best value found: ", f(best_position))

#I think all we need to do is define the bounds as a unit sphere and then map the rotations to their given coordinates. 

#1. if current_val >= threshold, save coordinate to array called tabu_list
#2. in each step, check if we are close to any coordinates in the tabu_list.
#3. if we are close to any coordinates in the tabu list, set current_val to 0. 