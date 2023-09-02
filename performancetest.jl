using LinearAlgebra, StaticArrays, Rotations, BenchmarkTools, Profile, PProf, Statistics

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
include("mrpkeeplist.jl")

# Resolution factor at each round
prevresolution = 2.0
resolution = 0.2

# Number of overlaps required to keep an mrp point as being relevant
keepcriterion = 20 
#keepcriterion = length(badreflectionindices)

allr_dest = allr[destindices]
allr_badrefl = allr[badreflectionindices]
centre = SVector{3,Float64}(0.0,0.0,0.0)

# @profile 

function benchmark_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,allr_dest,
    rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,centre,firstround)
    test_nums = 10
    println("starting tests")

    #println("Resolution: $resolution. Firstround: $firstround")
    benchmark_result = @benchmark make_mrpkeeplist($resolution, $keepcriterion, $allu, $allv, $allphi, $allcirclecentre, $allcirclenormal, $allcircleradius, 
    $allr_dest, $rdiff, $badreflloc_normed, $allr_badrefl, $deltamax, $prevresolution, $centre, $firstround) samples=test_nums
    save_bm_results(benchmark_result, resolution)
end
#save stats to a variable including metadata
#print when done

#final_text = ["Run, Metric, Time, \n"]
final_text = ["Resolution, Time, \n"]
function save_bm_results(benchmark_result, resolution)
    # print detailed results   
    run_num = length(final_text)
    # get the total time taken for all samples, in nanoseconds
    total_time_ns = benchmark_result.times

    # convert to seconds
    total_time_s = total_time_ns / 1e9
    result_text = "$(resolution), $(mean(total_time_s)),\n"
    #result_text *= "$(run_num), Minimum time, $(time(minimum(benchmark_result))),\n"
    #result_text *= "$(run_num), Mean time, $(time(mean(benchmark_result))),\n"
    #result_text *= "$(run_num), Maximum time, $(time(maximum(benchmark_result))),\n"
    
    push!(final_text, result_text)
        
end 


# Start with the initial resolution
resolution = 0.2

benchmark_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,
allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,2.0,centre,true)

# Update resolution to 0.05
prevresolution = resolution
resolution = 0.05

benchmark_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,
allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,centre,true)

# Update resolution to 0.01
prevresolution = resolution
resolution = 0.01

benchmark_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,
allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,centre,false)

# Update resolution to 0.001
prevresolution = resolution
resolution = 0.001

benchmark_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,
allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,centre,false)


print(join(final_text, ""))
