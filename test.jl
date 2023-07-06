#---
#NOTE TO SELF
# Run through Julia REPL in order to see plots
#----

using LinearAlgebra, StaticArrays, Rotations, BenchmarkTools, Profile, PProf

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

keptmrps = make_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,
    allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,centre,true)

# pprof()

keptmrps

# Eliminate symmetries
filter!(i -> i != [0.0,0.0,0.0],keptmrps)
filter!(i -> i != [0.0,1.0,0.0],keptmrps)
filter!(i -> i != [0.0,-1.0,0.0],keptmrps)

prevresolution = resolution
resolution = 0.25*resolution

newkeptmrps = Vector{SVector{3,Float64}}()
for i in keptmrps
    append!(newkeptmrps,make_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,i,false))
end

newkeptmrps

# Eliminate symmetries
# filter!(i -> i != [0.0,0.0,0.0],newkeptmrps)
# filter!(i -> i != [0.0,1.0,0.0],newkeptmrps)
# filter!(i -> i != [0.0,-1.0,0.0],newkeptmrps)


prevresolution = resolution
resolution = 0.2*resolution

newnewkeptmrps = Vector{SVector{3,Float64}}()
for i in newkeptmrps
    append!(newnewkeptmrps,make_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,i,false))
end

newnewkeptmrps

prevresolution = resolution
resolution = 0.1*resolution

newnewnewkeptmrps = Vector{SVector{3,Float64}}()
for i in newnewkeptmrps
    append!(newnewnewkeptmrps,make_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,i,false))
end

newnewnewkeptmrps

#=
prevresolution = resolution
resolution = 0.1*resolution

finalkeptmrps = Vector{SVector{3,Float64}}()
for i in newnewnewkeptmrps
    append!(finalkeptmrps,make_mrpkeeplist(resolution,keepcriterion,allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius,allr_dest,rdiff,badreflloc_normed,allr_badrefl,deltamax,prevresolution,i,false))
end
=#

finalkeptmrps = newnewnewkeptmrps #this is the variable we are looking at

println("The length of the array is: ", length(finalkeptmrps))
#=
xyzforplot = fill(0.0,length(finalkeptmrps),3)
for j in 1:length(finalkeptmrps)
    xyzforplot[j,:] = finalkeptmrps[j]
end

using PyPlot

figure()
scatter3D(xyzforplot[:,1],xyzforplot[:,2],xyzforplot[:,3])
display(gcf())
xyzforplot

figure()
scatter(xyzforplot[:,3],xyzforplot[:,1])
display(gcf())

# Find sources and destinations corresponding to a given MRP
# i = Int(ceil(rand()*length(finalkeptmrps)))
# testmrp = finalkeptmrps[i] 
# testmrp = SVector{3,Float64}(0.0,1.0,0.0)
testmrp = [1.0,0.0,0.0]
relcircles = find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,sqrt(3.0)/2.0*resolution,testmrp)
rellines = find_rellines(badreflloc_normed,sqrt(3.0)/2.0*resolution,testmrp)
relgreatcircles = find_relgreatcircles(badreflloc_normed,sqrt(3.0)/2.0*resolution,testmrp)
sourcelist = vcat(allhkl[sourceindices[cld.(relcircles,2)]],allhkl[badreflectionindices[rellines]],allhkl[badreflectionindices[relgreatcircles]])
destlist = vcat(allhkl[destindices[cld.(relcircles,2)]],allhkl[badreflectionindices[rellines]],allhkl[badreflectionindices[relgreatcircles]])
rdifflist = vcat(rdiff[cld.(relcircles,2)],fill(0.0,length(rellines)+length(relgreatcircles)))

destindices[cld.(relcircles,2)]

function showfulllist(list)
    show(IOContext(stdout, :limit => false), "text/plain", sourcelist)
    return nothing
end

showfulllist(sourcelist)
showfulllist(destlist)

sourcelistmodded = copy(sourcelist)

for i in sourcelistmodded
    i[2] = -i[2]
end

# NOTE: Not everything is just getting rotated from an item with the opposite j value, but most are.
anydifferences = sourcelistmodded - destlist
notthesame = findall(i -> i != [0, 0, 0], anydifferences)
sourcelist[notthesame]
destlist[notthesame]

showfulllist(sourcelist[notthesame])
showfulllist(destlist[notthesame])

####### TESTING THE PI ROTATION MORE GENERALLY
# What happens if I test the pi-rotation on all of the reflections?

badreflectionhkl = Vector{SVector{3,Int64}}()

(allloc,allr,allhkl,badreflectionindices) = make_allhkl(recipbasis,rmax,deltamax,badreflectionhkl)

allloc
badreflectionindices

# Create lists of indices for destinations and sources
(destindices,sourceindices) = findpairindices(badreflectionindices,allr,deltamax)

destindices

# Normalised locations (useful for computing and testing rotations)
allloc_normed = normalizelocs(allloc) 
badreflloc_normed = allloc_normed[badreflectionindices]

# List of differences in distance from the origin for paired points
rdiff = reflpair_rdiff(destindices,sourceindices,allr)

# Compute all of the vectors and parameters needed for constructing the circles
(allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius) = reflpair_circledetails(destindices,sourceindices,allloc_normed)

testmrp = SVector{3,Float64}(0.0,1.0,0.0)
relcircles = find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,sqrt(3.0)/2.0*resolution,testmrp)
rellines = find_rellines(badreflloc_normed,sqrt(3.0)/2.0*resolution,testmrp)
relgreatcircles = find_relgreatcircles(badreflloc_normed,sqrt(3.0)/2.0*resolution,testmrp)
sourcelist = vcat(allhkl[sourceindices[cld.(relcircles,2)]],allhkl[badreflectionindices[rellines]],allhkl[badreflectionindices[relgreatcircles]])
destlist = vcat(allhkl[destindices[cld.(relcircles,2)]],allhkl[badreflectionindices[rellines]],allhkl[badreflectionindices[relgreatcircles]])


=#
#= THINGS TO DO
* Set up more elegant system for iterative updating of MRPs
* Set up conversion of MRPs into axis/angle, rotation matrix, or other interpretation
* Set up system to identify number of twinning rotations for a given MRP (should be easy)
* Set up system to identify which pairs are involved for a given MRP (model system below)
* Set up File IO??
* DONE - Test with data from HKL file
=#