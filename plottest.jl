using LinearAlgebra, StaticArrays, Rotations, PyPlot

include("findpointpairs.jl")

# Basis vectors for reciprocal crystal lattice
astar = [1.0,0,0]
bstar = [0,1.0,0]
cstar = [0,0,1.0]

recipbasis = [astar bstar cstar] 

# Maximum distance from origin to be considered
rmax = 5.0

# List of hkl values for bad reflections to be investigated
badreflectionhkl = Vector{Vector{Int64}}()

# Distance (in reciprocal Angstroms) between real reflection and twinned reflection that could lead to confusion (a bad reflection)
deltamax = 0.1

# Run the makeallhkl function 
(allloc,allr,allhkl,badreflectionindices) = make_allhkl(recipbasis,rmax,deltamax)

# Create lists of indices for destinations and sources
(destindices,sourceindices) = findpairindices(badreflectionindices,allr,deltamax)

# Normalised locations (useful for computing and testing rotations)
allloc_normed = normalizelocs(allloc) 

# List of differences in distance from the origin for paired points
rdiff = reflpair_rdiff(destindices,sourceindices,allr)

# Compute all of the vectors and parameters needed for constructing the circles
(allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius) = reflpair_circledetails(destindices,sourceindices,allloc_normed)


# PLOT PLAYING BELOW HERE

usefulrangeforplot = -1.0:0.01:1

# Single plot

i = Int(ceil(rand()*length(allu)))

thetaforplot = allphi[i]*usefulrangeforplot 

xyzforplot = fill(0.0,length(thetaforplot),3)
for j in 1:length(thetaforplot)
    xyzforplot[j,:] = allcirclecentre[i] + sin(thetaforplot[j])*allu[i] + cos(thetaforplot[j])*allv[i]
end
xyzforplot


f = plot3D(xyzforplot[:,1],xyzforplot[:,2],xyzforplot[:,3])
display(gcf())



# Bundle of plots
figure()
maxi = 96 #length(allu)
for i in 1:maxi
    thetaforplot = allphi[i]*usefulrangeforplot 
    
    xyzforplot = fill(0.0,length(thetaforplot),3)
    for j in 1:length(thetaforplot)
        xyzforplot[j,:] = allcirclecentre[i] + sin(thetaforplot[j])*allu[i] + cos(thetaforplot[j])*allv[i]
    end
    
    f = plot3D(xyzforplot[:,1],xyzforplot[:,2],xyzforplot[:,3])
end

display(gcf())