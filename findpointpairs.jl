using LinearAlgebra, StaticArrays

# Process HKL values to remove duplicates and produce a set only in the main hemisphere
function removedupehkl(hkl_list)

    # Move everything to main hemisphere
    hkltemp = copy(hkl_list)
    mainhemisphere!.(hkltemp)

    # Keep only items that appear only once
    sort!(hkltemp)
    return unique!(hkltemp)

end

# Process HKL values to remove duplicates and produce a set only in the main hemisphere
function removedupehkl!(hkl_list)

    # Move everything to main hemisphere
    mainhemisphere!.(hkl_list)

    # Keep only items that appear only once
    sort!(hkl_list)
    unique!(hkl_list)

end

# Generate the HKL values and locations for reflections that could be involved in twinning
# a, b, c are basis vectors in reciprocal space
# rmax is max distance from the origin in reciprocal space for reflections to keep
# deltamax is the max distance of a reflection from a rotated reflection to create an overlaps
# badreflectionhkl is a list of hkl values for bad reflections 
function make_allhkl(a,b,c,rmax,deltamax,badreflectionhkl=Vector{Vector{Int64}}())

    return make_allhkl([a b c],rmax,deltamax,badreflectionhkl)

end

# Generate a vector of all of the vectors [x,y,z] where x, y, z are each taken from given lists in all possible combinations.
function cartesianprod(xlist,ylist,zlist)
    # Below is based on https://stackoverflow.com/questions/56120583/n-dimensional-cartesian-product-of-a-set-in-julia
    return vec(collect.(collect(Iterators.product(xlist,ylist,zlist))))
end


# Generate the HKL values and locations for reflections that could be involved in twinning
# recipbasis is the matrix whose columns are the basis vectors in reciprocal space
# rmax is max distance from the origin in reciprocal space for reflections to keep
# deltamax is the max distance of a reflection from a rotated reflection to create an overlaps
# badreflectionhkl is a list of hkl values for bad reflections 
function make_allhkl(
    recipbasis,
    rmax,
    deltamax,
    badreflectionhkl=Vector{Vector{Int64}}()
    )
    
    # Move any input bad reflection indices into the main hemisphere and remove duplicates
    removedupehkl!(badreflectionhkl)

    # Reciprocal lattice basis vector that is smallest in norm
    smallestnorm = minimum(norm.(eachcol(recipbasis)))
    
    # Maximum number of lattice basis vectors to consider in each direction (2.5 is a safety factor)
    maxhkl = Int((2.5*rmax) ÷ smallestnorm)
    
    # Create list of integers from -maxhkl to maxhkl
    intrangelist = collect(-maxhkl:maxhkl)

    allhkl = cartesianprod(intrangelist,intrangelist,intrangelist)
    removedupehkl!(allhkl)
    numhkl = length(allhkl)

    # Initialise locations, distances, and flags to indicate points within rmax of origin in main hemisphere
    allloc = fill(fill(0.0,3),numhkl)
    allr = fill(0.0,numhkl)
    inside_rmax = fill(false,numhkl)

    # Loop through allhkl and compute locations (allloc), distances from the origin (allr) and whether 
    # distance from origin exceeds rmax + deltamax. I feel like there must be a better way of doing this.
    for i in 1:length(allhkl)
        allloc[i] = recipbasis*allhkl[i]
        allr[i] = norm(allloc[i])
        inside_rmax[i] = allr[i] ≤ rmax + deltamax # NOTE: including +deltamax here to include potential sources that are outside rmax
    end

    # Test to ensure that any specified bad reflections have been found and throw an error if not
    badreflectionindices = indexin(badreflectionhkl,allhkl)
    for i in badreflectionindices
        if i === nothing
            error("Bad reflection indices not found in generated hkl indices. Increase rmax?")
        elseif allr[i] > rmax # NOTE: not including +deltamax here because potential targets must be within rmax
            error("Bad reflection indices found, but outside reflections being kept. Increase rmax?")
        end
    end

    # Filter allloc and allhkl down to the points within the sphere
    allloc = allloc[inside_rmax]
    allr = allr[inside_rmax]
    allhkl = allhkl[inside_rmax]

    # Sort allloc and allhkl in order of distance from the origin, removing origin in the process
    distorder = sortperm(allr)[2:end]
    allr = allr[distorder]
    allloc = allloc[distorder]
    allhkl = allhkl[distorder]
    

    # Find indices of bad reflections or specify that all points within rmax of the origin in upper hemisphere are bad reflections
    if length(badreflectionhkl) > 0
        # THIS SEEMS REALLY SLOW! (Maybe find a way of speeding it up?)
        badreflectionindices = indexin(badreflectionhkl,allhkl)
    else
        badreflectionindices = findall(allr .≤ rmax)
        removehemisphere!(badreflectionindices,allhkl)
    end

    return (allloc,allr,allhkl,badreflectionindices)

end

# Remove points in the lower hemisphere from a list of hkl indices
function removehemisphere!(hkllist)
    filter!(x->checkhemisphere(x),hkllist)
end

# Remove points in the lower hemisphere from a list of indices refering to a set of hkls
function removehemisphere!(indexlist,allhkl)
    indicestokeep = checkhemisphere.(allhkl)
    indexlist = filter!(x->indicestokeep[x],indexlist)
end

# For removing duplicate points, check whether the point is in the upper hemisphere (true) or lower hemisphere (false)
# with points on the plane being dealt with appropriately (looking at y and then x values)
function checkhemisphere(hkl)
    if hkl[3] > 0
        return true
    elseif hkl[3] < 0
        return false
    elseif hkl[2] > 0
        return true
    elseif hkl[2] < 0
        return false
    elseif hkl[1] > 0
        return true
    else
        return false
    end
end

# Move to main hemisphere
function mainhemisphere!(hkl)
    if !checkhemisphere(hkl)
        hkl.*=(-1)
    end
end


# Generate the lists of point indices associated with possible twinning
function findpairindices(badreflectionindices,allr,deltamax)
    
    # Lazy version without preallocation, using append. Was struggling to get a good estimate of array size.
    # Could also have made a vector of vectors and then crushed them down somehow?
    destindices = Vector{Int64}()
    sourceindices = Vector{Int64}()


    for i in badreflectionindices

        # Find all indices between with r between allr[i]-deltamax and allr[i]+deltamax (not including i)
        sourcesbefore = searchsortedfirst(allr,allr[i]-deltamax):(i-1)
        sourcesafter = (i+1):searchsortedlast(allr,allr[i]+deltamax)
        
        # Append the new destination and source indices on to the current version.
        append!(destindices,fill(i,length(sourcesbefore)))
        append!(destindices,fill(i,length(sourcesafter)))
        append!(sourceindices,sourcesbefore)
        append!(sourceindices,sourcesafter)


    end

    return (destindices,sourceindices)

end

function reflpair_rdiff(destindices,sourceindices,allr)
    
    # Number of pairs
    numpairs = length(destindices)
    if length(sourceindices) != numpairs
        error("Lists of destination indices and source indices must be the same")
    end

    # Differences in distance from origin
    return allr[sourceindices] - allr[destindices]

end

# Normalise a three-dimensional vector and return a static vector (useful for later calculations?)
function normalizelocs(locsarray::Vector)
    return normalize.(SVector{3,Float64}.(locsarray))
end

# Normalise a list of three-dimensional vectors and return a list of static vectors (useful for later calculations?)
function normalizelocs(locsarray::Array)
    numlocs = size(locsarray,2)
    normedlocs = Vector{SVector{3,Float64}}(undef,numlocs)
    for i in 1:numlocs
        normedlocs[i] = normalize(SVector{3,Float64}(locsarray[:,i]))
    end
    return normedlocs
end

#= Small functions that are useful for generating circle information below (based on trigonometry 
    and the relationship between the radius of the circle with the location of its centre =#
costocsc(x) = 1.0/sqrt(1.0 - x^2)
costocot(x) = x/sqrt(1.0 - x^2)
radiusfromcentre(circlecentre) = sqrt(1.0 + norm(circlecentre)^2)

#= Generate the (scaled) u vector, v vector, phi value, and list of circle information given 
    a list of indices for destinations/bad reflections (destindices)
    a list of indices for original true reflections (sourceindices)
    normalised (important!) reflection locations
=#
function reflpair_circledetails(destindices,sourceindices,loc_normed)
    
    # Number of pairs
    numpairs = length(destindices)
    if length(sourceindices) != numpairs
        error("Lists of destination indices and source indices must be the same")
    end

    # u vectors, v vectors, phi values, circle radii, normal vectors
    allu = Vector{SVector{3,Float64}}(undef,2*numpairs)
    allv = Vector{SVector{3,Float64}}(undef,2*numpairs)
    allcosphi = Vector{Float64}(undef,2*numpairs)
    for i in 1:numpairs
        
        # Temporaries for source reflection and destination (bad) reflection
        desttemp = loc_normed[destindices[i]]
        sourcetemp = loc_normed[sourceindices[i]]
        
        # Unit vectors for u and v (note important to consider that source could be opposite)
        allu[2i-1] = normalize(desttemp + sourcetemp)
        allu[2i] = normalize(desttemp - sourcetemp)
        allv[2i-1] = normalize(cross(allu[2i-1],desttemp))
        allv[2i] = normalize(cross(allu[2i],desttemp))

        # Dot product to get cos of half-angle between source and destination
        allcosphi[2i-1] = dot(allu[2i-1],desttemp)
        allcosphi[2i] = dot(allu[2i],desttemp)
    end

    # Half of the angle between destination and source
    allphi = acos.(allcosphi)

    # Centre of circle (using conversion from cos(phi) to cot(phi))
    allcirclecentre = -costocot.(allcosphi).*allv

    # Unit normal to circle
    allcirclenormal = cross.(allu,allv)

    # Radius of circle
    allcircleradius = costocsc.(allcosphi)

    # Renormalise u and v with radius so they are basis vectors of circle 
    allu .*= allcircleradius
    allv .*= allcircleradius

    return(allu,allv,allphi,allcirclecentre,allcirclenormal,allcircleradius)
end