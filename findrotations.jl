using Rotations

# Generate array of integer zeros for counting the number of circles that each point (in mrp space) is close to. 
# In the first round, this grid should go from -1 to 1 (the entire cube); 
# in later rounds it should go from -0.5 to 0.5 (a unit cube around each centre point)
# For the first round, this is achieved by treating the previous resolution as if it were 2.0
function make_mrparray(resolution,prevresolution=2.0)
    numstepsinhalf = Integer(cld(0.5*prevresolution,resolution))
    return fill(0,2*numstepsinhalf+1,2*numstepsinhalf+1,2*numstepsinhalf+1)
end

# Generate range based on resolutions that indicates the relative mrp values of points in the mrparray
function make_relrange(resolution,prevresolution=2.0)
    numstepsinhalf = Integer(cld(0.5*prevresolution,resolution))
    return -numstepsinhalf*resolution:resolution:numstepsinhalf*resolution
end

# Convert from location in MRP space to relative coordinates (without rounding) about a given centre point
function mrptocoord_unrounded(mrploc,resolution,prevresolution=2.0,centre=SVector{3,Float64}(0.0,0.0,0.0))
    return (mrploc-centre)/resolution .+ cld(0.5*prevresolution,resolution) .+ 1.0
end

# Convert from coordinates to an absolute location in MRP space based on relative range and centre
function coordtomrp(mrpcoord,relrange,centre=SVector{3,Float64}(0.0,0.0,0.0))
    SVector{3,Float64}(relrange[mrpcoord[1]],relrange[mrpcoord[2]],relrange[mrpcoord[3]]) + centre
end

# Calculate the appropriate tolerance in MRP space given a tolerance in distance space and the radius of the bad reflection
# NOTE: deltamax needs to be adjusted if some of the tolerance has been "used up" by the fact that the r values are different.
calc_mrptol(deltamax,r_refl) = deltamax/(2.0*r_refl)

# Calculate position in MRP space given circle details (theta,u,v,centre)
circleangletomrp(theta,u,v,centre) = centre + cos(theta)*v + sin(theta)*u

# Take unrounded coordinate for point in MRP space and convert to a block of coordinates within a given tolerance
#=  NOTE: The tolerance coordtol must already be in coordinate units not in MRP space units
    NOTE ALSO: THIS SEEMS TO BE THE CURRENT BOTTLENECK ON SPEED! 
    No immediate ideas on how to speed it up apart from maybe using eightwayround when coordtol < 1? (Introduced below) =#
function roundcoordtoblock(coord_unrounded,coordtol)
    lowerlimits_below = Int.(ceil.(coord_unrounded.-coordtol)) .- 1
    upperlimits = Int.(floor.(coord_unrounded.+coordtol))
    xdiff = upperlimits[1]-lowerlimits_below[1]
    ydiff = upperlimits[2]-lowerlimits_below[2]
    zdiff = upperlimits[3]-lowerlimits_below[3]
    numentries = xdiff*ydiff*zdiff
    allvecs = Vector{SVector{3,Int}}(undef,numentries)
    for xcount in 1:xdiff, ycount in 1:ydiff, zcount in 1:zdiff
        allvecs[zdiff*ydiff*(xcount-1) + zdiff*(ycount-1) + zcount] = 
        SVector{3,Int}(lowerlimits_below[1]+xcount, lowerlimits_below[2]+ycount, lowerlimits_below[3]+zcount)
    end
    return allvecs
end

# Take a given point and just round up and down
function eightwayround(unrounded)
    listofroundings = fill([0,0,0],8)
    listofroundings[1] = [floor(Int,unrounded[1]),floor(Int,unrounded[2]),floor(Int,unrounded[3])]
    listofroundings[2] = [floor(Int,unrounded[1]),floor(Int,unrounded[2]), ceil(Int,unrounded[3])]
    listofroundings[3] = [floor(Int,unrounded[1]), ceil(Int,unrounded[2]),floor(Int,unrounded[3])]
    listofroundings[4] = [floor(Int,unrounded[1]), ceil(Int,unrounded[2]), ceil(Int,unrounded[3])]
    listofroundings[5] = [ ceil(Int,unrounded[1]),floor(Int,unrounded[2]),floor(Int,unrounded[3])]
    listofroundings[6] = [ ceil(Int,unrounded[1]),floor(Int,unrounded[2]), ceil(Int,unrounded[3])]
    listofroundings[7] = [ ceil(Int,unrounded[1]), ceil(Int,unrounded[2]),floor(Int,unrounded[3])]
    listofroundings[8] = [ ceil(Int,unrounded[1]), ceil(Int,unrounded[2]), ceil(Int,unrounded[3])]
    return SVector{3,Int}.(listofroundings)
end

# Determine whether a given point is within a given disttol of a given circle (characterised by centre, normal, radius) in mrp space
function fullcircledistcheck(mrploc,circlecentre,circlenormal,circleradius,disttol)

    # Distance away from plane of circle (given that plane goes through origin)
    outofplane = dot(mrploc,circlenormal)

    # Vector within plane of circle
    inplane = mrploc - outofplane*circlenormal

    # Total distance from circle (including out of plane and in-plane parts), and test within disttol
    outofplane^2 + (norm(inplane-circlecentre)-circleradius)^2 ≤ disttol^2
end


# Calculate the distance of a given point from a given circle (characterised by centre, normal, radius) in mrp space
function fullcircledist(mrploc,circlecentre,circlenormal,circleradius)
    # Distance away from plane of circle (given that plane goes through origin)
    outofplane = dot(mrploc,circlenormal)

    # Vector within plane of circle
    inplane = mrploc - outofplane*circlenormal

    # Total distance from circle (including out of plane and in-plane parts), and test within disttol
    sqrt(dot(outofplane,outofplane) + (norm(inplane-circlecentre)-circleradius)^2)
end

# Increment the grid representing number of circle intersections in MRP space at a single point
function arrayincrement!(mrp3darray,vec)
    mrp3darray[vec[1],vec[2],vec[3]] += 1
    return nothing
end

# Increment the grid representing number of circle intersections in MRP space at a list of points
function manyarrayincrements!(mrp3darray,listofvecs)
    for i in listofvecs
        arrayincrement!(mrp3darray,i)
    end
    return nothing
end

function find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,disttol,centre)
    # Probably a way to rewrite this using findall?
    circlestokeep = Vector{Int64}()
    for i in 1:length(allcirclecentre)
        if fullcircledistcheck(centre,allcirclecentre[i],allcirclenormal[i],allcircleradius[i],disttol)
            append!(circlestokeep,i)
        end
    end
    return circlestokeep
end

function find_relgreatcircles(allcirclenormal,disttol,centre)
    # Probably a way to rewrite this using findall?
    circlestokeep = Vector{Int64}()
    for i in 1:length(allcirclenormal)
        if fullcircledistcheck(centre,SVector{3,Float64}(0.0,0.0,0.0),allcirclenormal[i],1.0,disttol)
            append!(circlestokeep,i)
        end
    end
    return circlestokeep
end


function coordtoint(coord,arraysize)
    return (coord[1]-1)*arraysize^2 + (coord[2]-1)*arraysize + coord[3]
end

function inttocoord(num,arraysize)
    return SVector{3,Int}(cld(num,arraysize^2),rem(cld(num,arraysize)-1,arraysize)+1,rem(num-1,arraysize)+1)
end

function find_rellines(badreflloc_normed,disttol,centre)
    # Probably a way to rewrite this using findall?
    linestokeep = Vector{Int64}()
    for i in 1:length(badreflloc_normed)
        if linedistcheck(centre,badreflloc_normed[i],disttol)
            append!(linestokeep,i)
        end
    end
    return linestokeep
end

function linedistcheck(mrploc,linevec,disttol)
    # Vector between point and line
    perptoline = mrploc - dot(mrploc,linevec)*linevec

    # Total distance from line, and test within disttol
    dot(perptoline,perptoline) ≤ disttol^2
end

function greatcirclebasis(circlenormal)
    if circlenormal[3] != 1.0
        uvec = SVector{3,Float64}(normalize([-circlenormal[2],circlenormal[1],0.0]))
        vvec = cross(circlenormal,uvec)
    else
        uvec = SVector{3,Float64}(1.0,0.0,0.0)
        vvec = SVector{3,Float64}(0.0,1.0,0.0)
    end
    return (uvec,vvec)
end

function manygreatcirclebasis(allcirclenormal)
    alluvec = Vector{SVector{3,Float64}}(undef,length(allcirclenormal))
    allvvec = Vector{SVector{3,Float64}}(undef,length(allcirclenormal))
    for i in 1:length(allcirclenormal)
        (alluvec[i],allvvec[i]) = greatcirclebasis(allcirclenormal[i])
    end
    return (alluvec,allvvec)
end