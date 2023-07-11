function findpoint(
    target_point,
    allu,
    allcirclecentre,
    allcirclenormal,
    allcircleradius,
    allr_dest,
    dist_from_point=0.001 #the distance from the point described where we check for overlaps.
    )
    centre = SVector{3,Float64}(0.0,0.0,0.0)
    circlestokeep = find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,1,centre)
    # Initialize the overlap count=#
    overlap_count = 0
    mrptol = dist_from_point
    overlap_count += checkcirclerotations(target_point, mrptol, allu, allcirclecentre, allcirclenormal, allcircleradius, circlestokeep)
    return overlap_count
end

function checkcirclerotations(target_point, mrptol, allu, allcirclecentre, allcirclenormal, allcircleradius, circlestokeep)
    # The number of circles to consider
    numcircles = length(circlestokeep)
    overlap_count = 0
    # For each circle, we will check if the target_point is within the circle
    for i in 1:numcircles
        
        # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
        mrptol = 0.1
        # Check if the target_point is within mrptol of the circle
        if circledistcheck(target_point, allcirclecentre[i], allcirclenormal[i], allcircleradius[i], mrptol)
            # If it is, increment the overlap count
            overlap_count += 1
        end
    end
    return overlap_count
end

# Determine whether a given point is within a given disttol of a given circle (characterised by centre, normal, radius) in mrp space
function circledistcheck(mrploc,circlecentre,circlenormal,circleradius,disttol)

    # Distance away from plane of circle (given that plane goes through origin)
    outofplane = dot(mrploc,circlenormal)

    # Vector within plane of circle
    inplane = mrploc - outofplane*circlenormal

    result = outofplane^2 + (norm(inplane-circlecentre)-circleradius)^2 ≤ disttol^2

    # Total distance from circle (including out of plane and in-plane parts), and test within disttol
    return result
end

#=    function find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,disttol,centre)
        # Probably a way to rewrite this using findall?
        circlestokeep = Vector{Int64}()
        for i in 1:length(allcirclecentre)
            if fullcircledistcheck(centre,allcirclecentre[i],allcirclenormal[i],allcircleradius[i],disttol)
                append!(circlestokeep,i)
            end
        end
        return circlestokeep
    end
    # Return the overlap count for the target point
    return overlap_count
end=#