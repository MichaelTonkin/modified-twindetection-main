function findpoint(
    target_point,
    allu,
    allcirclecentre,
    allcirclenormal,
    allcircleradius,
    allr_dest,
    badreflloc_normed,
    dist_from_point=0.001 #the distance from the point described where we check for overlaps.
    )
    centre = SVector{3,Float64}(0.0,0.0,0.0)
    circlestokeep = find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,1,centre)
    # Initialize the overlap count=#
    overlap_count = 0
    mrptol = dist_from_point
    overlap_count += checkcirclerotations(target_point, mrptol, allcirclecentre, allcirclenormal, allcircleradius, circlestokeep
    , dist_from_point)

    overlap_count += checklinerotations(target_point, dist_from_point, badreflloc_normed)

    return overlap_count
end

function checkcirclerotations(target_point, mrptol, allcirclecentre, allcirclenormal, allcircleradius,
    circlestokeep, dist_from_point)
    # The number of circles to consider
    numcircles = length(circlestokeep)
    overlap_count = 0
    # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
    mrptol = dist_from_point
    # For each circle, we will check if the target_point is within the circle
    for i in 1:numcircles       

        # Check if the target_point is within mrptol of the circle
        if circledistcheck(target_point, allcirclecentre[i], allcirclenormal[i], allcircleradius[i], mrptol)
            # If it is, increment the overlap count
            overlap_count += 1
        end
    end
    return overlap_count
end

function checklinerotations(target_point, dist_from_point, badreflloc_normed)
   
    #= Quick loop through all bad reflection points to see which are associated with lines to keep =#
    linestokeep = find_rellines(badreflloc_normed,1,centre)

    # Strip out only the lines to keep
    numlines = length(linestokeep)
    badrefl_lines = badreflloc_normed[linestokeep]

    overlap_count = 0
    # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
    mrptol = dist_from_point
    # For each circle, we will check if the target_point is within the circle
    for i in 1:numlines       

        # Check if the target_point is within mrptol of the circle
        if linedistcheck(target_point,badrefl_lines[i],mrptol)
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

function linedistcheck(mrploc,linevec,disttol)
    # Vector between point and line
    perptoline = mrploc - dot(mrploc,linevec)*linevec

    # Total distance from line, and test within disttol
    result = dot(perptoline,perptoline) ≤ disttol^2

    return result
end
