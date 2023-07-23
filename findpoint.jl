function findpoint(
    target_point,
    allcirclecentre,
    allcirclenormal,
    allcircleradius,
    badreflloc_normed,
    allr_dest,
    allr_badrefl,
    )
    centre = SVector{3,Float64}(0.0,0.0,0.0)
    circlestokeep = find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,1,centre)
    # Initialize the overlap count=#
    overlap_count = 0
    overlap_count += checkcirclerotations(target_point, allcirclecentre, allcirclenormal, allcircleradius, circlestokeep, allr_dest)

    overlap_count += checklinerotations(target_point, badreflloc_normed, allr_badrefl)
    overlap_count += checkgreatcirclerotations(target_point, badreflloc_normed, allr_badrefl)

    return overlap_count
end

function findpoint(target_point, options)

    allcirclecentre = options["allcirclecentre"]
    allcirclenormal = options["allcirclenormal"]
    allcircleradius = options["allcircleradius"]
    badreflloc_normed = options["badreflloc_normed"]
    allr_dest = options["allr_dest"]
    allr_badrefl = options["allr_badrefl"]

    centre = SVector{3,Float64}(0.0,0.0,0.0)
    circlestokeep = find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,1,centre)
    # Initialize the overlap count
    overlap_count = 0
    overlap_count += checkcirclerotations(target_point, allcirclecentre, allcirclenormal, allcircleradius, circlestokeep, allr_dest)

    overlap_count += checklinerotations(target_point, badreflloc_normed, allr_badrefl)
    overlap_count += checkgreatcirclerotations(target_point, badreflloc_normed, allr_badrefl)

    return overlap_count
end

function checkcirclerotations(target_point, allcirclecentre, allcirclenormal, allcircleradius,
    circlestokeep, allr_dest)
    # The number of circles to consider
    numcircles = length(circlestokeep)
    overlap_count = 0
    allr_dest_circles = allr_dest[cld.(circlestokeep,2)]
    # For each circle, we will check if the target_point is within the circle
    for i in 1:numcircles       
        
        # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
        mrptol = calc_mrptol(sqrt(deltamax^2-rdiff[cld(i,2)]^2),allr_dest_circles[i])
        coordtol = 0.5*sqrt(3.0) + mrptol
        
        # Check if the target_point is within mrptol of the circle
        if circledistcheck(target_point, allcirclecentre[i], allcirclenormal[i], allcircleradius[i],coordtol)
            # If it is, increment the overlap count
            overlap_count += 1
        end
    end
    return overlap_count
end

function checklinerotations(target_point, badreflloc_normed, allr_badrefl)
   
    #= Quick loop through all bad reflection points to see which are associated with lines to keep =#
    linestokeep = find_rellines(badreflloc_normed,1,centre)

    # Strip out only the lines to keep
    numlines = length(linestokeep)
    badrefl_lines = badreflloc_normed[linestokeep]
    allr_lines = allr_badrefl[linestokeep]

    overlap_count = 0

    # For each circle, we will check if the target_point is within the circle
    for i in 1:numlines   

        # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
        mrptol = calc_mrptol(deltamax,allr_lines[i])    
        coordtol = 0.5*sqrt(3.0) + mrptol
        # Check if the target_point is within mrptol of the circle
        if linedistcheck(target_point,badrefl_lines[i], coordtol)
            # If it is, increment the overlap count
            overlap_count += 1
        end
    end
    
    return overlap_count
end

function checkgreatcirclerotations(target_point, badreflloc_normed, allr_badrefl)

    #= Quick loop through all bad reflection points to see which are associated with great circles to keep =#
    greatcirclestokeep = find_relgreatcircles(badreflloc_normed, 1.0, centre)
    allr_greatcircle = allr_badrefl[greatcirclestokeep]

    # The number of circles to consider
    numgreatcircles = length(greatcirclestokeep)
    overlap_count = 0

    # Strip out only the great circless to keep
    allnormal_greatcircle = badreflloc_normed[greatcirclestokeep]

    # For each circle, we will check if the target_point is within the circle
    for i in 1:numgreatcircles       
        # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
        mrptol = calc_mrptol(deltamax,allr_greatcircle[i])
        coordtol = 0.5*sqrt(3.0) + mrptol
        # Check if the target_point is within mrptol of the circle
        if fullcircledistcheck(target_point, SVector{3,Float64}(0,0,0), allnormal_greatcircle[i], 1.0, coordtol)
            # If it is, increment the overlap count
            overlap_count += 1
        end
    end
    return overlap_count
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
