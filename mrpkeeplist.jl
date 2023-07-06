
# Find list of mrps with many overlapping circles
function make_mrpkeeplist(
    resolution,
    keepcriterion,
    allu,
    allv,
    allphi,
    allcirclecentre,
    allcirclenormal,
    allcircleradius,
    allr_dest,
    rdiff,
    badreflloc_normed,
    allr_badrefl,
    deltamax,
    prevresolution=2.0,
    centre=SVector{3,Float64}(0.0,0.0,0.0),
    firstround=true
    )

    #= SOMETHING TO FIX!!! - Using the prev coordtol rather than 
    prev resolution in determining which circles etc. are relevant...???=#

    # Generate integer array that will store the number of overlaps near each point in mrp space
    mrp3darray = make_mrparray(resolution,prevresolution)

    # Generate range of relative mrp values (useful for converting coord in mrp array to mrp value)
    relrange = make_relrange(resolution,prevresolution)

    ######## CIRCLES FROM PAIRS OF POINTS

    #= Quick loop through all overlaps to check which ones are relevant for checking (based on whether centre is within
    an appropriate amount of prevresolution of the circle) =#
    circlestokeep = find_relcircles(allcirclecentre,allcirclenormal,allcircleradius,sqrt(3.0)/2.0*prevresolution,centre)

    # Strip out only the circles to keep 
    numcircles = length(circlestokeep)
    allu = allu[circlestokeep]
    allv = allv[circlestokeep]
    allphi = allphi[circlestokeep]
    allcirclecentre = allcirclecentre[circlestokeep]
    allcirclenormal = allcirclenormal[circlestokeep]
    allcircleradius = allcircleradius[circlestokeep]
    allr_dest_circles = allr_dest[cld.(circlestokeep,2)]

    #TODO: these don't actually appear to be saving the variable we want to anything...
    mrp3darray = plot_circles(resolution, numcircles, deltamax, rdiff, allr_dest_circles, firstround, allphi, allcircleradius, centre, allu, allv, allcirclecentre, allcirclenormal, prevresolution, mrp3darray, relrange)
    plot_lines(badreflloc_normed, prevresolution, centre, allr_badrefl, deltamax, resolution, firstround, relrange, mrp3darray)
    process_all_great_circles(badreflloc_normed, allr_badrefl, deltamax, resolution, firstround, prevresolution, centre, relrange, mrp3darray)
    final_mrps = coordtomrp.(findall(x->x≥keepcriterion,mrp3darray),(make_relrange(resolution,prevresolution),),(centre,))
    
    # Find MRP values for pixels where value recorded in array is larger than keep criterion
    return final_mrps
end

function plot_circles(resolution, numcircles, deltamax, rdiff, allr_dest_circles, firstround, allphi, allcircleradius, centre, allu, allv, allcirclecentre, allcirclenormal, prevresolution, mrp3darray, relrange)
    
    currcirclerefs = Vector{Int64}()
    
    for i in 1:numcircles
        coordtol = get_mrp_tol_circ(resolution, deltamax, rdiff, allr_dest_circles, i)
        if firstround
            thetarange = gen_theta_range_using_tol(coordtol, resolution, allcircleradius, allphi, i)
        else
            thetarange = gen_theta_range_using_mid(coordtol, resolution, allcircleradius, allphi, centre, allu, allv, allcirclecentre, allcirclenormal, prevresolution, i)
        end
        currcirclerefs = filter_circle_points(thetarange, allu, allv, allcirclecentre, mrp3darray, resolution, prevresolution, centre, coordtol, allcirclenormal, allcircleradius, relrange, currcirclerefs, i)    
    end

return currcirclerefs
end

function get_mrp_tol_circ(resolution, deltamax, rdiff, allr_dest_circles, i)
    mrptol = calc_mrptol(sqrt(deltamax^2-rdiff[cld(i,2)]^2),allr_dest_circles[i])
    coordtol = 0.5*sqrt(3.0) + mrptol/resolution
    return coordtol
end

function gen_theta_range_using_tol(coordtol, resolution, allcircleradius, allphi, i)
    thetarange = make_relrange(coordtol*resolution/allcircleradius[i],2*allphi[i])
    return thetarange
end

function gen_theta_range_using_mid(coordtol, resolution, allcircleradius, allphi, centre, allu, allv, allcirclecentre, allcirclenormal, prevresolution, i)
    centrevec_incircleplane = centre - dot(centre,allcirclenormal[i])*allcirclenormal[i] - allcirclecentre[i]
    theta_mid = atan(dot(centrevec_incircleplane,allu[i]),dot(centrevec_incircleplane,allv[i]))
    thetarange = theta_mid .+ make_relrange(coordtol*resolution/allcircleradius[i],sqrt(3.0)*prevresolution)
    return thetarange
end


function filter_circle_points(thetarange, allu, allv, allcirclecentre, mrp3darray, resolution, prevresolution, centre, coordtol, allcirclenormal, allcircleradius, relrange, currcirclerefs, i)
    # Loop through range of theta values
    for theta in thetarange
        # Convert theta to MRP space location
        mrpfromtheta = circleangletomrp(theta,allu[i],allv[i],allcirclecentre[i])
        possiblecoords = construct_mrp_coords(mrpfromtheta, resolution, prevresolution, centre, coordtol)
    
        # Filter out coordinates outside the range of the grid
        filter_coords!(possiblecoords, mrp3darray)
    
        # MRP locations
        possiblemrps = coordtomrp.(possiblecoords,(relrange,),(centre,))
    
        possiblecoords, possiblemrps = filter_outside_sphere(possiblemrps, possiblecoords, resolution)
        possiblecoords = filter_far_from_circle(possiblemrps, possiblecoords, allcirclecentre[i], allcirclenormal[i], allcircleradius[i], coordtol, resolution)
            
        # Convert coordinates to integers
        coordrefs = coordtoint.(possiblecoords,length(relrange))

        # In integer representation, sort, append to the list for the current circle, and remove any repeats
        # sort!(coordrefs)
        append!(currcirclerefs,coordrefs)
        sort!(currcirclerefs)
        unique!(currcirclerefs)
        return currcirclerefs
    end

    # Convert to coordinates from integer representations
    currcirclecoords = inttocoord.(currcirclerefs,length(relrange))

    # Increment the MRP array at all the points indicated as being close to this circle
    manyarrayincrements!(mrp3darray,currcirclecoords)

end

function filter_coords!(possiblecoords::Array, mrp3darray::Array)
    for k in 1:3
        filter!(coord -> coord[k] ≤ size(mrp3darray,1), possiblecoords)
        filter!(coord -> coord[k] ≥ 1, possiblecoords)
    end
end

function filter_outside_sphere(possiblemrps, possiblecoords, resolution)
    # Filter out coordinates where the entire pixel lies outside the unit sphere
    insidesphere = norm.(possiblemrps) .≤ 1.0+resolution*sqrt(3.0)/2.0
    possiblecoords = possiblecoords[insidesphere]
    possiblemrps = possiblemrps[insidesphere]
    return possiblecoords, possiblemrps
end

function filter_far_from_circle(possiblemrps, possiblecoords, allcirclecentre, allcirclenormal, allcircleradius, coordtol, resolution)
    # Filter out coordinates where the entire pixel lies too far away from the circle
    closetocircle = fullcircledistcheck.(possiblemrps,(allcirclecentre,),(allcirclenormal,),allcircleradius,coordtol*resolution)
    possiblecoords = possiblecoords[closetocircle]
    return possiblecoords
end

function construct_mrp_coords(mrpfromtheta, resolution, prevresolution, centre, coordtol)
    # Construct block of MRP coords around this MRP space point
    #= NOTE: The sqrt(5)/2 factor comes from a worst case scenario of finding lattice points 
    within coordtol of a straight line that we're discretising.
    I'm still not confident of this! =#
    if coordtol < 2.0/sqrt(5.0)
        possiblecoords = eightwayround(mrptocoord_unrounded(mrpfromtheta,resolution,prevresolution,centre))
    else
        possiblecoords = roundcoordtoblock(mrptocoord_unrounded(mrpfromtheta,resolution,prevresolution,centre),coordtol*sqrt(5.0)/2.0)
    end
    return possiblecoords
end

    ######## LINES FROM INDIVIDUAL POINTS
function plot_lines(badreflloc_normed, prevresolution, centre, allr_badrefl, deltamax, resolution, firstround, relrange, mrp3darray)
    #= Quick loop through all bad reflection points to see which are associated with lines to keep =#
    linestokeep = find_rellines(badreflloc_normed, sqrt(3.0)/2.0*prevresolution, centre)
    
    # Strip out only the lines to keep
    numlines = length(linestokeep)
    badrefl_lines = badreflloc_normed[linestokeep]
    allr_lines = allr_badrefl[linestokeep]

    # Plot each line as points (within tolerance) in mrp3darray
    for i in 1:numlines
        # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
        coordtol = get_mrp_tol_line(deltamax, allr_lines, i, resolution)

        srange = get_s_range(firstround, coordtol, resolution, prevresolution, badrefl_lines, centre, i)
        currlinerefs = Vector{Int64}()

        for s in srange
            possiblecoords = construct_mrp_coords(coordtol, resolution, prevresolution, centre, s, badrefl_lines, i)
            filter_coords!(possiblecoords, mrp3darray)
            possiblemrps = coordtomrp.(possiblecoords, (relrange,), (centre,))
            filter_points(possiblecoords, possiblemrps, resolution, badrefl_lines, i, coordtol)
            coordrefs = coordtoint.(possiblecoords,length(relrange))
            append_unique_sorted!(currlinerefs, coordrefs)
        end
        currlinecoords = inttocoord.(currlinerefs,length(relrange))
        manyarrayincrements!(mrp3darray, currlinecoords)
    end
end

function get_mrp_tol_line(deltamax, allr_lines, i, resolution)
    # Tolerance in MRP "coordinate" space (i.e. zoomed out by 1/resolution) for pixels to be stored
    # Note that this includes sqrt(3)/2 (a ball around the centre point) as well as an additional component based on mrptol
    mrptol = calc_mrptol(deltamax,allr_lines[i])
    coordtol = 0.5*sqrt(3.0) + mrptol/resolution
    return coordtol
end

function get_s_range(firstround, coordtol, resolution, prevresolution, badrefl_lines, centre, i)
    if firstround 
        # In the first round, take a range of theta values from -pi to pi in steps of an MRP-space tolerance (based on coordtol)
        srange = make_relrange(coordtol*resolution, 2.0)
    else
        # In all other rounds, calculate an smid value based on the closest point between the circle and the centre.
        s_mid = dot(badrefl_lines[i], centre)
        srange = s_mid .+ make_relrange(coordtol*resolution, sqrt(3.0)*prevresolution)
    end
    return srange
end

function construct_mrp_coords(coordtol, resolution, prevresolution, centre, s, badrefl_lines, i)
    mrpfroms = s*badrefl_lines[i]
    if coordtol < 2.0/sqrt(5.0)
        possiblecoords = eightwayround(mrptocoord_unrounded(mrpfroms, resolution, prevresolution, centre))
    else
        possiblecoords = roundcoordtoblock(mrptocoord_unrounded(mrpfroms, resolution, prevresolution, centre), coordtol*sqrt(5.0)/2.0)
    end
    return possiblecoords
end

function filter_points(possiblecoords, possiblemrps, resolution, badrefl_lines, i, coordtol)
    insidesphere = norm.(possiblemrps) .≤ 1.0 + resolution*sqrt(3.0)/2.0
    possiblecoords = possiblecoords[insidesphere]
    possiblemrps = possiblemrps[insidesphere]

    closetoline = linedistcheck.(possiblemrps, (badrefl_lines[i],), coordtol*resolution)
    possiblecoords = possiblecoords[closetoline]
end

function append_unique_sorted!(currlinerefs, coordrefs)
    append!(currlinerefs, coordrefs)
    sort!(currlinerefs)
    unique!(currlinerefs)
end

function filter_coords!(possiblecoords, mrp3darray)
    for k in 1:3
        filter!(coord -> coord[k] ≤ size(mrp3darray,1), possiblecoords)
        filter!(coord -> coord[k] ≥ 1, possiblecoords)
    end
end


######## GREAT CIRCLES FROM INDIVIDUAL POINTS
function process_all_great_circles(badreflloc_normed, allr_badrefl, deltamax, resolution, firstround, prevresolution, centre, relrange, mrp3darray)
    greatcirclestokeep = find_relgreatcircles(badreflloc_normed, sqrt(3.0)/2.0*prevresolution, centre)

    numgreatcircles = length(greatcirclestokeep)
    allnormal_greatcircle = badreflloc_normed[greatcirclestokeep]
    allr_greatcircle = allr_badrefl[greatcirclestokeep]

    (allugreatcircle, allvgreatcircle) = manygreatcirclebasis(allnormal_greatcircle)

    for i in 1:numgreatcircles
        coordtol = get_mrp_tol_great(deltamax, allr_greatcircle, i, resolution)
        thetarange = get_theta_range(firstround, coordtol, resolution, prevresolution, allugreatcircle, allvgreatcircle, allnormal_greatcircle, centre, i)
        currcirclerefs = Vector{Int64}()

        for theta in thetarange
            possiblecoords = construct_great_circle_coords(coordtol, resolution, prevresolution, theta, allugreatcircle, allvgreatcircle, centre, i)
            filter_coords!(possiblecoords, mrp3darray)
            possiblemrps = coordtomrp.(possiblecoords, (relrange,), (centre,))
            filter_great_circle_points(possiblecoords, possiblemrps, resolution, allnormal_greatcircle, i, coordtol)
            coordrefs = coordtoint.(possiblecoords, length(relrange))
            append_unique_sorted!(currcirclerefs, coordrefs)
        end

        currcirclecoords = inttocoord.(currcirclerefs, length(relrange))
        manyarrayincrements!(mrp3darray, currcirclecoords)
    end
end

function get_mrp_tol_great(deltamax, allr_greatcircle, i, resolution)
    mrptol = calc_mrptol(deltamax, allr_greatcircle[i])
    coordtol = 0.5*sqrt(3.0) + mrptol/resolution
    return coordtol
end

function get_theta_range(firstround, coordtol, resolution, prevresolution, allugreatcircle, allvgreatcircle, allnormal_greatcircle, centre, i)
    if firstround
        thetarange = make_relrange(coordtol*resolution, 2*pi)
    else
        centrevec_incircleplane = centre - dot(centre, allnormal_greatcircle[i])*allnormal_greatcircle[i]
        theta_mid = atan(dot(centrevec_incircleplane, allugreatcircle[i]), dot(centrevec_incircleplane, allvgreatcircle[i]))
        thetarange = theta_mid .+ make_relrange(coordtol*resolution, sqrt(3.0)*prevresolution)
    end
    return thetarange
end

function construct_great_circle_coords(coordtol, resolution, prevresolution, theta, allugreatcircle, allvgreatcircle, centre, i)
    mrpfromtheta = circleangletomrp(theta, allugreatcircle[i], allvgreatcircle[i], SVector{3, Float64}(0.0, 0.0, 0.0))
    if coordtol < 2.0/sqrt(5.0)
        possiblecoords = eightwayround(mrptocoord_unrounded(mrpfromtheta, resolution, prevresolution, centre))
    else
        possiblecoords = roundcoordtoblock(mrptocoord_unrounded(mrpfromtheta, resolution, prevresolution, centre), coordtol*sqrt(5.0)/2.0)
    end
    return possiblecoords
end

function filter_great_circle_points(possiblecoords, possiblemrps, resolution, allnormal_greatcircle, i, coordtol)
    insidesphere = norm.(possiblemrps) .≤ 1.0+resolution*sqrt(3.0)/2.0
    possiblecoords = possiblecoords[insidesphere]
    possiblemrps = possiblemrps[insidesphere]

    closetocircle = fullcircledistcheck.(possiblemrps, (SVector{3, Float64}(0, 0, 0),), (allnormal_greatcircle[i],), 1.0, coordtol*resolution)
    possiblecoords = possiblecoords[closetocircle]
end
