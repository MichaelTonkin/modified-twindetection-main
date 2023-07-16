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

    # Plot each circle as points (within tolerance) in mrp3darray
    for i in 1:numcircles
        
        # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
        #deltamax is radius of the cloud around B.
        #rdiff is the difference in the distance from the origin of A and B
        mrptol = calc_mrptol(sqrt(deltamax^2-rdiff[cld(i,2)]^2),allr_dest_circles[i])

        # Tolerance in MRP "coordinate" space (i.e. zoomed out by 1/resolution) for voxels to be stored
        # Note that this includes sqrt(3)/2 (a ball around the centre point) as well as an additional component based on mrptol
        coordtol = 0.5*sqrt(3.0) + mrptol/resolution
        if firstround
            # In the first round, take a range of thetavalues from -phi to phi in steps of an MRP-space tolerance (based on coordtol)
            # divided by the circle radius to get a resolution in angle space
            thetarange = make_relrange(coordtol*resolution/allcircleradius[i],2*allphi[i])
        else
            # In all other rounds, calculate a thetamid value based on the closest point between the circle and the centre.

            # Find the part of the centre vector that lies in the plane of the circle relative to the centre of said circle
            centrevec_incircleplane = centre - dot(centre,allcirclenormal[i])*allcirclenormal[i] - allcirclecentre[i]

            # Convert this into an angle based on u and v
            theta_mid = atan(dot(centrevec_incircleplane,allu[i]),dot(centrevec_incircleplane,allv[i]))
                           
            # Create a range based on assuming a straight line with (worst case) length based on going straight through centre of voxel
            thetarange = theta_mid .+ make_relrange(coordtol*resolution/allcircleradius[i],sqrt(3.0)*prevresolution)
        end

        # Initialise list of coordinates associated with current circle
        #        currcirclecoords = Vector{SVector{3,Int64}}()
        currcirclerefs = Vector{Int64}()



        # Loop through range of theta values
        for theta in thetarange

            # Convert theta to MRP space location
            mrpfromtheta = circleangletomrp(theta,allu[i],allv[i],allcirclecentre[i])
            
            # Construct block of MRP coords around this MRP space point
            #= NOTE: The sqrt(5)/2 factor comes from a worst case scenario of finding lattice points 
            within coordtol of a straight line that we're discretising.
            I'm still not confident of this! =#
            if coordtol < 2.0/sqrt(5.0)
                possiblecoords = eightwayround(mrptocoord_unrounded(mrpfromtheta,resolution,prevresolution,centre))
            else
                possiblecoords = roundcoordtoblock(mrptocoord_unrounded(mrpfromtheta,resolution,prevresolution,centre),coordtol*sqrt(5.0)/2.0)
            end
            
            # Filter out coordinates outside the range of the grid
            for k in 1:3
                filter!(coord -> coord[k] ≤ size(mrp3darray,1),possiblecoords)
                filter!(coord -> coord[k] ≥ 1,possiblecoords)
            end

            # MRP locations
            possiblemrps = coordtomrp.(possiblecoords,(relrange,),(centre,))

            # Filter out coordinates where the entire voxel lies outside the unit sphere
            insidesphere = norm.(possiblemrps) .≤ 1.0+resolution*sqrt(3.0)/2.0
            possiblecoords = possiblecoords[insidesphere]
            possiblemrps = possiblemrps[insidesphere]

            # Filter out coordinates where the entire voxel lies too far away from the circle
            closetocircle = fullcircledistcheck.(possiblemrps,(allcirclecentre[i],),(allcirclenormal[i],),allcircleradius[i],coordtol*resolution)
            possiblecoords = possiblecoords[closetocircle]

            # Convert coordinates to integers
            coordrefs = coordtoint.(possiblecoords,length(relrange))

            # In integer representation, sort, append to the list for the current circle, and remove any repeats
            # sort!(coordrefs)
            append!(currcirclerefs,coordrefs)
            sort!(currcirclerefs)
            unique!(currcirclerefs)

        end

        # Convert to coordinates from integer representations
        currcirclecoords = inttocoord.(currcirclerefs,length(relrange))

        # Increment the MRP array at all the points indicated as being close to this circle
        manyarrayincrements!(mrp3darray,currcirclecoords)

    end

    ######## LINES FROM INDIVIDUAL POINTS

    #= Quick loop through all bad reflection points to see which are associated with lines to keep =#
    linestokeep = find_rellines(badreflloc_normed,sqrt(3.0)/2.0*prevresolution,centre)

    # Strip out only the lines to keep
    numlines = length(linestokeep)
    badrefl_lines = badreflloc_normed[linestokeep]
    allr_lines = allr_badrefl[linestokeep]

    # Plot each line as points (within tolerance) in mrp3darray
    for i in 1:numlines
        
        # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
        mrptol = calc_mrptol(deltamax,allr_lines[i])

        # Tolerance in MRP "coordinate" space (i.e. zoomed out by 1/resolution) for voxels to be stored
        # Note that this includes sqrt(3)/2 (a ball around the centre point) as well as an additional component based on mrptol
        coordtol = 0.5*sqrt(3.0) + mrptol/resolution
        if firstround
            # In the first round, take a range of theta values from -pi to pi in steps of an MRP-space tolerance (based on coordtol)
            srange = make_relrange(coordtol*resolution,2.0)
        else
            # In all other rounds, calculate an smid value based on the closest point between the circle and the centre.
            s_mid = dot(badrefl_lines[i],centre)

            # Create a range based on worst case length based on going straight through centre of voxel
            srange = s_mid .+ make_relrange(coordtol*resolution,sqrt(3.0)*prevresolution)
        end

        # Initialise list of coordinates associated with current circle
        currlinerefs = Vector{Int64}()

        # Loop through range of s values
        for s in srange

            # Convert s to MRP space location
            mrpfroms = s*badrefl_lines[i]
            
            # Construct block of MRP coords around this MRP space point
            #= NOTE: The sqrt(5)/2 factor comes from a worst case scenario of finding lattice points 
            within coordtol of a straight line that we're discretising.
            I'm still not confident of this! =#
            if coordtol < 2.0/sqrt(5.0)
                possiblecoords = eightwayround(mrptocoord_unrounded(mrpfroms,resolution,prevresolution,centre))
            else
                possiblecoords = roundcoordtoblock(mrptocoord_unrounded(mrpfroms,resolution,prevresolution,centre),coordtol*sqrt(5.0)/2.0)
            end
            
            # Filter out coordinates outside the range of the grid
            for k in 1:3
                filter!(coord -> coord[k] ≤ size(mrp3darray,1),possiblecoords)
                filter!(coord -> coord[k] ≥ 1,possiblecoords)
            end

            # MRP locations
            possiblemrps = coordtomrp.(possiblecoords,(relrange,),(centre,))

            # Filter out coordinates where the entire voxel lies outside the unit sphere
            insidesphere = norm.(possiblemrps) .≤ 1.0+resolution*sqrt(3.0)/2.0
            possiblecoords = possiblecoords[insidesphere]
            possiblemrps = possiblemrps[insidesphere]

            # Filter out coordinates where the entire voxel lies too far away from the line
            closetoline = linedistcheck.(possiblemrps,(badrefl_lines[i],),coordtol*resolution)
            possiblecoords = possiblecoords[closetoline]

            # Convert coordinates to integers
            coordrefs = coordtoint.(possiblecoords,length(relrange))

            # In integer representation, sort, append to the list for the current line, and remove any repeats
            # sort!(coordrefs)
            append!(currlinerefs,coordrefs)
            sort!(currlinerefs)
            unique!(currlinerefs)

        end

        # Convert to coordinates from integer representations
        currlinecoords = inttocoord.(currlinerefs,length(relrange))

        # Increment the MRP array at all the points indicated as being close to this circle
        manyarrayincrements!(mrp3darray,currlinecoords)

    end



    ######## GREAT CIRCLES FROM INDIVIDUAL POINTS

    #= Quick loop through all bad reflection points to see which are associated with great circles to keep =#
    greatcirclestokeep = find_relgreatcircles(badreflloc_normed,sqrt(3.0)/2.0*prevresolution,centre)

    # Strip out only the great circless to keep
    numgreatcircles = length(greatcirclestokeep)
    allnormal_greatcircle = badreflloc_normed[greatcirclestokeep]
    allr_greatcircle = allr_badrefl[greatcirclestokeep]

    # Find basis vectors for the great circles
    (allugreatcircle,allvgreatcircle) = manygreatcirclebasis(allnormal_greatcircle)

    # Plot each line as points (within tolerance) in mrp3darray
    for i in 1:numgreatcircles
        
        # Tolerance in MRP space. If a point in MRP space is within this distance of a circle then it is a possible twin rotation.
        mrptol = calc_mrptol(deltamax,allr_greatcircle[i])

        # Tolerance in MRP "coordinate" space (i.e. zoomed out by 1/resolution) for voxels to be stored
        # Note that this includes sqrt(3)/2 (a ball around the centre point) as well as an additional component based on mrptol
        coordtol = 0.5*sqrt(3.0) + mrptol/resolution
        if firstround
            # In the first round, take a range of thetavalues from -pi to pi in steps of an MRP-space tolerance (based on coordtol)
            thetarange = make_relrange(coordtol*resolution,2*pi)
        else
            # In all other rounds, calculate a thetamid value based on the closest point between the circle and the centre.

            # Find the part of the centre of voxel vector that lies in the plane of the great circle
            centrevec_incircleplane = centre - dot(centre,allnormal_greatcircle[i])*allnormal_greatcircle[i]

            # Convert this into an angle based on u and v
            theta_mid = atan(dot(centrevec_incircleplane,allugreatcircle[i]),dot(centrevec_incircleplane,allvgreatcircle[i]))
            
            # Create a range based on assuming a straight line with (worst case) length based on going straight through centre of voxel
            thetarange = theta_mid .+ make_relrange(coordtol*resolution,sqrt(3.0)*prevresolution)
        end

         # Initialise list of coordinates associated with current circle
        #        currcirclecoords = Vector{SVector{3,Int64}}()
        currcirclerefs = Vector{Int64}()

        # Loop through range of theta values
        for theta in thetarange

            # Convert theta to MRP space location
            mrpfromtheta = circleangletomrp(theta,allugreatcircle[i],allvgreatcircle[i],SVector{3,Float64}(0.0,0.0,0.0))
            
            # Construct block of MRP coords around this MRP space point
            #= NOTE: The sqrt(5)/2 factor comes from a worst case scenario of finding lattice points 
            within coordtol of a straight line that we're discretising.
            I'm still not confident of this! =#
            if coordtol < 2.0/sqrt(5.0)
                possiblecoords = eightwayround(mrptocoord_unrounded(mrpfromtheta,resolution,prevresolution,centre))
            else
                possiblecoords = roundcoordtoblock(mrptocoord_unrounded(mrpfromtheta,resolution,prevresolution,centre),coordtol*sqrt(5.0)/2.0)
            end
            
            # Filter out coordinates outside the range of the grid
            for k in 1:3
                filter!(coord -> coord[k] ≤ size(mrp3darray,1),possiblecoords)
                filter!(coord -> coord[k] ≥ 1,possiblecoords)
            end

            # MRP locations
            possiblemrps = coordtomrp.(possiblecoords,(relrange,),(centre,))

            # Filter out coordinates where the entire voxel lies outside the unit sphere
            insidesphere = norm.(possiblemrps) .≤ 1.0+resolution*sqrt(3.0)/2.0
            possiblecoords = possiblecoords[insidesphere]
            possiblemrps = possiblemrps[insidesphere]

            # Filter out coordinates where the entire voxel lies too far away from the circle
            closetocircle = fullcircledistcheck.(possiblemrps,(SVector{3,Float64}(0,0,0),),(allnormal_greatcircle[i],),1.0,coordtol*resolution)
            possiblecoords = possiblecoords[closetocircle]

            # Convert coordinates to integers
            coordrefs = coordtoint.(possiblecoords,length(relrange))

            # In integer representation, sort, append to the list for the current circle, and remove any repeats
            # sort!(coordrefs)
            append!(currcirclerefs,coordrefs)
            sort!(currcirclerefs)
            unique!(currcirclerefs)

        end

        # Convert to coordinates from integer representations
        currcirclecoords = inttocoord.(currcirclerefs,length(relrange))

        # Increment the MRP array at all the points indicated as being close to this circle
        manyarrayincrements!(mrp3darray,currcirclecoords)

    end



    # Find MRP values for voxels where value recorded in array is larger than keep criterion
    return coordtomrp.(findall(x->x≥keepcriterion,mrp3darray),(make_relrange(resolution,prevresolution),),(centre,))

end