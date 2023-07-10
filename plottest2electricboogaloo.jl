using Plots

# Find list of mrps with many overlapping circles
function plottest2electricboogaloo(
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
    
pyplot()

# Assuming allcirclecentre is an array of SVector{3,Float64}
circle_x = [coord[1] for coord in allcirclecentre]
circle_y = [coord[2] for coord in allcirclecentre]
circle_z = [coord[3] for coord in allcirclecentre]

scatter3d(circle_x, circle_y, circle_z, color=:blue, legend=false, markersize=5)  # Adjust markersize as needed

# To plot the unit sphere, we can use the surface function with spherical coordinates.
theta = 0:0.01:2π
phi = 0:0.01:π
x = cos.(theta) * sin.(phi)'
y = sin.(theta) * sin.(phi)'
z = ones(length(theta)) * cos.(phi)'

surface!(x, y, z)

title!("Circles in MRP Space")
xlabel!("X")
ylabel!("Y")
zlabel!("Z")

gui()
end