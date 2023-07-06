function findpoint(
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
    )
# Define the MRP parameters of the point of interest
mrp_x, mrp_y, mrp_z = 0, 0, 0  # these are just examples; use the actual values you're interested in

mrp3darray = make_mrparray(resolution,prevresolution)

# Create MRP vector
mrp_vector = SVector{3,Float64}(mrp_x, mrp_y, mrp_z)

# Convert MRP vector to voxel coordinate
voxel_coord = mrptocoord_unrounded(mrp_vector, resolution, prevresolution, centre)

# Round the coordinates to the nearest voxel grid point
rounded_voxel_coords = eightwayround(voxel_coord)

# Ensure coordinates are within the grid boundaries
for rounded_voxel_coord in rounded_voxel_coords
    if CartesianIndex(rounded_voxel_coord...) in axes(mrp3darray)
        # If within boundaries, retrieve the number of overlaps from the mrp3darray
        num_overlaps = mrp3darray[rounded_voxel_coord...]
        println("Number of overlapping rotations at the MRP point ", rounded_voxel_coord, " is: ", num_overlaps)
    else
        # If out of boundaries, handle accordingly
        println("MRP point ", rounded_voxel_coord, " is outside of the grid boundaries")
    end
end



end