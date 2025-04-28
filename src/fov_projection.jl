"""
    transform_shape(asteroid::SpiceAsteroid, from::String, to::String, et::Float64, abcorr::String, obs::String) -> shape_new

Transform the shape of an asteroid to a new reference frame.

# Arguments
- `asteroid` : Asteroid
- `from`     : Reference frame to transform from (asteroid-fixed frame)
- `to`       : Reference frame to transform to
- `et`       : Ephemeris time
- `abcorr`   : Aberration correction
- `obs`      : Observing body name

# Returns
- `shape_new` : Transformed shape model
"""
function transform_shape(asteroid::SpiceAsteroid, from::String, to::String, et::Float64, abcorr::String, obs::String)

    Rot = RotMatrix{3}(SPICE.pxform(from, to, et))
    target_pos = SVector{3}(SPICE.spkpos(asteroid.static.name, et, to, abcorr, obs)[1]) * 1000

    shape = asteroid.static.shape
    nodes = [Rot * node + target_pos for node in shape.nodes]
    faces = shape.faces
    
    face_centers  = [AsteroidThermoPhysicalModels.face_center(nodes[face]) for face in faces]
    face_normals  = [AsteroidThermoPhysicalModels.face_normal(nodes[face]) for face in faces]
    face_areas    = [AsteroidThermoPhysicalModels.face_area(nodes[face])   for face in faces]
    visiblefacets = shape.visiblefacets

    shape_new = ShapeModel(nodes, faces, face_centers, face_normals, face_areas, visiblefacets)
    
    return shape_new
end


"""
    focal_length(fov_angle::Float64, n_pixel::Int) -> f

Calculate a focal length of a camera.

# Arguments
- `fov_angle` : Field of view angle [deg]
- `n_pixel`   : Number of pixels
"""
focal_length(fov_angle::Real, n_pixel::Int) = n_pixel / (2 * tan(deg2rad(fov_angle) / 2))


"""
    project_point_fov(p::SVector{3, Float64}, fov_angles::Tuple{Float64, Float64}, img_size::Tuple{Int, Int}) -> Tuple{Int, Int}

Project a 3D point onto a 2D image plane using the field-of-view angles and image size.

# Arguments
- `p`          : 3D point in space (SVector{3, Float64})
- `fov_angles` : Tuple of field-of-view angles (width [deg], height [deg])
- `img_size`   : Tuple of image size (width [pixels], height [pixels])

# Returns
- `(u, v)`     : 2D pixel coordinates of the projected point (rounded to nearest integer)
"""
function project_point_to_fov(p::SVector{3, Float64}, fov_angles::Tuple{Float64, Float64}, img_size::Tuple{Int, Int})
    fov_x, fov_y = fov_angles
    width, height = img_size

    f_x = focal_length(fov_x, width)   # focal length in the x-direction [pixel]
    f_y = focal_length(fov_y, height)  # focal length in the y-direction [pixel]

    c_x = width  / 2  # x-coordinate of the principal point
    c_y = height / 2  # y-coordinate of the principal point

    u = round(Int, f_x * p[1] / p[3] + c_x)
    v = round(Int, f_y * p[2] / p[3] + c_y)
    
    return (u, v)
end


"""
    struct ProjectedFace

Struct for storing projected face information to a camera frame.

# Fields
- `u`          : pixel u-coordinate (column)
- `v`          : pixel v-coordinate (row)
- `z`          : depth (distance from camera)
- `face_index` : original face index
"""
struct ProjectedFace
    u::Int
    v::Int
    z::Float64
    face_index::Int
end


"""
    project_face_centers(shape::ShapeModel, fov_angles::Tuple{Float64, Float64}, img_size::Tuple{Int, Int}) -> projections

Project the face centers of a ShapeModel onto image coordinates.

# Arguments
- `shape`      : ShapeModel transformed into the camera frame
- `fov_angles` : Tuple of field-of-view angles (width [deg], height [deg])
- `img_size`   : Tuple of image size (width [pixels], height [pixels])

# Returns
- `projections` : Vector of `ProjectedFace` objects
"""
function project_face_centers(shape::ShapeModel, fov_angles::Tuple{Float64, Float64}, img_size::Tuple{Int, Int})
    width, height = img_size
    projected_faces = ProjectedFace[]

    for (i, center) in enumerate(shape.face_centers)
        z = center[3]
        z ≤ 0 && continue  # Skip if the face is behind the camera

        u, v = project_point_fov(center, fov_angles, img_size)

        if 1 ≤ u ≤ width && 1 ≤ v ≤ height
            push!(projected_faces, ProjectedFace(u, v, z, i))
        end
    end

    return projected_faces
end


"""
    map_temperature_to_image(projected_faces::Vector{ProjectedFace}, temperatures::Vector{Float64}, img_size::Tuple{Int, Int}, λ::Float64) -> Array{Float64,2}

Map temperatures to an image using the Planck function to calculate radiance.

# Arguments
- `projected_faces` : Vector of ProjectedFace objects with pixel coordinates and face indices
- `temperatures`    : Vector of temperatures [K] corresponding to face indices
- `img_size`        : Tuple of image size (width [pixels], height [pixels])
- `λ`               : Wavelength [m]

# Returns
- `img`             : 2D array of radiance values [W/m²/sr/μm]
"""
function map_temperature_to_image(projected_faces::Vector{ProjectedFace}, temperatures::Vector{Float64}, img_size::Tuple{Int, Int}, λ::Float64)
    width, height = img_size
    
    # Initialize image and z-buffer
    img = zeros(Float64, height, width)
    z_buffer = fill(Inf, height, width)
    
    # Physical constants
    h = 6.62607015e-34  # Planck constant [J s]
    c = 2.99792458e8    # Speed of light [m/s]
    k = 1.380649e-23    # Boltzmann constant [J/K]
    
    # Process each projected face
    for face in projected_faces
        u, v, z, face_index = face.u, face.v, face.z, face.face_index
        
        # Check if this face is closer to the camera than what's already in the z-buffer
        if z < z_buffer[v, u]
            # Get temperature for this face
            T = temperatures[face_index]
            
            # Calculate radiance using Planck function
            # B(λ,T) = (2hc²/λ⁵) / (exp(hc/λkT) - 1)
            numerator = 2.0 * h * c^2 / λ^5
            denominator = exp((h * c) / (λ * k * T)) - 1.0
            
            # Calculate spectral radiance [W/m²/sr/m]
            radiance = numerator / denominator
            
            # Convert to [W/m²/sr/μm]
            radiance_μm = radiance * 1.0e-6
            
            # Update image and z-buffer
            img[v, u] = radiance_μm
            z_buffer[v, u] = z
        end
    end
    
    return img
end
