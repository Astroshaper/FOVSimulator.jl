"""
    image_generation.jl

Functions for generating images based on camera FOV information.
Provides functionality for thermal image generation with bounding box optimization.
"""

"""
    focal_length(fov_angle::Float64, n_pixel::Int) -> f

Calculate a focal length of a camera.

# Arguments
- `fov_angle` : Field of view angle [deg]
- `n_pixel`   : Number of pixels
"""
focal_length(fov_angle::Real, n_pixel::Int) = n_pixel / (2 * tan(deg2rad(fov_angle) / 2))

"""
    generate_pixel_rays(cam::SpiceCamera) -> Matrix{Ray}

Generate rays for each pixel based on camera FOV information in the camera frame.

# Arguments
- `cam` : SpiceCamera object

# Returns
- `rays` : 2D array of rays corresponding to each pixel (same shape as the image)
"""
function generate_pixel_rays(cam::SpiceCamera)
    # Get camera parameters
    fov_angles = cam.static.fov_angles  # Field of view angles (width, height) in degrees
    img_size = cam.static.img_size      # Image size (width, height) in pixels
    width, height = img_size
    
    # Calculate focal lengths in pixels
    fov_x, fov_y = fov_angles
    f_x = focal_length(fov_x, width)   # focal length in the x-direction [pixel]
    f_y = focal_length(fov_y, height)  # focal length in the y-direction [pixel]
    
    # Calculate principal point (center of the image)
    c_x = width / 2
    c_y = height / 2
    
    # Initialize ray matrix
    rays = Matrix{Ray}(undef, height, width)
    
    # Generate rays for each pixel
    for v in 1:height
        for u in 1:width
            # Convert pixel coordinates to normalized image coordinates
            # (0,0) is at the center of the image
            x = (u - c_x - 0.5) / f_x  # -0.5 to get the center of the pixel
            y = (v - c_y - 0.5) / f_y
            
            # Create a ray in camera frame
            # Z-axis is the boresight direction
            ray_origin = @SVector zeros(3)
            ray_direction = normalize(SVector{3, Float64}(x, y, 1.0))
            
            # Create ray with camera position as origin and calculated direction
            rays[v, u] = Ray(ray_origin, ray_direction)
        end
    end
    
    return rays
end

"""
    generate_pixel_rays(cam::SpiceCamera, asteroid::SpiceAsteroid) -> Matrix{Ray}

Generate rays for each pixel based on camera FOV information and transform them to the asteroid-fixed frame.

# Arguments
- `cam`     : SpiceCamera object
- `asteroid`: SpiceAsteroid object

# Returns
- `rays`    : 2D array of rays in the asteroid-fixed frame corresponding to each pixel
"""
function generate_pixel_rays(cam::SpiceCamera, asteroid::SpiceAsteroid)
    # Check if camera and asteroid states are updated with the same ephemeris time
    if cam.state.et != asteroid.state.et
        error("The ephemeris time of the camera and asteroid must match. Update the camera and asteroid state before calling this function.")
    end

    # Check if camera and asteroid states use the same aberration correction flag
    if cam.state.abcorr != asteroid.state.abcorr
        error("The `abcorr` flag of the camera and asteroid must match. Update the camera and asteroid state before calling this function.")
    end

    # Get ephemeris time and aberration correction flag
    et = cam.state.et
    abcorr = cam.state.abcorr

    # Generate rays in camera frame
    rays = generate_pixel_rays(cam)
    
    # Get transformation from camera frame to asteroid-fixed frame
    from = cam.static.fov_frame        # Camera frame
    to = asteroid.static.asteroid_fixed_frame  # Asteroid-fixed frame
    
    # Calculate rotation matrix from camera frame to asteroid-fixed frame
    Rot = RotMatrix{3}(SPICE.pxform(from, to, et))
    
    # Get camera position in asteroid-fixed frame
    camera_pos_asteroid_frame = SVector{3}(SPICE.spkpos(cam.static.name, et, to, abcorr, asteroid.static.name)[1]) * 1000
    
    # Get dimensions of the ray matrix
    height, width = size(rays)
    
    # Transform each ray from camera frame to asteroid-fixed frame
    for v in 1:height
        for u in 1:width
            # Transform ray direction from camera frame to asteroid-fixed frame
            ray_direction = Rot * rays[v, u].direction
            
            # Update ray with new origin and direction
            rays[v, u] = Ray(camera_pos_asteroid_frame, ray_direction)
        end
    end
    
    return rays
end

"""
    generate_intersection_map(cam::SpiceCamera, asteroid::SpiceAsteroid) -> intersection_map::Matrix{RayShapeIntersectionResult}

Generate an intersection map by casting rays from the camera to the asteroid shape model.
Uses bounding box optimization for faster computation.

# Arguments
- `cam`     : SpiceCamera object
- `asteroid`: SpiceAsteroid object

# Returns
- `intersections` : 2D array of intersection results for each pixel
"""
function generate_intersection_map(cam::SpiceCamera, asteroid::SpiceAsteroid)
    # Generate rays for each pixel in the asteroid-fixed frame
    rays = generate_pixel_rays(cam, asteroid)
    
    # Get image dimensions
    height, width = size(rays)
    
    # Initialize array to store intersection results
    intersection_map = Matrix{RayShapeIntersectionResult}(undef, height, width)
    
    # Perform intersection test for each pixel
    for v in 1:height
        for u in 1:width
            intersection_map[v, u] = intersect_ray_shape(rays[v, u], asteroid)
        end
    end
    
    return intersection_map
end

"""
    generate_image_radiance(intersection_map::Matrix{RayShapeIntersectionResult}, emissivities::Vector{Float64}, temperatures::Vector{Float64}) -> image::Matrix{Float64}

Generate a thermal image based on intersection results and surface properties.

# Arguments
- `intersection_map` : 2D array of intersection results (from `generate_intersection_map`)
- `emissivities`     : Emissivity values (0-1) for each face, length should match the number of faces
- `temperatures`     : Temperature values [K] for each face, length should match the number of faces

# Returns
- `image`         : 2D array of radiance values [W/m²/sr] for each pixel
"""
function generate_image_radiance(intersection_map::Matrix{RayShapeIntersectionResult}, emissivities::Vector{Float64}, temperatures::Vector{Float64})
    # Get image dimensions
    height, width = size(intersection_map)
    
    # Initialize array to store radiance values
    image = zeros(Float64, height, width)
    
    # Calculate radiance for each pixel
    for v in 1:height
        for u in 1:width
            # Get intersection result
            intersection = intersection_map[v, u]
            
            # If intersection exists
            if intersection.hit
                # Get face index
                face_index = intersection.face_index
                
                # Get emissivity and temperature
                ε = emissivities[face_index]
                T = temperatures[face_index]
                
                # Calculate radiance using Stefan-Boltzmann law
                # E = ε·σT⁴ [W/m²] (total emitted power per unit area)
                # L = E/π [W/m²/sr] (radiance assuming Lambertian surface)
                radiance = ε * σ_SB * T^4 / π  # [W/m²/sr]
                
                # Store radiance
                image[v, u] = radiance
            else
                # No intersection (space)
                image[v, u] = 0.0
            end
        end
    end
    
    return image
end
