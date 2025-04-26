"""

"""
function transform_shape(asteroid::SpiceAsteroid, from::String, to::String, et::Float64, abcorr::String, obs::String)

    Rot = RotMatrix{3}(SPICE.pxform(from, to, et))
    obs_pos = SVector{3}(SPICE.spkpos(asteroid.static.name, et, to, abcorr, obs)[1]) * 1000

    nodes = [Rot * node + obs_pos for node in asteroid.static.shape.nodes]
    faces = asteroid.static.shape.faces
    
    face_centers = [AsteroidThermoPhysicalModels.face_center(nodes[face]) for face in faces]
    face_normals = [AsteroidThermoPhysicalModels.face_normal(nodes[face]) for face in faces]
    face_areas   = [AsteroidThermoPhysicalModels.face_area(nodes[face])   for face in faces]

    visiblefacets = asteroid.static.shape.visiblefacets

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

"""
function project_point_fov(p::SVector{3, Float64}, fov_angles::Tuple{Float64, Float64}, img_size::Tuple{Int, Int})
    fov_x, fov_y = fov_angles
    width, height = img_size

    f_x = focal_length(fov_x, width)   # focal length in the x-direction [pixel]
    f_y = focal_length(fov_y, height)  # focal length in the y-direction [pixel]

    c_x = width  / 2  # x-coordinate of the principal point
    c_y = height / 2  # y-coordinate of the principal point

    u = f_x * p[1] / p[3] + c_x
    v = f_y * p[2] / p[3] + c_y
    
    return (round(Int, u), round(Int, v))

    # Points behind the camera (p[3]<0) are not visible
end
