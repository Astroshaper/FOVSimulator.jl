"""
    mutable struct SpiceAsteroid

# Fields
- `_name_`   : Asteroid name
- `position` : Asteroid position
- `velocity` : Asteroid velocity
- `shape`    : Shape model
"""
mutable struct SpiceAsteroid
    _name_   ::String
    position ::SVector{3, Float64}
    velocity ::SVector{3, Float64}
    shape    ::ShapeModel
end


"""
"""
function SpiceAsteroid(_name_::String, shape::ShapeModel)
    asteroid = SpiceAsteroid(_name_, zeros(3), zeros(3), shape)

    return asteroid
end


function Base.show(io::IO, asteroid::SpiceAsteroid)
    msg =  "Asteroid parameters\n"
    msg *= "-------------------\n"
    msg *= "Asteroid name : $(asteroid._name_)\n"
    msg *= "Position        : $(asteroid.position)\n"
    msg *= "Velocity        : $(asteroid.velocity)\n"
    msg *= "-------------------\n"
    
    print(io, msg)
    println(asteroid.shape)
end


"""
    update!(asteroid::SpiceAsteroid, et::Float64, ref::String, abcorr::String, obs::String)

Update position and velocity vectors of asteroid at a reference frame `ref` and ephemeris time `et`.

# Arguments
- `asteroid` : Asteroid
- `et`       : Ephemeris time
- `ref`      : Target frame
- `abcorr`   : Aberration correction flag.
- `obs`      : Observing body name.
"""
function update!(asteroid::SpiceAsteroid, et::Float64, ref::String, abcorr::String, obs::String)

    state, _ = SPICE.spkezr(asteroid._name_, et, ref, abcorr, obs)
    asteroid.position = state[1:3] * 1000
    asteroid.velocity = state[4:6] * 1000

    return
end


"""

"""
function transform_shape(asteroid::SpiceAsteroid, from::String, to::String, et::Float64, abcorr::String, obs::String)

    Rot = RotMatrix{3}(SPICE.pxform(from, to, et))
    obs_pos = SVector{3}(SPICE.spkpos(asteroid._name_, et, to, abcorr, obs)[1]) * 1000

    nodes = [Rot * node + obs_pos for node in asteroid.shape.nodes]
    faces = asteroid.shape.faces
    
    face_centers = [AsteroidThermoPhysicalModels.face_center(nodes[face]) for face in faces]
    face_normals = [AsteroidThermoPhysicalModels.face_normal(nodes[face]) for face in faces]
    face_areas   = [AsteroidThermoPhysicalModels.face_area(nodes[face])   for face in faces]

    visiblefacets = asteroid.shape.visiblefacets

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
focal_length(fov_angle::Float64, n_pixel::Int) = n_pixel / (2 * tan(deg2rad(fov_angle) / 2))


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

    # カメラの背面(p[3]<0)は写らない
end