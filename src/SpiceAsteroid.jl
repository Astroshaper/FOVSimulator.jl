"""
    struct SpiceAsteroidStatic

A structure holding static asteroid information.

# Fields
- `name`  : Asteroid name
- `shape` : Shape model
"""
struct SpiceAsteroidStatic
    name::String
    shape::ShapeModel
end

"""
    mutable struct SpiceAsteroidState

A structure holding dynamic asteroid state.

# Fields
- `position` : Asteroid position
- `velocity` : Asteroid velocity
"""
mutable struct SpiceAsteroidState
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
end

"""
    struct SpiceAsteroid

An asteroid model based on SPICE kernels.

# Fields
- `static` : Static asteroid information
- `state`  : Dynamic asteroid state
"""
struct SpiceAsteroid
    static::SpiceAsteroidStatic
    state::SpiceAsteroidState
end


"""
    SpiceAsteroid(name::String, shape::ShapeModel)

Construct an asteroid model.

# Arguments
- `name`  : Asteroid name
- `shape` : Shape model
"""
function SpiceAsteroid(name::String, shape::ShapeModel)
    # Create static information structure
    static = SpiceAsteroidStatic(name, shape)
    
    # Initialize dynamic information
    position = @SVector zeros(3)
    velocity = @SVector zeros(3)
    
    # Create dynamic information structure
    state = SpiceAsteroidState(position, velocity)
    
    # Create asteroid object
    asteroid = SpiceAsteroid(static, state)
    
    return asteroid
end


function Base.show(io::IO, asteroid::SpiceAsteroid)
    msg =  "Asteroid parameters\n"
    msg *= "-------------------\n"
    msg *= "Asteroid name : $(asteroid.static.name)\n"
    msg *= "Position      : $(asteroid.state.position)\n"
    msg *= "Velocity      : $(asteroid.state.velocity)\n"
    msg *= "-------------------\n"
    
    print(io, msg)
    println(asteroid.static.shape)
end


"""
    update!(asteroid::SpiceAsteroid, et::Float64, ref::String, abcorr::String, obs::String)

Update asteroid state at the specified time and frame.

# Arguments
- `asteroid` : Asteroid
- `et`       : Ephemeris time
- `ref`      : Target frame
- `abcorr`   : Aberration correction flag
- `obs`      : Observing body name
"""
function update!(asteroid::SpiceAsteroid, et::Float64, ref::String, abcorr::String, obs::String)
    state, _ = SPICE.spkezr(asteroid.static.name, et, ref, abcorr, obs)
    asteroid.state.position = state[1:3] * 1000
    asteroid.state.velocity = state[4:6] * 1000

    return
end


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
