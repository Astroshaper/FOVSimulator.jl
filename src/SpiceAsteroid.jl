"""
    struct SpiceAsteroidStatic

A structure holding static asteroid information.

# Fields
- `name`                : Asteroid name
- `shape`               : Shape model
- `asteroid_fixed_frame`: Name of the asteroid-fixed reference frame
"""
struct SpiceAsteroidStatic
    name::String
    shape::ShapeModel
    asteroid_fixed_frame::String
end

"""
    mutable struct SpiceAsteroidState

A structure holding dynamic asteroid state.

# Fields
- `et`       : Ephemeris time of the last update
- `ref`      : Reference frame of the last update
- `abcorr`   : Aberration correction flag of the last update
- `obs`      : Observing body name of the last update
- `position` : Asteroid position
- `velocity` : Asteroid velocity
"""
mutable struct SpiceAsteroidState
    et::Float64
    ref::String
    abcorr::String
    obs::String
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
    SpiceAsteroid(name::String, shape::ShapeModel, asteroid_fixed_frame::String)

Construct an asteroid model.

# Arguments
- `name`                : Asteroid name
- `shape`               : Shape model
- `asteroid_fixed_frame`: Name of the asteroid-fixed reference frame
"""
function SpiceAsteroid(name::String, shape::ShapeModel, asteroid_fixed_frame::String)
    # Create static information structure
    static = SpiceAsteroidStatic(name, shape, asteroid_fixed_frame)
    
    # Initialize dynamic information
    et = 0.0
    ref = ""
    abcorr = ""
    obs = ""
    position = @SVector zeros(3)
    velocity = @SVector zeros(3)
    
    # Create dynamic information structure
    state = SpiceAsteroidState(et, ref, abcorr, obs, position, velocity)
    
    # Create asteroid object
    asteroid = SpiceAsteroid(static, state)
    
    return asteroid
end


function Base.show(io::IO, asteroid::SpiceAsteroid)
    msg =  "Asteroid parameters\n"
    msg *= "-------------------\n"
    msg *= "Asteroid name : $(asteroid.static.name)\n"
    if !isempty(asteroid.static.asteroid_fixed_frame)
        msg *= "Fixed frame   : $(asteroid.static.asteroid_fixed_frame)\n"
    end
    if !isempty(asteroid.state.ref)
        msg *= "Last update   : et = $(asteroid.state.et), ref = $(asteroid.state.ref), "
        msg *= "abcorr = $(asteroid.state.abcorr), obs = $(asteroid.state.obs)\n"
    end
    msg *= "Position      : $(asteroid.state.position)\n"
    msg *= "Velocity      : $(asteroid.state.velocity)\n"
    msg *= "-------------------\n"
    
    print(io, msg)
    println(io, asteroid.static.shape)
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
    # Update state parameters
    asteroid.state.et = et
    asteroid.state.ref = ref
    asteroid.state.abcorr = abcorr
    asteroid.state.obs = obs
    
    # Update position and velocity
    state, _ = SPICE.spkezr(asteroid.static.name, et, ref, abcorr, obs)
    asteroid.state.position = state[1:3] * 1000
    asteroid.state.velocity = state[4:6] * 1000

    return
end

"""
    intersect_ray_shape(ray::Ray, asteroid::SpiceAsteroid) -> RayShapeIntersectionResult

Perform ray-shape intersection test with a `SpiceAsteroid`.

# Arguments
- `ray`      : Ray at the asteroid-fixed frame
- `asteroid` : `SpiceAsteroid` object

# Returns
- `RayShapeIntersectionResult` object containing the intersection test result
"""
function intersect_ray_shape(ray::AsteroidShapeModels.Ray, asteroid::SpiceAsteroid)
    return AsteroidShapeModels.intersect_ray_shape(ray, asteroid.static.shape)
end
