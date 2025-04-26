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
