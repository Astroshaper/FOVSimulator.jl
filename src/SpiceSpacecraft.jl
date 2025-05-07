"""
    struct SpiceSpacecraftStatic

A structure holding static spacecraft information.

# Fields
- `name` : Spacecraft name
"""
struct SpiceSpacecraftStatic
    name::String
end

"""
    mutable struct SpiceSpacecraftState

A structure holding dynamic spacecraft state.

# Fields
- `et`         : Ephemeris time of the last update
- `ref`        : Reference frame of the last update
- `abcorr`     : Aberration correction flag of the last update
- `obs`        : Observing body name of the last update
- `position`    : Spacecraft position
- `velocity`    : Spacecraft velocity
- `instruments` : Onboard instruments (cameras, etc.)
"""
mutable struct SpiceSpacecraftState
    et::Float64
    ref::String
    abcorr::String
    obs::String
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
    instruments::Dict{String, SpiceCamera}
end

"""
    struct SpiceSpacecraft

A spacecraft model based on SPICE kernels.

# Fields
- `static` : Static spacecraft information
- `state`  : Dynamic spacecraft state
"""
struct SpiceSpacecraft
    static::SpiceSpacecraftStatic
    state::SpiceSpacecraftState
end


"""
    SpiceSpacecraft(name::String)

Construct a spacecraft model.

# Arguments
- `name` : Spacecraft name
"""
function SpiceSpacecraft(name::String)
    # Create static information structure
    static = SpiceSpacecraftStatic(name)
    
    # Initialize dynamic information
    et = 0.0
    ref = ""
    abcorr = ""
    obs = ""
    position = @SVector zeros(3)
    velocity = @SVector zeros(3)
    instruments = Dict{String, SpiceCamera}()
    
    # Create dynamic information structure
    state = SpiceSpacecraftState(et, ref, abcorr, obs, position, velocity, instruments)
    
    # Create spacecraft object
    spacecraft = SpiceSpacecraft(static, state)
    
    return spacecraft
end


function Base.show(io::IO, spacecraft::SpiceSpacecraft)
    msg =  "Spacecraft parameters\n"
    msg *= "---------------------\n"
    msg *= "Spacecraft name : $(spacecraft.static.name)\n"
    if !isempty(spacecraft.state.ref)
        msg *= "Last update    : et = $(spacecraft.state.et), ref = $(spacecraft.state.ref), "
        msg *= "abcorr = $(spacecraft.state.abcorr), obs = $(spacecraft.state.obs)\n"
    end
    msg *= "Position        : $(spacecraft.state.position)\n"
    msg *= "Velocity        : $(spacecraft.state.velocity)\n"
    msg *= "Instruments     : $(keys(spacecraft.state.instruments))\n"
    msg *= "---------------------\n"
    
    print(io, msg)
end


"""
    add_instrument!(spacecraft::SpiceSpacecraft, instrument::SpiceCamera)

Add an instrument to the spacecraft.

# Arguments
- `spacecraft` : Spacecraft
- `instrument` : Instrument to add (camera, etc.)
"""
function add_instrument!(spacecraft::SpiceSpacecraft, instrument::SpiceCamera)
    spacecraft.state.instruments[instrument.static.name] = instrument

    return
end

"""
    update!(spacecraft::SpiceSpacecraft, et::Float64, ref::String, abcorr::String, obs::String)

Update spacecraft state at the specified time and frame.

# Arguments
- `spacecraft` : Spacecraft
- `et`         : Ephemeris time
- `ref`        : Target frame
- `abcorr`     : Aberration correction flag
- `obs`        : Observing body name

c.f. https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html
"""
function update!(spacecraft::SpiceSpacecraft, et::Float64, ref::String, abcorr::String, obs::String)
    # Update state parameters
    spacecraft.state.et = et
    spacecraft.state.ref = ref
    spacecraft.state.abcorr = abcorr
    spacecraft.state.obs = obs
    
    # Update position and velocity
    state, _ = SPICE.spkezr(spacecraft.static.name, et, ref, abcorr, obs)
    spacecraft.state.position = state[1:3] * 1000
    spacecraft.state.velocity = state[4:6] * 1000

    # Update instruments
    for (_, instrument) in spacecraft.state.instruments
        update!(instrument, et, ref, abcorr, obs)
    end

    return
end
