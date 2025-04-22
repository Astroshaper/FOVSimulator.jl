"""
    mutable struct SpiceSpaceCraft

# Fields
- `_name_`      : Spacecraft name
- `position`    : Spacecraft position
- `velocity`    : Spacecraft velocity
- `instruments` : Instruments
"""
mutable struct SpiceSpaceCraft
    _name_       ::String

    position    ::SVector{3, Float64}
    velocity    ::SVector{3, Float64}
    instruments ::Dict{String, SpiceCamera}
end


function SpiceSpaceCraft(_name_::String)
    spacecraft = SpiceSpaceCraft(_name_, zeros(3), zeros(3), Dict{String, SpiceCamera}())

    return spacecraft
end


function Base.show(io::IO, spacecraft::SpiceSpaceCraft)
    msg =  "Spacecraft parameters\n"
    msg *= "---------------------\n"
    msg *= "Spacecraft name : $(spacecraft._name_)\n"
    msg *= "Position        : $(spacecraft.position)\n"
    msg *= "Velocity        : $(spacecraft.velocity)\n"
    msg *= "Instruments     : $(keys(spacecraft.instruments))\n"
    msg *= "-----------------\n"
    
    print(io, msg)
end


function add_instrument!(spacecraft::SpiceSpaceCraft, instrument::SpiceCamera)
    spacecraft.instruments[instrument._name_] = instrument

    return
end

"""
    update!(spacecraft::SpiceSpaceCraft, et::Float64, ref::String, abcorr::String, obs::String)

Upadate a spacecraft state vector at ephemeris time `et` and a reference frame `ref`.

# Arguments
- `spacecraft` : Target spacecraft.
- `et`         : Observer epoch.
- `ref`        : Reference frame of output state vector.
- `abcorr`     : Aberration correction flag.
- `obs`        : Observing body name.

c.f. https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html
"""
function update!(spacecraft::SpiceSpaceCraft, et::Float64, ref::String, abcorr::String, obs::String)
    state, _ = SPICE.spkezr(spacecraft._name_, et, ref, abcorr, obs)

    spacecraft.position = state[1:3] * 1000
    spacecraft.velocity = state[4:6] * 1000

    for (_, instrument) in spacecraft.instruments
        update!(instrument, et, ref, abcorr, obs)
    end

    return
end
