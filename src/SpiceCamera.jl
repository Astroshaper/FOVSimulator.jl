
"""
    mutable struct SpiceCamera

# Fields
- `_name_`       : Instrument name
- `_id_`         : Instrument ID
- `_fov_shape_`  : Instrument FOV shape as defined in the SPICE kernel.
- `_fov_frame_`  : Name of the frame in which FOV vectors are defined as defined in the SPICE kernel.
- `_boresight_`  : Boresight vector as defined in the SPICE kernel.
- `_fov_bounds_` : FOV boundary vectors as defined in the SPICE kernel.
- `boresight`    : Boresight vector at a frame/epoch.
- `fov_bounds`   : FOV boundary vectors at a frame/epoch.
- `position`     : Instrument position
- `velocity`     : Instrument velocity
"""
mutable struct SpiceCamera
    _name_ ::String
    _id_   ::Int

    _fov_shape_  ::String
    _fov_frame_  ::String
    _boresight_  ::SVector{3, Float64}
    _fov_bounds_ ::Vector{SVector{3, Float64}}

    boresight  ::SVector{3, Float64}
    fov_bounds ::Vector{SVector{3, Float64}}

    position ::SVector{3, Float64}
    velocity ::SVector{3, Float64}
end


function SpiceCamera(_name_::String, _id_::Int)
    _fov_shape_, _fov_frame_, _boresight_, _fov_bounds_ = SPICE.getfov(_id_)

    boresight = similar(_boresight_)
    fov_bounds = similar(_fov_bounds_)

    position = @SVector zeros(3)
    velocity = @SVector zeros(3)

    cam = SpiceCamera(_name_, _id_, _fov_shape_, _fov_frame_, _boresight_, _fov_bounds_, boresight, fov_bounds, position, velocity)

    return cam
end


function Base.show(io::IO, cam::SpiceCamera)
    msg =  "Camera parameters\n"
    msg *= "-----------------\n"
    msg *= "Instrument name      : $(cam._name_)\n"
    msg *= "Instrument ID        : $(cam._id_)\n"
    msg *= "FOV shape            : $(cam._fov_shape_)\n"
    msg *= "FOV reference frame  : $(cam._fov_frame_)\n"
    msg *= "Boresight vector     : $(cam._boresight_)\n"
    msg *= "FOV boundary vectors : \n"
    for v in cam._fov_bounds_
        msg *= "    $v\n"
    end
    msg *= "-----------------\n"
    msg *= "Current boresight vector     : $(cam.boresight)\n"
    msg *= "Current FOV boundary vectors : \n"
    for v in cam.fov_bounds
        msg *= "    $v\n"
    end
    msg *= "-----------------\n"
    msg *= "Camera position : $(cam.position)\n"
    msg *= "Camera velocity : $(cam.velocity)\n"
    msg *= "-----------------\n"
    
    print(io, msg)
end


"""
    update!(cam::SpiceCamera, target_frame::String, et::Float64)

Update a boresight vector and FOV boundary vectors at an ephemeris time `et` and reference frame `ref`.

# Arguments
- `cam`    : Camera
- `et`     : Ephemeris time
- `ref`    : Target frame
- `abcorr` : Aberration correction flag.
- `obs`    : Observing body name.
"""
function update!(cam::SpiceCamera, et::Float64, ref::String, abcorr::String, obs::String)

    Rot = RotMatrix{3}(SPICE.pxform(cam._fov_frame_, ref, et))

    cam.boresight = Rot * cam._boresight_

    for (i, v) in enumerate(cam._fov_bounds_)
        cam.fov_bounds[i] = Rot * v
    end

    state, _ = SPICE.spkezr(cam._name_, et, ref, abcorr, obs)
    cam.position = state[1:3] * 1000
    cam.velocity = state[4:6] * 1000

    return
end
