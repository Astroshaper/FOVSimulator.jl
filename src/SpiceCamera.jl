"""
    struct SpiceCameraStatic

A structure holding static camera information.

# Fields
- `name`       : Camera name
- `id`         : Camera ID
- `fov_shape`  : FOV shape ("RECTANGLE", "CIRCLE", "ELLIPSE", etc.)
- `fov_frame`  : FOV reference frame
- `boresight`  : Boresight vector in the FOV reference frame
- `bounds`     : FOV boundary vectors in the FOV reference frame
- `fov_angles` : Field of view angles (width, height) in degrees
- `img_size`   : Image size (width, height) in pixels
"""
struct SpiceCameraStatic
    name       ::String
    id         ::Int
    fov_shape  ::String
    fov_frame  ::String
    boresight  ::SVector{3, Float64}
    bounds     ::Vector{SVector{3, Float64}}
    fov_angles ::Tuple{Float64, Float64}
    img_size   ::Tuple{Int, Int}
end

"""
    mutable struct SpiceCameraState

A structure holding dynamic camera state.

# Fields
- `et`       : Ephemeris time of the last update
- `ref`      : Reference frame of the last update
- `abcorr`   : Aberration correction flag of the last update
- `obs`      : Observing body name of the last update
- `boresight` : Boresight vector in the current frame
- `bounds`    : FOV boundary vectors in the current frame
- `position`  : Camera position
- `velocity`  : Camera velocity
"""
mutable struct SpiceCameraState
    et::Float64
    ref::String
    abcorr::String
    obs::String
    boresight::SVector{3, Float64}
    bounds::Vector{SVector{3, Float64}}
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
end

"""
    struct SpiceCamera

A camera model based on SPICE kernels.

# Fields
- `static` : Static camera information
- `state`  : Dynamic camera state
"""
struct SpiceCamera
    static::SpiceCameraStatic
    state::SpiceCameraState
end


function SpiceCamera(name::String, id::Int, fov_angles::Tuple{Float64, Float64}, img_size::Tuple{Int, Int})
    # Get static information from SPICE kernel
    fov_shape, fov_frame, boresight, bounds = SPICE.getfov(id)
    
    # Create static information structure
    static = SpiceCameraStatic(name, id, fov_shape, fov_frame, boresight, bounds, fov_angles, img_size)
    
    # Initialize dynamic information
    et = 0.0
    ref = ""
    abcorr = ""
    obs = ""
    state_boresight = similar(boresight)
    state_bounds = similar(bounds)
    position = @SVector zeros(3)
    velocity = @SVector zeros(3)
    
    # Create dynamic information structure
    state = SpiceCameraState(et, ref, abcorr, obs, state_boresight, state_bounds, position, velocity)
    
    # Create camera object
    cam = SpiceCamera(static, state)
    
    return cam
end


function Base.show(io::IO, cam::SpiceCamera)
    msg =  "Camera parameters\n"
    msg *= "-----------------\n"
    msg *= "Instrument name      : $(cam.static.name)\n"
    msg *= "Instrument ID        : $(cam.static.id)\n"
    msg *= "FOV shape            : $(cam.static.fov_shape)\n"
    msg *= "FOV reference frame  : $(cam.static.fov_frame)\n"
    msg *= "FOV angles (h,v)     : $(cam.static.fov_angles) deg\n"
    msg *= "Image size (w,h)     : $(cam.static.img_size) px\n"
    msg *= "Boresight vector     : $(cam.static.boresight)\n"
    msg *= "FOV boundary vectors : \n"
    for v in cam.static.bounds
        msg *= "    $v\n"
    end
    msg *= "-----------------\n"
    if !isempty(cam.state.ref)
        msg *= "Last update   : et = $(cam.state.et), ref = $(cam.state.ref), "
        msg *= "abcorr = $(cam.state.abcorr), obs = $(cam.state.obs)\n"
    end
    msg *= "Current boresight vector     : $(cam.state.boresight)\n"
    msg *= "Current FOV boundary vectors : \n"
    for v in cam.state.bounds
        msg *= "    $v\n"
    end
    msg *= "-----------------\n"
    msg *= "Camera position : $(cam.state.position)\n"
    msg *= "Camera velocity : $(cam.state.velocity)\n"
    msg *= "-----------------\n"
    
    print(io, msg)
end


"""
    update!(cam::SpiceCamera, et::Float64, ref::String, abcorr::String, obs::String)

Update camera state at the specified time and frame.

# Arguments
- `cam`    : Camera
- `et`     : Ephemeris time
- `ref`    : Target frame
- `abcorr` : Aberration correction flag
- `obs`    : Observing body name
"""
function update!(cam::SpiceCamera, et::Float64, ref::String, abcorr::String, obs::String)
    # Update state parameters
    cam.state.et = et
    cam.state.ref = ref
    cam.state.abcorr = abcorr
    cam.state.obs = obs
    
    # Calculate rotation matrix
    Rot = RotMatrix{3}(SPICE.pxform(cam.static.fov_frame, ref, et))

    # Update boresight vector
    cam.state.boresight = Rot * cam.static.boresight

    # Update boundary vectors
    for (i, v) in enumerate(cam.static.bounds)
        cam.state.bounds[i] = Rot * v
    end

    # Update position and velocity
    state, _ = SPICE.spkezr(cam.static.name, et, ref, abcorr, obs)
    cam.state.position = state[1:3] * 1000
    cam.state.velocity = state[4:6] * 1000

    return
end
