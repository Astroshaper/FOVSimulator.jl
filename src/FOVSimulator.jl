module FOVSimulator

using LinearAlgebra
using Rotations
using StaticArrays

import SPICE
using AsteroidThermoPhysicalModels

include("SpiceCamera.jl")
export SpiceCameraStatic, SpiceCameraState, SpiceCamera, update!

include("SpiceSpacecraft.jl")
export SpiceSpacecraftStatic, SpiceSpacecraftState, SpiceSpacecraft, add_instrument!

include("SpiceAsteroid.jl")
export SpiceAsteroidStatic, SpiceAsteroidState, SpiceAsteroid, transform_shape

## Utility functions
angle_rad(v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real}) = acos(clamp(normalize(v1) ⋅ normalize(v2), -1.0, 1.0))
angle_deg(v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real}) = rad2deg(angle_rad(v1, v2))

angle_rad(v1::AbstractVector{<:AbstractVector{<:Real}}, v2::AbstractVector{<:AbstractVector{<:Real}}) = angle_rad.(v1, v2)
angle_deg(v1::AbstractVector{<:AbstractVector{<:Real}}, v2::AbstractVector{<:AbstractVector{<:Real}}) = angle_deg.(v1, v2)


"""
    solar_phase_angle(sun, target, observer) -> ∠STO

Calculate a sun-target-observer angle (phase angle).

# Arguments
- `sun`     : Sun position vector
- `target`  : Target position vector
- `observer`: Observer position vector

# Return
- ∠STO : Sun-target-observer angle (phase angle) [rad]
"""
solar_phase_angle(sun::AbstractVector{<:Real}, target::AbstractVector{<:Real}, observer::AbstractVector{<:Real}) = angle_rad(sun - target, observer - target)


"""
    solar_elongation_angle

Calculate a sun-observer-target angle (solar elongation angle).

# Arguments
- `sun`     : Sun position vector
- `observer`: Observer position vector
- `target`  : Target position vector

# Return
- ∠SOT : Sun-observer-target angle (solar elongation angle) [rad]
"""
solar_elongation_angle(sun::AbstractVector{<:Real}, observer::AbstractVector{<:Real}, target::AbstractVector{<:Real}) = angle_rad(sun - observer, target - observer)


export angle_rad, angle_deg
export solar_phase_angle, solar_elongation_angle

end # module FOVSimulator
