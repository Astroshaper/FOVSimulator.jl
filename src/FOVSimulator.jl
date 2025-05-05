module FOVSimulator

using LinearAlgebra
using Rotations
using StaticArrays

import SPICE
using AsteroidThermoPhysicalModels

include("SpiceCamera.jl")
include("SpiceSpacecraft.jl")
include("SpiceAsteroid.jl")
export SpiceCamera
export SpiceSpacecraft, add_instrument!
export SpiceAsteroid
export update!

include("util.jl")
export angle_rad, angle_deg
export solar_phase_angle, solar_elongation_angle

# include("fov_projection.jl")
# export simulate_image

include("ray_intersection.jl")
export Ray, intersect_ray_triangle, intersect_ray_shape

include("ray_generation.jl")
export generate_pixel_rays

end # module FOVSimulator
