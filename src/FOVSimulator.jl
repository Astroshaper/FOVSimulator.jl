module FOVSimulator

using LinearAlgebra
using Rotations
using StaticArrays

import SPICE
using AsteroidThermoPhysicalModels

# Physical constants
const σ_SB = 5.670374419e-8  # Stefan-Boltzmann constant [W/m²/K⁴]

include("util.jl")
export angle_rad, angle_deg
export solar_phase_angle, solar_elongation_angle

include("ray_intersection.jl")
export Ray, intersect_ray_triangle, intersect_ray_shape
export BoundingBox, compute_bounding_box, intersect_ray_bounding_box

include("SpiceCamera.jl")
include("SpiceSpacecraft.jl")
include("SpiceAsteroid.jl")
export SpiceCamera
export SpiceSpacecraft, add_instrument!
export SpiceAsteroid
export update!

include("image_generation.jl")
export generate_pixel_rays, generate_thermal_image

end # module FOVSimulator
