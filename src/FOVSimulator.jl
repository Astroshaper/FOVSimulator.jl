module FOVSimulator

using AsteroidShapeModels
using LinearAlgebra
using Rotations
using StaticArrays

import SPICE

# Physical constants
const σ_SB = 5.670374419e-8  # Stefan-Boltzmann constant [W/m²/K⁴]

include("SpiceCamera.jl")
include("SpiceSpacecraft.jl")
include("SpiceAsteroid.jl")
export SpiceCamera
export SpiceSpacecraft, add_instrument!
export SpiceAsteroid
export update!

include("image_generation.jl")
export generate_intersection_map, generate_image_temperature, generate_image_radiance

end # module FOVSimulator
