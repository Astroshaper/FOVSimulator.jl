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

include("util.jl")
export angle_rad, angle_deg
export solar_phase_angle, solar_elongation_angle

include("fov_projection.jl")
export ProjectedFace, project_face_centers, map_temperature_to_image

end # module FOVSimulator
