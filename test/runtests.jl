using AsteroidShapeModels
using CairoMakie
using Downloads
using FOVSimulator
using LinearAlgebra
using Rotations
using SPICE
using StaticArrays
using Test

import SPICE

include("HERA.jl")
include("ray_intersection_vs_DSK.jl")  # Ray intersection tests with SPICE comparison
include("TIRI_image_Didymos.jl")       # Didymos and Dimorphos thermal simulation tests
