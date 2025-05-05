using AsteroidThermoPhysicalModels
using FOVSimulator
using Downloads
using LinearAlgebra
using SPICE
using StaticArrays
using Test

@testset "Utility functions" begin
    @test angle_deg([1, 0, 0], [ 0, 1, 0]) ≈ 90
    @test angle_deg([1, 0, 0], [ 1, 0, 0]) ≈ 0
    @test angle_deg([1, 0, 0], [-1, 0, 0]) ≈ 180

    sun = [1, 0, 0]
    tgt = [0, 0, 0]
    obs = [0, 1, 0]

    @test solar_phase_angle(sun, tgt, obs) ≈ deg2rad(90)
    @test solar_elongation_angle(sun, obs, tgt) ≈ deg2rad(45)
end

include("HERA.jl")

include("ray_intersection_simple_cases.jl")  # Simple ray intersection tests
include("ray_intersection_vs_DSK.jl")        # Ray intersection tests with SPICE comparison
