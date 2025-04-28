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


@testset "Projection functions" begin
    @test FOVSimulator.focal_length(90, 100) ≈ 50

    p = SVector(1.0, 0.0, 1.0)
    fov = (90.0, 90.0)
    img = (100, 100)
    @test FOVSimulator.project_point_fov(p, fov, img) == (100, 50)
end


include("HERA/HERA.jl")
