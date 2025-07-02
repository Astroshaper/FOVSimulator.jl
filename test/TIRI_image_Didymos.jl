"""
    TIRI_image_Didymos.jl

Test code for thermal physics simulation of Didymos and Dimorphos,
and generation of thermal images with TIRI camera.
"""

@testset "Didymos & Dimorphos Thermal Image Simulation" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |   Test: Didymos & Dimorphos Thermal Image Simulation   |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    ##==== Download Files ====##

    # List of SPICE kernels from hera_crema_2_1.tm
    # https://spiftp.esac.esa.int/data/SPICE/HERA/kernels/mk/hera_crema_2_1.tm
    paths_kernel = [
        "ck/hera_sc_crema_2_1_LPO_241007_270303_f181203_v01.bc",
        # "ck/hera_juventas_sc_cruise_v02.bc",
        # "ck/hera_milani_sc_cruise_v02.bc",

        "fk/hera_v14.tf",
        "fk/hera_ops_v05.tf",
        "fk/hera_dsk_surfaces_v05.tf",
        "fk/estrack_v04.tf",
        # "fk/hera_milani_v05.tf",
        # "fk/hera_juventas_v06.tf",

        # "dsk/deimos_k005_tho_v02.bds",
        # "dsk/phobos_m003_gas_v01.bds",
        # "dsk/g_01165mm_spc_obj_didy_0000n00000_v003.bds",
        # "dsk/g_00243mm_spc_obj_dimo_0000n00000_v004.bds",

        # "ik/hera_afc_v05.ti",
        # "ik/hera_hsh_v02.ti",
        # "ik/hera_palt_v02.ti",
        "ik/hera_tiri_v03.ti",
        # "ik/hera_smc_v01.ti",
        # "ik/hera_str_v01.ti",
        # "ik/hera_juventas_jura_v00.ti",
        # "ik/hera_juventas_navcam_v00.ti",
        # "ik/hera_juventas_ccam_v00.ti",
        # "ik/hera_juventas_lidar_v00.ti",
        # "ik/hera_juventas_str_v01.ti",
        # "ik/hera_milani_aspect_v01.ti",
        # "ik/hera_milani_vista_v01.ti",
        # "ik/hera_milani_navcam_v01.ti",
        # "ik/hera_milani_lidar_v01.ti",
        # "ik/hera_milani_mlrh_v02.ti",
        # "ik/hera_milani_str_v00.ti",

        "lsk/naif0012.tls",

        "pck/pck00011.tpc",
        "pck/de-403-masses.tpc",
        "pck/hera_didymos_v06.tpc",

        "sclk/hera_fict_181203_v01.tsc",
        # "sclk/hera_juventas_fict_250123_v01.tsc",
        # "sclk/hera_milani_fict_250123_v01.tsc",

        "spk/de432s.bsp",
        "spk/estrack_v04.bsp",
        "spk/mar097_20160314_20300101.bsp",
        "spk/didymos_hor_000101_500101_v01.bsp",
        "spk/didymos_crema_2_1_ECP_PDP_DCP_261120_270515_v01.bsp",
        "spk/hera_dart_impact_site_v04.bsp",
        "spk/hera_sci_v01.bsp",
        "spk/hera_struct_v02.bsp",
        "spk/hera_sc_crema_2_1_LPO_241007_261202_v01.bsp",
        "spk/hera_sc_crema_2_1_ECP_PDP_DCP_261125_270303_v01.bsp",
        # "spk/hera_juventas_struct_v00.bsp",
        # "spk/hera_juventas_cruise_v03.bsp",
        # "spk/hera_milani_cruise_v03.bsp",
    ]

    paths_shape = [
        "phobos_m003_gas_v01.obj",
        "deimos_k005_tho_v02.obj",
        # "g_01165mm_spc_obj_didy_0000n00000_v003.obj",
        # "g_00243mm_spc_obj_dimo_0000n00000_v004.obj",
        "g_50677mm_rad_obj_didy_0000n00000_v001.obj",
        "g_01332mm_lgt_obj_dimo_000n00000_v001.obj",
    ]

    for path_kernel in paths_kernel
        url_kernel = "https://spiftp.esac.esa.int/data/SPICE/HERA/kernels/$(path_kernel)"
        filepath = joinpath("kernel", path_kernel)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_kernel, filepath)
    end

    for path_shape in paths_shape
        url_shape = "https://spiftp.esac.esa.int/data/SPICE/HERA/kernels/dsk/$(path_shape)"
        filepath = joinpath("shape", path_shape)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_shape, filepath)
    end

    ##==== Load data with SPICE ====##
    for path_kernel in paths_kernel
        filepath = joinpath("kernel", path_kernel)
        SPICE.furnsh(filepath)
    end
    
    # Thermophysical model parameters
    P₁ = SPICE.convrt(2.2593, "hours", "seconds")  # Didymos rotation period [seconds]
    P₂ = SPICE.convrt(11.93, "hours", "seconds")  # Dimorphos rotation period [seconds]
    
    # TPM calculation cycles and time steps
    n_cycle = 10           # Number of cycles to calculate
    n_step_in_cycle = 180  # Number of time steps per rotation
    
    # TPM calculation start and end times
    et_begin = SPICE.utc2et("2027-02-01T00:00:00")  # TPM calculation start time
    et_end   = et_begin + P₂ * n_cycle  # TPM calculation end time
    et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycle+1)
    
    # Prepare ephemeris data
    ephem = (
        time = collect(et_range),
        sun1 = [SVector{3}(SPICE.spkpos("SUN"      , et, "DIDYMOS_FIXED"  , "None", "DIDYMOS"  )[1]) * 1000 for et in et_range],
        sun2 = [SVector{3}(SPICE.spkpos("SUN"      , et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1]) * 1000 for et in et_range],
        sec  = [SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED"  , "None", "DIDYMOS"  )[1]) * 1000 for et in et_range],
        P2S  = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED"  , "DIMORPHOS_FIXED", et)) for et in et_range],
        S2P  = [RotMatrix{3}(SPICE.pxform("DIMORPHOS_FIXED", "DIDYMOS_FIXED"  , et)) for et in et_range],
    )
    
    # Load shape models
    path_shape1 = joinpath("shape", "g_50677mm_rad_obj_didy_0000n00000_v001.obj")
    path_shape2 = joinpath("shape", "g_01332mm_lgt_obj_dimo_000n00000_v001.obj")

    shape1 = AsteroidShapeModels.load_shape_obj(path_shape1; scale=1000, with_face_visibility=false)
    shape2 = AsteroidShapeModels.load_shape_obj(path_shape2; scale=1000, with_face_visibility=false)
    
    n_face1 = length(shape1.faces)  # Number of faces in Didymos
    n_face2 = length(shape2.faces)  # Number of faces in Dimorphos

    # Thermophysical parameters
    k  = 0.125   # Thermal conductivity [W/m/K]
    ρ  = 2170.0  # Density [kg/m³]
    Cₚ = 600.0   # Heat capacity [J/kg/K]
    
    # Thermal skin depth calculation (temporarily commented out due to dependency on AsteroidThermoPhysicalModels)
    # l₁ = thermal_skin_depth(P₁, k, ρ, Cₚ)  # Thermal skin depth of Didymos
    # l₂ = thermal_skin_depth(P₂, k, ρ, Cₚ)  # Thermal skin depth of Dimorphos
    # Γ = thermal_inertia(k, ρ, Cₚ)          # Thermal inertia [tiu]
    
    # Optical properties
    R_vis = 0.059  # Visible reflectance [-]
    R_ir  = 0.0    # Thermal infrared reflectance [-]
    ε     = 0.9    # Emissivity [-]
    
    # Heat conduction equation calculation parameters
    z_max = 0.6                 # Bottom depth of heat conduction equation [m]
    n_depth = 41                # Number of divisions in depth direction
    Δz = z_max / (n_depth - 1)  # Step size in depth direction [m]
    
    # # Thermophysical parameter settings
    # thermo_params1 = AsteroidThermoPhysicalModels.ThermoParams(P₁, l₁, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    # thermo_params2 = AsteroidThermoPhysicalModels.ThermoParams(P₂, l₂, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    
    # # TPM settings
    # stpm1 = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape1, thermo_params1;
    #     SELF_SHADOWING = true,
    #     SELF_HEATING   = true,
    #     SOLVER         = AsteroidThermoPhysicalModels.ExplicitEulerSolver(thermo_params1),
    #     BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
    #     BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    # )
    
    # stpm2 = AsteroidThermoPhysicalModels.SingleAsteroidTPM(shape2, thermo_params2;
    #     SELF_SHADOWING = true,
    #     SELF_HEATING   = true,
    #     SOLVER         = AsteroidThermoPhysicalModels.ExplicitEulerSolver(thermo_params2),
    #     BC_UPPER       = AsteroidThermoPhysicalModels.RadiationBoundaryCondition(),
    #     BC_LOWER       = AsteroidThermoPhysicalModels.InsulationBoundaryCondition(),
    # )
    
    # btpm = AsteroidThermoPhysicalModels.BinaryAsteroidTPM(stpm1, stpm2; MUTUAL_SHADOWING=true, MUTUAL_HEATING=false)
    # AsteroidThermoPhysicalModels.init_temperature!(btpm, 200.)
    
    # # Run TPM
    # times_to_save = ephem.time[end-n_step_in_cycle:end]  # Save temperatures for the final rotation period
    # face_ID_pri = [1, 2, 3, 4, 10]  # Face indices to save surface temperatures for primary asteroid
    # face_ID_sec = [1, 2, 3, 4, 20]  # Face indices to save surface temperatures for secondary asteroid
    
    # @info "Running thermophysical model..."
    # result = AsteroidThermoPhysicalModels.run_TPM!(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec; show_progress=false)
    # @info "Thermophysical model completed."
    
    # # TIRI image simulation
    # TIRI_ID    = -91200        # TIRI NAIF ID
    # fov_angles = (13.3, 10.0)  # TIRI field of view angles (width, height) [degrees]
    # img_size   = (1024, 768)   # TIRI image size (width, height) [pixels]
    
    # # Spacecraft and camera setup
    # hera = SpiceSpacecraft("HERA")
    # tiri = SpiceCamera("HERA_TIRI", TIRI_ID, fov_angles, img_size)
    # add_instrument!(hera, tiri)
    
    # # Asteroid setup
    # didymos = SpiceAsteroid("DIDYMOS", shape1, "DIDYMOS_FIXED")
    # dimorphos = SpiceAsteroid("DIMORPHOS", shape2, "DIMORPHOS_FIXED")
    
    # # Simulation time setup
    # et_idx = 1  # Index of ephemeris time to use
    # et     = result.pri.times_to_save[et_idx]
    # ref    = "HERA_TIRI"
    # abcorr = "NONE"
    # obs    = "HERA_TIRI"
    
    # # Update spacecraft, camera, and asteroid states
    # update!(hera     , et, ref, abcorr, obs)
    # update!(didymos  , et, ref, abcorr, obs)
    # update!(dimorphos, et, ref, abcorr, obs)
    
    # # Get surface temperatures
    # temperatures1 = result.pri.surface_temperature[:, et_idx]
    # temperatures2 = result.sec.surface_temperature[:, et_idx]
    
    # # Generate intersection maps
    # @info "Generating intersection maps..."
    # intersection_map1 = generate_intersection_map(hera.state.instruments["HERA_TIRI"], didymos)
    # intersection_map2 = generate_intersection_map(hera.state.instruments["HERA_TIRI"], dimorphos)
    
    # # Generate thermal image
    # @info "Generating thermal image..."
    # img = generate_image_temperature(
    #     intersection_map1, temperatures1,
    #     intersection_map2, temperatures2
    # )
    
    # # Image display function
    # function show_image(img::Array{Float64,2})
    #     fig = Figure()
    #     ax = Axis(fig[1, 1]; aspect = DataAspect())  # Create axis with 1:1 aspect ratio
        
    #     hm = image!(ax, img; colormap = :inferno)   # Display image
    #     Colorbar(fig[1, 2], hm, label = "Temperature [K]")  # Add colorbar
        
    #     return fig
    # end
    
    # # Display and save image
    # fig = show_image(img)
    # save("didymos_dimorphos_thermal_image.png", fig)
    
    # # Tests
    # @test size(img) == (img_size[2], img_size[1])  # Confirm correct image size
    # @test maximum(img) > 0.0  # Confirm that image contains temperature data
    
    # Clear SPICE kernels
    SPICE.kclear()
end
