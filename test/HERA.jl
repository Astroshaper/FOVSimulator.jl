@testset "HERA SPICE kernels" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |               Test: HERA SPICE kernels                 |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

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
        "g_01165mm_spc_obj_didy_0000n00000_v003.obj",
        "g_00243mm_spc_obj_dimo_0000n00000_v004.obj",
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

    ##==== Test SPICE integration for HERA_TIRI camera ====##
    @testset "HERA_TIRI SPICE integration" begin
        TIRI_ID    = -91200
        fov_angles = (13.3, 10.0)  # TIRI's FOV angles (width, height) in degrees
        img_size   = (1024, 768)   # TIRI's image size (width, height) in pixels
        cam = SpiceCamera("HERA_TIRI", TIRI_ID, fov_angles, img_size)
        
        ## Validate static parameters
        @test cam.static.name == "HERA_TIRI"
        @test cam.static.id == TIRI_ID
        @test cam.static.fov_shape == "RECTANGLE"
        @test cam.static.fov_frame == "HERA_TIRI"
        @test cam.static.boresight == SVector(0.0, 0.0, 1.0)  # Boresight vector at the camera reference frame
        @test length(cam.static.bounds) == 4                    # "RECTANGLE" FOV has 4 boundary vectors
        
        ## Update dynamic state at an epoch
        et = SPICE.utc2et("2025-06-01T00:00:00")
        update!(cam, et, "J2000", "LT+S", "HERA_SPACECRAFT")

        @test norm(cam.state.position) > 0  # Ensure a non-zero position vector is returned
        @test norm(cam.state.velocity) > 0  # Ensure a non-zero velocity vector is returned
        @test !isempty(cam.state.bounds)    # Ensure rotated boundary vectors are populated
    end

    ##==== Unit tests for custom structure: `SpiceAsteroid` ====##
    @testset "SpiceAsteroid basic struct" begin
        path_shape_deimos = joinpath("shape", paths_shape[2])
        shape_deimos = AsteroidShapeModels.load_shape_obj(path_shape_deimos; scale=1000, with_face_visibility=false, with_bvh=true)
        deimos = SpiceAsteroid("DEIMOS", shape_deimos, "DEIMOS_FIXED")

        @test deimos.static.name                 == "DEIMOS"
        @test deimos.static.shape                == shape_deimos
        @test deimos.static.asteroid_fixed_frame == "DEIMOS_FIXED"
        @test deimos.state.position == SVector(0.0, 0.0, 0.0)
        @test deimos.state.velocity == SVector(0.0, 0.0, 0.0)

        @test length(deimos.static.shape.nodes) == 2522
        @test length(deimos.static.shape.faces) == 5040
        
        str = sprint(show, deimos)
        @test occursin("Asteroid parameters", str)
        @test occursin("Fixed frame", str)
        @test occursin("DEIMOS_FIXED", str)
    end

    ##==== Unit tests for custom structure: `SpiceSpacecraft` ====##
    @testset "SpiceSpacecraft basic struct" begin
        hera = SpiceSpacecraft("HERA");
    
        @test hera.static.name == "HERA"
        @test hera.state.position == SVector(0.0, 0.0, 0.0)
        @test hera.state.velocity == SVector(0.0, 0.0, 0.0)
        @test isempty(hera.state.instruments)

        TIRI_ID = -91200
        fov_angles = (45.0, 45.0)  # Horizontal and vertical field of view angles (degrees)
        img_size = (256, 256)      # Image size (width, height in pixels)
        tiri = SpiceCamera("HERA_TIRI", TIRI_ID, fov_angles, img_size)
        add_instrument!(hera, tiri)

        @test length(hera.state.instruments) == 1
        @test haskey(hera.state.instruments, "HERA_TIRI")
        @test hera.state.instruments["HERA_TIRI"].static.name == "HERA_TIRI"
    end

    SPICE.kclear()
end
