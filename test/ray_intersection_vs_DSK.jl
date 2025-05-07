
@testset "Ray-Shape Intersection vs. SPICE/DSK" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |       Test: Ray-Shape Intersection vs. SPICE/DSK       |
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
        "dsk/g_01165mm_spc_obj_didy_0000n00000_v003.bds",
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
        # "phobos_m003_gas_v01.obj",
        # "deimos_k005_tho_v02.obj",
        "g_01165mm_spc_obj_didy_0000n00000_v003.obj",
        # "g_00243mm_spc_obj_dimo_0000n00000_v004.obj",
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

    #### FOVSimulator.jlで実装したレイと形状モデルの交差判定 ####

    obj_path = joinpath("shape", "g_01165mm_spc_obj_didy_0000n00000_v003.obj")
    shape = AsteroidThermoPhysicalModels.load_shape_obj(obj_path; scale=1000, find_visible_facets=false)
    println(shape)

    TIRI_ID    = -91200
    fov_angles = (13.3, 10.0)  # TIRI's FOV angles (width, height) in degrees
    img_size   = (1024, 768)   # TIRI's image size (width, height) in pixels
    cam = SpiceCamera("HERA_TIRI", TIRI_ID, fov_angles, img_size)

    utc = "2027-02-01T01:00:00"
    et = SPICE.utc2et(utc)
    
    ref    = "DIDYMOS_FIXED"
    abcorr = "NONE"  # 光行差補正なし
    # abcorr = "LT+S"  # 光行差補正あり（交差点が0.0012 m程度変わる）
    obs    = "DIDYMOS"

    update!(cam, et, ref, abcorr, obs)
    ray = Ray(cam.state.position, cam.state.boresight)
    
    # バウンディングボックスを計算
    bbox = compute_bounding_box(shape)
    
    intersection = intersect_ray_shape(ray, shape, bbox)  # 交差判定の結果
    
    # @show ray.origin
    # @show ray.direction
    # @show intersection.distance
    # @show intersection.point

    #### SPICEのsincpt関数を用いた交差判定 ####
    ## cf. https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/sincpt_c.html

    spoint, trgepc, srfvec = SPICE.sincpt(
        "DSK/UNPRIORITIZED",          # 形状モデルの種類（DSKカーネルが読み込まれている場合は自動的にDSKが使用される）
        obs,                          # 対象天体の名前
        et,                           # エフェメリス時刻
        ref,                          # 対象天体の座標系
        abcorr,                       # 光行差補正フラグ
        cam.static.name,              # 観測者の名前
        cam.static.fov_frame,         # 参照座標系
        collect(cam.static.boresight) # レイの方向ベクトル。ccallのエラーを防ぐため、Vector{Float64}型に変換してsincptに渡す
    )

    spoint *= 1000  # 単位をkmからmに変換
    srfvec *= 1000  # 単位をkmからmに変換

    # @show spoint  # Surface intercept point on the target body
    # @show trgepc  # Intercept epoch
    # @show srfvec  # Vector from observer to intercept point

    #### 交差判定結果の比較 ####

    diff = norm(spoint - intersection.point)  # 交差判定結果の差 [m]
    @test diff < 0.01  # 交差判定結果の差が許容誤差未満であることを確認

    println("Intersection point [m]")
    println("    ∘ FOVSimulator.jl : $(intersection.point)")
    println("    ∘ SPICE/DSK       : $spoint")
    println("    → Difference between them : $diff m")
    println()
    
    #### 実行時間の比較 ####
    
    println("Computation time")
    print("    ∘ intersect_ray_shape in FOVSimulator.jl :")
    @time intersect_ray_shape(ray, shape, bbox)
    print("    ∘ sincpt in SPICE.jl                     :")
    @time SPICE.sincpt("DSK/UNPRIORITIZED", obs, et, ref, abcorr, "HERA", "HERA_TIRI", collect(cam.static.boresight))
    println()

    SPICE.kclear()  # SPICEカーネルをアンロード
end
