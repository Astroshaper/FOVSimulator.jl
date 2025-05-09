"""
    TIRI_image_Didymos.jl

DidymosとDimorphosの熱物理シミュレーションを行い、TIRIカメラによる熱画像を生成するテストコード。
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
    
    # 熱物理モデルのパラメータ
    P₁ = SPICE.convrt(2.2593, "hours", "seconds")  # Didymosの自転周期 [秒]
    P₂ = SPICE.convrt(11.93, "hours", "seconds")  # Dimorphosの自転周期 [秒]
    
    # TPMの計算サイクル数と時間ステップ
    n_cycle = 10  # 計算するサイクル数
    n_step_in_cycle = 180  # 1回転あたりの時間ステップ数
    
    # TPMの計算開始時刻と終了時刻
    et_begin = SPICE.utc2et("2027-02-01T00:00:00")  # TPM計算開始時刻
    et_end   = et_begin + P₂ * n_cycle  # TPM計算終了時刻
    et_range = range(et_begin, et_end; length=n_step_in_cycle*n_cycle+1)
    
    # エフェメリスデータの準備
    ephem = (
        time = collect(et_range),
        sun1 = [SVector{3}(SPICE.spkpos("SUN"      , et, "DIDYMOS_FIXED"  , "None", "DIDYMOS"  )[1]) * 1000 for et in et_range],
        sun2 = [SVector{3}(SPICE.spkpos("SUN"      , et, "DIMORPHOS_FIXED", "None", "DIMORPHOS")[1]) * 1000 for et in et_range],
        sec  = [SVector{3}(SPICE.spkpos("DIMORPHOS", et, "DIDYMOS_FIXED"  , "None", "DIDYMOS"  )[1]) * 1000 for et in et_range],
        P2S  = [RotMatrix{3}(SPICE.pxform("DIDYMOS_FIXED"  , "DIMORPHOS_FIXED", et)) for et in et_range],
        S2P  = [RotMatrix{3}(SPICE.pxform("DIMORPHOS_FIXED", "DIDYMOS_FIXED"  , et)) for et in et_range],
    )
    
    # 形状モデルの読み込み
    path_shape1 = joinpath("shape", "g_50677mm_rad_obj_didy_0000n00000_v001.obj")
    path_shape2 = joinpath("shape", "g_01332mm_lgt_obj_dimo_000n00000_v001.obj")

    shape1 = AsteroidThermoPhysicalModels.load_shape_obj(path_shape1; scale=1000, find_visible_facets=false)
    shape2 = AsteroidThermoPhysicalModels.load_shape_obj(path_shape2; scale=1000, find_visible_facets=false)
    
    n_face1 = length(shape1.faces)  # Didymosの面の数
    n_face2 = length(shape2.faces)  # Dimorphosの面の数

    # 熱物理パラメータ
    k  = 0.125   # 熱伝導率 [W/m/K]
    ρ  = 2170.0  # 密度 [kg/m³]
    Cₚ = 600.0   # 熱容量 [J/kg/K]
    
    # 熱浸透深さの計算
    l₁ = AsteroidThermoPhysicalModels.thermal_skin_depth(P₁, k, ρ, Cₚ)  # Didymosの熱浸透深さ
    l₂ = AsteroidThermoPhysicalModels.thermal_skin_depth(P₂, k, ρ, Cₚ)  # Dimorphosの熱浸透深さ
    Γ = AsteroidThermoPhysicalModels.thermal_inertia(k, ρ, Cₚ)          # 熱慣性 [tiu]
    
    # 光学特性
    R_vis = 0.059  # 可視光の反射率 [-]
    R_ir  = 0.0    # 熱赤外の反射率 [-]
    ε     = 0.9    # 放射率 [-]
    
    # 熱伝導方程式の計算パラメータ
    z_max = 0.6   # 熱伝導方程式の下端深さ [m]
    n_depth = 41  # 深さ方向の分割数
    Δz = z_max / (n_depth - 1)  # 深さ方向のステップ幅 [m]
    
    # # 熱物理パラメータの設定
    # thermo_params1 = AsteroidThermoPhysicalModels.ThermoParams(P₁, l₁, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    # thermo_params2 = AsteroidThermoPhysicalModels.ThermoParams(P₂, l₂, Γ, R_vis, R_ir, ε, z_max, Δz, n_depth)
    
    # # TPMの設定
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
    
    # # TPMの実行
    # times_to_save = ephem.time[end-n_step_in_cycle:end]  # 最終回転周期の温度を保存
    # face_ID_pri = [1, 2, 3, 4, 10]  # 主小惑星の表面温度を保存する面のインデックス
    # face_ID_sec = [1, 2, 3, 4, 20]  # 副小惑星の表面温度を保存する面のインデックス
    
    # @info "Running thermophysical model..."
    # result = AsteroidThermoPhysicalModels.run_TPM!(btpm, ephem, times_to_save, face_ID_pri, face_ID_sec; show_progress=false)
    # @info "Thermophysical model completed."
    
    # # TIRI画像シミュレーション
    # TIRI_ID    = -91200        # TIRIのNAIF ID
    # fov_angles = (13.3, 10.0)  # TIRIの視野角（幅、高さ）[度]
    # img_size   = (1024, 768)   # TIRIの画像サイズ（幅、高さ）[ピクセル]
    
    # # スペースクラフトとカメラの設定
    # hera = SpiceSpacecraft("HERA")
    # tiri = SpiceCamera("HERA_TIRI", TIRI_ID, fov_angles, img_size)
    # add_instrument!(hera, tiri)
    
    # # 小惑星の設定
    # didymos = SpiceAsteroid("DIDYMOS", shape1, "DIDYMOS_FIXED")
    # dimorphos = SpiceAsteroid("DIMORPHOS", shape2, "DIMORPHOS_FIXED")
    
    # # シミュレーション時刻の設定
    # et_idx = 1  # 使用するエフェメリス時刻のインデックス
    # et     = result.pri.times_to_save[et_idx]
    # ref    = "HERA_TIRI"
    # abcorr = "NONE"
    # obs    = "HERA_TIRI"
    
    # # スペースクラフト、カメラ、小惑星の状態を更新
    # update!(hera     , et, ref, abcorr, obs)
    # update!(didymos  , et, ref, abcorr, obs)
    # update!(dimorphos, et, ref, abcorr, obs)
    
    # # 表面温度の取得
    # temperatures1 = result.pri.surface_temperature[:, et_idx]
    # temperatures2 = result.sec.surface_temperature[:, et_idx]
    
    # # 交差マップの生成
    # @info "Generating intersection maps..."
    # intersection_map1 = generate_intersection_map(hera.state.instruments["HERA_TIRI"], didymos)
    # intersection_map2 = generate_intersection_map(hera.state.instruments["HERA_TIRI"], dimorphos)
    
    # # 熱画像の生成
    # @info "Generating thermal image..."
    # img = generate_image_temperature(
    #     intersection_map1, temperatures1,
    #     intersection_map2, temperatures2
    # )
    
    # # 画像の表示関数
    # function show_image(img::Array{Float64,2})
    #     fig = Figure()
    #     ax = Axis(fig[1, 1]; aspect = DataAspect())  # 縦横比1:1の軸を作成
        
    #     hm = image!(ax, img; colormap = :inferno)   # 画像を表示
    #     Colorbar(fig[1, 2], hm, label = "Temperature [K]")  # カラーバーを追加
        
    #     return fig
    # end
    
    # # 画像の表示と保存
    # fig = show_image(img)
    # save("didymos_dimorphos_thermal_image.png", fig)
    
    # # テスト
    # @test size(img) == (img_size[2], img_size[1])  # 画像サイズが正しいか確認
    # @test maximum(img) > 0.0  # 画像に温度データが含まれているか確認
    
    # SPICEカーネルのクリア
    SPICE.kclear()
end
