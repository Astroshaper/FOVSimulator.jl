
"""
    intersect_ray_shape_spice(ray::Ray, target_name::String, target_frame::String, et::Float64, abcorr::String, observer::String) -> Tuple{Bool, Float64, SVector{3, Float64}}

SPICEの`sincpt`関数を用いてレイと形状モデルの交差判定を行う。

# 引数
- `ray`          : レイ
- `target_name`  : 対象天体の名前
- `target_frame` : 対象天体の座標系
- `et`           : エフェメリス時刻
- `abcorr`       : 光行差補正フラグ
- `observer`     : 観測者の名前

# 戻り値
- `hit`      : 交差があればtrue、なければfalse
- `distance` : レイの始点から交点までの距離
- `point`    : 交点の座標
"""
function intersect_ray_shape_spice(ray::Ray, target_name::String, target_frame::String, et::Float64, abcorr::String, observer::String)
    try
        # レイの方向ベクトルを単位ベクトルに変換
        direction = normalize(ray.direction)
        
        # SPICEの`sincpt`関数を呼び出し
        spoint, trgepc, srfvec, found = SPICE.sincpt(
            "ELLIPSOID",  # 形状モデルの種類（DSKカーネルが読み込まれている場合は自動的にDSKが使用される）
            target_name,  # 対象天体の名前
            et,           # エフェメリス時刻
            target_frame, # 対象天体の座標系
            abcorr,       # 光行差補正フラグ
            observer,     # 観測者の名前
            "J2000",      # 参照座標系
            ray.direction # レイの方向ベクトル
        )
        
        if found
            # 交点の座標
            point = SVector{3, Float64}(spoint)
            
            # レイの始点から交点までの距離
            distance = norm(point - ray.origin)
            
            return true, distance, point
        else
            return false, 0.0, @SVector zeros(3)
        end
    catch e
        # エラーが発生した場合
        @warn "SPICE sincpt error: $e"
        return false, 0.0, @SVector zeros(3)
    end
end

"""
    compare_intersection_methods(
        shape_model::AsteroidThermoPhysicalModels.ShapeModel,
        target_name::String,
        target_frame::String,
        rays::Vector{Ray},
        et::Float64,
        abcorr::String,
        observer::String,
        tolerance::Float64=1e-3
    ) -> Dict

SPICEの`sincpt`関数と自前実装（Möller–Trumbore法）による交差判定を比較する。

# 引数
- `shape_model`  : 形状モデル
- `target_name`  : 対象天体の名前
- `target_frame` : 対象天体の座標系
- `rays`         : レイのベクトル
- `et`           : エフェメリス時刻
- `abcorr`       : 光行差補正フラグ
- `observer`     : 観測者の名前
- `tolerance`    : 許容誤差（メートル）

# 戻り値
- 結果を格納した辞書
"""
function compare_intersection_methods(
    shape_model::AsteroidThermoPhysicalModels.ShapeModel,
    target_name::String,
    target_frame::String,
    rays::Vector{Ray},
    et::Float64,
    abcorr::String,
    observer::String,
    tolerance::Float64=1e-3
)
    results = Dict(
        "total_rays" => length(rays),
        "mt_hits" => 0,
        "spice_hits" => 0,
        "both_hits" => 0,
        "errors" => Float64[],
        "mt_time" => 0.0,
        "spice_time" => 0.0
    )
    
    for (i, ray) in enumerate(rays)
        # Möller–Trumbore法による交差判定
        mt_start = time()
        result = intersect_ray_shape(ray, shape_model)
        mt_time = time() - mt_start
        results["mt_time"] += mt_time
        
        mt_hit = result.hit
        mt_distance = result.distance
        mt_point = result.point
        mt_face = result.face_index
        
        # SPICEの`sincpt`関数による交差判定
        spice_start = time()
        spice_hit, spice_distance, spice_point = intersect_ray_shape_spice(ray, target_name, target_frame, et, abcorr, observer)
        spice_time = time() - spice_start
        results["spice_time"] += spice_time
        
        # 結果の集計
        if mt_hit
            results["mt_hits"] += 1
        end
        
        if spice_hit
            results["spice_hits"] += 1
        end
        
        # 両方の方法で交点が見つかった場合
        if mt_hit && spice_hit
            results["both_hits"] += 1
            
            # 交点間の距離（誤差）
            error = norm(mt_point - spice_point)
            push!(results["errors"], error)
            
            # 詳細な結果を表示
            println("Ray $i: MT hit = $mt_hit, SPICE hit = $spice_hit, Error = $error m")
            println("  MT point: $mt_point")
            println("  SPICE point: $spice_point")
            println("  MT time: $(mt_time*1000) ms, SPICE time: $(spice_time*1000) ms")
            
            # 誤差が許容範囲を超える場合
            if error > tolerance
                @warn "Ray $i: Error ($error m) exceeds tolerance ($tolerance m)"
            end
        elseif mt_hit != spice_hit
            # 片方の方法でのみ交点が見つかった場合
            println("Ray $i: MT hit = $mt_hit, SPICE hit = $spice_hit (Discrepancy)")
            if mt_hit
                println("  MT point: $mt_point")
            end
            if spice_hit
                println("  SPICE point: $spice_point")
            end
        end
    end
    
    # 平均誤差の計算
    if !isempty(results["errors"])
        results["mean_error"] = mean(results["errors"])
        results["max_error"] = maximum(results["errors"])
    else
        results["mean_error"] = NaN
        results["max_error"] = NaN
    end
    
    return results
end

@testset "Ray Intersection Tests" begin
    # テスト用の形状モデルを読み込む
    shape_file = joinpath("test", "shape", "deimos_k005_tho_v02.obj")
    shape_model = AsteroidThermoPhysicalModels.load_shape_obj(shape_file; scale=1000, find_visible_facets=false)
    
    # SPICEカーネルを読み込む
    spice_kernels = [
        joinpath("test", "kernel", "lsk", "naif0012.tls"),
        joinpath("test", "kernel", "pck", "pck00011.tpc"),
        joinpath("test", "kernel", "dsk", "deimos_k005_tho_v02.bds")
    ]
    
    for kernel in spice_kernels
        if isfile(kernel)
            SPICE.furnsh(kernel)
        else
            @warn "SPICE kernel not found: $kernel"
        end
    end
    
    # テスト用のパラメータ
    target_name = "DEIMOS"
    target_frame = "DEIMOS_FIXED"
    et = 0.0  # エフェメリス時刻（J2000）
    abcorr = "NONE"  # 光行差補正なし
    observer = "MARS"
    
    # テスト用のレイを生成
    camera_position = @SVector [10000.0, 0.0, 0.0]  # カメラの位置（デイモスから10km離れた位置）
    target_position = @SVector [0.0, 0.0, 0.0]      # デイモスの中心
    n_rays = 10                                      # レイの数
    spread_angle = 5.0                               # レイの広がり角（度）
    
    rays = generate_test_rays(camera_position, target_position, n_rays, spread_angle)
    
    # 交差判定の比較
    results = compare_intersection_methods(
        shape_model,
        target_name,
        target_frame,
        rays,
        et,
        abcorr,
        observer,
        1e-3  # 許容誤差（メートル）
    )
    
    # 結果の表示
    println("===== Ray Intersection Test Results =====")
    println("Total rays: $(results["total_rays"])")
    println("MT hits: $(results["mt_hits"])")
    println("SPICE hits: $(results["spice_hits"])")
    println("Both hits: $(results["both_hits"])")
    println("Mean error: $(results["mean_error"]) m")
    println("Max error: $(results["max_error"]) m")
    println("MT time: $(results["mt_time"]*1000) ms")
    println("SPICE time: $(results["spice_time"]*1000) ms")
    println("========================================")
    
    # テスト
    @test results["mt_hits"] > 0
    @test results["spice_hits"] > 0
    @test results["both_hits"] > 0
    @test !isnan(results["mean_error"])
    @test results["mean_error"] < 1e-3  # 平均誤差が1mm未満
    
    # SPICEカーネルをアンロード
    SPICE.kclear()
end
