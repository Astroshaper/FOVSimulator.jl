
@testset "Ray-Triangle Intersection" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |            Test: Ray-Triangle Intersection             |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)

    # 単一の三角形を定義
    # この三角形はxy平面上にある
    v1 = @SVector [0.0, 0.0, 0.0]  # 原点
    v2 = @SVector [1.0, 0.0, 0.0]  # x軸上の点
    v3 = @SVector [0.0, 1.0, 0.0]  # y軸上の点
    
    # テストケース1: 三角形に交差するレイ
    ray1 = Ray(@SVector([0.25, 0.25, 1.0]), @SVector([0.0, 0.0, -1.0]))
    result1 = intersect_ray_triangle(ray1, v1, v2, v3)
    
    @test result1.hit == true
    @test result1.distance ≈ 1.0
    @test result1.point ≈ @SVector([0.25, 0.25, 0.0])
    
    # テストケース2: 三角形に交差しないレイ（三角形の外側を通過）
    ray2 = Ray(@SVector([2.0, 2.0, 1.0]), @SVector([0.0, 0.0, -1.0]))
    result2 = intersect_ray_triangle(ray2, v1, v2, v3)
    
    @test result2.hit == false
    
    # テストケース3: 三角形に交差しないレイ（三角形と平行）
    ray3 = Ray(@SVector([0.5, 0.5, 1.0]), @SVector([1.0, 0.0, 0.0]))
    result3 = intersect_ray_triangle(ray3, v1, v2, v3)
    
    @test result3.hit == false
    
    # テストケース4: 三角形の頂点を通過するレイ
    ray4 = Ray(@SVector([0.0, 0.0, 1.0]), @SVector([0.0, 0.0, -1.0]))
    result4 = intersect_ray_triangle(ray4, v1, v2, v3)
    
    @test result4.hit == true
    @test result4.distance ≈ 1.0
    @test result4.point ≈ @SVector([0.0, 0.0, 0.0])
    
    # テストケース5: 三角形の辺上を通過するレイ
    ray5 = Ray(@SVector([0.5, 0.0, 1.0]), @SVector([0.0, 0.0, -1.0]))
    result5 = intersect_ray_triangle(ray5, v1, v2, v3)
    
    @test result5.hit == true
    @test result5.distance ≈ 1.0
    @test result5.point ≈ @SVector([0.5, 0.0, 0.0])
    
    # テストケース6: 三角形の裏側からのレイ（バックフェイスカリング）
    ray6 = Ray(@SVector([0.25, 0.25, -1.0]), @SVector([0.0, 0.0, 1.0]))
    result6 = intersect_ray_triangle(ray6, v1, v2, v3)
    
    @test result6.hit == true  # バックフェイスカリングを行わない場合
    @test result6.distance ≈ 1.0
    @test result6.point ≈ @SVector([0.25, 0.25, 0.0])
    
    # テストケース7: レイの始点が三角形上にある場合
    ray7 = Ray(@SVector([0.25, 0.25, 0.0]), @SVector([0.0, 0.0, -1.0]))
    result7 = intersect_ray_triangle(ray7, v1, v2, v3)
    
    @test result7.hit == false  # 数値誤差により、始点が三角形上にある場合は交差しないと判定される
    
    # テストケース8: レイの始点が三角形の裏側にある場合
    ray8 = Ray(@SVector([0.25, 0.25, -1.0]), @SVector([0.0, 0.0, -1.0]))
    result8 = intersect_ray_triangle(ray8, v1, v2, v3)
    
    @test result8.hit == false  # レイが三角形から遠ざかる方向に進む場合は交差しない
end


@testset "Ray-Shape Intersection" begin
    # 単一の三角形を定義
    # この三角形はxy平面上にある
    v1 = @SVector [0.0, 0.0, 0.0]  # 原点
    v2 = @SVector [1.0, 0.0, 0.0]  # x軸上の点
    v3 = @SVector [0.0, 1.0, 0.0]  # y軸上の点

    nodes = [v1, v2, v3]
    faces = [@SVector([1, 2, 3])]  # 頂点インデックスは1から始まる
    
    face_centers  = [AsteroidThermoPhysicalModels.face_center(nodes[face]) for face in faces]
    face_normals  = [AsteroidThermoPhysicalModels.face_normal(nodes[face]) for face in faces]
    face_areas    = [AsteroidThermoPhysicalModels.face_area(nodes[face])   for face in faces]
    visiblefacets = [AsteroidThermoPhysicalModels.VisibleFacet[] for _ in faces]

    shape = AsteroidThermoPhysicalModels.ShapeModel(nodes, faces, face_centers, face_normals, face_areas, visiblefacets)
    
    # レイを定義
    ray = Ray(@SVector([0.25, 0.25, 1.0]), @SVector([0.0, 0.0, -1.0]))
    
    # バウンディングボックスを計算
    bbox = compute_bounding_box(shape)
    
    # 交差判定
    result = intersect_ray_shape(ray, shape, bbox)
    
    # テスト
    @test result.hit == true
    @test result.distance ≈ 1.0
    @test result.point ≈ @SVector([0.25, 0.25, 0.0])
    @test result.face_index == 1
end
