"""
    ray_intersection.jl

小惑星の形状モデルに対するレイと表面の交差判定処理を実装する。
Möller–Trumbore法を用いたレイと三角形メッシュの交差判定を提供する。
"""


"""
    Ray

レイを表す構造体。

# フィールド
- `origin`    : レイの始点
- `direction` : レイの方向ベクトル（正規化済み）
"""
struct Ray
    origin::SVector{3, Float64}
    direction::SVector{3, Float64}
    
    function Ray(origin::AbstractVector{<:Real}, direction::AbstractVector{<:Real})
        # 方向ベクトルを正規化
        dir_normalized = normalize(direction)
        return new(SVector{3, Float64}(origin), SVector{3, Float64}(dir_normalized))
    end
end

"""
    Base.show(io::IO, ray::Ray)

Rayオブジェクトを表示するための関数。

# 引数
- `io`  : 出力ストリーム
- `ray` : 表示するRayオブジェクト
"""
function Base.show(io::IO, ray::Ray)
    msg = """
    Ray:
        ∘ origin    = $(ray.origin)
        ∘ direction = $(ray.direction)
    """
    print(io, msg)
end

"""
    RayTriangleIntersectionResult

レイと三角形の交差判定結果を表す構造体。

# フィールド
- `hit`      : 交差があればtrue、なければfalse
- `distance` : レイの始点から交点までの距離
- `point`    : 交点の座標
"""
struct RayTriangleIntersectionResult
    hit::Bool
    distance::Float64
    point::SVector{3, Float64}
end

# レイと三角形の交差なしの結果を表す定数
const NO_INTERSECTION_RAY_TRIANGLE = RayTriangleIntersectionResult(false, NaN, SVector(NaN, NaN, NaN))

"""
    Base.show(io::IO, result::RayTriangleIntersectionResult

RayTriangleIntersectionResultオブジェクトを表示するための関数。

# 引数
- `io`     : 出力ストリーム
- `result` : 表示するRayTriangleIntersectionResultオブジェクト
"""
function Base.show(io::IO, result::RayTriangleIntersectionResult)
    if result.hit
        msg = """
        Ray-Triangle Intersection:
            ∘ hit      = $(result.hit)
            ∘ distance = $(result.distance)
            ∘ point    = $(result.point)
        """
        print(io, msg)
    else
        msg = """
        Ray-Triangle Intersection:
            ∘ hit = $(result.hit)
        """
        print(io, msg)
    end
end

"""
    intersect_ray_triangle(ray::Ray, v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real}, v3::AbstractVector{<:Real}) -> RayTriangleIntersectionResult

Möller–Trumbore法を用いてレイと三角形の交差判定を行う。

# 引数
- `ray` : レイ
- `v1`  : 三角形の頂点1
- `v2`  : 三角形の頂点2
- `v3`  : 三角形の頂点3

# 戻り値
- 交差判定結果を格納した`RayTriangleIntersectionResult`オブジェクト
"""
function intersect_ray_triangle(ray::Ray, v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real}, v3::AbstractVector{<:Real})
    # エッジベクトル
    e1 = v2 - v1
    e2 = v3 - v1
    
    # レイの方向ベクトルとe2の外積
    p = cross(ray.direction, e2)
    
    # 行列式
    det = dot(e1, p)
    
    # 行列式がほぼ0の場合、レイは三角形と平行
    if abs(det) < 1e-8
        return NO_INTERSECTION_RAY_TRIANGLE
    end
    
    inv_det = 1.0 / det
    
    # レイの始点からv1へのベクトル
    t = ray.origin - v1
    
    # u座標の計算
    u = dot(t, p) * inv_det
    
    # u座標が三角形の範囲外
    if u < 0.0 || u > 1.0
        return NO_INTERSECTION_RAY_TRIANGLE
    end
    
    # tとe1の外積
    q = cross(t, e1)
    
    # v座標の計算
    v = dot(ray.direction, q) * inv_det
    
    # v座標が三角形の範囲外、またはu+vが1を超える
    if v < 0.0 || u + v > 1.0
        return NO_INTERSECTION_RAY_TRIANGLE
    end
    
    # 交点までの距離
    distance = dot(e2, q) * inv_det
    
    # 交点がレイの正の方向にある場合
    if distance > 0.0
        # 交点の座標
        point = ray.origin + distance * ray.direction
        return RayTriangleIntersectionResult(true, distance, point)
    end
    
    # 交点がレイの負の方向にある場合（レイの背後）
    return NO_INTERSECTION_RAY_TRIANGLE
end

"""
    RayShapeIntersectionResult

レイと形状モデルの交差判定結果を表す構造体。

# フィールド
- `hit`        : 交差があればtrue、なければfalse
- `distance`   : レイの始点から交点までの距離
- `point`      : 交点の座標
- `face_index` : 交差した面のインデックス
"""
struct RayShapeIntersectionResult
    hit::Bool
    distance::Float64
    point::SVector{3, Float64}
    face_index::Int
end

# レイと形状モデルの交差なしの結果を表す定数
const NO_INTERSECTION_RAY_SHAPE = RayShapeIntersectionResult(false, NaN, SVector(NaN, NaN, NaN), 0)

"""
    Base.show(io::IO, result::RayShapeIntersectionResult)

RayShapeIntersectionResultオブジェクトを表示するための関数。

# 引数
- `io`     : 出力ストリーム
- `result` : 表示するRayShapeIntersectionResultオブジェクト
"""
function Base.show(io::IO, result::RayShapeIntersectionResult)
    if result.hit
        msg = """
        Ray-Shape Intersection:
            ∘ hit        = $(result.hit)
            ∘ distance   = $(result.distance)
            ∘ point      = $(result.point)
            ∘ face_index = $(result.face_index)
        """
        print(io, msg)
    else
        msg = """
        Ray-Shape Intersection:
            ∘ hit = $(result.hit)
        """
        print(io, msg)
    end
end


"""
    intersect_ray_shape(ray::Ray, shape::AsteroidThermoPhysicalModels.ShapeModel) -> RayShapeIntersectionResult

Möller–Trumbore法を用いてレイと形状モデルの交差判定を行う。

# 引数
- `ray`   : レイ
- `shape` : 形状モデル

# 戻り値
- 交差判定結果を格納した`RayShapeIntersectionResult`オブジェクト
"""
function intersect_ray_shape(ray::Ray, shape::AsteroidThermoPhysicalModels.ShapeModel)
    min_distance   = Inf
    closest_point  = @SVector zeros(3)
    hit_face_index = 0
    hit_any        = false
    
    # 形状モデルの全ての面について交差判定
    for (i, face) in enumerate(shape.faces)

        ## Early-out 1: Backface culling (skip if ray is coming from the backside of the face)
        n̂ = shape.face_normals[i]
        dot(ray.direction, n̂) ≥ 0 && continue  # skip this face

        ## Early-out 2: Visibility check from observer
        c = shape.face_centers[i]
        dot(c - ray.origin, n̂) ≥ 0 && continue  # skip this face
        
        # 面を構成する頂点
        v1 = shape.nodes[face[1]]
        v2 = shape.nodes[face[2]]
        v3 = shape.nodes[face[3]]
        
        # レイと三角形の交差判定
        result = intersect_ray_triangle(ray, v1, v2, v3)
        
        # 交差があり、これまでの最小距離よりも近い場合
        if result.hit && result.distance < min_distance
            min_distance   = result.distance
            closest_point  = result.point
            hit_face_index = i
            hit_any        = true
        end
    end
    
    if hit_any
        return RayShapeIntersectionResult(true, min_distance, closest_point, hit_face_index)
    else
        return NO_INTERSECTION_RAY_SHAPE
    end
end
