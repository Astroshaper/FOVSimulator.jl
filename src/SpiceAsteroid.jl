"""
    struct SpiceAsteroidStatic

静的な小惑星情報を保持する構造体

# Fields
- `name`  : 小惑星名
- `shape` : 形状モデル
"""
struct SpiceAsteroidStatic
    name::String
    shape::ShapeModel
end

"""
    struct SpiceAsteroidState

動的な小惑星状態を保持する構造体

# Fields
- `position` : 小惑星の位置
- `velocity` : 小惑星の速度
"""
struct SpiceAsteroidState
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
end

"""
    mutable struct SpiceAsteroid

SPICEカーネルに基づく小惑星モデル

# Fields
- `static` : 静的な小惑星情報
- `state`  : 動的な小惑星状態
"""
mutable struct SpiceAsteroid
    static::SpiceAsteroidStatic
    state::SpiceAsteroidState
end


"""
    SpiceAsteroid(name::String, shape::ShapeModel)

小惑星モデルを構築する

# Arguments
- `name`  : 小惑星名
- `shape` : 形状モデル
"""
function SpiceAsteroid(name::String, shape::ShapeModel)
    # 静的情報の構造体を作成
    static = SpiceAsteroidStatic(name, shape)
    
    # 動的情報の初期化
    position = @SVector zeros(3)
    velocity = @SVector zeros(3)
    
    # 動的情報の構造体を作成
    state = SpiceAsteroidState(position, velocity)
    
    # 小惑星オブジェクトを作成
    asteroid = SpiceAsteroid(static, state)
    
    return asteroid
end


function Base.show(io::IO, asteroid::SpiceAsteroid)
    msg =  "Asteroid parameters\n"
    msg *= "-------------------\n"
    msg *= "Asteroid name : $(asteroid.static.name)\n"
    msg *= "Position      : $(asteroid.state.position)\n"
    msg *= "Velocity      : $(asteroid.state.velocity)\n"
    msg *= "-------------------\n"
    
    print(io, msg)
    println(asteroid.static.shape)
end


"""
    update!(asteroid::SpiceAsteroid, et::Float64, ref::String, abcorr::String, obs::String)

指定された時刻とフレームで小惑星の状態を更新する

# Arguments
- `asteroid` : 小惑星
- `et`       : 暦時間
- `ref`      : 目標フレーム
- `abcorr`   : 収差補正フラグ
- `obs`      : 観測者名
"""
function update!(asteroid::SpiceAsteroid, et::Float64, ref::String, abcorr::String, obs::String)
    state, _ = SPICE.spkezr(asteroid.static.name, et, ref, abcorr, obs)
    asteroid.state.position = state[1:3] * 1000
    asteroid.state.velocity = state[4:6] * 1000

    return
end


"""

"""
function transform_shape(asteroid::SpiceAsteroid, from::String, to::String, et::Float64, abcorr::String, obs::String)

    Rot = RotMatrix{3}(SPICE.pxform(from, to, et))
    obs_pos = SVector{3}(SPICE.spkpos(asteroid.static.name, et, to, abcorr, obs)[1]) * 1000

    nodes = [Rot * node + obs_pos for node in asteroid.static.shape.nodes]
    faces = asteroid.static.shape.faces
    
    face_centers = [AsteroidThermoPhysicalModels.face_center(nodes[face]) for face in faces]
    face_normals = [AsteroidThermoPhysicalModels.face_normal(nodes[face]) for face in faces]
    face_areas   = [AsteroidThermoPhysicalModels.face_area(nodes[face])   for face in faces]

    visiblefacets = asteroid.static.shape.visiblefacets

    shape_new = ShapeModel(nodes, faces, face_centers, face_normals, face_areas, visiblefacets)
    
    return shape_new
end


"""
    focal_length(fov_angle::Float64, n_pixel::Int) -> f

Calculate a focal length of a camera.

# Arguments
- `fov_angle` : Field of view angle [deg]
- `n_pixel`   : Number of pixels
"""
focal_length(fov_angle::Float64, n_pixel::Int) = n_pixel / (2 * tan(deg2rad(fov_angle) / 2))


"""

"""
function project_point_fov(p::SVector{3, Float64}, fov_angles::Tuple{Float64, Float64}, img_size::Tuple{Int, Int})
    fov_x, fov_y = fov_angles
    width, height = img_size

    f_x = focal_length(fov_x, width)   # focal length in the x-direction [pixel]
    f_y = focal_length(fov_y, height)  # focal length in the y-direction [pixel]

    c_x = width  / 2  # x-coordinate of the principal point
    c_y = height / 2  # y-coordinate of the principal point

    u = f_x * p[1] / p[3] + c_x
    v = f_y * p[2] / p[3] + c_y
    
    return (round(Int, u), round(Int, v))

    # カメラの背面(p[3]<0)は写らない
end
