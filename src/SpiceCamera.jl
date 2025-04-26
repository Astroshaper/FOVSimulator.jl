
"""
    struct SpiceCameraStatic

静的なカメラ情報を保持する構造体

# Fields
- `name`       : カメラ名
- `id`         : カメラID
- `fov_shape`  : FOV形状（"RECTANGLE", "CIRCLE", "ELLIPSE"など）
- `fov_frame`  : FOVの参照フレーム
- `boresight` : FOV参照フレームにおけるボアサイトベクトル
- `bounds`    : FOV参照フレームにおける境界ベクトル
"""
struct SpiceCameraStatic
    name::String
    id::Int
    fov_shape::String
    fov_frame::String
    boresight::SVector{3, Float64}
    bounds::Vector{SVector{3, Float64}}
end

"""
    struct SpiceCameraState

動的なカメラ状態を保持する構造体

# Fields
- `boresight` : 現在のフレームにおけるボアサイトベクトル
- `bounds`    : 現在のフレームにおける境界ベクトル
- `position`  : カメラの位置
- `velocity`  : カメラの速度
"""
struct SpiceCameraState
    boresight::SVector{3, Float64}
    bounds::Vector{SVector{3, Float64}}
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
end

"""
    mutable struct SpiceCamera

SPICEカーネルに基づくカメラモデル

# Fields
- `static` : 静的なカメラ情報
- `state`  : 動的なカメラ状態
"""
mutable struct SpiceCamera
    static::SpiceCameraStatic
    state::SpiceCameraState
end


function SpiceCamera(name::String, id::Int)
    # SPICEカーネルから静的情報を取得
    fov_shape, fov_frame, boresight, bounds = SPICE.getfov(id)
    
    # 静的情報の構造体を作成
    static = SpiceCameraStatic(name, id, fov_shape, fov_frame, boresight, bounds)
    
    # 動的情報の初期化
    boresight_state = similar(boresight)
    bounds_state = similar(bounds)
    position = @SVector zeros(3)
    velocity = @SVector zeros(3)
    
    # 動的情報の構造体を作成
    state = SpiceCameraState(boresight_state, bounds_state, position, velocity)
    
    # カメラオブジェクトを作成
    cam = SpiceCamera(static, state)
    
    return cam
end


function Base.show(io::IO, cam::SpiceCamera)
    msg =  "Camera parameters\n"
    msg *= "-----------------\n"
    msg *= "Instrument name      : $(cam.static.name)\n"
    msg *= "Instrument ID        : $(cam.static.id)\n"
    msg *= "FOV shape            : $(cam.static.fov_shape)\n"
    msg *= "FOV reference frame  : $(cam.static.fov_frame)\n"
    msg *= "Boresight vector     : $(cam.static.boresight)\n"
    msg *= "FOV boundary vectors : \n"
    for v in cam.static.bounds
        msg *= "    $v\n"
    end
    msg *= "-----------------\n"
    msg *= "Current boresight vector     : $(cam.state.boresight)\n"
    msg *= "Current FOV boundary vectors : \n"
    for v in cam.state.bounds
        msg *= "    $v\n"
    end
    msg *= "-----------------\n"
    msg *= "Camera position : $(cam.state.position)\n"
    msg *= "Camera velocity : $(cam.state.velocity)\n"
    msg *= "-----------------\n"
    
    print(io, msg)
end


"""
    update!(cam::SpiceCamera, et::Float64, ref::String, abcorr::String, obs::String)

指定された時刻とフレームでカメラの状態を更新する

# Arguments
- `cam`    : カメラ
- `et`     : 暦時間
- `ref`    : 目標フレーム
- `abcorr` : 収差補正フラグ
- `obs`    : 観測者名
"""
function update!(cam::SpiceCamera, et::Float64, ref::String, abcorr::String, obs::String)
    # 回転行列を計算
    Rot = RotMatrix{3}(SPICE.pxform(cam.static.fov_frame, ref, et))

    # ボアサイトベクトルを更新
    cam.state.boresight = Rot * cam.static.boresight

    # 境界ベクトルを更新
    for (i, v) in enumerate(cam.static.bounds)
        cam.state.bounds[i] = Rot * v
    end

    # 位置と速度を更新
    state, _ = SPICE.spkezr(cam.static.name, et, ref, abcorr, obs)
    cam.state.position = state[1:3] * 1000
    cam.state.velocity = state[4:6] * 1000

    return
end
