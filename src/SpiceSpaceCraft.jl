"""
    struct SpiceSpacecraftStatic

静的な宇宙機情報を保持する構造体

# Fields
- `name` : 宇宙機名
"""
struct SpiceSpacecraftStatic
    name::String
end

"""
    struct SpiceSpacecraftState

動的な宇宙機状態を保持する構造体

# Fields
- `position`    : 宇宙機の位置
- `velocity`    : 宇宙機の速度
- `instruments` : 搭載機器（カメラなど）
"""
struct SpiceSpacecraftState
    position::SVector{3, Float64}
    velocity::SVector{3, Float64}
    instruments::Dict{String, SpiceCamera}
end

"""
    mutable struct SpiceSpacecraft

SPICEカーネルに基づく宇宙機モデル

# Fields
- `static` : 静的な宇宙機情報
- `state`  : 動的な宇宙機状態
"""
mutable struct SpiceSpacecraft
    static::SpiceSpacecraftStatic
    state::SpiceSpacecraftState
end


"""
    SpiceSpacecraft(name::String)

宇宙機モデルを構築する

# Arguments
- `name` : 宇宙機名
"""
function SpiceSpacecraft(name::String)
    # 静的情報の構造体を作成
    static = SpiceSpacecraftStatic(name)
    
    # 動的情報の初期化
    position = @SVector zeros(3)
    velocity = @SVector zeros(3)
    instruments = Dict{String, SpiceCamera}()
    
    # 動的情報の構造体を作成
    state = SpiceSpacecraftState(position, velocity, instruments)
    
    # 宇宙機オブジェクトを作成
    spacecraft = SpiceSpacecraft(static, state)
    
    return spacecraft
end


function Base.show(io::IO, spacecraft::SpiceSpacecraft)
    msg =  "Spacecraft parameters\n"
    msg *= "---------------------\n"
    msg *= "Spacecraft name : $(spacecraft.static.name)\n"
    msg *= "Position        : $(spacecraft.state.position)\n"
    msg *= "Velocity        : $(spacecraft.state.velocity)\n"
    msg *= "Instruments     : $(keys(spacecraft.state.instruments))\n"
    msg *= "---------------------\n"
    
    print(io, msg)
end


"""
    add_instrument!(spacecraft::SpiceSpacecraft, instrument::SpiceCamera)

宇宙機に機器を追加する

# Arguments
- `spacecraft` : 宇宙機
- `instrument` : 追加する機器（カメラなど）
"""
function add_instrument!(spacecraft::SpiceSpacecraft, instrument::SpiceCamera)
    spacecraft.state.instruments[instrument.static.name] = instrument

    return
end

"""
    update!(spacecraft::SpiceSpacecraft, et::Float64, ref::String, abcorr::String, obs::String)

指定された時刻とフレームで宇宙機の状態を更新する

# Arguments
- `spacecraft` : 宇宙機
- `et`         : 暦時間
- `ref`        : 目標フレーム
- `abcorr`     : 収差補正フラグ
- `obs`        : 観測者名

c.f. https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html
"""
function update!(spacecraft::SpiceSpacecraft, et::Float64, ref::String, abcorr::String, obs::String)
    state, _ = SPICE.spkezr(spacecraft.static.name, et, ref, abcorr, obs)

    spacecraft.state.position = state[1:3] * 1000
    spacecraft.state.velocity = state[4:6] * 1000

    for (_, instrument) in spacecraft.state.instruments
        update!(instrument, et, ref, abcorr, obs)
    end

    return
end
