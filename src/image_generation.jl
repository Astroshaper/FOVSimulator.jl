"""
    image_generation.jl

カメラの視野情報からピクセルごとのレイを生成し、熱画像を生成する機能を提供する。
バウンディングボックスを用いた高速化機能も提供する。
"""

"""
    focal_length(fov_angle::Float64, n_pixel::Int) -> f

Calculate a focal length of a camera.

# Arguments
- `fov_angle` : Field of view angle [deg]
- `n_pixel`   : Number of pixels
"""
focal_length(fov_angle::Real, n_pixel::Int) = n_pixel / (2 * tan(deg2rad(fov_angle) / 2))

"""
    generate_pixel_rays(cam::SpiceCamera) -> Matrix{Ray}

SpiceCameraオブジェクトからカメラの視野情報を取得して、カメラ座標系において各ピクセルに対応するレイを生成する。

# 引数
- `cam` : SpiceCameraオブジェクト

# 戻り値
- `rays` : 各ピクセルに対応するレイの2次元配列（画像サイズと同じ形状）
"""
function generate_pixel_rays(cam::SpiceCamera)
    # Get camera parameters
    fov_angles = cam.static.fov_angles  # Field of view angles (width, height) in degrees
    img_size = cam.static.img_size      # Image size (width, height) in pixels
    width, height = img_size
    
    # Calculate focal lengths in pixels
    fov_x, fov_y = fov_angles
    f_x = focal_length(fov_x, width)   # focal length in the x-direction [pixel]
    f_y = focal_length(fov_y, height)  # focal length in the y-direction [pixel]
    
    # Calculate principal point (center of the image)
    c_x = width / 2
    c_y = height / 2
    
    # Initialize ray matrix
    rays = Matrix{Ray}(undef, height, width)
    
    # Generate rays for each pixel
    for v in 1:height
        for u in 1:width
            # Convert pixel coordinates to normalized image coordinates
            # (0,0) is at the center of the image
            x = (u - c_x - 0.5) / f_x  # -0.5 to get the center of the pixel
            y = (v - c_y - 0.5) / f_y
            
            # Create a ray in camera frame
            # Z-axis is the boresight direction
            ray_origin = @SVector zeros(3)
            ray_direction = normalize(SVector{3, Float64}(x, y, 1.0))
            
            # Create ray with camera position as origin and calculated direction
            rays[v, u] = Ray(ray_origin, ray_direction)
        end
    end
    
    return rays
end

"""
    generate_pixel_rays(cam::SpiceCamera, asteroid::SpiceAsteroid) -> Matrix{Ray}

SpiceCameraオブジェクトからカメラの視野情報を取得して、各ピクセルに対応するレイを生成し、
小惑星固定座標系に変換する。

# 引数
- `cam`     : SpiceCameraオブジェクト
- `asteroid`: SpiceAsteroidオブジェクト
- `et`      : エフェメリス時刻
- `abcorr`  : 光行差補正フラグ（例："LT+S"）

# 戻り値
- `rays`    : 小惑星固定座標系における各ピクセルに対応するレイの2次元配列（画像サイズと同じ形状）
"""
function generate_pixel_rays(cam::SpiceCamera, asteroid::SpiceAsteroid)

    if cam.state.et != asteroid.state.et
        error("The ephemeris time of the camera and asteroid must match. Update the camera and asteroid state before calling this function.")
    end

    if cam.state.abcorr != asteroid.state.abcorr
        error("The `abcorr` flag of the camera and asteroid must match. Update the camera and asteroid state before calling this function.")
    end

    et = cam.state.et
    abcorr = cam.state.abcorr

    # Generate rays in camera frame
    rays = generate_pixel_rays(cam)
    
    # Get transformation from camera frame to asteroid-fixed frame
    from = cam.static.fov_frame        # Camera frame
    to = asteroid.static.asteroid_fixed_frame  # Asteroid-fixed frame
    
    # Calculate rotation matrix from camera frame to asteroid-fixed frame
    Rot = RotMatrix{3}(SPICE.pxform(from, to, et))
    
    # Get camera position in asteroid-fixed frame
    camera_pos_asteroid_frame = SVector{3}(SPICE.spkpos(cam.static.name, et, to, abcorr, asteroid.static.name)[1]) * 1000
    
    # Get dimensions of the ray matrix
    height, width = size(rays)
    
    # Transform each ray from camera frame to asteroid-fixed frame
    for v in 1:height
        for u in 1:width
            # Transform ray direction from camera frame to asteroid-fixed frame
            ray_direction = Rot * rays[v, u].direction
            
            # Update ray with new origin and direction
            rays[v, u] = Ray(camera_pos_asteroid_frame, ray_direction)
        end
    end
    
    return rays
end

"""
    generate_thermal_image(rays::Matrix{Ray}, asteroid::SpiceAsteroid, emissivities::Vector{Float64}, temperatures::Vector{Float64}) -> Matrix{Float64}

`generate_pixel_rays`関数で生成したレイ、小惑星オブジェクト、各面に与えた温度と放射率をもとに、赤外線カメラの模擬画像を作成する。
小惑星オブジェクトに格納されているバウンディングボックスを使用して高速化する。

# 引数
- `rays`         : 各ピクセルに対応するレイの2次元配列（`generate_pixel_rays`関数で生成）
- `asteroid`     : 小惑星オブジェクト
- `emissivities` : 各面の放射率（0-1）、face数と同じ長さ
- `temperatures` : 各面の温度 [K]、face数と同じ長さ

# 戻り値
- `image`        : 画素ごとの放射輝度 [W/m²/sr] を格納した2次元配列
"""
function generate_thermal_image(rays::Matrix{Ray}, asteroid::SpiceAsteroid, emissivities::Vector{Float64}, temperatures::Vector{Float64})
    shape = asteroid.static.shape  # Get shape model from asteroid object
    
    # Check if the number of emissivities and temperatures matches the number of faces
    n_face = length(shape.faces)
    if length(emissivities) != n_face || length(temperatures) != n_face
        error("The number of emissivities and temperatures must match the number of faces in the shape model.")
    end
    
    # Initialize image with the dimensions of the ray matrix
    height, width = size(rays)
    image = zeros(Float64, height, width)
    
    # Process each pixel
    for v in 1:height
        for u in 1:width
            # Get ray for this pixel
            ray = rays[v, u]
            
            # Perform ray-shape intersection using the asteroid's bounding box
            intersection = intersect_ray_shape(ray, asteroid)
            
            # If intersection occurred
            if intersection.hit
                # Get emissivity and temperature for this face
                ε = emissivities[intersection.face_index]
                T = temperatures[intersection.face_index]
                
                # Calculate total radiance using Stefan-Boltzmann law with emissivity
                # E = ε·σT⁴ [W/m²] (total emitted power per unit area)
                # L = E/π [W/m²/sr] (radiance assuming Lambertian surface)
                radiance = ε * σ_SB * T^4 / π  # [W/m²/sr]
                
                image[v, u] = radiance
            else
                # If no intersection, set pixel to zero (or some other value)
                image[v, u] = 0.0
            end
        end
    end
    
    return image
end
