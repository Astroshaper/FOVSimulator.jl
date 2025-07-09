# FOVSimulator.jl

`FOVSimulator.jl` is a Julia package for simulating camera field of view (FOV) in space missions. It specializes in camera simulations for observing celestial bodies such as asteroids, providing high-precision position and attitude calculations using SPICE kernels and image generation capabilities through ray tracing.

## Features

- High-precision position and attitude calculations for spacecraft, cameras, and celestial bodies using SPICE kernels
- Ray tracing-based intersection detection for asteroid shape models
- Thermal image simulation in conjunction with thermophysical models
- Fast intersection detection algorithms using bounding boxes
- Support for multiple bodies (e.g., binary asteroid systems)

## Installation

You can install the package using Julia's package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/Astroshaper/FOVSimulator.jl")
```

Or in the package manager mode in Julia REPL (press `]` to enter):

```julia
] add https://github.com/Astroshaper/FOVSimulator.jl
```

## Dependencies

- [`LinearAlgebra.jl`](https://github.com/JuliaLang/LinearAlgebra.jl)
- [`Rotations.jl`](https://github.com/JuliaGeometry/Rotations.jl)
- [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl)
- [`SPICE.jl`](https://github.com/JuliaAstro/SPICE.jl)
- [`AsteroidShapeModels.jl`](https://github.com/Astroshaper/AsteroidShapeModels.jl)

## Basic Usage

### 1. Import Required Modules

```julia
using FOVSimulator
using SPICE
```

### 2. Load SPICE Kernels

```julia
# Load SPICE kernels
SPICE.furnsh("path/to/kernel/file.tm")
```

### 3. Set Up Camera

```julia
# Specify camera NAIF ID, FOV angles, and image size
camera_id = -12345  # Camera NAIF ID
fov_angles = (10.0, 8.0)  # FOV angles (width, height) [degrees]
img_size = (1024, 768)  # Image size (width, height) [pixels]

# Create camera object
camera = SpiceCamera("CAMERA_NAME", camera_id, fov_angles, img_size)
```

### 4. Set Up Spacecraft

```julia
# Create spacecraft object
spacecraft = SpiceSpacecraft("SPACECRAFT_NAME")

# Mount camera on spacecraft
add_instrument!(spacecraft, camera)
```

### 5. Set Up Asteroid

```julia
# Load shape model
# Load an asteroid shape model
# - `path/to/shape.obj` is the path to your OBJ file (mandatory)
# - `scale` : scale factor for the shape model (e.g., 1000 for km to m conversion)
# - `with_face_visibility` : whether to build face-to-face visibility graph for illumination checking and thermophysical modeling
# - `with_bvh` : whether to build BVH for ray tracing
shape = load_shape_obj("path/to/shape.obj"; scale=1000, with_face_visibility=false, with_bvh=true)

# Create asteroid object
asteroid = SpiceAsteroid("ASTEROID_NAME", shape, "ASTEROID_FIXED_FRAME")
```

### 6. Update States

```julia
# Set simulation time
et     = SPICE.str2et("2025-01-01T12:00:00")  # Ephemeris time
ref    = "J2000"                              # Reference frame
abcorr = "NONE"                               # Aberration correction flag
obs    = "EARTH"                              # Observer

# Update spacecraft and asteroid states
update!(spacecraft, et, ref, abcorr, obs)
update!(asteroid, et, ref, abcorr, obs)
```

### 7. Generate Intersection Map

```julia
# Perform intersection detection for each pixel
intersection_map = generate_intersection_map(spacecraft.state.instruments["CAMERA_NAME"], asteroid)
```

### 8. Generate Image

```julia
# Generate thermal image using temperature data
temperatures = fill(300.0, length(asteroid.static.shape.faces))  # Example: set all faces to 300K
image = generate_image_temperature(intersection_map, temperatures)

# Generate radiance image using emissivity and temperature data
emissivities = fill(0.9, length(asteroid.static.shape.faces))  # Example: set all faces to 0.9 emissivity
radiance_image = generate_image_radiance(intersection_map, emissivities, temperatures)
```

## Advanced Usage

### Binary Asteroid Simulation

```julia
# Set up primary and secondary asteroids
primary = SpiceAsteroid("PRIMARY", shape_primary, "PRIMARY_FIXED_FRAME")
secondary = SpiceAsteroid("SECONDARY", shape_secondary, "SECONDARY_FIXED_FRAME")

# Update states
update!(spacecraft, et, ref, abcorr, obs)
update!(primary, et, ref, abcorr, obs)
update!(secondary, et, ref, abcorr, obs)

# Generate intersection maps
intersection_map1 = generate_intersection_map(spacecraft.state.instruments["CAMERA_NAME"], primary)
intersection_map2 = generate_intersection_map(spacecraft.state.instruments["CAMERA_NAME"], secondary)

# Generate thermal image (considering both asteroids)
image = generate_image_temperature(
    intersection_map1, temperatures_primary,
    intersection_map2, temperatures_secondary
)
```

### Integration with Thermophysical Models

You can integrate with thermophysical modeling packages to calculate asteroid surface temperatures and use them for thermal image simulation. See the test file `test/TIRI_image_Didymos.jl` for a detailed example.

## Key Components

### `SpiceCamera`

Represents a camera model based on SPICE kernels. Manages static camera information (name, ID, FOV shape, etc.) and dynamic state (position, velocity, etc.).

### `SpiceSpacecraft`

Represents a spacecraft model based on SPICE kernels. Manages static spacecraft information and dynamic state, and can carry multiple instruments (cameras, etc.).

### `SpiceAsteroid`

Represents an asteroid model based on SPICE kernels. Manages static asteroid information (name, shape model, fixed frame, etc.) and dynamic state (position, velocity, etc.).

### Ray Tracing Components

The ray tracing functionality (Ray, BoundingBox, intersection detection functions) is now provided by the [`AsteroidShapeModels.jl`](https://github.com/Astroshaper/AsteroidShapeModels.jl) package. This includes:

- `Ray`: Represents a ray for ray tracing
- `BoundingBox`: Represents a bounding box for shape models
- `intersect_ray_triangle`: Detects intersection between a ray and a triangle
- `intersect_ray_shape`: Detects intersection between a ray and a shape model
- `intersect_ray_bounding_box`: Detects intersection between a ray and a bounding box

### Image Generation Functions

- `generate_pixel_rays`: Generates rays corresponding to each pixel of a camera
- `generate_intersection_map`: Generates intersection results for each pixel and an asteroid surface
- `generate_image_temperature`: Generates a thermal image from intersection results and surface temperatures
- `generate_image_radiance`: Generates a radiance image from intersection results, emissivities, and surface temperatures

## License

This package is released under the MIT License. See the [`LICENSE`](https://github.com/Astroshaper/FOVSimulator.jl/blob/main/LICENSE) file for details.
