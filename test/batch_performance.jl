"""
Test batch ray processing performance
"""

@testset "Batch Ray Processing Performance" begin
    msg = """\n
    ┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
    |          Test: Batch Ray Processing Performance        |
    ┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛
    """
    println(msg)
    
    # Download test files if needed
    paths_kernel = [
        "lsk/naif0012.tls",
        "pck/pck00011.tpc",
        "fk/hera_v14.tf",
        "ik/hera_tiri_v03.ti",
        "spk/de432s.bsp",
    ]
    
    path_shape = "shape/g_50677mm_rad_obj_didy_0000n00000_v001.obj"
    
    for path_kernel in paths_kernel
        url_kernel = "https://spiftp.esac.esa.int/data/SPICE/HERA/kernels/$(path_kernel)"
        filepath = joinpath(@__DIR__, "kernel", path_kernel)
        mkpath(dirname(filepath))
        isfile(filepath) || Downloads.download(url_kernel, filepath)
    end
    
    url_shape = "https://spiftp.esac.esa.int/data/SPICE/HERA/kernels/dsk/g_50677mm_rad_obj_didy_0000n00000_v001.obj"
    filepath_shape = joinpath(@__DIR__, path_shape)
    mkpath(dirname(filepath_shape))
    isfile(filepath_shape) || Downloads.download(url_shape, filepath_shape)
    
    # Load SPICE kernels
    for path_kernel in paths_kernel
        filepath = joinpath(@__DIR__, "kernel", path_kernel)
        SPICE.furnsh(filepath)
    end
    
    # Load shape model with BVH
    shape = AsteroidShapeModels.load_shape_obj(filepath_shape; scale=1000, with_face_visibility=false, with_bvh=true)
    
    # Create a simple test setup
    img_size = (1024, 768)   # Image size same as HERA/TIRI
    
    # Create rays directly in asteroid-fixed frame
    height, width = img_size
    rays = Matrix{Ray}(undef, height, width)
    
    # Camera at 1000 km distance looking at origin
    camera_pos = @SVector [1000e3, 0.0, 0.0]
    
    # Simple pinhole camera model
    fov_x = deg2rad(13.3)
    fov_y = deg2rad(10.0)
    
    for v in 1:height
        for u in 1:width
            # Normalized image coordinates (-1 to 1)
            x = 2.0 * (u - 1) / (width - 1) - 1.0
            y = 2.0 * (v - 1) / (height - 1) - 1.0
            
            # Ray direction
            dx = -1.0  # Looking towards origin
            dy = tan(fov_x/2) * x
            dz = tan(fov_y/2) * y
            
            ray_direction = normalize(@SVector [dx, dy, dz])
            rays[v, u] = Ray(camera_pos, ray_direction)
        end
    end
    
    println("\nShape model:")
    println("- Number of faces : $(length(shape.faces))")
    println("- Image size      : $(img_size)")
    println("- Total rays      : $(prod(img_size))")
    
    # Test 1: Loop-based approach (old method)
    function generate_intersection_map_loop(rays, shape)
        height, width = size(rays)
        intersection_map = Matrix{RayShapeIntersectionResult}(undef, height, width)
        
        for v in 1:height
            for u in 1:width
                intersection_map[v, u] = intersect_ray_shape(rays[v, u], shape)
            end
        end
        
        return intersection_map
    end
    
    # Test 2: Batch approach (new method)
    function generate_intersection_map_batch(rays, shape)
        return intersect_ray_shape(rays, shape)
    end
    
    # Warm up
    generate_intersection_map_loop(rays[1:10, 1:10], shape)
    generate_intersection_map_batch(rays[1:10, 1:10], shape)
    
    # Performance comparison
    println("\nPerformance comparison:")
    print("- Loop-based approach : ")
    t_loop = @elapsed result_loop = generate_intersection_map_loop(rays, shape)
    println("$(round(t_loop, digits=5)) seconds")
    
    print("- Batch approach      : ")
    t_batch = @elapsed result_batch = generate_intersection_map_batch(rays, shape)
    println("$(round(t_batch, digits=5)) seconds")
    
    speedup = t_loop / t_batch
    println("- Speedup             : $(round(speedup, digits=1))x")
    
    # Verify results are identical
    @test size(result_loop) == size(result_batch)
    for i in eachindex(result_loop)
        @test result_loop[i].hit == result_batch[i].hit
        if result_loop[i].hit
            @test result_loop[i].face_idx == result_batch[i].face_idx
            @test result_loop[i].distance ≈ result_batch[i].distance rtol=1e-10
        end
    end
    
    # Clean up
    SPICE.kclear()
end