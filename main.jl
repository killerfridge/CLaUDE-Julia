module CLaudE
include("atmosphere.jl")
include("constants.jl")
using .Atmos_Mod
using .Constants
using Plots
using BenchmarkTools

atmos = initialize_atmosphere(Int32(1500))
lats, lons, temps, esv, nsv = interpolate_pixels(atmos)

function plotting(atmos::Atmosphere, quiver_resample::Integer=4)
    

    lons_grid = range(0, τ, length=config.res)
    lats_grid = range(0, π, length=config.res)

    lons_gridded = lons_grid' .* ones(config.res)
    lats_gridded = ones(config.res) .* lats_grid
    atmos_temps = temps(lats_grid, lons_grid)
    atmos_esv = esv(lats_grid, lons_grid)
    atmos_nsv = nsv(lats_grid, lons_grid)

    plt = scatter(lons_gridded, lats_gridded, marker_z = atmos_temps, marker = (:rect, 10), markerstrokewidth = 0.1, aspect_ratio = 1, label = false)
    scatter!(lons, lats, markersize=2, markeralpha=0.7)
    quiver!(
        lons_gridded[1:quiver_resample:end, 1:quiver_resample:end],
        lats_gridded[1:quiver_resample:end, 1:quiver_resample:end],
        quiver=(
            atmos_esv[1:quiver_resample:end, 1:quiver_resample:end],
            atmos_nsv[1:quiver_resample:end, 1:quiver_resample:end]
            )
        )

    display(plt)
end
sun_lon = 0.0
#=for i in 1:100

    for pixel in atmos.pixels
        update_temp!(pixel)
        update_velocity!(pixel)
        advect!(pixel)
    end

    plotting(atmos, 4)

    sun_lon += constants.dt * τ / constants.day
    
end =#

# int_atmos = InterpolatedAtmosphere(atmos)


t = @benchmark update_atmos!(atmos)
println("Median")
println(median(t))
println("Maximum")
println(maximum(t))
println("Minimum")
println(minimum(t))

end