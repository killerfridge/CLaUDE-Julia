module CLaUDE
using Pkg
Pkg.activate("claude")
include("atmosphere.jl")
include("constants.jl")
using .Atmos_Mod
using .Constants
using Plots
using DelimitedFiles
using BenchmarkTools

atmos = initialize_atmosphere(Int32(2500))
# lats, lons, temps, esv, nsv = interpolate_pixels(atmos)

function plotting(atmos::Atmosphere, quiver_resample::Integer=4)
    
    lats, lons, temps, esv, nsv = interpolate_pixels(atmos)

    lons_grid = range(0, τ, length=config.res * 2)
    lats_grid = range(0, π, length=config.res)

    lons_gridded = lons_grid' .* ones(config.res)
    lats_gridded = ones(config.res) .* lats_grid
    atmos_temps = temps(lats_grid, lons_grid)
    atmos_ewv = esv(lats_grid, lons_grid)
    atmos_nsv = nsv(lats_grid, lons_grid)

    open("logs_lon.csv", "w") do io
        writedlm(io, lons_gridded[1:quiver_resample:end, 1:quiver_resample:end],)
    end
    open("logs_lats.csv", "w") do io
        writedlm(io, lats_gridded[1:quiver_resample:end, 1:quiver_resample:end],)
    end
    open("ewv.csv", "w") do io
        writedlm(io, atmos_ewv[1:quiver_resample:end, 1:quiver_resample:end],)
    end
    open("nsv.csv", "w") do io
        writedlm(io, atmos_nsv[1:quiver_resample:end, 1:quiver_resample:end],)
    end

    plt = scatter(
        lons_gridded,
        lats_gridded,
        marker_z = atmos_temps,
        marker = (:rect, 10),
        markerstrokewidth = 0.1,
        aspect_ratio = 1,
        label = false,
        )
    plot!(xlims=(0, τ), ylims=(0, π), showaxis=false, grid=false)
    scatter!(lons, lats, markersize=2, markeralpha=0.7, leg=false)

        #=quiver!(
            lons_gridded[1:quiver_resample:end, 1:quiver_resample:end],
            lats_gridded[1:quiver_resample:end, 1:quiver_resample:end],
            quiver=(
                atmos_ewv[1:quiver_resample:end, 1:quiver_resample:end],
                atmos_nsv[1:quiver_resample:end, 1:quiver_resample:end]
                )
            )=#

    return plt
end


function plotloop(atmos,loops, sun_lon=0.)

    anim = Animation()

    for i in 1:loops

        update_atmos!(atmos; fast=true, s=0.5, sun_lon=sun_lon)

        plotting(atmos, 4)
        frame(anim)

        sun_lon += constants.dt * τ / constants.day

    end

    gif(anim, "anim.gif", fps=32)

end


plotloop(atmos, 160)

end