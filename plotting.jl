module Plotting_Mod

using Plots
include("atmosphere.jl")
using .Atmos_Mod

export plotting

const τ = π * 2

function plotting(inter_atmos::Main.CLaUDE.Atmos_Mod.InterpolatedAtmosphere, quiver_resample::Integer=4)
    
    lats = inter_atmos.lats
    lons = inter_atmos.lons
    temps = inter_atmos.temps
    ewv = inter_atmos.ewv
    nsv = inter_atmos.nsv

    lons_grid = range(0, τ, length=config.res * 2)
    lats_grid = range(0, π, length=config.res)

    lons_gridded = lons_grid' .* ones(config.res)
    lats_gridded = ones(config.res) .* lats_grid
    atmos_temps = temps(lats_grid, lons_grid)
    atmos_ewv = ewv(lats_grid, lons_grid)
    atmos_nsv = nsv(lats_grid, lons_grid)

    plt = scatter(
        lons_gridded,
        lats_gridded,
        marker_z = atmos_temps,
        marker = (:rect, 10),
        markerstrokewidth = 0.0,
        aspect_ratio = 1,
        label = false,
        )
    scatter!(lons, lats, markersize=2, markeralpha=0.7, leg=false)
    plot!(xlims=(0, τ), ylims=(0, π), showaxis=false, grid=false, showgrid=false)
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

end