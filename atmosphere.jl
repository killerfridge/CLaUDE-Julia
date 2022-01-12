module Atmos_Mod
include("constants.jl")
using .Constants
using PyCall

interpolate = pyimport("scipy.interpolate")

export Pixel, update_temp!, update_velocity!, advect!, interpolate, push!, Atmosphere

mutable struct Pixel
    lat::Float64
    lon::Float64
    temp::Float64
    velocity_east_west::Float64
    velocity_north_south::Float64
    coriolis::Float64
    function Pixel(lat::Float64, lon::Float64)
        coriolis::Float64 = cos(lat) * 1e-5
        temp::Float64 = 270 + (20 * sin(lat))
        return new(lat, lon, temp, 0., 0., coriolis)
    end
end


mutable struct Atmosphere
    pixels::Vector{Pixel}
    name::String
    function Atmosphere()
        return new(Vector(), "Atmosphere")
    end
end


function update_temp!(pixel::Pixel, constants::Constant=constants, sun_lon::Float64=0.)

    pixel.temp = constants.dt * (
        constants.solar *
        (1 - constants.albedo) *
        max(0, sin(pixel.lat)) * 
        max(0, sin(lon - sun_lon)) - 
        (constants.boltzman) * pixel.temp^4
    ) / constants.heatcapacity

    return pixel

end

function update_velocity!(pixel::Pixel, constants::Constant=constants)
    # TODO: recreate the field_d_lat/field_d_lon functions from the original (0s as placeholders)
    pixel.velocity_east_west -= constants.dt * (
        pixel.velocity_east_west * 0 + 
        pixel.velocity_north_south * 0 +
        pixel.coriolis * pixel.velocity_north_south + 0
    )
    # TODO: recreate the field_d_lat/field_d_lon functions from the original (0s as placeholders)
    pixel.velocity_north_south -= constants.dt * (
        pixel.velocity_east_west * 0 + 
        pixel.velocity_north_south * 0 +
        pixel.coriolis * pixel.velocity_east_west + 0
    )

    return pixel
end

function advect!(pixel::Pixel, constants::Constant=constants)

    # TODO: recreate the field_d_lat/field_d_lon functions from the original

    pixel.temp -= constants.dt * (
        pixel.temp * 0 + 
        pixel.velocity_east_west * 0 + 
        pixel.temp * 0 +
        pixel.velocity_north_south * 0
    )

    return pixel

end

function Base.push!(s::Atmosphere, x::Pixel)
    Base.push!(s.pixels, x)
end


end