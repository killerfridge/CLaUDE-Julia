module Atmos_Mod
include("constants.jl")
include("fibonacci.jl")
using .Constants
using .Fibonacci_Mod
using PyCall

interpolate = pyimport("scipy.interpolate")

export Pixel, 
update_temp!,
update_velocity!, 
advect!,
interpolate,
push!, 
Atmosphere, 
initialize_atmosphere,
getindex, 
interpolate_pixels,
InterpolatedAtmosphere,
update_atmos!,
Constant,
constants,
config,
τ

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
    function Atmosphere()
        return new(Vector())
    end
end

mutable struct InterpolatedAtmosphere
    pixels::Vector{Pixel}
    lats::Vector{Float64}
    lons::Vector{Float64}
    temps
    ewv
    nsv
    function InterpolatedAtmosphere(x::Atmosphere, s::Number=4)
        return new(x.pixels, interpolate_pixels(x, s)...)
    end
end


function update_temp!(pixel::Pixel; constants::Constant=constants, sun_lon::Float64=0.)

    pixel.temp += constants.dt * (
        constants.solar *
        (1 - constants.albedo) *
        max(0, sin(pixel.lat)) * 
        max(0, sin(pixel.lon - sun_lon)) - 
        (constants.boltzman) * pixel.temp^4
    ) / constants.heatcapacity

    return pixel

end

function update_temp!(atmos::Atmosphere; constants::Constant=constants, sun_lon::Float64=0.)

    for pixel in atmos.pixels
        update_temp!(pixel; constants=constants, sun_lon=sun_lon)
    end

    return atmos
end

# A couple of helper functions for the update velocity/update temp/advect functions
field_d_lon(field, lat, lon, constants=constants) = field(lat, lon, dtheta=1)[1]/constants.radius

field_d_lat(field, lat, lon, constants=constants) = field(lat, lon, dphi=0)[1]/(constants.radius*sin(lat))

function update_velocity!(pixel::Pixel, ewv::PyCall.PyObject, nsv::PyCall.PyObject, temps::PyCall.PyObject, constants::Constant=constants)
    # TODO: recreate the field_d_lat/field_d_lon functions from the original (0s as placeholders)
    pixel.velocity_east_west -= constants.dt * (
        pixel.velocity_east_west * field_d_lon(ewv, pixel.lat, pixel.lon, constants) + 
        pixel.velocity_north_south * field_d_lat(ewv, pixel.lat, pixel.lon, constants) +
        pixel.coriolis * pixel.velocity_north_south + field_d_lat(temps, pixel.lat, pixel.lon, constants)
    )
    # TODO: recreate the field_d_lat/field_d_lon functions from the original (0s as placeholders)
    pixel.velocity_north_south -= constants.dt * (
        pixel.velocity_east_west * field_d_lon(nsv, pixel.lat, pixel.lon, constants) + 
        pixel.velocity_north_south * field_d_lat(nsv, pixel.lat, pixel.lon, constants)+
        pixel.coriolis * pixel.velocity_east_west + field_d_lat(temps, pixel.lat, pixel.lon, constants)
    )

    return pixel
end

function update_velocity!(atmos::Atmosphere, constants::Constant=constants)
    inter_atmos = InterpolatedAtmosphere(atmos)
    for pixel in atmos.pixels
        update_velocity!(pixel, inter_atmos.ewv, inter_atmos.nsv, inter_atmos.temps, constants)
    end

    return atmos
end

function update_velocity!(atmos::Atmosphere, inter_atmos::InterpolatedAtmosphere, constants::Constant=constants)
    for pixel in atmos.pixels
        update_velocity!(pixel, inter_atmos.ewv, inter_atmos.nsv, inter_atmos.temps, constants)
    end

    return atmos
end


function advect!(pixel::Pixel, ewv, nsv, temps, constants::Constant=constants)

    # TODO: recreate the field_d_lat/field_d_lon functions from the original

    pixel.temp -= constants.dt * (
        pixel.temp * field_d_lon(ewv, pixel.lat, pixel.lon, constants) + 
        pixel.velocity_east_west * field_d_lon(temps, pixel.lat, pixel.lon, constants) + 
        pixel.temp * field_d_lat(nsv, pixel.lat, pixel.lon, constants) +
        pixel.velocity_north_south * field_d_lat(temps, pixel.lat, pixel.lon, constants)
    )

    return pixel

end

function advect!(atmos::Atmosphere, constants::Constant=constants)
    inter_atmos = InterpolatedAtmosphere(atmos)
    for pixel in atmos.pixels
        advect!(pixel, inter_atmos.ewv, inter_atmos.nsv, inter_atmos.temps, constants)
    end

    return atmos
end

function advect!(atmos::Atmosphere, inter_atmos::InterpolatedAtmosphere, constants::Constant=constants)

    for pixel in atmos.pixels
        advect!(pixel, inter_atmos.ewv, inter_atmos.nsv, inter_atmos.temps, constants)
    end

    return atmos
end

function update_atmos!(atmos::Atmosphere; constants::Constant=constants, sun_lon::Float64=0., s::Number=4)

    update_temp!(atmos; constants=constants, sun_lon=sun_lon)
    update_velocity!(atmos, constants)
    advect!(atmos, constants)

end

function update_atmos!(atmos::Atmosphere; constants::Constant=constants, sun_lon::Float64=0., fast::Bool=true, s::Number=4)

    # inter_atmos = InterpolatedAtmosphere(atmos, s)
    update_temp!(atmos; constants=constants, sun_lon=sun_lon)
    # update_velocity!(atmos, inter_atmos, constants)
    # advect!(atmos, inter_atmos, constants)
    return InterpolatedAtmosphere(atmos, s)

end


function Base.push!(s::Atmosphere, x::Pixel)
    Base.push!(s.pixels, x)
end

function initialize_atmosphere(samples::Integer)

    atmosphere = Atmosphere()

    points = fibonacci_sphere(Int32(samples))
    for point in points
        push!(atmosphere, Pixel(point[2], point[1]))
    end

    return atmosphere

end

function Base.getindex(a::Pixel, i::Symbol)

    return getfield(a, i)

end

function interpolate_pixels(atmos::Atmosphere, s::Number=4)
    lats = [pixel.lat for pixel in atmos.pixels]
    lons = [pixel.lon for pixel in atmos.pixels]
    temp_list = [pixel.temp for pixel in atmos.pixels]
    ewv_list = [pixel.velocity_east_west for pixel in atmos.pixels]
    nsv_list = [pixel.velocity_north_south for pixel in atmos.pixels]

    temps = interpolate.SmoothSphereBivariateSpline(lats, lons, temp_list, s=s)
    ewv = interpolate.SmoothSphereBivariateSpline(lats, lons, ewv_list, s=s)
    nsv = interpolate.SmoothSphereBivariateSpline(lats, lons, nsv_list, s=s)

    return lats, lons, temps, ewv, nsv

end

end