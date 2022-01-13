module Fibonacci_Mod

export fibonacci_sphere

function fibonacci_sphere(samples::Int32)::Vector{Tuple{Float64, Float64}}
    points::Vector{Tuple{Float64, Float64}} = []
    ϕ::Float64 = π * (3. - √5)
    for i = 0:samples
        y::Float64 = 1 - (i / samples) * 2
        radius::Float64 = √(1 - y^2)
        golden_angle::Float64 = ϕ * i

        x::Float64 = cos(golden_angle) * radius
        z::Float64 = sin(golden_angle) * radius

        θ::Float64 = acos(z)

        φ::Float64 = atan(y, x) + π
        if 0.05 * π < θ < 0.95 * π
            push!(points, (φ, θ))
        end
    end
    return points
end

end