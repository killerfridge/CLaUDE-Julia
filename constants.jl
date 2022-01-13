module Constants

export constants, Constant, config, τ

struct Constant
    solar::Int16
    heatcapacity::Int32
    albedo::Float32
    radius::Int32
    dt::Int32
    day::Int32
    year::Float32
    boltzman::Float64

    function Constant()
        return new(
            1370,
            1e5, 
            0.3, 
            6.4e6,
            18*60, 
            60 * 60 * 24,
            365.25 * 60 * 60 * 24,
            5.67e-8
            )
    end
end

struct Config
    res::Int32
    function Config()
        return new(75)
    end
end

constants = Constant()
config = Config()
τ = 2 * π

end