import Base.@kwdef

@kwdef struct EarthProperties
    mu::Real
    j_2::Real
    r::Real
    irradiance::Real
end

@kwdef struct SolarPanel{T<:Real}
    normal::Vector{T}
    efficiency::T
    area::T
end

# assumedly a transmitter, maybe rename as such
@kwdef struct Antenna{T<:Real}
    normal::Vector{T}
    pattern::Function
    power::T
end

@kwdef struct PowerProperties{T<:Real}
    capacity::T
    consumption::T      # later, make this a lookup by Mode
    solarpanels::Vector{SolarPanel}
end

@kwdef struct DataProperties{T<:Real}
    capacity::T
    transmit::T
    antennas::Vector{Antenna}
end

@kwdef struct MassProperties{T<:Real}
    mass::T
    inertia::Matrix{T}
end

@kwdef struct Mode
    name::String
    entry::Function
    power::Real
    data::Real
end

@kwdef struct SpacecraftProperties
    name::String
    power::PowerProperties
    # data::DataProperties
    mass::MassProperties
    # modes::Vector{Mode}
end

@kwdef struct Mission
    name::String
    # modes::Vector{Mode}
    spacecraft::SpacecraftProperties
    # targets::Vector{FrameFixedTarget}
end

#  make all these 
@kwdef struct LEOSimulation
    earth::EarthProperties
    mission::Mission
    tspan::Vector{<:Real}
    dt::Real
    initstate::Vector{<:Real}
end
