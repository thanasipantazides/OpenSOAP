import Base.@kwdef

@kwdef struct Component
    name::String
    id::Int64
    power::PowerProperties
    data::DataProperties
    mass::MassProperties
    procurement::ProcurementProperties
end

@kwdef mutable struct PowerProperties{T<:Real}
    value::T
    capacity::T      # should be BOL capacity. Track DoD with a state variable.
    consumption::T
    efficiency::T
    panels::Vector{SolarPanel{T}}
end

@kwdef mutable struct SolarPanel{T<:Real}
    normal::Vector{T}
    area::T
    efficiency::T
end

@kwdef mutable struct DataProperties
    value<:Real
    capacity<:Real
    production<:Real
end

@kwdef mutable struct ProcurementProperties
    vendor::String
    manufacturer::String
    part_number::String
    cost<:Real
end

@kwdef mutable struct MassProperties{T<:Real}
    inertia::Matrix{T}
    bbox::Matrix{T}
    model::String
end