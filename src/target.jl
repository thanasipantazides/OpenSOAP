using SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies, SatelliteToolboxGeomagneticField
using LinearAlgebra

const igrf_deg = 13
# global _igrf_P = Matrix{Float64}(undef, igrf_deg + 1, igrf_deg + 1)
# global _igrf_dP = similar(igrf_P)

abstract type AbstractTarget end

struct SunTarget<:AbstractTarget
    name::String
    iers_eops
end
struct GroundTarget<:AbstractTarget
    name::String
    lla::Vector{<:Real}
    direction::Vector{<:Real}
    cone::Real
    iers_eops
end
struct CelestialTarget<:AbstractTarget
    name::String
    eci::Vector{<:Real}
    iers_eops
end
struct MagneticTarget<:AbstractTarget
    name::String
    latitude_masks::Vector{Pair{Real,Real}}
    iers_eops
end

function position_eci(target::SunTarget, t_jd_s::Real)
    return SatelliteToolboxCelestialBodies.sun_position_mod(t_jd_s/3600/24)
end
function position_eci(target::CelestialTarget, t_jd_s::Real)
    return target.eci
end
function position_eci(target::GroundTarget, t_jd_s::Real)
    dcm = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/3600/24, target.iers_eops)
    return dcm*geodetic_to_ecef(target.lla[1], target.lla[2], target.lla[3])
end
function position_eci(target::MagneticTarget, position_eci::Vector{<:Real}, t_jd_s::Real)
    dcm = r_eci_to_ecef(J2000(),ITRF(), t_jd_s/3600/24, target.iers_eops)
    pos_lla = ecef_to_geodetic(dcm*position_eci)
    field_ned = igrf(t_jd_s/3600/24/365.25 - 4712, pos_lla[3], pos_lla[1], pos_lla[2], Val(:geodetic), max_degree=13, show_warnings=false)
    return Vector(dcm'*ned_to_ecef(field_ned, pos_lla..., translate=false))
end
function position_ecef(target::SunTarget, t_jd_s::Real)
    r_eci_to_ecef(J2000(), ITRF(), t_jd_s/3600/24, target.iers_eops)*position_eci(target, t_jd_s)
end
function position_ecef(target::GroundTarget, t_jd_s::Real)
    return Vector(geodetic_to_ecef(target.lla[1], target.lla[2], target.lla[3]))
end
function position_ecef(target::MagneticTarget, position_ecef::Vector{<:Real}, t_jd_s::Real)
    pos_lla = ecef_to_geodetic(position_ecef)
    field_ned = igrf(t_jd_s/3600/24/365.25 - 4712, pos_lla[3], pos_lla[1], pos_lla[2], Val(:geodetic), max_degree=13, show_warnings=false)
    return Vector(ned_to_ecef(field_ned, pos_lla..., translate=false))
end
function position_lla(position_eci::Vector{<:Real}, t_jd_s::Real, eops)
    dcm = r_eci_to_ecef(J2000(), ITRF(), t_jd_s/3600/24, eops)
    pos_lla = ecef_to_geodetic(dcm*position_eci)
    return Vector([pos_lla...])
end

function Base.isequal(a::SunTarget, b::SunTarget)
    return a.name == b.name
end
function Base.hash(a::SunTarget)
    return hash(a.name)
end
function Base.isequal(a::GroundTarget, b::GroundTarget)
    return a.name == b.name && all(a.lla == b.lla) && all(a.direction == b.direction) && a.cone == b.cone
end
function Base.hash(a::GroundTarget)
    return hash(a.name, hash(a.lla..., hash(a.direction..., hash(a.cone))))
end

function visibility(target::GroundTarget, time::Real, pos_I::AbstractVector{<:Real})
    target_I = position_eci(target, time)
    visibility = (pos_I - target_I)'*target_I / norm(pos_I - target_I) / norm(target_I) > cos(target.cone)
    return visibility
end

function visibility(target::SunTarget, time::Real, pos_I::AbstractVector{<:Real})
    r_E = EARTH_EQUATORIAL_RADIUS
    sun_I = position_eci(target, time)
    visibility = sun_I'*pos_I / norm(sun_I) / norm(pos_I) > -sqrt(max(0, 1 - r_E^2 / norm(pos_I)^2))
    return visibility
end

function visibility(target::MagneticTarget, time::Real, pos_I::AbstractVector{<:Real})
    dcm = r_eci_to_ecef(J2000(), ITRF(), time/3600/24, target.iers_eops)
    pos_lla = ecef_to_geodetic(dcm*pos_I)
    for p in target.latitude_masks
        if p.first <= pos_lla[1] <= p.second
            return true
        end
    end
    return false
end

@doc raw"""
    visibility_history(target, soln)

Compute when the target (one of `SunTarget` or `GroundTarget`) is visible from the ECI position vector captured at each time `soln["state"][1:3, i]`

"""
function visibility_history(target::T, soln::Dict) where T<:AbstractTarget
    visibility = zeros(length(soln["time"]))
    if typeof(target) == SunTarget
        r_E = EARTH_EQUATORIAL_RADIUS
        for i in 1:length(soln["time"])
            sun_I = position_eci(target, soln["time"][i])
            sat_I = soln["state"][1:3,i]
            visibility[i] = sun_I'*sat_I / norm(sun_I) / norm(sat_I) > -sqrt(1 - r_E^2 / norm(sat_I)^2)
        end
    elseif typeof(target) == GroundTarget
        # println("in visibility_history, name: ", target.name)
        for i in 1:length(soln["time"])
            target_I = position_eci(target, soln["time"][i])
            sat_I = soln["state"][1:3,i]
            visibility[i] = (sat_I - target_I)'*target_I / norm(sat_I - target_I) / norm(target_I) > cos(target.cone)
        end
    elseif typeof(target) == MagneticTarget 
        for i in 1:length(soln["time"])
            visibility[i] = true
        end
    else
        error("Unsupported target type!")
    end
    return visibility
end

function visibility_history(target::T, times::Vector{<:Real}, states::State{<:Real}) where T<: AbstractTarget
    visibility = zeros(length(times))
    if typeof(target) == SunTarget
        r_E = EARTH_EQUATORIAL_RADIUS
        for i in 1:length(times)
            sun_I = position_eci(target, times[i])
            sat_I = states[i].position
            visibility[i] = sun_I'*sat_I / norm(sun_I) / norm(sat_I) > -sqrt(1 - r_E^2 / norm(sat_I)^2)
        end
    elseif typeof(target) == GroundTarget
        # println("in visibility_history, name: ", target.name)
        for i in 1:length(times)
            target_I = position_eci(target, times[i])
            sat_I = states[i].position
            visibility[i] = (sat_I - target_I)'*target_I / norm(sat_I - target_I) / norm(target_I) > cos(target.cone)
        end
    else
        error("Unsupported target type!")
    end
    return visibility
end

