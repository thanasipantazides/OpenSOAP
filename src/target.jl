using SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using LinearAlgebra

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
function position_ecef(target::SunTarget, t_jd_s::Real)
    r_eci_to_ecef(J2000(), ITRF(), t_jd_s/3600/24, target.iers_eops)*position_eci(target, t_jd_s)
end
function position_ecef(target::GroundTarget, t_jd_s::Real)
    return Vector(geodetic_to_ecef(target.lla[1], target.lla[2], target.lla[3]))
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
    else
        error("Unsupported target type!")
    end
    return visibility
end