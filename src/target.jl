using SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using LinearAlgebra

@doc raw"""
    FrameFixedTarget

A structure holding configuration information for a target attached to a reference frame.
"""
@kwdef struct FrameFixedTarget
    name::String
    id::Int64
    "fixed position of the target in frame"
    position::Vector{<:Real}    # position (fixed) of object in :frame
    "axis of target visibility cone"
    direction::Vector{<:Real}   # cone axis, in :frame
    "half-angle of target visibility cone, in radians"
    cone::Real                  # in radians
    "reference frame of "
    frame::Symbol               # ECEF or ECI or ICRS or LLA
    iers_eops
end

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

function Base.isequal(a::FrameFixedTarget, b::FrameFixedTarget)
    return a.id == b.id && a.name == b.name && all(a.position == b.position) && all(a.direction == b.direction) && a.cone == b.cone && a.frame == b.frame
end
function Base.hash(a::FrameFixedTarget)
    return hash(a.name, hash(a.id, hash(a.position..., hash(a.direction..., hash(a.cone)))))
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
        println("in visibility_history, name: ", target.name)
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

# function visibility_history(target::FrameFixedTarget, soln::Dict)
#     if target.frame == :ECEF
#         visibility = zeros(length(soln["time"]))
#         for i in 1:length(soln["time"])
#             # ECEF to ECI
#             C_IF = Matrix(r_ecef_to_eci(ITRF(), J2000(), soln["time"][i]/24/3600, target.iers_eops))
                
#             # sat position
#             sat_I = soln["state"][1:3,i]
            
#             # groundstations position
#             # this is a hack---need to check that the target actually is in ECEF.
#             target_I = C_IF*target.position
#             visibility[i] = (sat_I - target_I)'*target_I / norm(sat_I - target_I) / norm(target_I) > cos(target.cone)
#         end
#         return visibility
#     elseif target.frame == :ICRS
#         visibility = zeros(length(soln["time"]))
#         r_E = EARTH_EQUATORIAL_RADIUS

#         for i in 1:length(soln["time"])
#             sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(soln["time"][i]/3600/24)
#             sat_I = soln["state"][1:3,i]
#             visibility[i] = sun_I'*sat_I / norm(sun_I) / norm(sat_I) > -sqrt(1 - r_E^2 / norm(sat_I)^2)
#         end
#         return visibility
#     else
#         error("unusable frame time!")
#     end
# end