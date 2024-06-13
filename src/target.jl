using SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using LinearAlgebra

@kwdef struct FrameFixedTarget
    name::String
    id::Int64
    position::Vector{<:Real}    # position (fixed) of object in :frame
    direction::Vector{<:Real}   # cone axis, in :frame
    cone::Real                  # in radians
    frame::Symbol               # ECEF or ECI or ICRS or LLA
    iers_eops
end

function position(target::FrameFixedTarget, frame::Symbol, t_s::Real)::Vector{>:Real}
    if target.frame == frame
        return target.position
    end
    
    t_jd = t_s / 3600 / 24
    if target.frame == :ECEF
        if frame == :ECEF
            return target.position
        elseif frame == :ECI
            dcm = r_ecef_to_eci(ITRF(), J2000(), t_jd, frame.iers_eops)
            return dcm*target.position
        elseif frame == :LLA
            return collect(ecef_to_geodetic(target.position))
        else
            error("unknown frame type!")
        end

    elseif target.frame == :ECI
        if frame == :ECEF
            dcm = r_eci_to_ecef(J2000(), ITRF(), t_jd, frame.iers_eops)
            return dcm*target.position
        elseif frame == :ECI
            return target.position
        elseif frame == :LLA
            dcm = r_eci_to_ecef(J2000(), ITRF(), t_jd, frame.iers_eops)
            return collect(ecef_to_geodetic(dcm*target.position))
        else
            error("unknown frame type!")
        end
    elseif target.frame == :LLA
        if frame == :ECEF
            return geodetic_to_ecef(target.position[1], target.position[2], target.position[3])
        elseif frame == :ECI
            dcm = r_ecef_to_eci(ITRF(), J2000(), t_jd, frame.iers_eops)
            return dcm*geodetic_to_ecef(target.position[1], target.position[2], target.position[3])
        elseif frame == :LLA
            return target.position
        else
            error("unknown frame type!")
        end
    else
        error("unknown frame type!")
    end
end

function Base.isequal(a::FrameFixedTarget, b::FrameFixedTarget)
    return a.id == b.id && a.name == b.name && all(a.position == b.position) && all(a.direction == b.direction) && a.cone == b.cone && a.frame == b.frame
end

function Base.hash(a::FrameFixedTarget)
    return hash(a.name, hash(a.id, hash(a.position..., hash(a.direction..., hash(a.cone)))))
end

function visibility_history(target::FrameFixedTarget, soln::Dict)
    if target.frame == :ECEF
        visibility = zeros(length(soln["time"]))
        for i in 1:length(soln["time"])
            # ECEF to ECI
            C_IF = Matrix(r_ecef_to_eci(ITRF(), J2000(), soln["time"][i]/24/3600, target.iers_eops))
                
            # sat position
            sat_I = soln["state"][1:3,i]
            
            # groundstations position
            # this is a hack---need to check that the target actually is in ECEF.
            target_I = C_IF*target.position
            visibility[i] = (sat_I - target_I)'*target_I / norm(sat_I - target_I) / norm(target_I) > cos(target.cone)
        end
        return visibility
    elseif target.frame == :ICRS
        visibility = zeros(length(soln["time"]))
        r_E = EARTH_EQUATORIAL_RADIUS

        for i in 1:length(soln["time"])
            sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(soln["time"][i]/3600/24)
            sat_I = soln["state"][1:3,i]
            visibility[i] = sun_I'*sat_I / norm(sun_I) / norm(sat_I) > -sqrt(1 - r_E^2 / norm(sat_I)^2)
        end
        return visibility
    else
        error("unusable frame time!")
    end
end