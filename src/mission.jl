import Base.@kwdef
import SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using LinearAlgebra

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
    production::T
    transmit::T
    # antennas::Vector{Antenna}
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
    data::DataProperties
    mass::MassProperties
    # modes::Vector{Mode}
end

@kwdef struct Mission
    name::String
    # modes::Vector{Mode}
    spacecraft::SpacecraftProperties
    targets::Vector{<:AbstractTarget}
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

function can_see_sun(state::Vector{<:Real}, t_jd_s::Real, params)
    sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t_jd_s/3600/24)
    pos_I = state[1:3]
    can_see_sun = sun_I'*pos_I/norm(sun_I)/norm(pos_I) > -sqrt(1 - (params.earth.r / norm(pos_I))^2)
    return can_see_sun
end

function can_see_groundtarget(state::Vector{<:Real}, t_jd_s::Real, target::GroundTarget, params)
    gnd_I = position_eci(target, t_jd_s)
    pos_I = state[1:3]
    can_see_gnd = (pos_I - gnd_I)'*gnd_I / norm(pos_I - gnd_I) / norm(gnd_I) > cos(target.cone)
    return can_see_gnd
end

function mission_stats(soln::Dict, target_histories::Dict, params)

    printstyled("\nDisplaying mission statistics\n", bold=true)
    period = 2*π*sqrt(norm(soln["state"][1:3,1])^3/SatelliteToolboxBase.GM_EARTH)

    println("Mean observing time fraction for target:")
    for (key, val) in target_histories
        change = diff(val)
        starts = findall(i->(i > 0), change)
        stops = findall(i->(i < 0), change)
        if val[1] == 1
            starts = [1; starts]
        end 
        if length(starts) == 0 && length(stops) == 0
            # println("found no contacts for key!")
            continue
        end
        durations = [soln["time"][stops[i]] - soln["time"][starts[i]] for i in 1:min(length(starts), length(stops))]
        println("\t- ", key, ":\t", Statistics.mean(durations)/period)
    end

    in_science = zeros(length(soln["time"]))
    in_sun_no_science = zeros(length(soln["time"]))
    power = zeros(length(soln["time"]))
    power_counter = 0
    for i in 1:length(soln["time"])
        pos_I = soln["state"][1:3,i]
        pos_F = SatelliteToolboxTransformations.r_eci_to_ecef(J2000(), ITRF(), soln["time"][i]/3600/24, params.mission.targets[1].iers_eops)*pos_I
        lla = SatelliteToolboxTransformations.ecef_to_geocentric(pos_F)[:]
        latitude = lla[1]
        if abs(latitude) >= 35*pi/180 && abs(latitude) <= 70*pi/180
            in_science[i] = true
        else
            if target_histories["sun"][i] == 1.0
                in_sun_no_science[i] = true
            end
        end
        if i > 1
            if soln["state"][19,i] - soln["state"][19,i - 1] != 0
                power[i] = (soln["state"][19,i] - soln["state"][19,i - 1])/(soln["time"][i] - soln["time"][i - 1])
                power_counter += 1
            end
        end
    end
    println("Mean time fraction in science latitudes:\t", sum(in_science)/length(in_science))
    println("Mean time fraction in sun w/o science:\t", sum(in_sun_no_science)/length(in_sun_no_science))
    println("Mean power:\t", sum(power)/power_counter)
    println()
end