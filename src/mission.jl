import Base.@kwdef
import SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxCelestialBodies
using LinearAlgebra, Statistics, Printf

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

    println("Mean observing fraction per orbit for target:")

    total_gs_passes = 0
    suntarget_i = findfirst(target -> typeof(target) === SunTarget, params.mission.targets)
    suntarget = params.mission.targets[suntarget_i]
    println(suntarget)
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
        # println("\t- ", key, ":\t", Statistics.mean(durations)/period)
        @printf "\t- %s:\t%.3f\n" key Statistics.mean(durations)/period
        
        if key != suntarget.name
            total_gs_passes += min(length(starts), length(stops))
        end
    end

    in_safe = zeros(length(soln["time"]))
    in_power = zeros(length(soln["time"]))
    # in_sun_no_science = zeros(length(soln["time"]))
    in_downlink = zeros(length(soln["time"]))
    in_science = zeros(length(soln["time"]))
    power = zeros(length(soln["time"]))
    data = zeros(length(soln["time"]))
    power_counter = 0
    data_counter = 0

    for i in eachindex(soln["time"])
        if soln["state"][21,i] == 1
            in_safe[i] = true
        elseif soln["state"][21,i] == 2
            in_power[i] = true
        elseif soln["state"][21,i] == 3
            in_downlink[i] = true
        elseif soln["state"][21,i] == 4
            in_science[i] = true
        else
            println("got unknown state!")
        end
        if i > 1
            if soln["state"][19,i] - soln["state"][19,i - 1] != 0
                power[i] = (soln["state"][19,i] - soln["state"][19,i - 1])/(soln["time"][i] - soln["time"][i - 1])
                power_counter += 1
            end
            if soln["state"][20,i] - soln["state"][20,i - 1] != 0
                data[i] = (soln["state"][20,i] - soln["state"][20,i - 1])/(soln["time"][i] - soln["time"][i - 1])
                data_counter += 1
            end
        end
    end

    printstyled("\nMean time fraction per state:\n", bold=true)
    # println("\tscience:\t",     sum(in_science)/length(in_science))
    # println("\tpower:\t\t",       sum(in_power)/length(in_power))
    # println("\tdownlink:\t",    sum(in_downlink)/length(in_downlink))
    # println("\tsafe:\t\t",        sum(in_safe)/length(in_safe))
    # println("Mean net power:\t", sum(power)/power_counter)
    # println("Mean net data:\t", sum(data)/data_counter)

    @printf "\tscience:\t%.3f\n"      sum(in_science)/length(in_science)
    @printf "\tpower:\t\t%.3f\n"      sum(in_power)/length(in_power)
    @printf "\tdownlink:\t%.3f\n"     sum(in_downlink)/length(in_downlink)
    @printf "\tsafe:\t\t%.3f\n"       sum(in_safe)/length(in_safe)
    @printf "Net power:\t\t%.3f W\n" sum(power)/power_counter
    @printf "Net data:\t\t%.3f kbps\n"  sum(data)/data_counter/1e3
    @printf "Ground passes per day:\t%.6f\n"  total_gs_passes/(soln["time"][end] - soln["time"][1])*3600*24

    println()
end