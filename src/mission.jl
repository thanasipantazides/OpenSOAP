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

struct ReactionWheelProperties{T<:Real}
    wheel_axes::Matrix{T}
    momentum::Vector{T}     # limit, Nms
    torque::Vector{T}       # limit, Nm
    momentum_env::Polyhedron{T}
    torque_env::Polyhedron{T}
    
    ReactionWheelProperties(wheel_axes::Matrix{S}, momenta::Vector{S}, torques::Vector{S}) where S<:Real = begin
        nwheels = length(wheel_axes[1,:])
        m_normals = zeros(nwheels*(nwheels - 1),3)
        m_distances = zeros(nwheels*(nwheels - 1))
        t_normals = zeros(nwheels*(nwheels - 1),3)
        t_distances = zeros(nwheels*(nwheels - 1))
        v = 1
        for i in 1:nwheels
            for j in 1:nwheels
                if i == j
                    continue
                end
                # note: these `n` define the normals to bounding planes 
                n = cross(wheel_axes[:,i])*wheel_axes[:,j]
                n = n/LinearAlgebra.norm(n)
                
                sgns = [sign(wheel_axes[:,k]'*n) for k in 1:nwheels]
                
                m_distances[v] = (wheel_axes*(momenta.*sgns))'*n
                m_normals[v,:] = n
                t_distances[v] = (wheel_axes*(torques.*sgns))'*n
                t_normals[v,:] = n
                v += 1
            end
        end 

        mP = Polyhedron(m_normals, m_distances)
        tP = Polyhedron(t_normals, t_distances)
        return new{S}(wheel_axes, momenta, torques, mP, tP)
    end
end

@kwdef struct AttitudeProperties
    wheels::ReactionWheelProperties
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
    attitude::AttitudeProperties
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

function mission_stats(soln::Dict, target_histories::Dict, params)::Dict{String, Float64}

    # printstyled("\nDisplaying mission statistics\n", bold=true)
    period = 2*π*sqrt(norm(soln["state"][1:3,1])^3/SatelliteToolboxBase.GM_EARTH)

    # println("Mean observing fraction per orbit for target:")

    output = "Mean observing fraction per orbit for target:\n"

    total_gs_passes = 0
    suntarget_i = findfirst(target -> typeof(target) === SunTarget, params.mission.targets)
    suntarget = params.mission.targets[suntarget_i]

    for (key, val) in target_histories
        change = diff(val)
        starts = findall(i -> (i > 0), change)
        stops = findall(i -> (i < 0), change)
        
        if val[1] == 1
            starts = [1; starts]
        end
        
        if length(starts) == 0 && length(stops) == 0
            continue
        end
    
        durations = [soln["time"][stops[i]] - soln["time"][starts[i]] for i in 1:min(length(starts), length(stops))]
    
        # Save durations in the timings dictionary (or whatever you called it)
        timings[key] = durations
    
        if key != suntarget.name
            total_gs_passes += min(length(starts), length(stops))
        end
    end
    
    # After you've collected all durations, now print out the averages safely:
    for (key, durations) in timings
        if isempty(durations)
            output *= @sprintf "\t- %s:\t(no data)\n" key
        else
            output *= @sprintf "\t- %s:\t%.3f\n" key Statistics.mean(durations) / period
        end
    end
    
    in_safe = zeros(length(soln["time"]))
    in_power = zeros(length(soln["time"]))
    in_downlink = zeros(length(soln["time"]))
    in_science = zeros(length(soln["time"]))
    power_in_downlink = zeros(length(soln["time"]))
    power_in_science = zeros(length(soln["time"]))
    power_in_safe = zeros(length(soln["time"]))
    power = zeros(length(soln["time"]))
    data = zeros(length(soln["time"]))
    slew_angle = zeros(length(soln["time"]))
    mode_angle = zeros(length(soln["time"]))
    maneuver_counter = 0
    mode_counter = 0
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
            if abs(soln["state"][19,i] - soln["state"][19,i - 1]) > 1e-3
                power[i] = (soln["state"][19,i] - soln["state"][19,i - 1])/(soln["time"][i] - soln["time"][i - 1])
                power_counter += 1
                if in_science[i] == 1
                    power_in_science[i] = (soln["state"][19,i] - soln["state"][19,i - 1])/(soln["time"][i] - soln["time"][i - 1])
                end
                if in_downlink[i] == 1
                    power_in_downlink[i] = (soln["state"][19,i] - soln["state"][19,i - 1])/(soln["time"][i] - soln["time"][i - 1])
                end
                if in_safe[i] == 1
                    power_in_safe[i] = (soln["state"][19,i] - soln["state"][19,i - 1])/(soln["time"][i] - soln["time"][i - 1])
                end
            end
            if abs(soln["state"][20,i] - soln["state"][20,i - 1]) > 1e-3
                data[i] = (soln["state"][20,i] - soln["state"][20,i - 1])/(soln["time"][i] - soln["time"][i - 1])
                data_counter += 1
            end
            last_C = LinearAlgebra.reshape(soln["state"][10:18, i - 1], 3, 3)
            this_C = LinearAlgebra.reshape(soln["state"][10:18, i], 3, 3)
            angle = axisangle(last_C*this_C')[2]
            if abs(angle) > 0.1*pi/180
                slew_angle[i] = angle
                maneuver_counter += 1
            end
            if soln["state"][21, i] != soln["state"][21, i - 1]
                mode_angle[i] = angle
                mode_counter += 1
            end
        end
    end

    # printstyled("\nMean time fraction per state:\n", bold=true)
    output *= "\nMean time fraction per state:\n"

    output *= @sprintf "\tscience:\t%.3f\n"     sum(in_science)/length(in_science)
    output *= @sprintf "\tpower:\t\t%.3f\n"     sum(in_power)/length(in_power)
    output *= @sprintf "\tdownlink:\t%.3f\n"    sum(in_downlink)/length(in_downlink)
    output *= @sprintf "\tsafe:\t\t%.3f\n"      sum(in_safe)/length(in_safe)
    output *= @sprintf "Net power:\t\t%.3f W\n" sum(power)/power_counter
    output *= @sprintf "Net data:\t\t%.3f kbps\n"  sum(data)/data_counter/1e3
    output *= @sprintf "Ground passes per day:\t%.6f\n"  total_gs_passes/(soln["time"][end] - soln["time"][1])*3600*24

    # @printf "Net power:\t%0.3f W\n" sum(power)/power_counter
    # @printf "Net power:\t%0.3f W\n" Statistics.mean(power)
    output *= @sprintf "Net power in science mode:\t%0.3f W\n" sum(power_in_science)/sum.([abs.(power_in_science) .> 1e-5])[1]
    output *= @sprintf "Net power in downlink mode:\t%0.3f W\n" sum(power_in_downlink)/sum.([abs.(power_in_downlink) .> 1e-5])[1]
    output *= @sprintf "Net power in safe mode:\t%0.3f W\n" sum(power_in_safe)/sum.([abs.(power_in_safe) .> 1e-5])[1]
    output *= @sprintf "Mean slew angle:\t%0.3f º\n" 180*sum(slew_angle)/maneuver_counter/pi
    output *= @sprintf "Mean conop slew:\t%0.3f º\n" 180*sum(mode_angle)/mode_counter/pi
    output *= @sprintf "Mean mode changes per orbit:\t%0.3f\n" mode_counter*period/(soln["time"][end] - soln["time"][1])
    open(joinpath("cases","sim_results.txt"), "w") do file
        Base.write(file, output)
    end
    println("wrote mission statistics to cases/sim_results.txt")

    return_val = Dict(
        "state_mean_frac_science"=>sum(in_science)/length(in_science),
        "state_mean_frac_power"=>sum(in_power)/length(in_power),
        "state_mean_frac_downlink"=>sum(in_downlink)/length(in_downlink),
        "state_mean_frac_safe"=>sum(in_safe)/length(in_safe)
    )

    return return_val
end