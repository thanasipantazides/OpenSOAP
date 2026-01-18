using LinearAlgebra
import SatelliteToolboxBase, SatelliteToolboxTransformations
using StaticArrays
using ProgressMeter

import Base: +, *

"""
    integrate_system(dynamics, initial, tspan, dt, params)

Integrate the dynamical system (specified by `dynamics!()`) over `tspan`, in timesteps `dt`.

Return a `Dict` with two `String` keys:
- `"time"`, the time history over which dynamics were integrated
- `"state"`, the state history of the system. This is a vector with the following entries:

| Index | Value             | Reference frame |
|-------|-------------------|-----------------|
| 1:3   | position of the spacecraft | ECI  |
| 4:6   | velocity of the spacecraft | ECI  |
| 7:9   | angular velocity of the spacecraft | ECI relative to Body |
| 10:18 | vectorized direction cosine matrix | ECI relative to Body |
| 19    | onboard battery level | — |
| 20    | onboard data storage level | — |
| 21    | state                      | — |

"""
function integrate_system(dynamics!::Function, initial::Vector{<:Real}, tspan::Vector{<:Real}, dt::Real, params)
    time = tspan[1]:dt:tspan[2]
    soln = Dict("time" => time, "state" => zeros(length(initial), length(time)))
    soln["state"][:, 1] = initial
    first = true
    @showprogress desc = "Total integration...\t" for ti in eachindex(time)
        if first
            first = false
            continue
        end
        temp = zeros(size(initial))
        rk4_step!(dynamics!, temp, soln["state"][:, ti-1], time[ti], dt, params)
        soln["state"][:, ti] = temp
    end

    return soln
end

function integrate_system(dynamics!::Function, initial::State{<:Real}, tspan::Vector{<:Real}, dt::Real, maneuver::Maneuver, params)
    time = tspan[1]:dt:tspan[2]
    soln = Dict("time" => time, "state" => Vector{State{<:Real}}(undef, length(time)))
    soln["state"][1] = initial

    first = true
    for ti in eachindex(time)
        if first
            first = false
            continue
        end
        # temp = State{Float64}()
        # println(ti, ": ", soln["state"][ti - 1])

        soln["state"][ti] = rk4_step_attitude!(dynamics!, soln["state"][ti-1], time[ti], maneuver, initial, params)
        # rk4_step!(dynamics!, temp, soln["state"][:,ti - 1], time[ti], dt, params)
    end

    return soln
end

@doc raw"""
    rk4_step!(dynamics!, x_new, x, t, dt, params)

Integrate the dynamical system (specified by `dynamics!()`) by one timestep `dt`, using the 4th order Runge-Kutta method.

!!! note
    This is not a pure RK4 implementation!
    The state fields 10:18 (attitude), 19 (power), and 20 (data) will be reprojected to their manifolds—states 10:18 will be projected to ``\mathrm{SO}(3)`` after each integration step, and states 19 and 20 will be clipped to the battery capacity and storage capacity, respectively (`params.mission.spacecraft.power.capacity` and `params.mission.spacecraft.data.capacity`).
"""
function rk4_step!(dynamics!::Function, x_new::Vector{<:Real}, x::Vector{<:Real}, t::Real, dt::Real, params)
    # Update spacecraft mode
    do_impax_conop!(x, t, params)

    k1 = zeros(size(x))
    k2 = zeros(size(x))
    k3 = zeros(size(x))
    k4 = zeros(size(x))
    dynamics!(k1, x, t, params)
    dynamics!(k2, x + dt / 2 * k1, t + dt / 2, params)
    dynamics!(k3, x + dt / 2 * k2, t + dt / 2, params)
    dynamics!(k4, x + dt * k3, t + dt, params)

    x_new[1:length(x)] = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # Project attitude state back to SO(3)
    (U, S, V) = svd(reshape(x_new[10:18], (3, 3)))
    x_new[10:18] = vec(V * diagm([1.0; 1.0; det(V) * det(U)]) * U')

    # Project onboard data and power back to acceptable ranges
    x_new[19] = max(min(x_new[19], params.mission.spacecraft.power.capacity), 0)
    x_new[20] = max(min(x_new[20], params.mission.spacecraft.data.capacity), 0)

end
# function rk4_step!(dynamics!::Function, x_new::Vector{<:Real}, x::Vector{<:Real}, t::Real, dt::Real, params)

function rk4_step_attitude!(dynamics!::Function, x::State{<:Real}, t::Real, maneuver::Maneuver, initial::State{<:Real}, params)
    # k1 = State{Float64}()
    # k2 = State{Float64}()
    # k3 = State{Float64}()
    # k4 = State{Float64}()

    dt = params.dt
    k1 = dynamics!(x, t, maneuver, initial, params)
    k2 = dynamics!(x + dt / 2 * k1, t + dt / 2, maneuver, initial, params)
    k3 = dynamics!(x + dt / 2 * k2, t + dt / 2, maneuver, initial, params)
    k4 = dynamics!(x + dt * k3, t + dt, maneuver, initial, params)

    x_new = x + params.dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

    # Project attitude state back to SO(3)
    (U, S, V) = svd(x_new.attitude)
    x_new.attitude = U * diagm([1.0; 1.0; det(V) * det(U)]) * V'

    # Project onboard data and power back to acceptable ranges
    x_new.battery = max(min(x_new.battery, params.mission.spacecraft.power.capacity), 0)
    x_new.storage = max(min(x_new.storage, params.mission.spacecraft.data.capacity), 0)

    return x_new
end

function dynamics_proto!(dx::Vector{<:Real}, x::Vector{<:Real}, t::Real, params)

end

@doc raw"""
    dynamics_orbit!(dx, x, t, params)

Implements dynamical equations for 20-DOF model of spacecraft: orbit dynamics, attitude dynamics, power, and data.
"""
function dynamics_orbit!(dx::Vector{<:Real}, x::Vector{<:Real}, t::Real, params)
    # state vector is:
    #   position            3
    #   velocity            3
    #   angular velocity    3
    #   attitude            9
    #   battery level       1
    #   disk storage        1
    #   operating state     1

    z = [0; 0; 1]
    dx[1:3] = x[4:6]
    dx[4:6] = -params.earth.mu / (norm(x[1:3])^3) .* x[1:3] + 3 * params.earth.mu * params.earth.j_2 * params.earth.r^2 / 2 / norm(x[1:3])^5 * ((5 / norm(x[1:3])^2 * (z' * x[1:3])^2 - 1) * x[1:3] - 2 * (z' * x[1:3]) * z)

    C_BI = reshape(x[10:18], (3, 3))

    dx[7:9] = params.mission.spacecraft.mass.inertia \ (-cross(x[7:9]) * params.mission.spacecraft.mass.inertia * x[7:9])
    dx[10:18] = vec(-cross(x[7:9]) * C_BI)

    # note: power/viewing calculations duplicate code from target/visibility_history()! Factor better.
    # accumulate input power for all solar panels
    sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t / 3600 / 24)
    sun_I_unit = sun_I / norm(sun_I)
    total_power = 0
    if can_see_sun(x, t, params)
        for panel in params.mission.spacecraft.power.solarpanels
            cosang = panel.normal' * C_BI * sun_I_unit
            if cosang > 0
                total_power += params.earth.irradiance * panel.efficiency * panel.area * cosang

                # println(cosang, " ", params.earth.irradiance*panel.efficiency*panel.area*cosang)
            end
        end
    end
    # lookup_power_consumption = [5.96; 6.09; 8.38; 22.73]
    lookup_data_generation = [3.54e4; 3.54e4; 3.54e4; 3.54e4]
    lookup_power_consumption = [5.96; 6.09; 8.38; 17.23]
    # dx[19] = total_power - params.mission.spacecraft.power.consumption
    dx[19] = total_power - lookup_power_consumption[Int64(x[21])]

    # accumulate downlink data for all ground station contacts
    total_data = 0
    for target in params.mission.targets
        if typeof(target) == GroundTarget
            if can_see_groundtarget(x, t, target, params)
                total_data += params.mission.spacecraft.data.transmit
                break
            end
        end
    end
    # dx[20] = params.mission.spacecraft.data.production - total_data
    dx[20] = lookup_data_generation[Int64(x[21])] - total_data
end

function do_impax_conop!(x::Vector{<:Real}, t::Real, params)

    pos_I = x[1:3]
    pos_F = SatelliteToolboxTransformations.r_eci_to_ecef(J2000(), ITRF(), t / 3600 / 24, params.mission.targets[1].iers_eops) * pos_I
    lla = SatelliteToolboxTransformations.ecef_to_geocentric(pos_F)[:]
    latitude = lla[1]
    C_BI = reshape(x[10:18], (3, 3))
    new_C_BI = C_BI

    tbd_low_battery_threshold = params.mission.spacecraft.power.capacity * 0.2
    tbd_low_battery_exit_threshold = params.mission.spacecraft.power.capacity * 0.25

    if x[21] == 1 && x[19] <= tbd_low_battery_exit_threshold
        # do charging: sun-point -z-axis
        sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t / 3600 / 24)
        new_C_BI = r_min_arc([0; 0; -1], C_BI * Vector(sun_I / norm(sun_I)))' * C_BI

        # set safe state flag
        x[21] = 1
        x[10:18] = vec(new_C_BI)
        x[7:9] = [0; 0; 0]
        return
    end
    txflag = false
    for target in params.mission.targets
        if typeof(target) == GroundTarget
            if can_see_groundtarget(x, t, target, params)
                txflag = true
                # set downlink state flag
                x[21] = 3
                gnd_I = pos_I - position_eci(target, t)
                new_C_BI = r_min_arc([0; 0; -1], C_BI * Vector(gnd_I / norm(gnd_I)))' * C_BI
                # just target the first groundstation you find
                break
            end
        end
    end
    if !txflag
        if abs(latitude) >= 35 * pi / 180 && abs(latitude) <= 70 * pi / 180

            # do science: nadir-point z-axis, interrupt for comms
            nadir_I = -pos_I / norm(pos_I)
            new_C_BI = r_min_arc([0; 0; 1], C_BI * nadir_I)' * C_BI
            # set science state flag
            x[21] = 4

        elseif can_see_sun(x, t, params)
            # do charging: sun-point -z-axis
            sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t / 3600 / 24)
            new_C_BI = r_min_arc([0; 0; -1], C_BI * Vector(sun_I / norm(sun_I)))' * C_BI
            # println([0;0;-1]'*new_C_BI'*Vector(sun_I/norm(sun_I)))

            # set charging state flag
            x[21] = 2
        else
            x[21] = 1
        end
    end

    if x[19] <= tbd_low_battery_threshold
        # do charging: sun-point -z-axis
        sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t / 3600 / 24)
        new_C_BI = r_min_arc([0; 0; -1], C_BI * Vector(sun_I / norm(sun_I)))' * C_BI

        # set safe state flag
        x[21] = 1
    end

    x[10:18] = vec(new_C_BI)
    x[7:9] = [0; 0; 0]
end

@doc raw"""
    Unimplemented!

(effectively). Will soon be embedded in a direct collocation-based optimal control algorithm.
"""
function dynamics_attitude!(x::State{<:Real}, t::Real, maneuver::Maneuver, initial::State{<:Real}, params::LEOSimulation)
    # todo: should extend to Vector{Maneuver}, ordered in time.

    # u_B = params.
    # I_B = params.mission.spacecraft.mass.inertia
    # C_BfI = maneuver.C
    # C_BiI = x.attitude

    # C_BfBi = C_BfI * C_BiI'
    # # C_BfBi = C_BiI * C_BfI'
    # ax = uncross(C_BfBi - C_BfBi')
    # # ax = rand(3)
    # # ax = ax/norm(ax)
    # torque_B = zeros(3)

    # if maneuver.tf - maneuver.dt < t < maneuver.tf - maneuver.dt / 2
    #     # println("maneuver: accelerating")
    #     # first half of maneuver: accelerating
    #     # println("axis: ", ax)
    #     torque_B = project(ax, params.mission.spacecraft.attitude.wheels.torque_env)
    #     # some thoughts:
    #         # Maneuver currently just prescribes start and end attitudes.
    #         # Add a function "control_axis" that takes a maneuver and params, and spits out rotation axis.
    #         # Explore: max amplitude bang-bang maneuvers for random start and end attitudes, separated by 180º. Get lower time bound.
    #         # need a decent inertia estimate to make this realistic.
    #         # syntax like `do_maneuver(maneuver::Maneuver, control::Controller)` might be expressive.
    # elseif maneuver.tf - maneuver.dt / 2 < t < maneuver.tf
    #     # println("maneuver: braking")
    #     # second half of maneuver: braking
    #     torque_B = -project(ax, params.mission.spacecraft.attitude.wheels.torque_env)
    # end

    # if t > maneuver.tf
    #     # if any(sum(x.attitude'*maneuver.C) != 3)
    #     #     println("maneuver done, failed to converge! Off by ", x.attitude .- maneuver.C)
    #     # end
    # end

    # # torque_B = 0 .*torque_B
    # # println(ax)
    # # out of maneuver
    # dx = State{Float64}()
    # dx.attitude = -cross(x.angular_velocity)*x.attitude
    # dx.angular_velocity = I_B \ (torque_B - cross(x.angular_velocity)*I_B*x.angular_velocity)

    # return dx

    I_B = params.mission.spacecraft.mass.inertia
    C_BfI = maneuver.C          # final attitude
    C_BiI = initial.attitude    # initial attitude
    C_BI = x.attitude    # current attitude

    # C_BiB = C_BiI * C_BI'
    # C_BfBi = C_BfI * C_BiI'
    # C_BfBi = C_BfI * C_BiI'
    # C_BfBi = C_BiI * C_BfI'
    axf, angf = axisangle(C_BfI)
    axi, angi = axisangle(C_BiI)
    ax, ang = axisangle(C_BI)
    # ax = rand(3)
    # ax = ax/norm(ax)
    torque_B = zeros(3)

    crossover = (angi + angf) / 2
    inmaneuver = (t >= maneuver.tf - maneuver.dt) && (t < maneuver.tf)
    # stopping_angvel = norm(x.angular_velocity) <= 1e-3
    # stopping_att = (ang - angf)^2 <= 1e-2 && (norm(ax - axf)) <= 1e-2

    if inmaneuver # restrict maneuvering to allotted time
        # if !stopping_angvel && !stopping_att
        if ang <= crossover
            torque_B = project(ax, params.mission.spacecraft.attitude.wheels.torque_env)
        elseif ang >= crossover
            torque_B = -project(ax, params.mission.spacecraft.attitude.wheels.torque_env)
        end
        # end
    end

    # torque_B = 0 .*torque_B
    # println(ax)
    # out of maneuver
    dx = State{Float64}()
    dx.attitude = -cross(x.angular_velocity) * x.attitude
    dx.angular_velocity = I_B \ (torque_B - cross(x.angular_velocity) * I_B * x.angular_velocity)

    return dx
end

function state_to_dict(time::Vector{S}, states::Vector{State{S}}) where S<:Real
    if length(states) != length(time) || length(states) == 0
        error("Arguments must have same, nonzero length")
    end
    soln = Dict("time" => time, "state" => zeros(length(states[1]), length(time)))
    for (k, state) in enumerate(states)
        soln["state"][:, k] = [
            state.position;
            state.velocity;
            state.angular_velocity;
            state.attitude[:];
            state.battery;
            state.storage;
            state.mode
        ]
    end
    return soln
end

function orbit_dynamics!(x_dot::State{S}, x::State{S}, t::S, dt::S, params) where {S<:Real}
    z = @SVector [0; 0; 1]
    x_dot.position = x.velocity
    x_dot.velocity = -params.earth.mu / (norm(x.position)^3) .* x.position + 3 * params.earth.mu * params.earth.j_2 * params.earth.r^2 / 2 / norm(x.position)^5 * ((5 / norm(x.position)^2 * (z' * x.position)^2 - 1) * x.position - 2 * (z' * x.position) * z)
end

function simulate_orbit!(sim::LEOSimulation, sim_config::Union{Nothing,Dict{String,String}}, times::Vector{S}, states::Vector{State{S}}, targets::Matrix{S}) where {S<:Real}

    if length(states) != length(times)
        states = cat(states, Vector{State{S}}(undef, length(times) - length(states)))
    end

    dt = 0.0
    k1 = State{S}()
    k2 = State{S}()
    k3 = State{S}()
    k4 = State{S}()
    @showprogress desc = "Propagating orbit...\t" for k in eachindex(times)
        if k == 1
            for (s, target) in enumerate(sim.mission.targets)
                targets[s, k] = visibility(target, times[k], states[k].position)
            end
            continue
        end

        dt = times[k] - times[k-1]
        orbit_dynamics!(k1, states[k-1], times[k-1], dt, sim)
        orbit_dynamics!(k2, states[k-1] + dt / 2 * k1, times[k-1] + dt / 2, dt, sim)
        orbit_dynamics!(k3, states[k-1] + dt / 2 * k2, times[k-1] + dt / 2, dt, sim)
        orbit_dynamics!(k4, states[k-1] + dt * k3, times[k-1] + dt, dt, sim)

        states[k] = states[k-1] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)

        for (s, target) in enumerate(sim.mission.targets)
            targets[s, k] = visibility(target, times[k], states[k].position)
        end
    end
end

function simulate_conop!(sim::LEOSimulation, sim_config::Union{Nothing,Dict{String,String}}, times::Vector{S}, states::Vector{State{S}}, targets::Matrix{S}, target_choice::Vector{Int64}, reference_directions::Matrix{S}) where {S<:Real}
    n_orbit = length(times)
    k = 2
    active_target = 0
    # states[k].mode will track the mode

    # usage:
    #   states[k].mode will give you a body vector
    #   reference_direction[:, k] will give you an inertial direction to point it
    #   target_selection[k] just might be helpful. Esp. visualization.

    transitions = hcat(zeros(length(sim.mission.targets)), diff(targets, dims=2))

    t_sun = findfirst(x -> isa(x, SunTarget), sim.mission.targets)
    t_mag = findfirst(x -> isa(x, MagneticTarget), sim.mission.targets)
    t_gs = findall(x -> isa(x, GroundTarget), sim.mission.targets)

    progbar = Progress(n_orbit, desc="Running target selection...\t")
    while k < length(times)
        update!(progbar, k)
        diffmask = targets[:, k] .> targets[:, k-1]

        if any(targets[t_gs, k] .!= 0)     # any groundstation is visible
            for (p, target) in enumerate(sim.mission.targets)
                if targets[p, k] == 1 && isa(target, GroundTarget)
                    active_target = p
                    break
                end
            end
            # target_set = findall(x->x == 1, targets[:,k])
            # active_target = findfirst(x->isa(x, GroundTarget), sim.mission.targets[target_set])
            # println(active_target)
            # active_target = findfirst(x->x == 1, targets[:, k]) # default to the first groundstation
            endtime_rel = findfirst(x -> x == 0, targets[active_target, k:end])

            if isnothing(endtime_rel)
                # sim ends before target goes out of view.
                # save data for this target
                k += 1
                continue
            end
            endtime = endtime_rel + k - 1
            duration = endtime - k

            # if there is another target with a longer duration that starts during this duration, set it as active target and continue.
            options = Vector{Pair{Int64,Int64}}(undef, 0) # record the (duration, index) for alternative targets
            for (p, target) in enumerate(sim.mission.targets)
                if isa(target, GroundTarget)
                    if p == active_target
                        continue
                    end
                    if any(transitions[p, k:endtime] .== 1)
                        # this other target p starts a window during active_target's window.
                        #   find out if that window has better duration.
                        other_endtime_rel = findfirst(x -> x == 0, targets[p, k:end])
                        if ~isnothing(other_endtime_rel)
                            other_endtime = other_endtime_rel + k - 1
                            other_duration = other_endtime - k
                            if other_duration > duration
                                push!(options, Pair(other_duration, p))
                            end
                        end
                    end
                end
            end
            if length(options) > 0
                best_duration, best_index = findmax(x -> x.first, options)
                if best_duration <= 0
                    throw("not incrementing (best)!")
                end
                active_target = best_index # update the target
                [target_choice[j] = active_target for j in k:best_duration+k-1]
                [states[j].mode = downlink::Modes for j in k:best_duration+k-1]
                # reference_directions[:,k:best_duration + k - 1] = ...
                [reference_directions[:, j] = position_eci(sim.mission.targets[active_target], times[j]) - states[j].position for j in k:best_duration+k-1]
                k += best_duration # skip to the end of the window for this target
                continue
            else
                if duration <= 0
                    println(targets[active_target, k:k+3])
                    throw("not incrementing (base)!")
                end
                [target_choice[j] = active_target for j in k:duration+k-1]
                [states[j].mode = downlink::Modes for j in k:duration+k-1]
                # reference_directions[:,k:best_duration + k - 1] = ...
                [reference_directions[:, j] = position_eci(sim.mission.targets[active_target], times[j]) - states[j].position for j in k:duration+k-1]
                k += duration
                continue
            end

        elseif targets[t_mag, k] != 0   # in science latitude
            # do science pointing
            active_target = t_mag
            target_choice[k] = active_target
            states[k].mode = science::Modes
            ref_dir = position_eci(sim.mission.targets[active_target], Vector(states[k].position), times[k])
            # if we are in southern magnetic hemisphere i.e. field lines point away from earth, want to actually point antiparallel to the lines.
            if ref_dir'*states[k].position > 0
                ref_dir = -ref_dir
            end
            reference_directions[:, k] = ref_dir

        elseif targets[t_sun, k] != 0    # sun is visible
            # do sun pointing
            active_target = t_sun
            target_choice[k] = active_target
            states[k].mode = charging::Modes
            reference_directions[:, k] = position_eci(sim.mission.targets[active_target], times[k])

        else
            states[k].mode = idle::Modes
            if k > 1
                # remain at the previous reference direction
                reference_directions[:, k] = reference_directions[:, k - 1]
            else
                reference_directions[:, k] = [0.0,0.0,1.0]
            end
        end
        k += 1
    end

    finish!(progbar)
end

# function assign_attitude_reference!(sim::LEOSimulation, sim_config::Union{Nothing,Dict{String,String}}, times::Vector{S}, states::Vector{State{S}}, targets::Matrix{S}, target_choice::Vector{Int64}, reference_directions::Matrix{S}, reference_attitude::Dict{Int64, Tuple{Vector{S}, Vector{State{S}}}}) where {S<:Real}
#     n_time = length(times)
#     coarse_dt = sim.dt

#     minute_k = Int(round(60 / coarse_dt))

#     mode_val = [state.mode for state in states]
#     mode_change_k = findall(x -> x != 0, diff(target_choice))
#     priority = Dict(
#         0 => downlink::Modes,
#         1 => science::Modes,
#         2 => charging::Modes,
#         3 => idle::Modes,
#         9 => safe::Modes,
#     )

#     @showprogress desc = "Assigning reference knot points...\t" for k in eachindex(times)
#         if k == length(times)
#             break
#         end
#         start_maneuver_k = 1
#         stop_maneuver_k = 1
#         # select mode changes
#         if target_choice[k+1] != target_choice[k]
#             if states[k+1].mode == Int(downlink::Modes)
#                 # transitioning to a groundstation
#                 if states[k].mode == Int(downlink::Modes)
#                     # split the difference in the maneuver
#                     halfminute_k = minute_k ÷ 2
#                     start_maneuver_k = max(1, k - halfminute_k)
#                     stop_maneuver_k = min(n_time, k + halfminute_k)
#                 else
#                     # consume the current mode
#                     start_maneuver_k = max(1, k - minute_k + 1)
#                     stop_maneuver_k = k
#                 end
#             elseif states[k+1].mode == Int(science::Modes)
#                 if states[k].mode == Int(downlink::Modes)
#                     # priority downlink
#                     start_maneuver_k = k
#                     stop_maneuver_k = min(n_time, k + minute_k)
#                 else
#                     # priority science
#                     start_maneuver_k = max(1, k - minute_k)
#                     stop_maneuver_k = k
#                 end
#             elseif states[k+1].mode == Int(charging::Modes)
#                 if states[k].mode == Int(downlink::Modes) || states[k].mode == Int(science::Modes)
#                     # priority downlink or science
#                     start_maneuver_k = k
#                     stop_maneuver_k = min(n_time, k + minute_k)
#                 else
#                     # priority power
#                     start_maneuver_k = max(1, k - minute_k)
#                     stop_maneuver_k = k
#                 end
#             else
#                 continue
#             end
#             [states[j].mode = Int(pointing::Modes) for j in start_maneuver_k:stop_maneuver_k]
#         end
#     end

#     # now, cover the entire state history again, selecting regions which are 
#     #   1. pointing mode (interpolate start to end)
#     #   2. science mode (track w/ interpolation)
#     #   3. downlink mode (track w/ interpolation)
#     #
#     # and generate attitude reference for all of them. Bridge these dynamic pointing 
#     # modes by assuming constant inertial pointing for
#     #   a. charging mode
#     #   b. idle mode

#     # fine_dt = 0.22 / sqrt(norm(Iinv * m_max))
#     fine_dt = 0.1
#     dt_mul = Int(round(coarse_dt/fine_dt))

#     @warn "using hardcoded body vectors for slew! Replace later."
#     dir_for_mode_B = Dict(
#         downlink::Modes=>[0.0;0.0;1.0],
#         science::Modes=>[0.0;0.0;1.0],
#         charging::Modes=>[0.0;0.0;-1.0],
#     )

#     snap_initial = false # for future use, snap initial attitude value
#     C_BI0 = states[1].attitude
#     dynamic_attitude_modes = [science::Modes, downlink::Modes, pointing::Modes]
    
#     # store fine time sequence (finetime, Vector{State} pair) in a dictionary keyed by the maneuver start time.
#     maneuver_reference = Dict{Int64, Tuple{Vector{S}, Vector{State{S}}}}()
#     progbar = Progress(n_time, desc="Interpolating reference attitude...\t")
#     k = 2
#     # C_BI0 = states[k].attitude
#     while k <= n_time
#         update!(progbar, k)
        
#         # find start point of maneuver (1, 2, 3 in above list)
#         #   find stop point of that maneuver (when target changes)
#         # rotinterp the maneuver endpoints (in attitude) with fine timestep to get attitude reference
#         #   then log map the attitude diff to get angular rate reference
#         # the other state elements can be directly interpolated (vectors).
#         # push the (finetime, State reference) pair into the dict with key k.
         
#         if Modes(states[k].mode) == pointing::Modes
#             this_mode_start_k = k
#             next_mode_start_k = findnext(x->x.mode != states[k].mode, states, k)
#             if isnothing(next_mode_start_k) break end
#             this_mode_final_k = next_mode_start_k - 1
            
#             # C_BI0 = I(3)
#             # the attitude once this maneuver is done:
#             println(Modes(states[next_mode_start_k].mode))
#             R_man = r_min_arc(reference_directions[:, next_mode_start_k], C_BI0'*dir_for_mode_B[Modes(states[next_mode_start_k].mode)])
#             println(tr(R_man'*R_man) + sum(R_man'*R_man - R_man*R_man'))
#             println(det(R_man))
#             C_BInext = R_man*C_BI0
            
#             # n_fine = dt_mul*(this_mode_final_k - k)
#             n_fine = 1*(this_mode_final_k - k)
#             time_fine = times[k]:fine_dt:times[this_mode_final_k]
            
#             this_state = states[k]
#             this_state.attitude = C_BI0
#             this_mode_final_state = states[this_mode_final_k]
#             this_mode_final_state.attitude = C_BInext
#             states_fine = interp(this_state, this_mode_final_state, n_fine)
            
#             reference_attitude[k] = (time_fine, states_fine)
#             k = next_mode_start_k
#             C_BI0 = C_BInext
            
#         elseif Modes(states[k].mode) in dynamic_attitude_modes
#             next_mode_start_k = findnext(x->x.mode != states[k].mode, states, k)
#             this_mode_final_k = next_mode_start_k - 1
            
            
#             # for j in k:this_mode_final_k-1
#             #     # can use k here rather than j, should be same for all j:
#             #     R_man = r_min_arc(C_BI0*reference_directions[:, j], dir_for_mode_B[Modes(states[j].mode)])
#             #     C_Bnext = C_BI0*R_man'
                
#             #     n_fine = dt_mul
#             #     time_fine = times[j]:fine_dt:times[j+1]
                
#             #     this_state = states[j]
#             #     this_state.attitude = C_BI0
#             #     this_mode_final_state = states[j+1]
#             #     this_mode_final_state.attitude = C_Bnext
#             #     states_fine = interp(this_state, this_mode_final_state, n_fine)
                
#             #     reference_attitude[j] = (time_fine, states_fine)
                
#             #     C_BI0 = C_Bnext
#             # end
#             # k = next_mode_start_k
#         end
        
#         # # entering a maneuver
#         # if Modes(states[k].mode) == pointing::Modes && Modes(states[k - 1].mode) != pointing::Modes
#         #     # find the endpoint of the maneuver. will act based on what is happening at the endpoint.
#         #     slew_stop_k = findnext(x->x.mode != Int(pointing::Modes), states, k)

#         #     if Modes(states[slew_stop_k].mode) == idle::Modes
#         #         @warn "found a slew into idle mode?"
#         #         k = slew_stop_k
#         #         continue
#         #     end

#         #     if Modes(states[slew_stop_k].mode) in dynamic_attitude_modes    
#         #         # println(k, ": ", Modes(states[k].mode), " -> ", Modes(states[stop_k].mode))
                
#         #         # slew to the start attitude of the next mode
#         #         C_AB = r_min_arc(reference_directions[:, k], dir_for_mode_B[Modes(states[slew_stop_k].mode)])
#         #         # C_AI = C_BI0*C_AB'
#         #         C_AI = C_AB

#         #         n_fine = dt_mul*(slew_stop_k - 1 - k)
#         #         time_fine = times[k]:fine_dt:times[slew_stop_k - 1]
#         #         s0 = states[k]
#         #         s0.attitude = C_BI0
#         #         sf = states[slew_stop_k - 1]
#         #         sf.attitude = C_AI
#         #         states_fine = interp(s0, sf, n_fine)
#         #         # states_fine = interp(sf, s0, n_fine)
#         #         reference_attitude[k] = (time_fine, states_fine)

#         #         # todo: now, also slew through the maneuver.
#         #         this_mode = states[slew_stop_k].mode
#         #         next_mode_k = findnext(x->x.mode != this_mode, states, slew_stop_k)
#         #         for s in dt_mul*(next_mode_k - slew_stop_k)
#         #             # todo: step through interpolate.
#         #         end
#         #         k = slew_stop_k
#         #         C_BI0 = C_AI


#         #         # note on r_min_arc usage:
#         #         #   pass in two vectors in the same frame (the initial frame). Then multiply the result by the initial frame, and get the final frame.  

#         #     # elseif Modes(states[k].mode) in dynamic_attitude_modes
#         #         # min arc slew between reference directions at each timestep during this maneuver (before next mode change) (making knot points)
#         #         #   then rotinterp the knot points.
#         #     # elseif Modes(states[k].mode) == charging::Modes
#         #         # attitude is min arc slew between last mode and mean sun direction of this mode
#         #     else
#         #         # hold last 
                
#         #     end
#         # end
#         k += 1
#     end
# end

function assign_attitude_reference!(sim::LEOSimulation, 
        sim_config::Union{Nothing,Dict{String,String}}, 
        times::Vector{S}, 
        states::Vector{State{S}},
        targets::Matrix{S}, 
        target_choice::Vector{Int64}, 
        reference_directions::Matrix{S}, 
        reference_attitude::Dict{Int64, Tuple{Vector{S}, Vector{State{S}}}}
    ) where {S<:Real}

    # for now, ignore fine time. reference attitude will be interpolated at each state value.
    
    minute_k = Int(round(60 / sim.dt))
    nt = length(times)
    
    track_modes = [downlink::Modes, science::Modes, charging::Modes]
    
    @warn "using hardcoded body vectors for slew! Replace later."
    dir_for_mode_B = Dict(
        downlink::Modes=>[0.0;0.0;1.0],
        science::Modes=>[0.0;0.0;1.0],
        charging::Modes=>[0.0;0.0;-1.0],
    )
    
    @showprogress desc = "Assigning reference knot points...\t" for k in eachindex(times)
        if k == length(times)
            break
        end
        start_maneuver_k = 1
        stop_maneuver_k = 1
        # select mode changes
        if target_choice[k+1] != target_choice[k]
            if states[k+1].mode == downlink::Modes
                # transitioning to a groundstation
                if states[k].mode == downlink::Modes
                    # split the difference in the maneuver
                    halfminute_k = minute_k ÷ 2
                    start_maneuver_k = max(1, k - halfminute_k)
                    stop_maneuver_k = min(nt, k + halfminute_k)
                else
                    # consume the current mode
                    start_maneuver_k = max(1, k - minute_k + 1)
                    stop_maneuver_k = k
                end
            elseif states[k+1].mode == science::Modes
                if states[k].mode == downlink::Modes
                    # priority downlink
                    start_maneuver_k = k
                    stop_maneuver_k = min(nt, k + minute_k)
                else
                    # priority science
                    start_maneuver_k = max(1, k - minute_k)
                    stop_maneuver_k = k
                end
            elseif states[k+1].mode == charging::Modes
                if states[k].mode == downlink::Modes || states[k].mode == science::Modes
                    # priority downlink or science
                    start_maneuver_k = k
                    stop_maneuver_k = min(nt, k + minute_k)
                else
                    # priority power
                    start_maneuver_k = max(1, k - minute_k)
                    stop_maneuver_k = k
                end
            else
                continue
            end
            [states[j].mode = pointing::Modes for j in start_maneuver_k:stop_maneuver_k]
        end
    end
    
    flat_modes = [states[k].mode for k in eachindex(times)]
    
    k = 1
    first = true
    while k < nt
        if first
            # todo
            first = false
            k += 1
            continue
        end
        
        # for every pointing, just interpolate between first and last attitude. 
        # for downlink and science, min arc slew between timesteps. Should be able to do that everywhere, actually.
        if states[k].mode == pointing::Modes
            C_BInext = I(3)
            k_next = findnext(flat_modes .!= flat_modes[k], k)
            if isnothing(k_next) # sim ends during this maneuver
                @warn "odd, finishing sim while pointing."
                C_BInext = states[k - 1].attitude
                k += 1
                continue
            end
            ref_norm = reference_directions[:,k_next] / norm(reference_directions[:,k_next])
            # C_BInext = states[k-1].attitude*r_min_arc(states[k-1].attitude*dir_for_mode_B[states[k_next].mode], ref_norm)
            # C_BInext = r_min_arc(states[k-1].attitude*dir_for_mode_B[states[k_next].mode], ref_norm)*states[k-1].attitude
            C_BInext = states[k-1].attitude*r_min_arc(dir_for_mode_B[states[k_next].mode], states[k-1].attitude'*ref_norm)
            
            seq = rotinterp(states[k-1].attitude, C_BInext, k_next - k)
            # todo: this is an apparent bug in rotinterp: result is relative to left arg, 
            # i.e. non-identity initial attitudes will be interpolated relative to the first attitude.
            for j in k:k_next - 1
                states[j].attitude = seq[:,:,j - k + 1]
                # states[j].attitude = states[k-1].attitude*seq[:,:,j - k + 1]
            end
            
            k = k_next
        elseif states[k].mode in track_modes
            ref_norm = reference_directions[:,k] / norm(reference_directions[:,k])
            # states[k].attitude = r_min_arc(states[k-1].attitude*dir_for_mode_B[states[k].mode], ref_norm)
            # states[k].attitude = states[k-1].attitude*r_min_arc(states[k-1].attitude'*ref_norm, dir_for_mode_B[states[k].mode])'
            states[k].attitude = states[k-1].attitude*r_min_arc(dir_for_mode_B[states[k].mode], states[k-1].attitude'*ref_norm)
            k += 1
        else
            states[k].attitude = states[k-1].attitude
            k += 1
        end
    end
end

function simulate_power_data!(sim::LEOSimulation, 
            sim_config::Union{Nothing,Dict{String,String}}, 
            times::Vector{S}, 
            states::Vector{State{S}},
            targets::Matrix{S}, 
            target_choice::Vector{Int64}
    ) where S<:Real

    for k in eachindex(times)
        if k == 1 
            continue 
        end
        
        # check modal power consumption, data production.
        # check visibility (eclipse + viewing) of sun for power production
        # check visibility (eclipse + viewing) of groundstations for data consumption
        #   add arg to record data sent to each groundstation
        
        
    end
end

function simulate(sim::LEOSimulation, sim_config::Union{Nothing,Dict{String,String}})

    times = Vector{Float64}(sim.tspan[1]:sim.dt:sim.tspan[2])
    n_orbit = length(times)

    # allocate:
    target_visibilities = Matrix{Float64}(undef, length(sim.mission.targets), n_orbit)
    target_choice = zeros(Int64, n_orbit) .- 4
    reference_directions = Matrix{Float64}(undef, 3, n_orbit) # pointing direction at each time
    reference_attitude = Array{Float64,3}(undef, 3, 3, n_orbit) # attitude reference at each time
    fine_reference_attitude = Dict{Int64, Tuple{Vector{Float64}, Vector{State{Float64}}}}()
    states = Vector{State{Float64}}(undef, n_orbit)
    # initial condition:
    states[1] = State{Float64}(sim)
    
    println("running.")

    println("simulating ", ((sim.tspan[2] - sim.tspan[1]) / 24 / 3600), " days")
    # propagate orbit
    simulate_orbit!(sim, sim_config, times, states, target_visibilities)
    # assign target based on conop
    simulate_conop!(sim, sim_config, times, states, target_visibilities, target_choice, reference_directions)
    # rectify attitude reference with pointings
    assign_attitude_reference!(sim, sim_config, times, states, target_visibilities, target_choice, reference_directions, fine_reference_attitude)

    simulate_power_data!(sim, sim_config, times, states, target_visibilities, target_choice)
    # fudging this for later usage:
    # [state.attitude = diagm([1.0; 1.0; 1.0]) for state in states]

    # run conop logic for targeting
    #   note: this is a priori---pointing logic can be interrupted by tumbling or low power triggers
    # flow:
    #   1. visibility masks and position -> mode
    #   2. mode -> target
    #   3. target and last target and -> attitude reference
    #   4. (justification -> attitude reference timing)
    #
    # attitude control
    #   along a fine grid, with steps defined by the agilitoid
    #   heuristic_control_stepsize = 0.22/sqrt(norm(Iinv*m_max))
    #   before computing, warn that "performing n maneuvers, each takes ~150 seconds to compute"
    #   ... and show progress bar (`using ProgressMeter`)
    #
    # power, data, flows
    #
    # return state trajectory, time, target visibility mask
    #
    return times, states, target_visibilities, target_choice, reference_directions, fine_reference_attitude
end
