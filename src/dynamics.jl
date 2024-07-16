using LinearAlgebra
import SatelliteToolboxBase, SatelliteToolboxTransformations

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
    soln = Dict("time"=>time, "state"=>zeros(length(initial), length(time)))
    soln["state"][:,1] = initial
    first = true
    for ti in eachindex(time)
        if first
            first = false
            continue
        end
        temp = zeros(size(initial))
        rk4_step!(dynamics!, temp, soln["state"][:,ti - 1], time[ti], dt, params)
        soln["state"][:,ti] = temp
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
    dynamics!(k2, x + dt/2*k1, t + dt/2, params)
    dynamics!(k3, x + dt/2*k2, t + dt/2, params)
    dynamics!(k4, x + dt*k3, t + dt, params)

    x_new[1:length(x)] = x + dt/6*(k1 + 2*k2 + 2*k3 + k4)

    # Project attitude state back to SO(3)
    (U, S, V) = svd(reshape(x_new[10:18], (3,3)))
    x_new[10:18] = vec(V*diagm([1.0; 1.0; det(V)*det(U)])*U')
    
    # Project onboard data and power back to acceptable ranges
    x_new[19] = max(min(x_new[19], params.mission.spacecraft.power.capacity), 0)
    x_new[20] = max(min(x_new[20], params.mission.spacecraft.data.capacity), 0)
    
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

    z = [0;0;1]
    dx[1:3] = x[4:6]
    dx[4:6] = -params.earth.mu/(norm(x[1:3])^3) .* x[1:3] + 3*params.earth.mu*params.earth.j_2*params.earth.r^2/2/norm(x[1:3])^5*((5/norm(x[1:3])^2*(z'*x[1:3])^2 - 1)*x[1:3] - 2*(z'*x[1:3])*z)
    
    C_BI = reshape(x[10:18], (3,3))

    dx[7:9] = params.mission.spacecraft.mass.inertia \ (-cross(x[7:9])*params.mission.spacecraft.mass.inertia*x[7:9])
    dx[10:18] = vec(-cross(x[7:9])*C_BI)

    # note: power/viewing calculations duplicate code from target/visibility_history()! Factor better.
    # accumulate input power for all solar panels
    sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t/3600/24)
    sun_I_unit = sun_I / norm(sun_I)
    total_power = 0
    if can_see_sun(x, t, params)
        for panel in params.mission.spacecraft.power.solarpanels
            cosang = panel.normal'*C_BI*sun_I_unit
            if cosang > 0
                total_power += params.earth.irradiance*panel.efficiency*panel.area*cosang
                
                # println(cosang, " ", params.earth.irradiance*panel.efficiency*panel.area*cosang)
            end
        end
    end
    lookup_power_consumption = [5.96; 6.09; 8.38; 22.73]
    lookup_data_generation = [1e2; 1e3; 1e4; 4e6]
    # lookup_consumption = [5.96; 6.09; 8.38; 17.23]
    # dx[19] = total_power - params.mission.spacecraft.power.consumption
    dx[19] = total_power - lookup_power_consumption[Int64(x[21])]

    # accumulate downlink data for all ground station contacts
    total_data = 0
    for target in params.mission.targets
        if typeof(target) == GroundTarget
            if can_see_groundtarget(x, t, target, params)
                total_data += params.mission.spacecraft.data.transmit
            end
        end
    end
    # dx[20] = params.mission.spacecraft.data.production - total_data
    dx[20] = lookup_data_generation[Int64(x[21])] - total_data
end

function do_impax_conop!(x::Vector{<:Real}, t::Real, params)

    pos_I = x[1:3]
    pos_F = SatelliteToolboxTransformations.r_eci_to_ecef(J2000(), ITRF(), t/3600/24, params.mission.targets[1].iers_eops)*pos_I
    lla = SatelliteToolboxTransformations.ecef_to_geocentric(pos_F)[:]
    latitude = lla[1]
    C_BI = reshape(x[10:18], (3,3))
    new_C_BI = C_BI

    tbd_low_battery_threshold = params.mission.spacecraft.power.capacity*0.2
    tbd_low_battery_exit_threshold = params.mission.spacecraft.power.capacity*0.25
    
    if x[21] == 1 && x[19] <= tbd_low_battery_exit_threshold
        # do charging: sun-point -z-axis
        sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t/3600/24)
        new_C_BI = r_min_arc([0;0;-1], C_BI*Vector(sun_I/norm(sun_I)))'*C_BI

        # set safe state flag
        x[21] = 1
        x[10:18] = vec(new_C_BI)
        x[7:9] = [0;0;0]
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
                new_C_BI = r_min_arc([0;0;-1], C_BI*Vector(gnd_I/norm(gnd_I)))'*C_BI
                # just target the first groundstation you find
                break
            end
        end
    end
    if !txflag
        if abs(latitude) >= 35*pi/180 && abs(latitude) <= 70*pi/180
            
            # do science: nadir-point z-axis, interrupt for comms
            nadir_I = -pos_I/norm(pos_I)
            new_C_BI = r_min_arc([0;0;1], C_BI*nadir_I)'*C_BI
            # set science state flag
            x[21] = 4

        elseif can_see_sun(x, t, params)
            # do charging: sun-point -z-axis
            sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t/3600/24)
            new_C_BI = r_min_arc([0;0;-1], C_BI*Vector(sun_I/norm(sun_I)))'*C_BI
            # println([0;0;-1]'*new_C_BI'*Vector(sun_I/norm(sun_I)))
            
            # set charging state flag
            x[21] = 2
        else
            x[21] = 1
        end
    end

    if x[19] <= tbd_low_battery_threshold
        # do charging: sun-point -z-axis
        sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t/3600/24)
        new_C_BI = r_min_arc([0;0;-1], C_BI*Vector(sun_I/norm(sun_I)))'*C_BI

        # set safe state flag
        x[21] = 1
    end

    x[10:18] = vec(new_C_BI)
    x[7:9] = [0;0;0]
end