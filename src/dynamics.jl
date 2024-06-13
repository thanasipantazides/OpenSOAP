using LinearAlgebra
import SatelliteToolboxBase

"""
    integrate_system(dynamics, initial, tspan, dt, params)

Integrate the dynamical system (specified by `dynamics!()`) over `tspan`, in timesteps `dt`.
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
        # println("ti: ",ti, ", ti - 1: ", ti - 1)
        # println("outer(last):",soln["state"][:,ti-1])
        temp = zeros(size(initial))
        rk4_step!(dynamics!, temp, soln["state"][:,ti - 1], time[ti], dt, params)
        soln["state"][:,ti] = temp
        # println("outer(this):",soln["state"][:,ti])
    end

    return soln
end

function rk4_step!(dynamics!::Function, x_new::Vector{<:Real}, x::Vector{<:Real}, t::Real, dt::Real, params)
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
    x_new[19] = min(x_new[19], params.mission.spacecraft.power.capacity)

    # Update spacecraft mode
end

function dynamics_proto!(dx::Vector{<:Real}, x::Vector{<:Real}, t::Real, params)

end

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
    sun_I = SatelliteToolboxCelestialBodies.sun_position_mod(t/3600/24)
    sun_I_unit = sun_I / norm(sun_I)
    total_power = 0
    can_see_sun = sun_I'*x[1:3]/norm(sun_I)/norm(x[1:3]) > -sqrt(1 - (params.earth.r / norm(x[1:3]))^2)
    if can_see_sun
        for panel in params.mission.spacecraft.power.solarpanels
            cosang = panel.normal'*C_BI*sun_I_unit
            if cosang > 0
                total_power += params.earth.irradiance*panel.efficiency*panel.area*cosang
            end
        end
    end
    dx[19] = total_power - params.mission.spacecraft.power.consumption
end