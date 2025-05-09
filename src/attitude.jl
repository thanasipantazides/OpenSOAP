using StaticArrays

struct Maneuver{T}
    tf::T   # time by which to complete the maneuver
    dt::T   # duration of maneuver
    C::Matrix{T}   # attitude at end of maneuver
end

# function bb_control(x::State{<:Real},  t::Real, maneuver::Maneuver, params::LEOSimulation)

    
# end

function backout_free_torque(soln, params)
    n = length(soln["time"])
    times = soln["time"]
    states = soln["state"]
    torques = zeros(3,n)
    for (i, state) in enumerate(states)
        if i == 1
            continue
        end
        # approximate the time derivative of velocity with backwards finite difference
        dw = (state.angular_velocity .- states[i - 1].angular_velocity)./(times[i] - times[i - 1])
        I_B = params.mission.spacecraft.mass.inertia
        torques[:,i] = I_B*dw + OpenSOAP.cross(state.angular_velocity)*I_B*state.angular_velocity
    end
    return torques
end