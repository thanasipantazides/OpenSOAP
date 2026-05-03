using LinearAlgebra
import Dates
import SatelliteToolboxTransformations, SatelliteToolboxBase, SatelliteToolboxCelestialBodies

function dynamics_orbit(x::PositionState, dt::Float64)::PositionState
    dv = -SatelliteToolboxBase.GM_EARTH / (norm(x.position_ECI)^3) .* x.position_ECI
    dx = x.velocity_ECI
    return PositionState(x.elapsed_time, dx, dv)
end

function dynamics_orbit(r_ECI::Vec3d, v_ECI::Vec3d, pert_F_ECI::Vec3d, dt::Float64, params)::Tuple{Vec3d, Vec3d}
    J2 = zeros(3)
    if params["do_J2"]
        z = [0.0;0.0;1.0]
        J2 = 3 * SatelliteToolboxBase.GM_EARTH * SatelliteToolboxBase.EGM_2008_J2 * SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS^2 / 2 / norm(r_ECI)^5 * ((5 / norm(r_ECI)^2 * dot(z,r_ECI)^2 - 1) * r_ECI - 2 * dot(z,r_ECI) * z)
    end
    dr = v_ECI
    dv = -SatelliteToolboxBase.GM_EARTH / (norm(r_ECI)^3) .* r_ECI + J2 + pert_F_ECI
    return (dr, dv)
end

function push_orbit(x::PositionState, dt::Float64)::PositionState
    # RK4
    k1 = dynamics_orbit(x, dt)
    k2 = dynamics_orbit(x + dt/2*k1, dt)
    k3 = dynamics_orbit(x + dt/2*k2, dt)
    k4 = dynamics_orbit(x + dt*k3, dt)
    res = x + dt/6*(k1 + 2*(k2 + k3) + k4)
    # res = x + dt*k1
    res.elapsed_time = x.elapsed_time + dt
    return res
end

function euler_rigid_dynamics(x::AttitudeState, dt::Float64, params)::AttitudeState
    dw = params["I_Body_inv"]*(-cross(x.angular_velocity_ECI_Body, params["I_Body"]*x.angular_velocity_ECI_Body))
    return AttitudeState(x.elapsed_time, dw, x.attitude_ECI_Body)
end

function euler_rigid_dynamics(w_Body::Vec3d, pert_M_Body::Vec3d, dt::Float64, params)::Vec3d
    dw = params["I_Body_inv"]*(pert_M_Body - cross(w_Body, params["I_Body"]*w_Body))
    return dw
end

function push_attitude(x::AttitudeState, dt::Float64, params)::AttitudeState
    # RK4
    k1 = euler_rigid_dynamics(x, dt, params)
    k2 = euler_rigid_dynamics(x + dt/2*k1, dt, params)
    k3 = euler_rigid_dynamics(x + dt/2*k2, dt, params)
    k4 = euler_rigid_dynamics(x + dt*k3, dt, params)
    res = x + dt/6*(k1 + 2*(k2 + k3) + k4)
    # res = x + dt*k1
    res.elapsed_time = x.elapsed_time + dt
    res.attitude_ECI_Body = x.attitude_ECI_Body*exp(cross(x.angular_velocity_ECI_Body)*dt)
    return res
end

function push_satellite(x::SatelliteState, dt::Float64, params)::SatelliteState
    newx = x
    
    ext_M_Body = Vec3d(0.0)
    ext_F_ECI = Vec3d(0.0)
    
    # RK4
    # wk1 = euler_rigid_dynamics(x.angular_velocity_ECI_Body,             ext_M_Body, dt, params)
    # wk2 = euler_rigid_dynamics(x.angular_velocity_ECI_Body + dt/2*wk1,  ext_M_Body, dt, params)
    # wk3 = euler_rigid_dynamics(x.angular_velocity_ECI_Body + dt/2*wk2,  ext_M_Body, dt, params)
    # wk4 = euler_rigid_dynamics(x.angular_velocity_ECI_Body + dt*wk3,    ext_M_Body, dt, params)
    # newx.angular_velocity_ECI_Body = x.angular_velocity_ECI_Body + dt/6*(wk1 + 2*(wk2 + wk3) + wk4)
    
    wr1, wv1 = dynamics_orbit(x.position_ECI,               x.velocity_ECI,                 ext_F_ECI, dt, params)
    wr2, wv2 = dynamics_orbit(x.position_ECI + dt/2*wr1,    x.velocity_ECI + dt/2*wv1,      ext_F_ECI, dt, params)
    wr3, wv3 = dynamics_orbit(x.position_ECI + dt/2*wr2,    x.velocity_ECI + dt/2*wv2,      ext_F_ECI, dt, params)
    wr4, wv4 = dynamics_orbit(x.position_ECI + dt*wr3,      x.velocity_ECI + dt*wv3,        ext_F_ECI, dt, params)
    newx.position_ECI = x.position_ECI + dt/6*(wr1 + 2*(wr2 + wr3) + wr4)
    newx.velocity_ECI = x.velocity_ECI + dt/6*(wv1 + 2*(wv2 + wv3) + wv4)
    
    # newx.attitude_ECI_Body = x.attitude_ECI_Body*exp(cross(x.angular_velocity_ECI_Body)*dt)
    
    newx.elapsed_time = x.elapsed_time + dt
    
    return newx
end 

function push_earth(x::EarthState, dt::Float64, t::Dates.DateTime, params)::EarthState
    x.elapsed_time += dt
    x.attitude_ECI_ECEF = SatelliteToolboxTransformations.r_eci_to_ecef(SatelliteToolboxTransformations.J2000(), SatelliteToolboxTransformations.ITRF(), SatelliteToolboxBase.date_to_jd(t), eop)
    return x
end

function push_sun(sun::SunState, x::SatelliteState, dt::Float64, t::Dates.DateTime, params)::SunState
    sun.elapsed_time += dt
    sun.position_ECI = SatelliteToolboxCelestialBodies.sun_position_mod(t)
    # todo: replace hardcode Earth radius
    sun.visible = sun.position_ECI'*x.position_ECI / norm(sun.position_ECI) / norm(x.position_ECI) > -sqrt(max(0, 1 - (SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS)^2 / norm(x.position_ECI)^2))
    return sun
end

function push_targets(targets::Vector{GroundState}, x::SatelliteState, dt::Float64, t::Dates.DateTime, params)::Vector{GroundState}
    min_priority = UInt16(0xffff)
    for target in targets
        target.elapsed_time += dt
        target.position_ECI = SatelliteToolboxTransformations.r_ecef_to_eci(SatelliteToolboxTransformations.ITRF(), SatelliteToolboxTransformations.J2000(), SatelliteToolboxBase.date_to_jd(t), eop)*SatelliteToolboxTransformations.geodetic_to_ecef(target.position_LLA...)
        # todo: replace hardcoded mask with params dict
        if (x.position_ECI - target.position_ECI)'*target.position_ECI / norm(x.position_ECI - target.position_ECI) / norm(target.position_ECI) > cos(pi/2 - 20*pi/180)
            target.visible = true
        else
            target.visible = false
        end
    end
    return targets
end

function mode_lookup(choice::T, params)::Union{ModeConfig, Nothing} where T<:Union{AbstractTarget, Nothing}
    k = keys(params["modes"])
    all_modes = [params["modes"][k] for k in keys(params["modes"])]
    mode = filter(m -> m.target_type === typeof(choice), all_modes)
    if length(mode) == 1
        return mode[1]
    else
        throw(DomainError)
    end

end

function clamp_attitude_align!(sat::SatelliteState, sun::SunState, earth::EarthState, targets::Vector{GroundState}, mode::ModeConfig, dt::Float64, params; secondary::Vec3d=Vec3d(0.0,0.0,1.0))
    from_Body = normalize(sat.attitude_ECI_Body*mode.direction_Body)
    to_ECI = normalize(sat.target_ECI - sat.position_ECI)
    # to_ECI = normalize(sat.position_ECI - sat.target_ECI)
    
    if isnan(norm(from_Body)) || isnan(norm(to_ECI)) # cases where either vector is NaN or zero
        # C_Body_Body
        # sat.attitude_ECI_Body = (r_min_arc(from_Body, to_ECI)*sat.attitude_ECI_Body)'
        return
    end
    C_diff = r_min_arc(to_ECI, from_Body)
    C_next = sat.attitude_ECI_Body'*C_diff
    
    # to clamp immediately to the target:
    # sat.attitude_ECI_Body = C_next
    
    # to slew at a fixed rate to the target:
    rot_axis_Body = normalize(uncross(log(C_diff)))
    sat.angular_velocity_ECI_Body = params["max_angular_rate"]*dt*rot_axis_Body
    C_push = exp(cross(sat.angular_velocity_ECI_Body))
    sat.attitude_ECI_Body = sat.attitude_ECI_Body'*C_push
end

function set_mode!(sat::SatelliteState, sun::SunState, earth::EarthState, targets::Vector{GroundState}, dt::Float64, t::Dates.DateTime, params)
    candidates = [sun; targets...]
    visibles = filter(c -> c.visible, candidates) # select the visible ones only
    if length(visibles) > 0
        choice = visibles[findmin(c -> c.priority, visibles)[2]] # select the lowest (highest) priority one only
        mode_conf = mode_lookup(choice, params)
        if !isnothing(mode_conf)
            sat.mode = mode_conf.id
            sat.target_ECI = choice.position_ECI
            clamp_attitude_align!(sat, sun, earth, targets, mode_conf, dt, params)
        end
        return
    end
    mode_conf = mode_lookup(nothing, params)
    sat.mode = mode_conf.id
    sat.target_ECI = Vec3d(NaN)
end

function set_power_data!(sat::SatelliteState, sun::SunState, earth::EarthState, targets::Vector{GroundState}, dt::Float64, t::Dates.DateTime, params)
    mode_conf = params["modes"][sat.mode]
    power_out = mode_conf.power_consumption
    data_in = mode_conf.data_production
    sun_cos = normalize(sun.position_ECI)'*sat.attitude_ECI_Body*normalize(mode_conf.direction_Body)
    power_in = sun.visible*params["irradiance"]*params["solar_panel_area"]*params["solar_panel_efficiency"]*max(0.0, sun_cos)
    
    # todo: replace with TargetConfig-specific downlink rate
    data_out = any([target.visible for target in targets])*params["downlink_rate"]
    
    # cases where sun vector is not set, or is behind the solar panels
    if isnan(power_in) || power_in < 0.0
        power_in = 0.0
    end
    
    # update battery level and rectify
    sat.battery_level = min(params["battery_max"], max(0.0, sat.battery_level + dt*(power_in - power_out)))
    sat.storage_level = min(params["storage_max"], max(0.0, sat.storage_level + dt*(data_in - data_out)))
end

function push_sim!(sat::SatelliteState, sun::SunState, earth::EarthState, targets::Vector{GroundState}, dt::Float64, t::Dates.DateTime, params)
    sat = push_satellite(sat, dt, params)
    sun = push_sun(sun, sat, dt, t, params)
    earth = push_earth(earth, dt, t, params)
    targets = push_targets(targets, sat, dt, t, params)
    
    set_mode!(sat, sun, earth, targets, dt, t, params)
    
    set_power_data!(sat, sun, earth, targets, dt, t, params)
end
