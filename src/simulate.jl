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

function step_orbit(x::PositionState, dt::Float64)::PositionState
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

function euler_rigid_dynamics(angular_velocity_ECI_Body::Vec3d, control_M_Body::Vec3d, pert_M_Body::Vec3d, dt::Float64, sat_config::SatelliteConfig)::Vec3d
    dw = sat_config.inertia_inv_Body*(pert_M_Body + control_M_Body - cross(angular_velocity_ECI_Body, sat_config.inertia_Body*angular_velocity_ECI_Body))
    return dw
end

function step_attitude!(x::SatelliteState, control_M_Body::Vec3d, pert_M_Body::Vec3d, dt::Float64, sat_config::SatelliteConfig, params)
    k1 = euler_rigid_dynamics(x.angular_velocity_ECI_Body, control_M_Body, pert_M_Body, dt, sat_config)
    k2 = euler_rigid_dynamics(x.angular_velocity_ECI_Body .+ dt/2*k1, control_M_Body, pert_M_Body, dt, sat_config)
    k3 = euler_rigid_dynamics(x.angular_velocity_ECI_Body .+ dt/2*k2, control_M_Body, pert_M_Body, dt, sat_config)
    k4 = euler_rigid_dynamics(x.angular_velocity_ECI_Body .+ dt*k3, control_M_Body, pert_M_Body, dt, sat_config)

    x.angular_velocity_ECI_Body = x.angular_velocity_ECI_Body + dt/6*(k1 + 2*(k2 + k3) + k4)
    x.attitude_ECI_Body = x.attitude_ECI_Body*exp(cross(x.angular_velocity_ECI_Body)*dt)

    x.net_moment_Body = control_M_Body + pert_M_Body
end

function step_attitude(x::AttitudeState, dt::Float64, params)::AttitudeState
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

function step_satellite(x::SatelliteState, dt::Float64, params)::SatelliteState
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

function step!(earth::EarthState, sat::SatelliteState, dt::Float64, t::Dates.DateTime, config::TargetConfig, params)
    earth.elapsed_time += dt
    earth.attitude_ECI_ECEF = SatelliteToolboxTransformations.r_eci_to_ecef(SatelliteToolboxTransformations.J2000(), SatelliteToolboxTransformations.ITRF(), SatelliteToolboxBase.date_to_jd(t), eop)
end

function step!(sun::SunState, sat::SatelliteState, dt::Float64, t::Dates.DateTime, config::TargetConfig, params)
    sun.elapsed_time += dt
    sun.position_ECI = SatelliteToolboxCelestialBodies.sun_position_mod(t)
    sun.visible = sun.position_ECI'*sat.position_ECI / norm(sun.position_ECI) / norm(sat.position_ECI) > -sqrt(max(0, 1 - (SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS)^2 / norm(sat.position_ECI)^2))
end

function step!(gs::GroundState, sat::SatelliteState, dt::Float64, t::Dates.DateTime, config::TargetConfig, params)
    gs.elapsed_time += dt
    gs.position_ECI = SatelliteToolboxTransformations.r_ecef_to_eci(SatelliteToolboxTransformations.ITRF(), SatelliteToolboxTransformations.J2000(), SatelliteToolboxBase.date_to_jd(t), eop)*SatelliteToolboxTransformations.geodetic_to_ecef(gs.position_LLA...)
    # todo: replace hardcoded mask with params dict
    if (sat.position_ECI - gs.position_ECI)'*gs.position_ECI / norm(sat.position_ECI - gs.position_ECI) / norm(gs.position_ECI) > cos(config.position_cone)
        gs.visible = true
    else
        gs.visible = false
    end
end

function step_earth(x::EarthState, dt::Float64, t::Dates.DateTime, params)::EarthState
    x.elapsed_time += dt
    x.attitude_ECI_ECEF = SatelliteToolboxTransformations.r_eci_to_ecef(SatelliteToolboxTransformations.J2000(), SatelliteToolboxTransformations.ITRF(), SatelliteToolboxBase.date_to_jd(t), eop)
    return x
end

function step_sun(sun::SunState, x::SatelliteState, dt::Float64, t::Dates.DateTime, params)::SunState
    sun.elapsed_time += dt
    sun.position_ECI = SatelliteToolboxCelestialBodies.sun_position_mod(t)
    # todo: replace hardcode Earth radius
    sun.visible = sun.position_ECI'*x.position_ECI / norm(sun.position_ECI) / norm(x.position_ECI) > -sqrt(max(0, 1 - (SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS)^2 / norm(x.position_ECI)^2))
    return sun
end

function step_targets(targets::Vector{GroundState}, x::SatelliteState, dt::Float64, t::Dates.DateTime, params)::Vector{GroundState}
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

function clamp_attitude_align!(sat::SatelliteState, mode::ModeConfig, dt::Float64, t::Dates.DateTime, target::AbstractTarget, sat_config::SatelliteConfig, params; exog::Union{Nothing, PerturbationMessage}=nothing, secondary::Vec3d=Vec3d(0.0,0.0,1.0))
    
    # target = targets[findfirst(t -> t.id == sat.target, targets)]
    from_Body = normalize(sat.attitude_ECI_Body*mode.direction_Body)
    to_ECI = normalize(target.position_ECI - sat.position_ECI)
    # to_ECI = normalize(sat.position_ECI - sat.target_ECI)
    
    if isnan(norm(from_Body)) || isnan(norm(to_ECI)) # cases where either vector is NaN or zero
        # C_Body_Body
        # sat.attitude_ECI_Body = (r_min_arc(from_Body, to_ECI)*sat.attitude_ECI_Body)'
        return
    end
    C_err = r_min_arc(to_ECI, from_Body)

    inner_attitude_dynamic!(sat, C_err, dt, t, sat_config, params; style=:pd, exog=exog)
end

function inner_attitude_dynamic!(sat::SatelliteState, C_err::Mat3d, dt::Float64, t::Dates.DateTime, sat_config::SatelliteConfig, params; style=:clamp, exog::Union{Nothing, PerturbationMessage}=nothing)

    m_pert_Body =  Vec3d(0.0)
    if !isnothing(exog)
        m_pert_Body = exog.moment_Body
    end
    # snap instantly to target attitude
    if style == :clamp
        C_next = sat.attitude_ECI_Body'*C_err
        sat.attitude_ECI_Body = C_next
        sat.angular_velocity_ECI_Body = Vec3d(0.0)

    # spherical linear interpolation of attitude
    elseif style == :slerp
        rot_axis_Body = normalize(uncross(log(C_err)))
        sat.angular_velocity_ECI_Body = sat_config.angular_rate_max*dt*rot_axis_Body
        C_push = exp(cross(sat.angular_velocity_ECI_Body))
        sat.attitude_ECI_Body = sat.attitude_ECI_Body'*C_push

    # slerp, but applied to the control input.
    elseif style == :uslerp
        rot_axis_Body = normalize(uncross(log(C_err)))
        kp = 1e-3 
        kd = 1e-3
        # u_B = kp*sat_config.angular_rate_max*dt*rot_axis_Body
        # u_B = kp*sat_config.angular_rate_max*dt*sat.attitude_ECI_Body'*rot_axis_Body

        # direct feedback, i.e. drive attitude to identity.
        # ax = uncross(log(sat.attitude_ECI_Body))
        # ud = Mat3d([0 -1 0; 1 0 0; 0 0 1])
        pax = sat.attitude_ECI_Body'*uncross(log(sat.attitude_ECI_Body*C_err))
        if norm(pax) == 0
            pax = Vec3d(0)
        else
            pax = normalize(pax)
        end
        u_B = -kp*sat_config.angular_rate_max*dt*pax - kd*sat.angular_velocity_ECI_Body
        step_attitude!(sat, Vec3d(u_B), m_pert_Body , dt, sat_config, params)
        
    # proportional-derivative control
    elseif style == :pd
        # todo: move these to params or SatelliteConfig. 
        kp = 5e-4
        kd = 1e-3
        # presume that we want zero angular velocity
        angular_velocity_ECI_Body_goal = Vec3d(0.0)
        angular_velocity_err = sat.angular_velocity_ECI_Body - angular_velocity_ECI_Body_goal
        # u_B = -kp*uncross(log(sat.attitude_ECI_Body'*C_err)) - kd*(angular_velocity_err)
        pax = sat.attitude_ECI_Body*uncross(log(sat.attitude_ECI_Body'*C_err))
        # if norm(pax) == 0
        #     pax = Vec3d(0)
        # else
        #     pax = normalize(pax)
        # end
        u_B = -kp*pax - kd*angular_velocity_err
        # println(norm(u_B), ", ", norm(sat.angular_velocity_ECI_Body))
        # u_B = Vec3d(0.0)

        step_attitude!(sat, Vec3d(u_B), m_pert_Body, dt, sat_config, params)

        if abs(residualSO3(sat.attitude_ECI_Body)) > 1.0
            throw(BoundsError())
        end
        # println(" > ", norm(sat.angular_velocity_ECI_Body))
        # println(" > ", residualSO3(sat.attitude_ECI_Body))
    else
        throw(ArgumentError(msg="invalid style value "*style))
    end
end

# todo: change targets::Vector{GroundState} to targets::Vector{AbstractTarget} so it includes SunTarget.
function set_mode!(sat::SatelliteState, targets::Vector{AbstractTarget}, target_configs::Vector{TargetConfig}, dt::Float64, t::Dates.DateTime, modes::Vector{ModeConfig}, mode_table::Matrix{UInt8}, sat_config::SatelliteConfig, params; exog::Union{Nothing, PerturbationMessage})
    
    visibility = map(c -> c.visible, targets)
    feasible = mode_table'*visibility
    best_priority = typemax(modes[1].priority)
    modechoice = nothing
    mode_k = 0
    for (m,mode) in enumerate(modes)
        if feasible[m] > 0 && mode.priority < best_priority
            modechoice = mode
            best_priority = mode.priority
            mode_k = m # gives the column index (mode index) in mode_table we have selected
        end
    end
    
    targetchoice = nothing
    if isnothing(modechoice) # due to lack of feasible targets
        idlemode = modes[findfirst(m -> length(m.target_ids) == 0, modes)]
        sat.mode = idlemode.id
        sat.target = typemax(IDType)
    else
        sat.mode = modechoice.id
        # find the target actually used for this mode, so we can point at it:
        targetchoice = targets[findfirst(t -> t != 0, visibility .* mode_table[:,mode_k])]
        sat.target = targetchoice.id
        
        clamp_attitude_align!(sat, modechoice, dt, t, targetchoice, sat_config, params; exog)
    end
    set_power_data!(sat, targets, target_configs, dt, t, modes, targetchoice, sat_config, params)
end

function set_power_data!(sat::SatelliteState, targets::Vector{AbstractTarget}, target_configs::Vector{TargetConfig}, dt::Float64, t::Dates.DateTime, modes::Vector{ModeConfig}, target::Union{Nothing, AbstractTarget}, sat_config::SatelliteConfig, params)
    sun = targets[findfirst(x -> isa(x, SunState), targets)]
    
    mode_conf = filter(m -> m.id == sat.mode, modes)[1] # todo: add arg modes::Dict{IDType, ModeConfig} to use for lookup.
    # mode_conf = params["modes"][sat.mode] 
    power_out = mode_conf.power_consumption
    data_in = mode_conf.data_production
    sun_cos = normalize(sun.position_ECI)'*sat.attitude_ECI_Body*normalize(mode_conf.direction_Body)
    # todo: pull irradiance from sim config
    power_in = 1360*sun.visible*sat_config.power_solar_panel_efficiency*sat_config.power_solar_panel_area*max(0.0, sun_cos)
    
    # todo: replace with TargetConfig-specific downlink rate
    data_out = 0.0
    if isa(target, GroundState)
        conf = find_config(target, target_configs)
        data_out = conf.data_consumption
    end
    
    # cases where sun vector is not set, or is behind the solar panels
    if isnan(power_in) || power_in < 0.0
        power_in = 0.0
    end
    
    sat.net_power = power_in - power_out
    sat.net_data = data_in - data_out
    # update battery level and rectify
    sat.battery_level = min(sat_config.power_battery_max, max(0.0, sat.battery_level + dt*(sat.net_power)))
    sat.storage_level = min(sat_config.data_storage_max, max(0.0, sat.storage_level + dt*(sat.net_data)))
end

function step_sim!(sat::SatelliteState, earth::EarthState, targets::Vector{AbstractTarget}, dt::Float64, t::Dates.DateTime, target_configs::Vector{TargetConfig}, modes::Vector{ModeConfig}, mode_table::Matrix{UInt8}, sat_config::SatelliteConfig, params; exog::Union{Nothing, PerturbationMessage}=nothing)
    sat = step_satellite(sat, dt, params)
    # sun = step_sun(sun, sat, dt, t, params)
    earth = step_earth(earth, dt, t, params)
    # targets = step_targets(targets, sat, dt, t, params)
    
    for (target, target_conf) in zip(targets, target_configs)
        step!(target, sat, dt, t, target_conf, params)
    end
    
    set_mode!(sat, targets, target_configs, dt, t, modes, mode_table, sat_config, params; exog=exog)
end
