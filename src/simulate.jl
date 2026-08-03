using LinearAlgebra
import Dates
import SatelliteToolboxTransformations,
    SatelliteToolboxBase, SatelliteToolboxCelestialBodies

function dynamics_orbit(x::PositionState, dt::Float64)::PositionState
    dv = -SatelliteToolboxBase.GM_EARTH / (norm(x.position_ECI)^3) .* x.position_ECI
    dx = x.velocity_ECI
    return PositionState(x.elapsed_time, dx, dv)
end

function dynamics_orbit(
    r_ECI::Vec3d,
    v_ECI::Vec3d,
    pert_F_ECI::Vec3d,
    dt::Float64,
    params,
)::Tuple{Vec3d,Vec3d}
    J2 = zeros(3)
    if params["do_J2"]
        z = [0.0; 0.0; 1.0]
        J2 =
            3 *
            SatelliteToolboxBase.GM_EARTH *
            SatelliteToolboxBase.EGM_2008_J2 *
            SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS^2 / 2 / norm(r_ECI)^5 *
            ((5 / norm(r_ECI)^2 * dot(z, r_ECI)^2 - 1) * r_ECI - 2 * dot(z, r_ECI) * z)
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

function euler_rigid_dynamics(
    angular_velocity_ECI_Body::Vec3d,
    control_M_Body::Vec3d,
    pert_M_Body::Vec3d,
    dt::Float64,
    sat_config::SatelliteConfig,
)::Vec3d
    dw =
        sat_config.inertia_inv_Body*(
            pert_M_Body + control_M_Body - cross(
                angular_velocity_ECI_Body,
                sat_config.inertia_Body*angular_velocity_ECI_Body,
            )
        )
    return dw
end

function step_attitude!(
    x::SatelliteState,
    control_M_Body::Vec3d,
    pert_M_Body::Vec3d,
    dt::Float64,
    sat_config::SatelliteConfig,
    params,
)
    k1 = euler_rigid_dynamics(
        x.angular_velocity_ECI_Body,
        control_M_Body,
        pert_M_Body,
        dt,
        sat_config,
    )
    k2 = euler_rigid_dynamics(
        x.angular_velocity_ECI_Body .+ dt/2*k1,
        control_M_Body,
        pert_M_Body,
        dt,
        sat_config,
    )
    k3 = euler_rigid_dynamics(
        x.angular_velocity_ECI_Body .+ dt/2*k2,
        control_M_Body,
        pert_M_Body,
        dt,
        sat_config,
    )
    k4 = euler_rigid_dynamics(
        x.angular_velocity_ECI_Body .+ dt*k3,
        control_M_Body,
        pert_M_Body,
        dt,
        sat_config,
    )

    x.angular_velocity_ECI_Body = x.angular_velocity_ECI_Body + dt/6*(k1 + 2*(k2 + k3) + k4)
    x.attitude_ECI_Body = x.attitude_ECI_Body*exp(cross(x.angular_velocity_ECI_Body)*dt)

    x.net_moment_Body = control_M_Body + pert_M_Body
end

function step_satellite(x::SatelliteState, dt::Float64, params)::SatelliteState
    newx = x

    ext_M_Body = Vec3d(0.0)
    ext_F_ECI = Vec3d(0.0)

    wr1, wv1 = dynamics_orbit(x.position_ECI, x.velocity_ECI, ext_F_ECI, dt, params)
    wr2, wv2 = dynamics_orbit(
        x.position_ECI + dt/2*wr1,
        x.velocity_ECI + dt/2*wv1,
        ext_F_ECI,
        dt,
        params,
    )
    wr3, wv3 = dynamics_orbit(
        x.position_ECI + dt/2*wr2,
        x.velocity_ECI + dt/2*wv2,
        ext_F_ECI,
        dt,
        params,
    )
    wr4, wv4 = dynamics_orbit(
        x.position_ECI + dt*wr3,
        x.velocity_ECI + dt*wv3,
        ext_F_ECI,
        dt,
        params,
    )
    newx.position_ECI = x.position_ECI + dt/6*(wr1 + 2*(wr2 + wr3) + wr4)
    newx.velocity_ECI = x.velocity_ECI + dt/6*(wv1 + 2*(wv2 + wv3) + wv4)

    # newx.attitude_ECI_Body = x.attitude_ECI_Body*exp(cross(x.angular_velocity_ECI_Body)*dt)

    newx.elapsed_time = x.elapsed_time + dt

    return newx
end

function step!(
    earth::EarthState,
    earth_config::EarthConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    modes::Vector{ModeConfig},
    dt::Float64,
    t::Dates.DateTime,
    params,
)
    earth.elapsed_time += dt
    earth.attitude_ECI_ECEF = SatelliteToolboxTransformations.r_eci_to_ecef(
        SatelliteToolboxTransformations.J2000(),
        SatelliteToolboxTransformations.ITRF(),
        SatelliteToolboxBase.date_to_jd(t),
        eop,
    )
end

function step!(
    sun::SunState,
    sun_config::TargetConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    modes::Vector{ModeConfig},
    dt::Float64,
    t::Dates.DateTime,
    params,
)
    sun.elapsed_time += dt
    sun.position_ECI = SatelliteToolboxCelestialBodies.sun_position_mod(t)
    set_visibility!(sun, sun_config, sat, sat_config, t, modes)
end

function step!(
    gs::GroundState,
    gs_config::TargetConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    modes::Vector{ModeConfig},
    dt::Float64,
    t::Dates.DateTime,
    params,
)
    gs.elapsed_time += dt
    gs.position_ECI =
        SatelliteToolboxTransformations.r_ecef_to_eci(
            SatelliteToolboxTransformations.ITRF(),
            SatelliteToolboxTransformations.J2000(),
            SatelliteToolboxBase.date_to_jd(t),
            eop,
        )*SatelliteToolboxTransformations.geodetic_to_ecef(gs.position_LLA...)
    set_visibility!(gs, gs_config, sat, sat_config, t, modes)
end

function step!(
    mag::MagneticFieldState,
    mag_config::MagneticFieldConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    modes::Vector{ModeConfig},
    dt::Float64,
    t::Dates.DateTime,
    params,
)
    mag.elapsed_time += dt
    C_ECEF_ECI = SatelliteToolboxTransformations.r_eci_to_ecef(
        SatelliteToolboxTransformations.J2000(),
        SatelliteToolboxTransformations.ITRF(),
        SatelliteToolboxBase.date_to_jd(t),
        eop,
    )
    pos_LLA = SatelliteToolboxTransformations.ecef_to_geodetic(C_ECEF_ECI*sat.position_ECI)
    field_ECI = Vec3d(
        C_ECEF_ECI'*SatelliteToolboxTransformations.ned_to_ecef(
            SatelliteToolboxGeomagneticField.igrf(
                DateFormats.yeardecimal(t),
                pos_LLA[3],
                pos_LLA[1],
                pos_LLA[2],
                Val(:geodetic);
                max_degree = Int64(mag_config.model_order),
                show_warnings = false,
            ),
            pos_LLA...,
            translate = false,
        ),
    )

    if mag_config.normalization == 0x00
        # no op
    elseif mag_config.normalization == 0x01 # nadir
        if sat.position_ECI'*field_ECI > 0
            field_ECI = -field_ECI
        end
    elseif mag_config.normalization == 0x02 # zenith
        if sat.position_ECI'*field_ECI < 0
            field_ECI = -field_ECI
        end
    else
        @warn "unhandled MagneticFieldConfig.normalization value $(mag_config.normalization)!\n"
    end

    mag.direction_ECI = field_ECI
    set_visibility!(mag, mag_config, sat, sat_config, t, modes)
end

function clamp_attitude_align!(
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    target::Union{Nothing,<:AbstractTarget},
    mode::Union{Nothing,ModeConfig},
    dt::Float64,
    t::Dates.DateTime,
    params;
    exog::Union{Nothing,PerturbationMessage} = nothing,
    secondary::Vec3d = Vec3d(0.0, 0.0, 1.0),
)


    # choose:
    #   - :pd for proportional-derivative feedback
    #   - :slerp for spherical linear interpolation
    #   - :snap to instantly snap attitude to target
    attitude_style = :pd

    if isnothing(mode)
        inner_attitude_dynamic!(
            sat,
            sat_config,
            sat.attitude_ECI_Body,
            dt,
            t,
            params;
            style = :pd,
            exog = exog,
        )
        return
    end

    from_Body = normalize(sat.attitude_ECI_Body*mode.direction_Body)
    to_ECI = Vec3d(0.0)
    if isnothing(target)
        to_ECI = from_Body
    else
        println(target)
        to_ECI = normalize(reference_direction(target) - sat.position_ECI)
    end

    if isnan(norm(from_Body)) || isnan(norm(to_ECI)) # cases where either vector is NaN or zero
        return
    end

    # C_err = r_min_arc(to_ECI, from_Body)
    C_target = Mat3d(r_min_arc(to_ECI, mode.direction_Body))


    inner_attitude_dynamic!(
        sat,
        sat_config,
        C_target,
        dt,
        t,
        params;
        style = attitude_style,
        exog = exog,
    )
end

function inner_attitude_dynamic!(
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    C_target::Mat3d,
    dt::Float64,
    t::Dates.DateTime,
    params;
    style = :clamp,
    exog::Union{Nothing,PerturbationMessage} = nothing,
)

    m_pert_Body = Vec3d(0.0) # todo: add real perturbation models
    if !isnothing(exog)
        m_pert_Body = exog.moment_Body
    end
    # snap instantly to target attitude
    if style == :snap
        sat.attitude_ECI_Body = C_target'
        sat.angular_velocity_ECI_Body = Vec3d(0.0)

        # spherical linear interpolation of attitude
    elseif style == :slerp
        # todo: this isn't quite right.

        rot_axis_Body = uncross(log(C_target*sat.attitude_ECI_Body))
        sat.angular_velocity_ECI_Body = sat_config.angular_rate_max*dt*rot_axis_Body
        C_push = exp(cross(sat.angular_velocity_ECI_Body))
        sat.attitude_ECI_Body = C_push'*sat.attitude_ECI_Body

        # proportional-derivative control
    elseif style == :pd
        # todo: move these to params or SatelliteConfig. 
        kp = 1e-3
        kd = 5e-3
        # presume that we want zero angular velocity
        angular_velocity_ECI_Body_goal = Vec3d(0.0)
        angular_velocity_err =
            sat.angular_velocity_ECI_Body - angular_velocity_ECI_Body_goal
        pax = uncross(log(C_target*sat.attitude_ECI_Body))

        u_B = -kp*pax - kd*angular_velocity_err
        # todo: limit magnitude of u_B/project to agilitoid
        step_attitude!(sat, Vec3d(u_B), m_pert_Body, dt, sat_config, params)

        if abs(residualSO3(sat.attitude_ECI_Body)) > 1.0
            throw(BoundsError("attitude fell out of SO(3)"))
        end
    else
        throw(ArgumentError("invalid style value "*string(style)))
    end
end

# todo: change targets::Vector{GroundState} to targets::Vector{AbstractTarget} so it includes SunTarget.

function set_mode!(
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    target_states::IDDict{<:AbstractTarget},
    target_configs::IDDict{<:AbstractConfig},
    constraints::IDDict{<:AbstractConstraint},
    modes::Vector{ModeConfig},
    dt::Float64,
    t::Dates.DateTime,
    params;
    exog::Union{Nothing,PerturbationMessage} = nothing,
)

    visible_mode_ids = Set{IDType}()
    allowed_modes = ModeConfig[]
    for mode in modes
        # list all the visible targets for this mode
        for tc in mode.target_ids
            # if *any* of these are visible, can do the mode.
            if target_states[target_configs[tc].dynamic_id].visible
                push!(visible_mode_ids, mode.id)
            end
        end

        satisfied = true
        for c in mode.constraint_ids
            # all constraints must be satisfied for the mode to run. 
            # But if there are no constraints on the mode, this loop won't run at all!
            satisfied = satisfied && check_constraint(sat, sat_config, t, constraints[c])
        end
        satisfied && push!(allowed_modes, mode)
    end

    # filter the allowed (constraint-satisfied) modes by those that have visible targets
    filter!(m->m.id in visible_mode_ids, allowed_modes)

    # pick the mode with highest priority
    modechoice = nothing
    if length(allowed_modes) > 0
        res = findmin(m -> m.priority, allowed_modes)
        if !isnothing(res)
            modechoice = allowed_modes[res[2]]
        end
    end

    # pick the target for this mode
    targetchoice = nothing
    if isnothing(modechoice) # due to lack of feasible targets
        idlemode = modes[findfirst(m -> length(m.target_ids) == 0, modes)]
        sat.mode = idlemode.id
        sat.target = typemax(IDType)
    else
        sat.mode = modechoice.id
        # find the target actually used for this mode, so we can point at it:
        for tc in modechoice.target_ids
            # if any of these are visible, can do the mode. This is finding the first.
            if target_states[target_configs[tc].dynamic_id].visible
                # note: if constraints/targets usage changes such that constraints depend on targets, this logic will need to be revised.
                sat.target = target_configs[tc].dynamic_id
                targetchoice = target_states[sat.target]
                break
            end
        end
    end
    clamp_attitude_align!(
        sat,
        sat_config,
        targetchoice,
        modechoice,
        dt,
        t,
        params;
        exog = exog,
    )
    set_power_data!(sat, sat_config, target_states, target_configs, modes, dt, t, params)
end



function set_power_data!(
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    target_states::IDDict{<:AbstractTarget},
    target_configs::IDDict{<:AbstractConfig},
    modes::Vector{ModeConfig},
    dt::Float64,
    t::Dates.DateTime,
    params,
)
    sun = nothing
    for ts in values(target_states)
        if isa(ts, SunState)
            sun = ts
        end
    end
    # sun = findfirst(x -> isa(x, SunState), values(target_states)...)
    # todo: add arg modes::Dict{IDType, ModeConfig} to use for lookup.
    mode_conf = filter(m -> m.id == sat.mode, modes)[1]
    # mode_conf = params["modes"][sat.mode] 
    power_out = mode_conf.power_consumption
    data_in = mode_conf.data_production
    sun_cos =
        normalize(sun.position_ECI)'*sat.attitude_ECI_Body*normalize(
            mode_conf.direction_Body,
        )
    # todo: pull irradiance from sim config
    power_in =
        1360*sun.visible*sat_config.power_solar_panel_efficiency*sat_config.power_solar_panel_area*max(
            0.0,
            sun_cos,
        )

    data_out = 0.0
    if sat.target != typemax(IDType)
        if isa(target_states[sat.target], GroundState)
            conf = find_config(target_states[sat.target], target_configs)
            data_out = conf.data_consumption
        end
    end
    # if isa(target_states, GroundState)
    # end

    # cases where sun vector is not set, or is behind the solar panels
    if isnan(power_in) || power_in < 0.0
        power_in = 0.0
    end

    sat.net_power = power_in - power_out
    sat.net_data = data_in - data_out
    # update battery level and rectify
    sat.battery_level =
        min(sat_config.power_battery_max, max(0.0, sat.battery_level + dt*(sat.net_power)))
    sat.storage_level =
        min(sat_config.data_storage_max, max(0.0, sat.storage_level + dt*(sat.net_data)))
end

function step_sim!(
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    target_states::IDDict{<:AbstractState},
    target_configs::IDDict{<:AbstractConfig},
    constraints::IDDict{<:AbstractConstraint},
    modes::Vector{ModeConfig},
    dt::Float64,
    time::Dates.DateTime,
    params;
    exog::Union{Nothing,PerturbationMessage} = nothing,
)
    # todo: rename this to "step_orbit" or something consistent. 
    # Most of the SatelliteState doesn't change until set_mode later on.
    sat = step_satellite(sat, dt, params)
    # NEW signature
    for tconf in values(target_configs)
        # extract the target state for this config
        tstate = target_states[tconf.dynamic_id]
        step!(tstate, tconf, sat, sat_config, modes, dt, time, params)
    end

    set_mode!(
        sat,
        sat_config,
        target_states,
        target_configs,
        constraints,
        modes,
        dt,
        time,
        params;
        exog = exog,
    )
end
