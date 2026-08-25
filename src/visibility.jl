
# note: this should be called back in simulation.jl mode  logic, *not* in set_visibility!()
function check_constraints(
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    t::Dates.DateTime,
    modes::Vector{ModeConfig},
    constraints::IDDict{<:AbstractConstraint},
)
    # lookup sat.mode in modes
    mode = findfirst(x->x.id == sat.mode, modes)
    if isnothing(mode)
        return true
    end

    res = true
    for cid in mode.constraint_ids
        # lookup the constraint in the mode
        # check this constraint
        res = res && check_constraint(sat, sat_config, t, constraints[cid])
    end

    return res
end

# define this for every new constraint type
function check_constraint(
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    t::Dates.DateTime,
    constraint::LLAConstraint,
)
    C_ECEF_ECI = SatelliteToolboxTransformations.r_eci_to_ecef(
        SatelliteToolboxTransformations.J2000(),
        SatelliteToolboxTransformations.ITRF(),
        SatelliteToolboxBase.date_to_jd(t),
        eop,
    )
    pos_LLA = SatelliteToolboxTransformations.ecef_to_geodetic(C_ECEF_ECI*sat.position_ECI)
    return all([
        min(constraint.lat...) <= pos_LLA[1] <= max(constraint.lat...),
        min(constraint.lon...) <= pos_LLA[2] <= max(constraint.lon...),
        min(constraint.alt...) <= pos_LLA[3] <= max(constraint.alt...),
    ])
end

function set_visibility!(
    sun::SunState,
    target_config::AbstractConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    t::Dates.DateTime,
)
    sun.visible =
        sun.position_ECI'*sat.position_ECI / norm(sun.position_ECI) /
        norm(sat.position_ECI) >
        -sqrt(
            max(
                0,
                1 -
                (SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS)^2 / norm(sat.position_ECI)^2,
            ),
        )
end

function set_visibility!(
    gs::GroundState,
    target_config::AbstractConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    t::Dates.DateTime,
)
    gs.visible =
        (sat.position_ECI - gs.position_ECI)'*gs.position_ECI /
        norm(sat.position_ECI - gs.position_ECI) / norm(gs.position_ECI) >
        cos(target_config.position_cone)
end

function set_visibility!(
    mag::MagneticFieldState,
    target_config::AbstractConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    t::Dates.DateTime,
)
    mag.visible = true
end
