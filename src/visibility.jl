function set_visibility!(sat::SatelliteState, x::SunState, config::TargetConfig)
    x.visible =
        x.position_ECI'*sat.position_ECI / norm(x.position_ECI) / norm(sat.position_ECI) >
        -sqrt(
            max(
                0,
                1 -
                (SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS)^2 / norm(sat.position_ECI)^2,
            ),
        )
end

function set_visibility!(sat::SatelliteState, x::GroundState, config::TargetConfig)
    x.visible =
        (sat.position_ECI - x.position_ECI)'*x.position_ECI /
        norm(sat.position_ECI - x.position_ECI) / norm(x.position_ECI) >
        cos(config.position_cone)
end
