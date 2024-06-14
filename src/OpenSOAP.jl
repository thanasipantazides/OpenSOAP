module OpenSOAP

include("base.jl")
include("target.jl")
include("mission.jl")

include("orbit.jl")
include("dynamics.jl")
include("visualize.jl")

export LEOSimulation, EarthProperties, Mission, SpacecraftProperties, PowerProperties, DataProperties, MassProperties, SolarPanel, Antenna, AbstractTarget, SunTarget, GroundTarget#, FrameFixedTarget
export cross, r_min_arc
export ECI, ECEF, ICRS, visibility_history, position_eci, position_ecef
export integrate_system, dynamics_orbit!
export plot!, plot_earth!, plot_spacecraft!, plot_targets!, plot_detail!, plot_visibilities!, plot_power!, load_earth_texture_to_ecef

greet() = print("Hello World!")

end # module OpenSOAP
