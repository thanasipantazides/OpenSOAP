module OpenSOAP

include("base.jl")
include("io.jl")
include("polyhedra.jl")
include("target.jl")
include("mission.jl")

include("orbit.jl")
include("attitude.jl")
include("dynamics.jl")
include("mission.jl")
include("control.jl")
include("visualize.jl")

export LEOSimulation, EarthProperties, Mission, SpacecraftProperties, PowerProperties, DataProperties, MassProperties, SolarPanel, Antenna, AbstractTarget, SunTarget, GroundTarget, MagneticTarget, State, Polyhedron, Maneuver #, FrameFixedTarget
export load_mission, load_groundstations, load_spacecraft, load_simulation
export cross, uncross, r_min_arc, r_random, r_euler3, r_euler2, r_euler1, ang321, axis, axisangle, project, rotinterp
export point_between
export can_see_sun, can_see_gnd
export ECI, ECEF, ICRS, visibility_history, position_eci, position_ecef
export simulate
export integrate_system, dynamics_orbit!, dynamics_attitude!, backout_free_torque
export plot!, polyplot!, plot_earth!, plot_spacecraft!, plot_targets!, plot_detail!, plot_visibilities!, plot_power!, plot_data!, plot_mode!, plot_frames!, plot_meta!, plot_state!, plot_magnetic_field!, obs_last_orbit_eci, lookat_orbit!, hud_conop!, plot_earth_static!, plot_targets_static!, load_earth_texture_to_ecef
export mission_stats, write

export Quaternion

greet() = print("Hello World!")

end # module OpenSOAP
