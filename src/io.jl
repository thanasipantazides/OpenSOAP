import Makie
using GeometryBasics
# for now, just serve a default mode list.
function load_mode_config(path::AbstractString)
    modelist = [
        ModeConfig(0x0101, "idle",      Nothing,    Makie.RGBAf(42/255, 133/255, 255/255), 5.0, 1e3, Vec3d(1.0,0.0,0.0)),
        ModeConfig(0x0102, "charging",  SunState,   Makie.RGBAf(255/255, 201/255, 74/255), 5.0, 1e6, Vec3d(-1.0,0.0,0.0)),
        ModeConfig(0x0103, "downlink",  TargetState, Makie.RGBAf(157/255, 226/255, 107/255), 23.0, 10e3, Vec3d(0.0,1.0,0.0)),
        # ModeConfig(0x0104, "science",   Nothing,    Makie.RGBAf(206/255, 155/255, 255/255), 20.0, 10e3, Vec3d(1.0,0.0,0.0))
    ]
    return modelist
end