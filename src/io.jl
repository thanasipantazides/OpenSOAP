import YAML, SatelliteToolboxTransformations

function load_mission(file::String, type::Type)

    data = YAML.load_file(file)
    println("\nloading ", length(data), " target records...")

    if type === GroundTarget
        eops = SatelliteToolboxTransformations.fetch_iers_eop()
        result = Vector{GroundTarget}(undef,length(data))
        for i in 1:length(result)
            this_gs = data[i]
            result[i] = GroundTarget(
                this_gs["name"],
                [this_gs["latitude"]*π/180; this_gs["longitude"]*π/180; this_gs["altitude"]],
                [0;0;1],
                (90 - this_gs["elevation"])*π/180,
                eops
            )
        end

        return result
    end
end