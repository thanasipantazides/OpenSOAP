using Test
using OpenSOAP

@testset "ser/des roundtrip tests" begin
    testjsonpath = joinpath("..", "config", "example.jsonc")
    if !isfile(testjsonpath)
        error("file $testjsonpath is not a file!")
    end
    sim, sat, sat_config, sim_env = load_config(Dict(load_jsonc(testjsonpath)))

    target_states = filter(p -> p isa AbstractTarget, sim_env)
    target_configs = filter(p -> p isa AbstractConfig && !(p isa ModeConfig), sim_env)
    constraints = filter(p -> p isa AbstractConstraint, sim_env)
    modes = filter(p -> p isa ModeConfig, sim_env)

    @testset "SatelliteState and Config" begin
        @test mut_struct_eq(sat, des(SatelliteState, ser(sat)))
        @test mut_struct_eq(sat_config, des(SatelliteConfig, ser(sat_config)))
    end

    @testset "TargetState and Config" begin
        for k in keys(target_states)
            @test mut_struct_eq(
                target_states[k],
                des(typeof(target_states[k]), ser(target_states[k])),
            )
        end
        for k in keys(target_configs)
            @test mut_struct_eq(
                target_configs[k],
                des(typeof(target_configs[k]), ser(target_configs[k])),
            )
        end
    end

    @testset "Constraints" begin
        for k in keys(constraints)
            @test mut_struct_eq(
                constraints[k],
                des(typeof(constraints[k]), ser(constraints[k])),
            )
        end
    end

    @testset "ModeConfig" begin
        for k in keys(modes)
            @test mut_struct_eq(modes[k], des(typeof(modes[k]), ser(modes[k])))
        end
    end

    whole_dict = merge(
        target_states,
        target_configs,
        constraints,
        Dict{IDType,ModeConfig}(m.id => m for m in modes),
        Dict{IDType,NetworkMessage}(sat.id => sat, sat_config.id => sat_config),
    )
    @testset "Packetization" begin
        all_flags = UInt16(0xabcd)
        for k in keys(whole_dict)
            out = packetize(whole_dict[k], all_flags)
            type, flags, len = behead(out[1:8])
            @test type == typeof(whole_dict[k]) && all_flags == flags

            @test mut_struct_eq(des(type, out[9:end]), whole_dict[k])
        end
    end

end
