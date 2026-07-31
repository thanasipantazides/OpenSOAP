using Test
using OpenSOAP

@testset "ID consistency" begin
    testjsonpath = joinpath("..", "config", "example_bool.jsonc")
    if isfile(testjsonpath)
        sim, sat, sat_config, target_states, target_configs, constraints, modes =
            load_config(Dict(load_jsonc(testjsonpath)))
        @test check_ids(sat, sat_config, target_states, target_configs, constraints, modes)
    end
end
