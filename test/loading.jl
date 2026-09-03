using Test
using OpenSOAP

@testset "ID consistency" begin
    testjsonpath = joinpath("..", "config", "example.jsonc")
    if isfile(testjsonpath)
        sim, sat, sat_config, sim_env = load_config(Dict(load_jsonc(testjsonpath)))
        @test check_ids(sat, sat_config, sim_env)
    end
end
