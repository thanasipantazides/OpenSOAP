using Distributed
if nprocs() < 3
    addprocs(3-nprocs())
end

println("running with $(nprocs()) workers")

@everywhere import OpenSOAP

function main()
    fconfig = joinpath("config", "example.jsonc")

    sim, sat, sat_config, targets, target_configs, constraints, modes =
        OpenSOAP.load_config(Dict(OpenSOAP.load_jsonc(fconfig)))

    println("launching core...")
    fcore = @spawnat workers()[1] OpenSOAP.run(
        sim,
        sat,
        sat_config,
        targets,
        target_configs,
        constraints,
        modes,
    )
    sleep(3)
    println("launching monitor...")
    fmon = @spawnat workers()[2] OpenSOAP.monitor(config_path = fconfig)

    # fm = fetch(fcore)

    sleep(1)

    return targets, target_configs
end
