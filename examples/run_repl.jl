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


    conn_repl =
        Threads.@spawn OpenSOAP.setup_server(OpenSOAP.SOAP_HOST, OpenSOAP.SOAP_REPL_PORT)

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

    connected = fetch(conn_repl)
    sleep(1)

    return connected, sim, sat, sat_config, targets, target_configs, constraints, modes
end
