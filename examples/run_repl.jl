using Distributed
using Observables

if nprocs() < 3
    addprocs(3-nprocs())
end

println("running with $(nprocs()) workers")

@everywhere import OpenSOAP

function main()
    fconfig = joinpath("config", "example.jsonc")

    sim, sat, sim_environment = OpenSOAP.load_config(Dict(OpenSOAP.load_jsonc(fconfig)))


    conn_repl =
        Threads.@spawn OpenSOAP.setup_server(OpenSOAP.SOAP_HOST, OpenSOAP.SOAP_REPL_PORT)

    println("launching core...")
    fcore = @spawnat workers()[1] OpenSOAP.run(sim, sat, sim_environment)
    sleep(3)
    println("launching monitor...")
    fmon = @spawnat workers()[2] OpenSOAP.monitor(config_path = fconfig)

    connected_sock = fetch(conn_repl)

    # proof of concept REPL update to a TargetState:
    #   - return to REPL from this main() function, get obs_targets
    #   - lookup in obs_targets[][id]
    #   - overwrite obs_targets[][id] or obs_targets[][id].field
    #   - call notify(obs_targets) and see results!

    targets2 = deepcopy(targets)
    obs_targets = Observable(targets2)
    on(obs_targets) do ot
        println("called back")
        println("keys: ", keys(ot))
        for k in keys(ot)
            if !(k in keys(targets)) || !OpenSOAP.mut_struct_eq(ot[k], targets[k])
                println("updating $k")
                packet = OpenSOAP.packetize(ot[k], 0x0000)
                OpenSOAP.write_transport(connected_sock, packet)
            end
        end
    end

    sleep(1)

    return connected_sock,
    sim,
    sat,
    sat_config,
    obs_targets,
    target_configs,
    constraints,
    modes
end
