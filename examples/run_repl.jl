using Distributed
using Observables

if nprocs() < 3
    addprocs(3-nprocs())
end

println("running with $(nprocs()) workers")

@everywhere import OpenSOAP

function main()
    fconfig = joinpath("config", "example.jsonc")

    sim, sat, sat_config, sim_environment =
        OpenSOAP.load_config(Dict(OpenSOAP.load_jsonc(fconfig)))


    conn_repl =
        Threads.@spawn OpenSOAP.setup_server(OpenSOAP.SOAP_HOST, OpenSOAP.SOAP_REPL_PORT)

    println("launching core...")
    fcore = @spawnat workers()[1] OpenSOAP.run(sim, sat, sat_config, sim_environment)
    sleep(3)
    println("launching monitor...")
    fmon = @spawnat workers()[2] OpenSOAP.monitor(config_path = fconfig)

    connected_sock = fetch(conn_repl)

    # proof of concept REPL update to a TargetState:
    #   - return to REPL from this main() function, get obs_targets
    #   - lookup in obs_targets[][id]
    #   - overwrite obs_targets[][id] or obs_targets[][id].field
    #   - call notify(obs_targets) and see results!

    sim_env = deepcopy(sim_environment)
    obs_env = Observable(sim_env)
    on(obs_env) do oe
        println("called back")
        for k in keys(oe)
            if !(k in keys(sim_environment)) ||
               !OpenSOAP.mut_struct_eq(oe[k], sim_environment[k])
                println("updating $k")
                packet = OpenSOAP.packetize(oe[k], 0x0000)
                OpenSOAP.write_transport(connected_sock, packet)
            end
        end
    end

    Threads.@spawn begin
        fetch(fcore)
        close(connected_sock.sock)
        sleep(1)
        # this has proven helpful in really closing sockets which 
        # may not be fully drained. Allows re-running quickly without
        # running into TIME_WAIT issues with EADDRINUSE.
        @everywhere GC.gc()
    end

    return connected_sock, sim, sat, sat_config, obs_env
end
