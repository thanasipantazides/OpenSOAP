using OpenSOAP
import Dates


function core_rate()
    testjsonpath = joinpath("config", "example.jsonc")
    if !isfile(testjsonpath)
        @error("bad config file path")
    end

    sim, sat, sat_config, target_states, target_configs, constraints, modes =
        load_config(Dict(load_jsonc(testjsonpath)))

    time = sim.start_time
    @time while sim.step_count < 10^5
        OpenSOAP.step_sim!(
            sat,
            sat_config,
            target_states,
            target_configs,
            constraints,
            modes,
            sim.time_step,
            time,
            Dict("do_J2"=>true);
        )

        # update timestep
        time = sim.start_time + Dates.Millisecond(1000*sat.elapsed_time)

        # serialize and send data
        sat_data = packetize(sat, 0x0000, sim.step_count)

        for target in values(target_states)
            t_data = packetize(target, 0x0000, sim.step_count)
        end

        sim.step_count += 1
    end
end
