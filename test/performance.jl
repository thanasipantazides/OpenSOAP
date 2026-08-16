using OpenSOAP
import Dates


function core_rate()
    testjsonpath = joinpath("config", "example.jsonc")
    if !isfile(testjsonpath)
        error("bad config file path")
    end

    println("loading config...")
    sim, sat, sat_config, target_states, target_configs, constraints, modes =
        load_config(Dict(load_jsonc(testjsonpath)))

    println("testing simulation speed...")
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

# e0eaa5b: 10^5 steps in 10.5 s
# 874c264: 10^5 steps in 11.1 s (includes calls to atmospheric model NRLMSISE00 each step)
