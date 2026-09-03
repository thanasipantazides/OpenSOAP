# global const unixsock"tname = "127.0.0.1"


function run(;
    config_path::AbstractString = joinpath("config", "example.jsonc"),
    do_repl = true,
)
    sim, sat, sat_config, sim_environment = load_config(Dict(load_jsonc(config_path)))
    return run(sim, sat, sat_config, sim_environment; do_repl)
end

function run(
    sim::SimConfig,
    sat::SatelliteState,
    sat_config::SatelliteConfig,
    sim_environment::IDDict{<:NetworkMessage};
    do_repl = true,
)# config_path::AbstractString = joinpath("config", "example.jsonc"))

    # sim, sat, sat_config, targets, target_configs, constraints, modes =
    #     load_config(Dict(load_jsonc(config_path)))

    # return targets, target_configs

    # sim.start_time = Dates.DateTime(2026,3,20,1,2,3)
    # sim.time_step = 1.0
    # sim.step_count = UInt64(0)

    pwidth = 12

    inertia_B = diagm([5, 10, 13])*1e-2

    params = Dict{String,Any}(
        "I_Body" => inertia_B,
        "I_Body_inv" => inv(inertia_B),
        "max_angular_rate" => 2e-2, # rad/s
        "battery_max" => 84*3600.0,
        "storage_max" => 8e9,
        "irradiance" => 1360.0,
        "solar_panel_area" => 0.2*0.3*2,
        "solar_panel_efficiency" => 0.29,
        "downlink_rate" => 17e6,
        "do_J2" => true,
    )

    target_states = filter(s -> typeof(s.second) <: AbstractState, sim_environment)
    target_configs = filter(s -> typeof(s.second) <: AbstractConfig, sim_environment)
    println(collect(values(target_configs)))

    # todo: replace this @show call with something that takes in a full sim_environement dict
    # println((modes, targets, target_configs, constraints))

    if do_repl
        println("connecting to REPL...")
        sock_repl = setup_client(SOAP_HOST, SOAP_REPL_PORT)
        println("connected!")
    end

    # udp 
    # sock =
    #     setup_server(SOAP_HOST, SOAP_CORE_PORT, SOAP_HOST, SOAP_MON_PORT)

    # tcp 
    println("waiting for monitor connection...")
    sock_mon = setup_server(SOAP_HOST, SOAP_CORE_PORT)
    println("connected!")

    # unix
    # sock = setup_server(SOAP_UNIX_SOCK)

    sim_env_lock = ReentrantLock()

    play = Ref(UInt8(0x01))
    do_quit = Ref(UInt8(0x00))
    playrate = Ref(UInt8(0x01))
    perturbation = Ref(PerturbationMessage(Vec3d(0.0), 0.0, Vec3d(0.0), 0.0))
    # perturb_moment = Ref(Vec3d(0.0))
    # perturb_moment_count = Ref(UInt32(0x0000))
    packlen = 1
    headbuff = zeros(UInt8, 8)
    buff = zeros(UInt8, packlen)
    send_targets_per_sat_update = Int(10)

    # handle received commands
    # todo: factor this out of the run() function.


    Threads.@spawn while do_quit[] != 0x01
        if !isreadable(sock_mon.sock)
            do_quit[] = 0x01
        end

        ret = read_transport(sock_mon)
        type = ret[1]
        cmd = nothing
        len = 0
        flags = nothing
        if type<:ControlMessage
            len = ret[2]
            flags = ret[3]
            cmd = ret[4]
        else
            continue
        end

        if type === PlayMessage
            play[] = cmd.message
            println("> playing: ", play[])
        end
        if type === RateMessage
            playrate[] = cmd.message
            println("> resting: ", playrate[])
        end
        if type === PerturbationMessage
            perturbation[] = cmd
            println("> perturbing: ", cmd)
        end
        if type === QuitMessage
            do_quit[] = cmd.message
        end
    end

    # look for updates from REPL
    do_repl && Threads.@spawn while do_quit[] != 0x01
        if !isreadable(sock_repl.sock)
            do_quit[] = 0x01
        end

        ret = read_transport(sock_repl)
        type = ret[1]
        len = ret[2]
        flags = ret[3]
        field = ret[4]
        println("< received $type from REPL")

        if type <: AbstractConfig || type <: AbstractConstraint || type <: AbstractState
            sim_environment[field.id] = field
            println("< changed $(field.id) to: $field")
            if field isa AbstractConfig || field isa AbstractConstraint
                # forward it to monitor for update. Otherwise, if it is an AbstractState, it will be updated during the integration step.
                newval = packetize(field, 0x0002)
                write_transport(sock_mon, newval)
                println("> core pushed config update to monitor")
            end
        elseif field isa AskMessage
            if field.message in keys(sim_environment)
                resp = packetize(sim_environment[field.message], 0x0003)
                write_transport(sock_repl, resp)
            end
        end
    end

    time = sim.start_time

    while do_quit[] != 0x01 #isopen(sock)
        try
            if play[] == 0x01
                # update all dynamic systems

                step_sim!(
                    sat,
                    sat_config,
                    sim_environment,
                    sim.time_step,
                    time,
                    params;
                    exog = perturbation[],
                )
                # update timestep
                time = sim.start_time + Dates.Millisecond(1000*sat.elapsed_time)

                # serialize and send data
                sat_data = packetize(sat, 0x0000)
                write_transport(sock_mon, sat_data)

                # don't need to update Earth/Sun/etc positions as often as spacecraft state---save bandwidth.
                if sim.step_count % send_targets_per_sat_update == 0
                    for k in keys(sim_environment)
                        if sim_environment[k] isa AbstractTarget
                            t_data = packetize(sim_environment[k], 0x0000)
                            write_transport(sock_mon, t_data)
                        end
                    end
                end

                # control playback rate
                if playrate[] != 0xff
                    if playrate[] == 0x00
                        sleep(0.01)
                    elseif sim.step_count % playrate[] == 0
                        sleep(0.01)
                    end
                end

                if perturbation[].moment_duration > 0 &&
                   norm(perturbation[].moment_Body) != 0.0
                    perturbation[].moment_duration -= 1
                else
                    perturbation[].moment_Body = Vec3d(0.0)
                end
            else
                sleep(0.1)
            end

            sim.step_count += 1
        catch e
            if typeof(e) <: InterruptException
                println("exiting...")
                close(sock_mon.sock)
                return -1
            else
                rethrow(e)
                continue
            end
        end
        # print(" "*worm(counter%pwidth, pwidth)*" \r")
    end

    # close(sock.sock)
end

# same as run(), but no socket communication. For testing.
function run_free(; config_path::AbstractString = joinpath("config", "example.jsonc"))

    sim, sat, sat_config, sim_environment = load_config(Dict(load_jsonc(config_path)))

    inertia_B = diagm([5, 10, 13])*1e-2

    params = Dict{String,Any}(
        "I_Body" => inertia_B,
        "I_Body_inv" => inv(inertia_B),
        "max_angular_rate" => 2e-2, # rad/s
        "battery_max" => 84*3600.0,
        "storage_max" => 8e9,
        "irradiance" => 1360.0,
        "solar_panel_area" => 0.2*0.3*2,
        "solar_panel_efficiency" => 0.29,
        "downlink_rate" => 17e6,
        "do_J2" => true,
    )

    earth_state = EarthState(
        IDType(0x00f0),
        0,
        SatelliteToolboxTransformations.r_eci_to_ecef(
            SatelliteToolboxTransformations.J2000(),
            SatelliteToolboxTransformations.ITRF(),
            SatelliteToolboxBase.date_to_jd(sim.start_time),
            eop,
        ),
    )

    targets = filter(c -> c.second isa AbstractTarget, sim_environment)

    packlen = 1

    time = sim.start_time

    # run simulation
    while true
        # update all dynamic systems
        step_sim!(sat, sat_config, sim_environment, sim.time_step, time, params;)
        # update timestep
        time = sim.start_time + Dates.Millisecond(1000*sat.elapsed_time)

        # serialize data
        sat_data = packetize(sat, 0x0000)
        earth_data = packetize(earth_state, 0x0000)

        for target in values(targets)
            t_data = packetize(target, 0x0000)
        end


        sim.step_count += 1
    end
end
