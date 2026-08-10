# global const unixsock"tname = "127.0.0.1"

function run(; config_path::AbstractString = joinpath("config", "example.jsonc"))

    sim, sat, sat_config, targets, target_configs, constraints, modes =
        load_config(Dict(load_jsonc(config_path)))

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

    println(collect(values(target_configs)))

    println((modes, targets, target_configs, constraints))

    # udp 
    # sock =
    #     setup_server(SOAP_HOST, SOAP_CORE_PORT, SOAP_HOST, SOAP_MON_PORT)

    # tcp 
    sock = setup_server(SOAP_HOST, SOAP_CORE_PORT)
    println("connected.")

    # unix
    # sock = setup_server(SOAP_UNIX_SOCK)

    play = Ref(UInt8(0x01))
    do_quit = Ref(UInt8(0x00))
    playrate = Ref(UInt8(0x01))
    perturbation = Ref(PerturbationMessage(Vec3d(0.0), 0.0, Vec3d(0.0), 0.0))
    # perturb_moment = Ref(Vec3d(0.0))
    # perturb_moment_count = Ref(UInt32(0x0000))
    packlen = 1
    headbuff = zeros(UInt8, 8)
    buff = zeros(UInt8, packlen)

    # handle received commands
    # todo: factor this out of the run() function.


    @async while do_quit[] != 0x01
        if !isreadable(sock.sock)
            do_quit[] = 0x01
        end

        ret = read_transport(sock)
        type = ret[1]
        cmd = nothing
        len = 0
        flags = nothing
        if type<:ControlMessage
            flags = ret[2]
            len = ret[3]
            cmd = ret[4]
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

    time = sim.start_time

    while do_quit[] != 0x01 #isopen(sock)
        try
            if play[] == 0x01
                # update all dynamic systems
                step_sim!(
                    sat,
                    sat_config,
                    targets,
                    target_configs,
                    constraints,
                    modes,
                    sim.time_step,
                    time,
                    params;
                    exog = perturbation[],
                )

                # update timestep
                time = sim.start_time + Dates.Millisecond(1000*sat.elapsed_time)

                # serialize and send data
                sat_data = packetize(sat, 0x0000, sim.step_count)
                write_transport(sock, sat_data)

                for target in values(targets)
                    t_data = packetize(target, 0x0000, sim.step_count)
                    write_transport(sock, t_data)
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
                close(sock.sock)
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

    sim, sat, sat_config, targets, target_configs, constraints, modes =
        load_config(Dict(load_jsonc(config_path)))

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

    packlen = 1

    time = sim.start_time

    # run simulation
    while true
        # update all dynamic systems
        step_sim!(
            sat,
            sat_config,
            targets,
            target_configs,
            constraints,
            modes,
            sim.time_step,
            time,
            params;
        )
        # update timestep
        time = sim.start_time + Dates.Millisecond(1000*sat.elapsed_time)

        # serialize data
        sat_data = packetize(sat, 0x0000, sim.step_count)
        earth_data = packetize(earth_state, 0x0000, sim.step_count)

        for target in values(targets)
            t_data = packetize(target, 0x0000, sim.step_count)
        end


        sim.step_count += 1
    end
end
