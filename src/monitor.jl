using GLMakie
import Makie
import Printf, FileIO

@enum CamState free body_focused target_focused nadir

function basis()
    return [Vec3d(1.0, 0.0, 0.0), Vec3d(0.0, 1.0, 0.0), Vec3d(0.0, 0.0, 1.0)]
end

function load_earth_texture_to_ecef(path::String)
    texture = FileIO.load(path)
    W_uv = length(texture[1, :])
    texture = circshift(texture, (0, round(W_uv)/2))
    return texture
end

function get_sphere_mesh()
    n = 100
    sphere = Tesselation(Sphere(Point3f(0, 0, 0), 6371e3), n)
    m = uv_normal_mesh(sphere)
    return m
end

function get_body_mesh(scale::Float64 = 0.25)
    bodyscale = scale*6371e3
    box = Rect3d(bodyscale*Point3f(-0.15, -0.1, -0.05), bodyscale*Vec3f(0.3, 0.2, 0.1))
    ps = coordinates(box)
    fs = faces(box)
    ns = face_normals(ps, fs)
    m = GeometryBasics.mesh(ps, fs, normal = ns)
    return m
end

function format_clock(met::Float64, start_time::Dates.DateTime)
    days = Int(round(met÷(3600*24)))
    hours = Int(round(met - days*3600*24))÷3600
    minutes = Int(round(met - days*3600*24 - hours*3600))÷60
    seconds = met - days*3600*24 - hours*3600 - minutes*60

    utc = start_time + Dates.Second(met)
    utc_s = Dates.format(utc, "dd u YYYY HH:MM:SS.sss")

    return Printf.@sprintf "MET: %3u day      %2u:%2u:%06.3f\nUTC:  %s" days hours minutes seconds utc_s
end

function disable_makie_cam_keyboard!(cam::Makie.Camera3D)
    # println("disabling Makie keyboard camera control!")
    exclude = [:reposition_button, :reset, :rotation_button, :translation_button]
    for item in keys(cam.controls)
        if !(item in exclude)
            cam.controls[item] = Observable{Any}(false)
        end
    end

    # remember to call update_cam!(lscene.scene, cam) after to apply changes
end

function handle_new_config(path::AbstractString)
    d = Dict(load_jsonc(path))
    res = load_config(d)
    return res
end

function get_maps()
    maps_files = readdir(joinpath(asset_path, "maps"))
    maps = [joinpath(asset_path, "maps", m) for m in maps_files if !isdir(m)]
    return maps
end

function monitor(; config_path::AbstractString = joinpath("config", "example.jsonc"))

    textures = get_maps()
    fdi = findfirst(f -> occursin("map_diffuse", f), textures)
    texturek = isnothing(fdi) ? 1 : fdi
    texture = Observable(load_earth_texture_to_ecef(textures[texturek]))

    earth_mesh_v = get_sphere_mesh()
    body_mesh_val = get_body_mesh()

    sim, sat, sat_config, target_states, target_configs, constraints, modes =
        handle_new_config(config_path)

    helpstring = """
        Commands:
        - h:\tshow help
        - v:\tverbose
        - q:\tquit
        - e:\tchange earth texture
        - f:\tshow/hide reference frames
        - l:\tshow/hide labels
        - p:\tchange camera projection
        - c:\tchange camera focus
        - k:\tperturb attitude
        - spacebar:\t\tplay/pause
        - left/right:\tslower/faster
        """

    eci_scale = 1.4
    ecef_scale = 1.2
    body_scale = 0.1
    angular_velocity_scale = 30
    moment_scale = 0.2

    idlecolor = RGBAf(42/255, 133/255, 255/255)
    gscolor = RGBAf(157/255, 226/255, 107/255)
    suncolor = RGBAf(255/255, 201/255, 74/255)

    pointbufflen = 1000
    obbufflen = 2*10^4
    start_time = sim.start_time

    pos_ECI = Observable([Point3d(NaN) for k = 1:pointbufflen])
    ob_time = Observable([start_time + Dates.Second(k) for k = 1:obbufflen])
    batt = Observable([NaN for k = 1:obbufflen])
    stor = Observable([NaN for k = 1:obbufflen])
    power = Observable([NaN for k = 1:obbufflen])
    data = Observable([NaN for k = 1:obbufflen])
    rSO3 = Observable([NaN for k = 1:obbufflen])
    target_dir_ECI = Observable([Point3d(NaN) for k = 1:2])
    target_rel_ECI = Observable(Vec3d(NaN))
    q_ECEF_ECI = Observable(Makie.Quaternion(0.0, 0.0, 0.0, 1.0))
    pos_sun_ECI = Observable(Vec3d(0.0))
    moment_ECI = Observable([Vec3d(NaN) for k = 1:2])
    angular_velocity_ECI = Observable([Vec3d(NaN) for k = 1:2])
    met = Observable(Float64(0.0))
    met_label = lift(met) do l
        return format_clock(l, start_time)
    end

    body_frame_visible = Observable(false)
    axes_label_visible = Observable(false)
    gs_cones_visible = Observable(false)

    cmd_label = Observable("")

    GLMakie.activate!(title = "Hello, Greetings, and Welcome.")
    GLMakie.set_theme!(Makie.Theme(fonts = (; regular = "Menlo")))

    al = AmbientLight(RGBf(0.4, 0.4, 0.4))
    dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), pos_sun_ECI[])
    fig = GLMakie.Figure(size = (730, 580))
    display(fig)
    ax = LScene(
        fig[1:5, 1:3],
        show_axis = false,
        scenekw = (lights = [al, dl],),
        tellwidth = false,
        tellheight = false,
    )

    cam = cameracontrols(ax.scene)

    ax_batt = Axis(fig[1, 4], ylabel = "battery [%]", tellwidth = false)
    ylims!(ax_batt, 0.0, 102)
    ax_power = Axis(fig[2, 4], ylabel = "power [W]", tellwidth = false)
    ax_stor = Axis(fig[3, 4], ylabel = "onboard\ndata [%]", tellwidth = false)
    ylims!(ax_stor, 0.0, 102)
    ax_data = Axis(fig[4, 4], ylabel = "datarate\n[Mbps]", tellwidth = false)
    ax_dbug = Axis(fig[5, 4], ylabel = "SO(3) res.", tellwidth = false)
    # hidedecorations!(ob)
    hidespines!(ax_batt)
    hidespines!(ax_power)
    hidespines!(ax_stor)
    hidespines!(ax_data)
    hidespines!(ax_dbug)
    linkxaxes!(ax_batt, ax_power, ax_stor, ax_data, ax_dbug)

    tailcolor = Observable(suncolor)
    debuginfo = Observable(rich("."))
    configfilename = Observable("default")
    debugvisible = Observable(true)

    body_frame_v =
        body_scale*6371e3*[
            0.0 0.0 0.0;
            1.0 0.0 0.0;
            NaN NaN NaN;
            0.0 0.0 0.0;
            0.0 1.0 0.0;
            NaN NaN NaN;
            0.0 0.0 0.0;
            0.0 0.0 1.0
        ]

    inertial_frame_v = 10*eci_scale*body_frame_v
    fixed_frame_v = 10*ecef_scale*body_frame_v
    body_frame_color = [fill(:red, 3); fill(:green, 3); fill(:blue, 2)]
    inertial_frame_color = [fill(:red, 3); fill(:green, 3); fill(:blue, 2)]

    lines!(ax, pos_ECI, color = tailcolor, linewidth = 2)
    lines!(ax, target_dir_ECI, color = tailcolor, linewidth = 1)
    lines!(ax, moment_ECI, color = :black, linewidth = 2)
    lines!(ax, angular_velocity_ECI, color = :grey, linewidth = 1)
    earth_mesh = GLMakie.mesh!(
        ax,
        earth_mesh_v,
        color = texture,
        diffuse = Vec3f(0.6),
        specular = Vec3f(0.1),
        # invert_normals = true,
    )
    body_frame = GLMakie.lines!(
        ax,
        body_frame_v[:, 1],
        body_frame_v[:, 2],
        body_frame_v[:, 3],
        color = body_frame_color,
        visible = body_frame_visible,
    )
    body_mesh = GLMakie.mesh!(
        ax,
        body_mesh_val,
        diffuse = Vec3f(0.7),
        specular = Vec3f(0.3),
        color = :grey,
    )
    inertial_frame = GLMakie.lines!(
        ax,
        inertial_frame_v[:, 1],
        inertial_frame_v[:, 2],
        inertial_frame_v[:, 3],
        color = inertial_frame_color,
        visible = body_frame_visible,
    )
    fixed_frame = GLMakie.lines!(
        ax,
        fixed_frame_v[:, 1],
        fixed_frame_v[:, 2],
        fixed_frame_v[:, 3],
        color = inertial_frame_color,
        visible = body_frame_visible,
    )
    # C_ECI_ECEF0 = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/3600/24, eops)
    gs_pts = Observable([Point3d(NaN)])
    gs_label_pts = Observable([Point3d(NaN)])
    gs_labels = Observable(["."])
    gs_col = Observable([RGBAf(0.0, 0.0, 0.0, 0.0)])
    gs_dict = Dict{UInt16,GroundState}()

    eci_label_pts = 1.1*inertial_frame_v
    eci_labels =
        [rich("X", subscript("I")), rich("Y", subscript("I")), rich("Z", subscript("I"))]
    ecef_label_pts = 1.1*fixed_frame_v
    ecef_labels =
        [rich("X", subscript("F")), rich("Y", subscript("F")), rich("Z", subscript("F"))]

    gs_scatter = GLMakie.meshscatter!(
        ax,
        gs_pts,
        color = gs_col,
        markersize = 0.02*6371e3,
        alpha = 1.0,
    )

    gs_cone_len = 1.1*norm(sat.position_ECI) - 6371e3
    gs_cone_scaling = gs_cone_len / 6371e3
    gs_src_pts = [
        (1 + gs_cone_scaling)*target_states[k].position_ECI for
        k in keys(target_states) if target_states[k] isa GroundState
    ]
    # gs_src_pts = Observable([
    #     (1 + gs_cone_scaling)*target_states[k].position_ECI for
    #     k in keys(target_states) if target_states[k] isa GroundState
    # ])
    cone = Tesselation(
        Cone(Point3d(0.0), gs_cone_len*Point3d(0.0, 0.0, 1.0), gs_cone_len*2),
        40,
    )

    # note: this meshscatter approach intrinsically assumes all groundstations have 
    # the same elevation mask. Unless we get creative with nonuniform scale!().
    gs_src_rots = [
        dcm_to_quat(
            r_min_arc(Point3d(0.0, 0.0, -1.0), normalize(target_states[k].position_ECI))',
        ) for k in keys(target_states) if target_states[k] isa GroundState
    ]

    gs_cones = GLMakie.meshscatter!(
        ax,
        gs_src_pts,
        marker = cone,
        markersize = 1,
        alpha = 0.2,
        color = :blue,
        rotation = gs_src_rots,
        visible = gs_cones_visible,
    )

    gs_names = GLMakie.text!(
        ax,
        gs_label_pts,
        text = gs_labels,
        fontsize = 10.0,
        align = (:center, :center),
        visible = axes_label_visible,
    )
    eci_names = GLMakie.text!(
        ax,
        6371e3*eci_scale*1.1*basis(),
        text = eci_labels,
        fontsize = 10.0,
        align = (:center, :center),
        visible = axes_label_visible,
    )
    ecef_names = GLMakie.text!(
        ax,
        6371e3*ecef_scale*1.1*basis(),
        text = ecef_labels,
        fontsize = 10.0,
        align = (:center, :center),
        visible = axes_label_visible,
    )

    clock = text!(ax, 0, 1, text = met_label, align = (:left, :top), space = :relative)
    cmd = text!(
        ax,
        0,
        0,
        text = cmd_label,
        align = (:left, :bottom),
        space = :relative,
        fontsize = 10.0,
    )
    dbg = text!(
        ax,
        1,
        0,
        text = debuginfo,
        align = (:right, :bottom),
        space = :relative,
        fontsize = 10.0,
        # glowcolor=tailcolor,
        # glowwidth=1.0,
        visible = debugvisible,
    )

    batt_lines = lines!(ax_batt, ob_time, batt, color = :black, linewidth = 1)
    stor_lines = lines!(ax_stor, ob_time, stor, color = :black, linewidth = 1)
    power_lines = lines!(ax_power, ob_time, power, color = :black, linewidth = 1)
    data_lines = lines!(ax_data, ob_time, data, color = :black, linewidth = 1)
    rSO3_lines = lines!(ax_dbug, ob_time, rSO3, color = :black, linewidth = 1)

    nrate = 100_000
    secondly_debug = 2.0
    bytecount = 0
    last_met = 0.0
    lastratetime = time()

    packlen = length(packetize(SatelliteState(), 0x0000, UInt64(0)))
    packlen = 1
    headbuff = zeros(UInt8, 8)

    play = Ref(UInt8(0x01))
    do_quit = Ref(false)
    sleep_size = Ref(UInt8(0x01))
    last_proj = Ref(UInt8(0x00))
    cam_state = Ref(UInt8(0x00))
    sleepstep = 8

    rx_buff = zeros(UInt8, SOAP_MAX_BUFF_LEN)
    # play = 0x01

    disable_makie_cam_keyboard!(cameracontrols(ax.scene))

    # handle file drag/drop
    on(events(ax).dropped_files) do drop
        for file in drop
            if occursin(".json", file)
                configfilename[] = file
                notify(configfilename)

                res = handle_new_config(configfilename[])
                println(res)
                return
            end
        end
    end

    # camera controls
    on(pos_ECI) do pos
        if CamState(cam_state[]) == body_focused
            cam.lookat[] = pos[end]
            update_cam!(ax.scene, cam)
        elseif CamState(cam_state[]) == nadir
            cam.lookat[] = pos[end]
            cam.eyeposition[] = 2*pos[end]
            update_cam!(ax.scene, cam)
        elseif CamState(cam_state[]) == target_focused
            cam.lookat[] = pos[end]
            cam.eyeposition[] =
                (pos[end] + normalize(target_rel_ECI[])) -
                6371e3*normalize(target_rel_ECI[])
            update_cam!(ax.scene, cam)
        elseif CamState(cam_state[]) == free
            # no op 
        else
            @warn "undefined camera state!"
        end
    end

    # I/O:

    # udp:
    # sock =
    #     setup_client(SOAP_HOST, SOAP_MON_PORT, SOAP_HOST, SOAP_CORE_PORT)

    # tcp
    sock = setup_client(SOAP_HOST, SOAP_CORE_PORT)

    # unix
    # sock = setup_client(SOAP_UNIX_SOCK)

    # handle keyboard input:
    on(events(ax).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat
            if event.key == Keyboard.space
                if play[] == 0x01
                    play[] = 0x00
                    cmd_label[] = "||"
                elseif play[] == 0x00
                    play[] = 0x01
                    cmd_label[] = ">"
                else
                    println("oh no!!!")
                end
                writeval = packetize(PlayMessage(UInt8(play[])), 0x0000, UInt64(1))
                write_transport(sock, writeval)
                println("> sent command ", writeval, " to server")
                notify(cmd_label)
            end
            if event.key == Keyboard.left && play[] == 0x01
                sleep_size[] = max(0, sleep_size[] - sleepstep)
                cmd_label[] = "<<"
                writeval = packetize(RateMessage(UInt8(sleep_size[])), 0x0000, UInt64(2))
                write_transport(sock, writeval)
                println("> sent command ", writeval, " to server")
                notify(cmd_label)
            end
            if event.key == Keyboard.right && play[] == 0x01
                sleep_size[] = min(0xff, sleep_size[] + sleepstep)
                cmd_label[] = ">>"
                writeval = packetize(RateMessage(UInt8(sleep_size[])), 0x0000, UInt64(2))
                write_transport(sock, writeval)
                println("> sent command ", writeval, " to server")
                notify(cmd_label)
            end
            if event.key == Keyboard.k && play[] == 0x01
                # todo: configure perturbation magnitude and duration elsewhere
                pert = PerturbationMessage(Vec3d(1e-2*rand(3)), 5, Vec3d(0.0), 1)
                writeval = packetize(pert, 0x0000, UInt64(0))
                cmd_label[] = "kicked!"
                write_transport(sock, writeval)
                notify(cmd_label)
            end
            if event.key == Keyboard.escape
                cmd_label[] = ""
                notify(cmd_label)
            end
            if event.key == Keyboard.h
                cmd_label[] = helpstring
                notify(cmd_label)
            end
            if event.key == Keyboard.e
                texturek = mod(texturek, length(textures)) + 1
                texture[] = load_earth_texture_to_ecef(textures[texturek])
                notify(texture)
                # return Consume(true)
            end
            if event.key == Keyboard.l
                if body_frame_visible[]
                    axes_label_visible[] = !axes_label_visible[]
                    notify(axes_label_visible)
                end
            end
            if event.key == Keyboard.p
                if last_proj[] == 0x00
                    cam.settings.projectiontype=Makie.Orthographic
                    update_cam!(ax.scene, cam)
                    last_proj[] = 0x01
                    cmd_label[] = "projection: orthographic"
                    notify(cmd_label)
                elseif last_proj[] == 0x01
                    cam.settings.projectiontype=Makie.Perspective
                    update_cam!(ax.scene, cam)
                    last_proj[] = 0x00
                    cmd_label[] = "projection: perspective"
                    notify(cmd_label)
                end
            end
            if event.key == Keyboard.c
                cam_state[] = (cam_state[] + 1) % length(instances(CamState))
                cmd_label[] = "camera: $(string(CamState(cam_state[])))"
                notify(cmd_label)
            end
            if event.key == Keyboard.v
                debugvisible[] = !debugvisible[]
                notify(debugvisible)
            end
            if event.key == Keyboard.f
                body_frame_visible[] = !body_frame_visible[]
                notify(body_frame_visible)
                if !body_frame_visible[] && axes_label_visible[]
                    axes_label_visible[] = false
                    notify(axes_label_visible)
                end
            end
            if event.key == Keyboard.q
                do_quit[] = true
                cmd_label[] = "exiting..."
                writeval = packetize(QuitMessage(0x01), 0x0000, UInt64(2))
                write_transport(sock, writeval)
                println("> sent command ", writeval, " to server")
                notify(cmd_label)
            end
        end
    end

    while !do_quit[]
        try

            ret = read_transport(sock; buff = rx_buff)
            type = ret[1]
            count = 0
            flags = 0
            len = 0
            simdata = UInt8[]
            if type<:AbstractState
                len = ret[2]
                flags = ret[3]
                count = ret[4]
                simdata = ret[5]
            end

            if typeof(simdata) === PositionState # unused
                circshift!(pos_ECI[], -1)
                pos_ECI[][end] = Point3d(simdata.position_ECI)
                notify(pos_ECI)

                # update mission clock
                met[] = simdata.elapsed_time
                notify(met)

            elseif typeof(simdata) === AttitudeState # unused
            # GLMakie.rotate!(body_mesh, dcm_to_quat(simdata.attitude_ECI_Body'))
            # GLMakie.translate!(body_mesh, pos_ECI[][end])

            elseif typeof(simdata) === MagneticFieldState
                target_states[simdata.id] = simdata

            elseif typeof(simdata) === EarthState
                q_ECEF_ECI[] = dcm_to_quat(simdata.attitude_ECI_ECEF)
                notify(q_ECEF_ECI)
                GLMakie.rotate!(earth_mesh, q_ECEF_ECI[])
                GLMakie.rotate!(fixed_frame, q_ECEF_ECI[])
                GLMakie.rotate!(ecef_names, q_ECEF_ECI[])
                gs_cones_visible[] && GLMakie.rotate!(gs_cones, q_ECEF_ECI[])

            elseif typeof(simdata) === SunState
                GLMakie.set_directional_light!(
                    ax,
                    color = RGBf(243/255, 241/255, 218/255),
                    direction = -Vec3d(simdata.position_ECI),
                )
                target_states[simdata.id] = simdata

            elseif typeof(simdata) === GroundState
                gs_dict[simdata.id] = simdata
                target_states[simdata.id] = simdata

            elseif typeof(simdata) === SatelliteState
                # update mission clock
                met[] = simdata.elapsed_time
                notify(met)

                # update position
                circshift!(pos_ECI[], -1)
                pos_ECI[][end] = Point3d(simdata.position_ECI)
                notify(pos_ECI)

                circshift!(batt[], -1)
                batt[][end] = simdata.battery_level / sat_config.power_battery_max * 100
                notify(batt)
                circshift!(stor[], -1)
                stor[][end] = simdata.storage_level / sat_config.data_storage_max * 100
                notify(stor)
                circshift!(power[], -1)
                power[][end] = simdata.net_power
                notify(power)
                circshift!(data[], -1)
                data[][end] = simdata.net_data / 1e6
                notify(data)
                circshift!(rSO3[], -1)
                rSO3[][end] = residualSO3(simdata.attitude_ECI_Body')
                notify(rSO3)
                circshift!(ob_time[], -1)
                ob_time[][end] = start_time + Dates.Millisecond(1000*simdata.elapsed_time)
                notify(ob_time)

                # update attitude
                q = dcm_to_quat(simdata.attitude_ECI_Body')
                GLMakie.rotate!(body_mesh, q)
                GLMakie.translate!(body_mesh, pos_ECI[][end])

                GLMakie.rotate!(body_frame, q)
                GLMakie.translate!(body_frame, pos_ECI[][end])

                # update color and target direction
                mode_conf = modes[findfirst(x -> x.id == simdata.mode, modes)]
                target_rel_ECI[] = simdata.attitude_ECI_Body*mode_conf.direction_Body
                notify(target_rel_ECI)

                # mode_conf = mode_table[simdata.mode]
                tailcolor[] = mode_conf.color

                if simdata.target == typemax(IDType)
                    target_dir_ECI[] = [Point3d(NaN), Point3d(NaN)]
                elseif isa(target_states[simdata.target], SunState)
                    target_dir_ECI[] = [
                        simdata.position_ECI,
                        simdata.position_ECI +
                        0.33*6371e3*normalize(
                            target_states[simdata.target].position_ECI -
                            simdata.position_ECI,
                        ),
                    ]
                elseif isa(target_states[simdata.target], GroundState)
                    target_dir_ECI[] =
                        [simdata.position_ECI, target_states[simdata.target].position_ECI]

                elseif isa(target_states[simdata.target], MagneticFieldState)
                    target_dir_ECI[] = [
                        simdata.position_ECI,
                        0.33*6371e3*normalize(target_states[simdata.target].direction_ECI) + simdata.position_ECI,
                    ]
                else
                    println("unidentified target type!")
                end

                angular_velocity_ECI[] = [
                    simdata.position_ECI,
                    simdata.position_ECI +
                    angular_velocity_scale*6371e3*simdata.angular_velocity_ECI_Body,
                ]
                moment_ECI[] = [
                    simdata.position_ECI,
                    simdata.position_ECI +
                    moment_scale*6371e3*normalize(
                        simdata.attitude_ECI_Body*simdata.net_moment_Body,
                    ),
                ]

                notify(tailcolor)
                notify(target_dir_ECI)
                notify(angular_velocity_ECI)
                notify(moment_ECI)

                target_name =
                    (simdata.target == typemax(IDType)) ? "nothing" :
                    find_config(target_states[simdata.target], target_configs).name

                debuginfo[] = rich(
                    "mode: ",
                    rich(
                        String(Printf.@sprintf "%20s" mode_conf.name),
                        "\n",
                        color = tailcolor[],
                    ),
                    "target: ",
                    rich(String(Printf.@sprintf "%20s" string(target_name)), "\n"),
                    "config: ",
                    rich(String(Printf.@sprintf "%20s" configfilename[]), "\n"),
                    "time gain: ",
                    rich(
                        String(
                            Printf.@sprintf "%18u:1" Int64(
                                round((met[] - last_met)/(time() - lastratetime)),
                            )
                        ),
                    ),
                )

                # debuginfo[] = rich("target: ", "nothin")
                notify(debuginfo)
            else
                println("unknown message type $(typeof(simdata))")
            end

            # update all groundstations
            gs_pts[] = [gs_dict[key].position_ECI for key in keys(gs_dict)]
            notify(gs_pts)
            gs_label_pts[] = [1.1*gs_dict[key].position_ECI for key in keys(gs_dict)]
            notify(gs_label_pts)
            gs_labels[] = [
                string(find_config(gs_dict[key], target_configs).name) for
                key in keys(gs_dict)
            ]
            notify(gs_labels)
            gs_col[] = [gs_dict[key].visible ? gscolor : idlecolor for key in keys(gs_dict)]
            notify(gs_col)

            bytecount += len + 8

            if time() - lastratetime > secondly_debug
                # Printf.@printf "rate: %0.3f MB/s\n" 1e-6*bytecount/(time() - lastratetime)
                # Printf.@printf "time gain: %u:1 \n" Int64(
                #     round((met[] - last_met)/(time() - lastratetime)),
                # )
                last_met = met[]
                lastratetime = time()
                bytecount = 0
                GLMakie.reset_limits!(ax_power)
                GLMakie.reset_limits!(ax_data)
            end
        catch e
            if typeof(e) <: InterruptException
                println("exiting...")
                GLMakie.closeall()
                close(sock.sock)
                return -1
            else
                rethrow(e)
                flush(sock)
                continue
            end
        end
    end

    close(sock.sock)
    cmd_label[] = "goodbye"
    notify(cmd_label)
    GLMakie.closeall()
    return 1
end

function stringify(something::AbstractArray)::String
    res = ""
    sep = ""
    for (k, el) in enumerate(something)
        if k > 1
            sep = ", "
        else
            sep = ""
        end
        res = res * sep * String(Printf.@sprintf "%.4e" something[k])
    end
    return res
end

mutable struct CSVLine
    utc::Dates.DateTime
    met::Float64
    t::Float64
    count::UInt64
    pos::Vec3d
    vel::Vec3d
    att::Mat3d
    batt::Float64
    stor::Float64
    mode::IDType
    target::IDType
end
function CSVLine()
    return CSVLine(
        Dates.DateTime(0),
        0.0,
        0.0,
        0x0,
        Vec3d(NaN),
        Vec3d(NaN),
        Mat3d(fill(NaN, (3, 3))),
        0.0,
        0.0,
        typemax(IDType),
        typemin(IDType),
    )
end
function update!(line::CSVLine, sat::SatelliteState)
    line.met = sat.elapsed_time
    line.pos = sat.position_ECI
    line.vel = sat.velocity_ECI
    line.att = sat.attitude_ECI_Body
    line.batt = sat.battery_level
    line.stor = sat.storage_level
    line.mode = sat.mode
    line.target = sat.target
end
function iscomplete(line::CSVLine)
    emptyl = CSVLine()
    if line.utc != emptyl.utc &&
       line.met != emptyl.met &&
       line.t != emptyl.t &&
       line.count != emptyl.count &&
       any(line.pos .!= emptyl.pos) &&
       any(line.vel .!= emptyl.vel) &&
       any(line.att .!= emptyl.att) &&
       line.batt != emptyl.batt &&
       line.stor != emptyl.stor &&
       line.mode != emptyl.mode &&
       line.target != emptyl.target
        return true
    else
        return false
    end
end

function stringify(line::CSVLine)
    out = "\n"
    out *= Dates.format(line.utc, "dd u YYYY HH:MM:SS.sss")
    out *= ", " * String(Printf.@sprintf "%.4e" line.met)
    out *= ", " * String(Printf.@sprintf "%.4e" line.t)
    out *= ", " * String(Printf.@sprintf "%u" line.count)
    out *= ", " * stringify(line.pos)
    out *= ", " * stringify(line.vel)
    out *= ", " * stringify(line.att)
    out *= ", " * String(Printf.@sprintf "%.4e" line.batt)
    out *= ", " * String(Printf.@sprintf "%.4e" line.stor)
    out *= ", " * String(Printf.@sprintf "%u" line.mode)
    out *= ", " * String(Printf.@sprintf "%u" line.target)
    return out
end

function write_csv(
    outfile::AbstractString;
    infile::AbstractString = joinpath("config", "example.jsonc"),
)

    sock = setup_client(SOAP_HOST, SOAP_CORE_PORT)

    outf = open(outfile, "a")

    sim, sat, sat_config, target_states, target_configs, constraints, modes =
        handle_new_config(infile)


    start_time = sim.start_time

    println("socket: ", sock)

    nrate = 100_000
    secondly_debug = 2.0
    bytecount = 0
    last_met = 0.0
    lastratetime = time()

    packlen = length(packetize(SatelliteState(), 0x0000, UInt64(0)))
    packlen = 1
    headbuff = zeros(UInt8, 8)
    buff = zeros(UInt8, packlen)

    play = Ref(UInt8(0x01))
    sleep = Ref(UInt8(0xff))
    last_proj = Ref(UInt8(0x00))
    sleepstep = 8

    header = "t [UTC], t [MET], t [s], counter [], pos_ECI_x [m], pos_ECI_y [m], pos_ECI_z [m], vel_ECI_x [m], vel_ECI_y [m], vel_ECI_z [m], att_ECI_Body_11 [], attitude_ECI_Body_12 [], attitude_ECI_Body_13 [], attitude_ECI_Body_21 [], attitude_ECI_Body_22 [], attitude_ECI_Body_23 [], attitude_ECI_Body_31 [], attitude_ECI_Body_32 [], attitude_ECI_Body_33 [], battery [J], data storage [b], mode, target"
    last_count = 0

    write(outf, header)

    this_line = CSVLine()

    rx_buff = zeros(UInt8, SOAP_MAX_BUFF_LEN)
    met = 0.0
    while !eof(sock.sock)
        try

            ret = read_transport(sock; buff = rx_buff)
            type = ret[1]
            count = 0
            flags = 0
            len = 0
            simdata = UInt8[]
            if type<:AbstractState
                len = ret[2]
                flags = ret[3]
                count = ret[4]
                simdata = ret[5]
            end

            if typeof(simdata) === EarthState
                # q_ECEF_ECI[] = dcm_to_quat(simdata.attitude_ECI_ECEF)
                # notify(q_ECEF_ECI)
                # GLMakie.rotate!(earth_mesh, q_ECEF_ECI[])

            elseif typeof(simdata) === SunState
                target_states[simdata.id] = simdata

            elseif typeof(simdata) === GroundState
                # gs_dict[simdata.id] = simdata
                target_states[simdata.id] = simdata
                # met[] = simdata.elapsed_time

            elseif typeof(simdata) === SatelliteState
                # update mission clock
                update!(this_line, simdata)
                met = simdata.elapsed_time
                this_line.count = this_count
                this_line.t = simdata.elapsed_time
                this_line.utc = start_time + Dates.Millisecond(1000*simdata.elapsed_time)
            else
                println("unknown message type")
            end

            # handle writing the CSV
            if iscomplete(this_line)
                # write the line to the file
                write(outf, stringify(this_line))
                # clear this_line
                this_line = CSVLine()
            end

            bytecount += len + 8

            if time() - lastratetime > secondly_debug
                Printf.@printf "\nstep, MET: %u, %.2f days\n" this_count met/3600/24
                Printf.@printf "rate:      %0.3f MB/s\n" 1e-6*bytecount/(
                    time() - lastratetime
                )
                Printf.@printf "time gain: %u:1 \n" Int64(
                    round((met - last_met)/(time() - lastratetime)),
                )
                Printf.@printf "file size: %u MB\n" Int64(round(stat(outfile).size/1e6))
                last_met = met
                lastratetime = time()
                bytecount = 0
            end
        catch e
            if typeof(e) <: InterruptException
                println("exiting...")
                close(sock)
                close(outf)
                break
            else
                rethrow(e)
                flush(sock)
                continue
            end
        end
    end
end
