using GLMakie, GeometryBasics, LinearAlgebra
import Makie
import SatelliteToolboxTransformations
import Sockets, Serialization
import Printf, FileIO

function load_earth_texture_to_ecef(path::String)
    texture = FileIO.load(path)
    texture = texture'
    W_uv = length(texture[:,1])
    texture = circshift(texture, (round(W_uv)/2, 0))
    return texture
end

function get_sphere_mesh(texture)
    nθ = length(texture[:,1])
    nφ = length(texture[1,:])
    θ = range(0, stop=2π, length=nθ)
    φ = range(0, stop=π, length=nφ)
    mesh_X = zeros(nθ, nφ)
    mesh_Y = zeros(nθ, nφ)
    mesh_Z = zeros(nθ, nφ)
    for i in 1:nθ
        for j in 1:nφ
            # a point in the mesh
            p = SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS*[cos(θ[i])*sin(φ[j]); sin(θ[i])*sin(φ[j]); cos(φ[j])]
            mesh_X[i,j] = p[1]
            mesh_Y[i,j] = p[2]
            mesh_Z[i,j] = p[3]
        end
    end
    return mesh_X, mesh_Y, mesh_Z
end

function test_texture()
    texture = load_earth_texture_to_ecef("/Users/thanasi/Documents/OpenSOAP/assets/map_diffuse.png")
    GLMakie.activate!()
    nθ = length(texture[:,1])
    nφ = length(texture[1,:])
    θ = range(0, stop=2π, length=nθ)
    φ = range(0, stop=π, length=nφ)
    mesh_X = zeros(nθ, nφ)
    mesh_Y = zeros(nθ, nφ)
    mesh_Z = zeros(nθ, nφ)
    for i in 1:nθ
        for j in 1:nφ
            # a point in the mesh
            p = SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS*[cos(θ[i])*sin(φ[j]); sin(θ[i])*sin(φ[j]); cos(φ[j])]
            mesh_X[i,j] = p[1]
            mesh_Y[i,j] = p[2]
            mesh_Z[i,j] = p[3]
        end
    end
    f,ax,pl = GLMakie.mesh(mesh_X,mesh_Y,mesh_Z,color=texture, invert_normals=true)
    display(f)
    return f
end

function get_body_mesh(scale::Float64=0.25)
    bodyscale = scale*6371e3
    box = Rect3d(bodyscale*Point3f(-0.15,-0.1,-0.05), bodyscale*Vec3f(0.3,0.2,0.1))
    ps = coordinates(box)
    fs = faces(box)
    ns = face_normals(ps, fs)
    m = GeometryBasics.mesh(ps, fs, normal=ns)
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
    println("disabling Makie keyboard camera control!")
    exclude = [:reposition_button, :reset, :rotation_button, :translation_button]
    for item in keys(cam.controls)
        if !(item in exclude)
            cam.controls[item] = Observable{Any}(false)
        end
    end
    
    # remember to call update_cam!(lscene.scene, cam) after to apply changes
end

function show()
    # server = Sockets.listen(unixsockname)
    sock = Sockets.UDPSocket()
    # bind(sock, unixsockname)
    sock = Sockets.connect(unixsockname)
    textureprefix = joinpath("/", "Users", "thanasi", "Documents", "OpenSOAP", "assets")
    textures = ["map_diffuse.png", "map_pol1.png", "map_pol2.png", "map_marble.png", "map_bathy.png", "map_bw.png", "map_snow.jpeg", "map_veggie.jpeg", "map_tissot.jpg", "map_cities.png"]
    texturek = 1
    texture = Observable(load_earth_texture_to_ecef(joinpath(textureprefix, textures[texturek])))
    X_ECI, Y_ECI, Z_ECI = get_sphere_mesh(texture[])
    
    tmp_groundstations = [
        Dict{String,Any}("name"=>"Svalbard",     "operator"=>"KSAT", "latitude"=>78.22, "longitude"=>15.65, "altitude"=>1718.0, "elevation"=>10.0, "viewable"=>false),
        Dict{String,Any}("name"=>"Punta Arenas", "operator"=>"AWS",  "latitude"=>-53.1666667, "longitude"=>-70.933333, "altitude"=>34.0, "elevation"=>10.0, "viewable"=>false)
    ]
    
    modes_table = mode_table(load_mode_config(""))
    
    pointbufflen = 1000
    obbufflen = 10^4
    start_time = Dates.DateTime(2026,3,20,1,2,3)
    pos_ECI = Observable([Point3d(NaN) for k in 1:pointbufflen])
    ob_time = Observable([start_time + Dates.Second(k) for k in 1:obbufflen])
    batt = Observable([NaN for k in 1:obbufflen])
    data = Observable([NaN for k in 1:obbufflen])
    target_dir_ECI = Observable([Point3d(NaN) for k in 1:2])
    q_ECEF_ECI = Observable(Makie.Quaternion(0.0,0.0,0.0,1.0))
    pos_sun_ECI = Observable(Vec3d(0.0))
    met = Observable(Float64(0.0))
    met_label = lift(met) do l
        return format_clock(l, start_time)
    end
    
    body_frame_visible = Observable(false)
    
    cmd_label = Observable("")
    helpstring = "Commands:\n- h:\tshow help\n- v:\tverbose\n- e:\tchange earth texture\n- f:\tshow/hide Body frame\n- p:\tchange camera projection\n- spacebar:\t\tplay/pause\n- left/right:\tslower/faster\n"
    
    body_mesh_val = get_body_mesh()
    
    GLMakie.activate!(title="Hello Smallsat")
    GLMakie.set_theme!(font="OCR-B")
    al = AmbientLight(RGBf(0.4, 0.4, 0.4))
    dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), pos_sun_ECI[])
    fig = GLMakie.Figure(size=(730, 580))
    display(fig)
    ax = LScene(fig[1:5,1:3], show_axis=false, scenekw=(lights=[al,dl],), tellwidth=false)
    ax_batt = Axis(fig[1,4], ylabel="battery", ylabelfont="OCR-B", tellwidth=false)
    ylims!(ax_batt, 0.0, 102)
    ax_data = Axis(fig[2,4], ylabel="onboard\ndata", ylabelfont="OCR-B", tellwidth=false)
    ylims!(ax_data, 0.0, 102)
    # hidedecorations!(ob)
    hidespines!(ax_batt)
    hidespines!(ax_data)
    # xlims!(ax, [-6e3, 6e3])
    # ylims!(ax, [-6e3, 6e3])
    # zlims!(ax, [-6e3, 6e3])
    
    idlecolor = RGBAf(42/255, 133/255, 255/255)
    gscolor = RGBAf(157/255, 226/255, 107/255)
    suncolor = RGBAf(255/255, 201/255, 74/255)
    
    tmp_modecolor = Dict{UInt16, RGBAf}(
        0x0003 => gscolor,
        0x0002 => suncolor,
        0x0001 => idlecolor
    )
    tmp_modename = Dict{UInt16, String}(
        0x0003 => "groundstation",
        0x0002 => "sun",
        0x0001 => "idle"
    )
    
    tailcolor = Observable(suncolor)
    debuginfo = Observable(rich("."))
    debugvisible = Observable(true) # todo: create a ControlMessage to enable/disable this
    
    body_frame_v = 0.1*6371e3*[ 0.0 0.0 0.0;
                                1.0 0.0 0.0;
                                NaN NaN NaN;
                                0.0 0.0 0.0;
                                0.0 1.0 0.0;
                                NaN NaN NaN;
                                0.0 0.0 0.0;
                                0.0 0.0 1.0]
    body_frame_color = [fill(:red, 3); fill(:green, 3); fill(:blue, 2)]
    
    lines!(ax, pos_ECI, color=tailcolor, linewidth=2)
    lines!(ax, target_dir_ECI, color=tailcolor, linewidth=1)
    earth_mesh = GLMakie.surface!(
        ax, 
        X_ECI, Y_ECI, Z_ECI,
        color=texture,
        diffuse=Vec3f(0.6),
        specular=Vec3f(0.1),
        invert_normals=true
    )
    body_frame = GLMakie.lines!(
        ax,
        body_frame_v[:,1], body_frame_v[:,2], body_frame_v[:,3],
        color=body_frame_color,
        visible=body_frame_visible
    )
    body_mesh = GLMakie.mesh!(
        ax,
        body_mesh_val,
        diffuse=Vec3f(0.7),
        specular=Vec3f(0.3),
        color=:grey
    )
    # C_ECI_ECEF0 = r_ecef_to_eci(ITRF(), J2000(), t_jd_s/3600/24, eops)
    gs_pts = Observable([Point3d(NaN)])
    gs_label_pts = Observable([Point3d(NaN)])
    gs_labels = Observable(["."])
    gs_col = Observable([RGBAf(0.0,0.0,0.0,0.0)])
    gs_dict = Dict{UInt16, TargetState}()
    
    gs_scatter = GLMakie.meshscatter!(
        ax,
        gs_pts,
        color=gs_col,
        markersize=0.03*6371e3
    )
    
    gs_names = GLMakie.text!(
        ax, 
        gs_label_pts,
        text=gs_labels,
        font="OCR-B",
        fontsize=10.0,
        align=(:center, :center)
    )
    clock = text!(
        ax, 
        0,1, 
        text=met_label, 
        align=(:left, :top), 
        space=:relative,
        font="OCR-B"
    )
    cmd = text!(
        ax,
        0,0,
        text=cmd_label,
        align=(:left, :bottom),
        space=:relative,
        font="OCR-B",
        fontsize=10.0
    )
    dbg = text!(
        ax,
        1,0,
        text=debuginfo,
        align=(:right,:bottom),
        space=:relative,
        font="OCR-B",
        fontsize=10.0,
        # glowcolor=tailcolor,
        # glowwidth=1.0,
        visible=debugvisible
    )
    
    batt_lines = lines!(
        ax_batt,
        ob_time,
        batt,
        color=:black,
        linewidth=1
    )
    data_lines = lines!(
        ax_data,
        ob_time,
        data,
        color=:black,
        linewidth=1
    )
    
    # println("server: ", server)
    # sock = Sockets.accept(server)
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
    # play = 0x01
    
    disable_makie_cam_keyboard!(cameracontrols(ax.scene))
    
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
                write(sock, writeval)
                println("> sent command ", writeval, " to server")
                notify(cmd_label)
            end
            if event.key == Keyboard.left && play[] == 0x01
                sleep[] = max(0, sleep[] - sleepstep)
                writeval = packetize(RateMessage(UInt8(sleep[])), 0x0000, UInt64(2))
                cmd_label[] = "<<"
                write(sock, writeval)
                println("> sent command ", writeval, " to server")
                notify(cmd_label)
            end
            if event.key == Keyboard.right && play[] == 0x01
                sleep[] = min(0xff, sleep[] + sleepstep)
                writeval = packetize(RateMessage(UInt8(sleep[])), 0x0000, UInt64(2))
                cmd_label[] = ">>"
                write(sock, writeval)
                println("> sent command ", writeval, " to server")
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
                texture[] = load_earth_texture_to_ecef(joinpath(textureprefix, textures[texturek]))
                notify(texture)
                # return Consume(true)
            end
            if event.key == Keyboard.p
                if last_proj[] == 0x00
                    cam = cameracontrols(ax.scene)
                    cam.settings.projectiontype=Makie.Orthographic
                    update_cam!(ax.scene, cam) 
                    last_proj[] = 0x01
                    cmd_label[] = "projection: orthographic"
                    notify(cmd_label)
                elseif last_proj[] == 0x01
                    cam = cameracontrols(ax.scene)
                    cam.settings.projectiontype=Makie.Perspective
                    update_cam!(ax.scene, cam)
                    last_proj[] = 0x00
                    cmd_label[] = "projection: perspective"
                    notify(cmd_label)
                end
            end
            if event.key == Keyboard.v
                debugvisible[] = !debugvisible[]
                notify(debugvisible)
            end
            if event.key == Keyboard.f
                body_frame_visible[] = !body_frame_visible[]
                notify(body_frame_visible)
            end
        end
    end
    
    while true
        try
            nbhead = readbytes!(sock, headbuff)
            type, flags, len = behead(headbuff)
            
            buff = zeros(UInt8, len)
            nbdata = readbytes!(sock, buff)
            count, simdata = des(type, buff)
            
            if typeof(simdata) === PositionState # unused
                circshift!(pos_ECI[], -1)
                pos_ECI[][end] = Point3d(simdata.position_ECI)
                notify(pos_ECI)
                
                # update mission clock
                met[] = simdata.elapsed_time
                notify(met)
                
            elseif typeof(simdata) === AttitudeState # unused
                GLMakie.rotate!(body_mesh, dcm_to_quat(simdata.attitude_ECI_Body'))
                GLMakie.translate!(body_mesh, pos_ECI[][end])
                
            elseif typeof(simdata) === EarthState
                q_ECEF_ECI[] = dcm_to_quat(simdata.attitude_ECI_ECEF')
                notify(q_ECEF_ECI)
                GLMakie.rotate!(earth_mesh, q_ECEF_ECI[])
                
            elseif typeof(simdata) === SunState
                GLMakie.set_directional_light!(
                    ax, 
                    color=RGBf(243/255, 241/255, 218/255), 
                    direction=-Vec3d(simdata.position_ECI)
                )
                
            elseif typeof(simdata) === TargetState
                gs_dict[simdata.id] = simdata
                
            elseif typeof(simdata) === SatelliteState
                # update mission clock
                met[] = simdata.elapsed_time
                notify(met)
                
                # update position
                circshift!(pos_ECI[], -1)
                pos_ECI[][end] = Point3d(simdata.position_ECI)
                notify(pos_ECI)
                
                circshift!(batt[], -1)
                # todo: replace with lookup of battery capacity
                batt[][end] = simdata.battery_level / 84/36.0 
                notify(batt)
                
                circshift!(data[], -1)
                # todo: replace with lookup of storage capacity
                data[][end] = simdata.storage_level / 8e9 * 100
                notify(data)
                
                circshift!(ob_time[], -1)
                ob_time[][end] = start_time + Dates.Millisecond(1000*simdata.elapsed_time)
                notify(ob_time)
                
                # update attitude
                GLMakie.rotate!(body_mesh, dcm_to_quat(simdata.attitude_ECI_Body'))
                GLMakie.translate!(body_mesh, pos_ECI[][end])
                
                GLMakie.rotate!(body_frame, dcm_to_quat(simdata.attitude_ECI_Body'))
                GLMakie.translate!(body_frame, pos_ECI[][end])
                
                # update color and target direction
                mode_conf = modes_table[simdata.mode]
                # tailcolor[] = tmp_modecolor[simdata.mode]
                tailcolor[] = mode_conf.color
                # tailcolor[] = :black
                if mode_conf.target_type === TargetState
                    target_dir_ECI[] = [simdata.position_ECI, simdata.target_ECI]
                elseif mode_conf.target_type === SunState # sun
                    target_dir_ECI[] = [simdata.position_ECI, simdata.position_ECI + 0.5*6371e3*normalize(simdata.target_ECI - simdata.position_ECI)]
                elseif mode_conf.target_type === Nothing # idle
                    target_dir_ECI[] = [Point3d(NaN), Point3d(NaN)]
                else
                    println("unidentified target type: ", typeof(mode_conf.target_type))
                end
                notify(tailcolor)
                notify(target_dir_ECI)
                
                # debuginfo[] = rich("target: ", rich(tmp_modename[simdata.mode], color=tmp_modecolor[simdata.mode]))
                debuginfo[] = rich(
                    "mode:   ", rich(String(Printf.@sprintf "%20s" mode_conf.name), "\n",  color=tailcolor[]),
                    "target: ", rich(String(Printf.@sprintf "%20s" string(mode_conf.target_type)))
                )
                
                # debuginfo[] = rich("target: ", "nothin")
                notify(debuginfo)
            else
                println("unknown message type")
            end
            
            # update all groundstations
            gs_pts[] = [gs_dict[key].position_ECI for key in keys(gs_dict)]
            notify(gs_pts)
            gs_label_pts[] = [1.1*gs_dict[key].position_ECI for key in keys(gs_dict)]
            notify(gs_label_pts)
            gs_labels[] = [string(gs_dict[key].id) for key in keys(gs_dict)]
            notify(gs_labels)
            gs_col[] = [gs_dict[key].visible ? gscolor : idlecolor for key in keys(gs_dict)]
            notify(gs_col)
            
            bytecount += len + 8
            
            @async GLMakie.reset_limits!(ax_batt)
            @async GLMakie.reset_limits!(ax_data)
            
            if time() - lastratetime > secondly_debug
                Printf.@printf "rate: %0.3f MB/s\n" 1e-6*bytecount/(time() - lastratetime)
                Printf.@printf "time gain: %u:1 \n" Int64(round((met[] - last_met)/(time() - lastratetime)))
                last_met = met[]
                lastratetime = time()
                bytecount = 0
            end
        catch e
            if typeof(e) <: InterruptException
                println("exiting...")
                close(sock)
                break
            else
                rethrow(e)
                flush(sock)
                continue
            end
        end
    end
end