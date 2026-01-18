using OpenSOAP, LinearAlgebra
import SatelliteToolboxBase, SatelliteToolboxTransformations, SatelliteToolboxGeomagneticField, SatelliteToolboxPropagators
using GLMakie
import CairoMakie

import Statistics

import DelimitedFiles
import DataFrames

import Dates, DateFormats

function orbinterp(time::Dates.DateTime, source_time::Vector{Dates.DateTime}, source_pos::Vector{Point3f})
    # next = searchsortedfirst(Dates.datetime2epochms(time), Dates.datetime2epochms.(source_time))
    next = searchsortedfirst(source_time, time)

    if next == length(source_time)
        # extrapolate end
    elseif next == 1
        # extrapolate start

    else
        # interpolate
        pos_start = source_pos[next - 1]
        pos_end = source_pos[next]
        time_start = source_time[next - 1]
        time_end = source_time[next]

        res = pos_start + (pos_end - pos_start) .* ((time - time_start) / (time_end - time_start))
    end

    return res
end



@doc raw"""
    mag_err()

This is a test of on-orbit magnetometer error.
"""
function mag_err()

    println("loading transformations...")
    eops = SatelliteToolboxTransformations.fetch_iers_eop()

    println("loading source data...")
    att_data, att_head = DelimitedFiles.readdlm("/Users/thanasi/Documents/IMPAX/ADCS/B-field/reference/2025-06-28-10-30-MainEst+RawMTMs/2025-06-28-10-30-MainEst.csv", ',', header=true)
    mag_data, mag_head = DelimitedFiles.readdlm("/Users/thanasi/Documents/IMPAX/ADCS/B-field/reference/2025-06-28-10-30-MainEst+RawMTMs/2025-06-28-10-30-RawMTMs.csv", ',', header=true)
    pos_data, pos_head = DelimitedFiles.readdlm("/Users/thanasi/Documents/IMPAX/ADCS/B-field/reference/obc_get_GNSS_BESTXYZ_Data_2025-06-27_12_09_15(in).csv", ',', header=true)

    att_df = DataFrames.DataFrame(att_data, vec(att_head))
    mag_df = DataFrames.DataFrame(mag_data, vec(mag_head))
    pos_df = DataFrames.DataFrame(pos_data, vec(pos_head))
    pos_df = sort!(pos_df, [:timestamp])

    att_n = length(att_df[:,"Remote Date/Time"])
    mag_n = length(mag_df[:,"Remote Date/Time"])
    pos_n = length(pos_df[:,"timestamp"])

    datetime_format = Dates.dateformat"mm/dd/yyyy HH:MM:SS p.s"
    att_time_raw = Vector{Dates.DateTime}(undef, att_n)
    mag_time_raw = Vector{Dates.DateTime}(undef, mag_n)
    pos_time_raw = Vector{Dates.DateTime}(undef, pos_n)
    
    for (k,t_att) in enumerate(att_df[:, "Remote Date/Time"])
        att_time_raw[k] = Dates.DateTime(t_att, datetime_format)
    end
    for (k,t_mag) in enumerate(mag_df[:, "Remote Date/Time"])
        mag_time_raw[k] = Dates.DateTime(t_mag, datetime_format)
    end
    for k in 1:pos_n
        pos_time_raw[k] = gpstime_to_utc(pos_df[k, "week"], pos_df[k, "seconds"]; warn=false)
    end

    # --- Position -------------------------------------------------------------------------------
    # now using the specific ECEF frame: ITRF
    # now using the specific ECI frame: TEME

    pos_itrf = Point3f[(pos_df[k, "px"], pos_df[k, "py"], pos_df[k, "pz"]) for k in 1:pos_n]
    vel_itrf = Point3f[(pos_df[k, "vx"], pos_df[k, "vy"], pos_df[k, "vz"]) for k in 1:pos_n]

    pos_interp_itrf = [orbinterp(t, pos_time_raw, pos_itrf) for t in att_time_raw]
    vel_interp_itrf = [orbinterp(t, pos_time_raw, vel_itrf) for t in att_time_raw]
    pos_interp_n = mag_n

    pos_interp_teme = Point3f[(0.0, 0.0, 0.0) for k in 1:pos_interp_n]
    vel_interp_teme = Point3f[(0.0, 0.0, 0.0) for k in 1:pos_interp_n]
    field_igrf_teme = Point3f[(0.0, 0.0, 0.0) for k in 1:pos_interp_n]

    for k in 1:pos_interp_n
        jd = SatelliteToolboxBase.date_to_jd(mag_time_raw[k])
        C_IF = SatelliteToolboxTransformations.r_ecef_to_eci(SatelliteToolboxTransformations.ITRF(), SatelliteToolboxTransformations.TEME(), jd, eops)
        pos_interp_teme[k] = C_IF*pos_interp_itrf[k]
        vel_interp_teme[k] = C_IF*vel_interp_itrf[k]

        decyear = DateFormats.yeardecimal(mag_time_raw[k])
        pos_lla = SatelliteToolboxTransformations.ecef_to_geodetic(pos_interp_itrf[k]) # NOTE: is this valid enough for an ITRF ECEF frame?
        field_ned = SatelliteToolboxGeomagneticField.igrf(decyear, pos_lla[3], pos_lla[1], pos_lla[2], Val(:geodetic), max_degree=13, show_warnings=false)
        # also scale field to µT instead of nT
        field_igrf_teme[k] = Vector(C_IF*SatelliteToolboxTransformations.ned_to_ecef(field_ned, pos_lla..., translate=false)) ./ 1e3
    end
    
    max_f = max(norm.(field_igrf_teme)...)
    field_norm_igrf_teme = [f ./ max_f for f in field_igrf_teme]

    pos_radius = [norm(p) for p in pos_interp_teme]

    # --- SGP4 propagator ----------------------------------------------------------------------------
    pos_jds = SatelliteToolboxBase.date_to_jd.(pos_time_raw)

    pos_teme = Point3f[(0.0,0.0,0.0) for k in 1:pos_n]
    vel_teme = Point3f[(0.0,0.0,0.0) for k in 1:pos_n]
    for k in 1:pos_n
        C_IF = SatelliteToolboxTransformations.r_ecef_to_eci(SatelliteToolboxTransformations.ITRF(), SatelliteToolboxTransformations.TEME(), pos_jds[k], eops)
        pos_teme[k] = C_IF*pos_itrf[k]
        vel_teme[k] = C_IF*vel_itrf[k]
    end
    # tle, P = SatelliteToolboxPropagators.Propagators.fit_mean_elements(Val(:SGP4), pos_jds, pos_teme, vel_teme)

    # --- Attitude -----------------------------------------------------------------------------------

    quat_raw = zeros(4, att_n)
    for k in 1:att_n
        quat_raw[:, k] = [att_df[k, "Estimated ORC quaternion Q0"]; att_df[k, "Estimated ORC quaternion Q1"]; att_df[k, "Estimated ORC quaternion Q2"]; att_df[k, "Estimated ORC quaternion Q3"]]
    end

    ypr_conv = zeros(3, att_n)
    ypr_raw =  zeros(3, att_n)
    for k in 1:att_n
        ypr_conv[:,k] = ang321(Matrix(r_from_quat(quat_raw[:, k])'))*180/pi
        ypr_raw[:,k] = [att_df[k, "Estimated yaw angle (deg)"]; att_df[k, "Estimated pitch angle (deg)"]; att_df[k, "Estimated roll angle (deg)"]]
    end

    # --- Construct ORC ------------------------------------------------------------------------------

    # inertial to ORC:
    C_OI = zeros(3,3,att_n)
    for k in 1:att_n
        x_OI = vel_interp_teme[k]
        x_OI = x_OI / norm(x_OI)
        z_OI = -pos_interp_teme[k]
        z_OI = z_OI / norm(z_OI)
        y_OI = LinearAlgebra.cross(z_OI, x_OI)
        C_OI[:,:,k] = [x_OI y_OI z_OI]'
    end

    # body to ORC
    C_OB = zeros(3,3,att_n)
    for k in 1:att_n
        C_OB[:,:,k] = r_from_quat(quat_raw[:,k])
    end

    # then TEME can be found from Body with <TEME vec> = C_OI'*C_OB*<body vec>


    # --- Magnetometer -------------------------------------------------------------------------------

    mag0_val_raw = zeros(3, mag_n)
    mag1_val_raw = zeros(3, mag_n)

    C_mag01 = [0 -1 0;
               -1 0 0;
               0 0 -1]

    C_BM = [1 0 0;
            0 1 0;
            0 0 1]
    for k in 1:mag_n
        mag0_val_raw[:, k] = [mag_df[k, "MAG0 raw vector X component (uT)"]; mag_df[k, "MAG0 raw vector Y component (uT)"]; mag_df[k, "MAG0 raw vector Z component (uT)"]]
        mag1_val_raw[:, k] = C_mag01*[mag_df[k, "MAG1 raw vector X component (uT)"]; mag_df[k, "MAG1 raw vector Y component (uT)"]; mag_df[k, "MAG1 raw vector Z component (uT)"]]
    end

    mean_mag_x_diff = Statistics.mean(mag0_val_raw[1,:] .- mag1_val_raw[1,:])
    mean_mag_y_diff = Statistics.mean(mag0_val_raw[2,:] .- mag1_val_raw[2,:])
    mean_mag_z_diff = Statistics.mean(mag0_val_raw[3,:] .- mag1_val_raw[3,:])
    println("mean MAG X difference: ", mean_mag_x_diff)
    println("mean MAG Y difference: ", mean_mag_y_diff)
    println("mean MAG Z difference: ", mean_mag_z_diff)

    mag1_corr = -[mean_mag_x_diff; mean_mag_y_diff; mean_mag_z_diff]

    do_correct_mag1 = false
    if do_correct_mag1
        for k in 1:mag_n
            mag1_val_raw[:,k] = mag1_val_raw[:,k] .- mag1_corr
        end
    end

    field_meas0_teme = Point3f[(0.0,0.0,0.0) for k in 1:mag_n]
    field_meas1_teme = Point3f[(0.0,0.0,0.0) for k in 1:mag_n]
    for k in 1:mag_n
        field_meas0_teme[k] = C_OI[:,:,k]'*C_OB[:,:,k]*C_BM*mag0_val_raw[:,k]
        field_meas1_teme[k] = C_OI[:,:,k]'*C_OB[:,:,k]*C_BM*mag1_val_raw[:,k]
    end

    # --- Wahba problem, static --------------------------------------------------------------------

    # convert the field value to a matrix and scale it to the magnetometer reading range:
    field_igrf_rect_teme = zeros(3,length(field_igrf_teme))
    for (k, p) in enumerate(field_igrf_teme)
        field_igrf_rect_teme[:,k] = [p[1]; p[2]; p[3]]
    end

    # determine mask for south pole GPS exclusion (don't include these measurements in the Wahba solution):
    dif = cat(0.0, diff(pos_radius), dims=1)

    skips = (dif .> 5e3) .|| (dif .< -5e3)
    pos_col = [s ? :red : :blue for s in skips]
    skiprange = [findfirst(skips), findlast(skips)]
    println("diffs:")
    println(skiprange[1])
    println(skiprange[2])

    weights = ones(mag_n)
    weights[skiprange[1]:skiprange[2]] .= 0

    C_BI = wahba(field_igrf_rect_teme, mag0_val_raw; weights=weights)
    # C_BI = wahba(mag0_val_raw, field_igrf_rect_teme)
    println("Computed C_BI from global Wahba solution:")
    display(C_BI)
    println("det, trace: ")
    println(det(C_BI))
    println(tr(C_BI'*C_BI))

    # --- Correct magnetometer -------------------------------------------------------------------

    do_rotate_mag_body = true
    if do_rotate_mag_body
        for k in 1:mag_n
            mag0_val_raw[:,k] = r_from_quat(quat_raw[:,k])*mag0_val_raw[:,k]
            mag1_val_raw[:,k] = r_from_quat(quat_raw[:,k])*mag1_val_raw[:,k]
        end
    end

    println(r_from_quat(quat_raw[:,1]))

    mag01_ang = Vector{Float64}(undef, mag_n)
    for k in 1:mag_n
        mag01_ang[k] = acos(dot(mag0_val_raw[:,k], mag1_val_raw[:,k])/norm(mag0_val_raw[:,k])/norm(mag1_val_raw[:,k]))*180/pi
    end

    wahba_static_igrf_body = zeros(3, mag_n)
    for k in 1:mag_n
        wahba_static_igrf_body[:, k] = C_BI*[field_igrf_teme[k][1]; field_igrf_teme[k][2]; field_igrf_teme[k][3]]
    end

    # --- Wahba problem, dynamic -----------------------------------------------------------------

    # approach: take a n step wide MAG0/MAG1 measurement and compare each with IGRF.
    C_dBI = zeros(3,3,mag_n)
    ypr_est = zeros(3, mag_n)
    wahba_dynamic_igrf_body = zeros(3, mag_n)
    
    window = 3
    
    k = 1
    while k <= mag_n - window
        if any(weights[k:k+window] .!= 1)
            k += 1
            continue
        end
        
        C_dBI[:,:,k] = wahba(field_igrf_rect_teme[:, k:k+window], mag0_val_raw[:, k:k+window]; weights=weights[k:k+window])

        wahba_dynamic_igrf_body[:,k] = C_dBI[:,:,k] * [field_igrf_teme[k][1]; field_igrf_teme[k][2]; field_igrf_teme[k][3]]
        ypr_est[:,k] = ang321(C_dBI[:,:,k])*180/pi
        k += 1
    end


    # --- Plotting -------------------------------------------------------------------------------
    # --------------- Magnetometer ---------------------------------------------------------------
    width_scale = 0.005
    length_scale = 0.0018
    # igrf_scale = 1e3
    r_E = 6371e3

    # GLMakie.activate!(title="OpenSOAP")
    CairoMakie.activate!()
    fig_mag_raw = Figure(size=(800, 800), backgroundcolor=:transparent) # create a new figure

    layout_mag = GridLayout(fig_mag_raw[1, 1])
    
    mag_x_ax = Axis(layout_mag[1,1], title="Raw Magnetometer", ylabel="MAG X [µT]")
    mag_y_ax = Axis(layout_mag[2,1], ylabel="MAG Y [µT]")
    mag_z_ax = Axis(layout_mag[3,1], ylabel="MAG Z [µT]")
    mag_scale_ax = Axis(layout_mag[4,1], ylabel="Field strength [µT]")

    datelabeler(ds) = begin
        res = Vector{String}(undef, length(ds))
        for (k,d) in enumerate(ds)
            res[k] = String(Dates.hour(d)) * ":" * String(Dates.minute(d)) * ":" * String(Dates.second(d))
        end
        return res
    end

    mag_diff_ax = Axis(layout_mag[5,1], xlabel="Time", ylabel="MAG0 vs. MAG1 [deg]", xtickformat = "HH:MM:SS")
    linkxaxes!(mag_x_ax, mag_y_ax, mag_z_ax, mag_scale_ax, mag_diff_ax)
    hidexdecorations!(mag_x_ax)
    hidexdecorations!(mag_y_ax)
    hidexdecorations!(mag_z_ax)
    hidexdecorations!(mag_scale_ax)
    
    plot!(mag_diff_ax, mag_time_raw, mag01_ang, color=:black)
    mag0_leg = lines!(mag_x_ax, mag_time_raw, mag0_val_raw[1,:], linewidth=2)
    mag1_leg = lines!(mag_x_ax, mag_time_raw, mag1_val_raw[1,:], linewidth=2)
    igrf_leg = lines!(mag_x_ax, mag_time_raw, [field_igrf_teme[k][1] for k in 1:mag_n], linewidth=2)
    wahba_static_leg = lines!(mag_x_ax, mag_time_raw, wahba_static_igrf_body[1,:], linewidth=2)
    wahba_dynamic_leg = lines!(mag_x_ax, mag_time_raw, wahba_dynamic_igrf_body[1,:], linewidth=2)
    lines!(mag_y_ax, mag_time_raw, mag0_val_raw[2,:], linewidth=2)
    lines!(mag_y_ax, mag_time_raw, mag1_val_raw[2,:], linewidth=2)
    lines!(mag_y_ax, mag_time_raw, [field_igrf_teme[k][2] for k in 1:mag_n], linewidth=2)
    lines!(mag_y_ax, mag_time_raw, wahba_static_igrf_body[2,:], linewidth=2)
    lines!(mag_y_ax, mag_time_raw, wahba_dynamic_igrf_body[2,:], linewidth=2)
    lines!(mag_z_ax, mag_time_raw, mag0_val_raw[3,:], linewidth=2)
    lines!(mag_z_ax, mag_time_raw, mag1_val_raw[3,:], linewidth=2)
    lines!(mag_z_ax, mag_time_raw, [field_igrf_teme[k][3] for k in 1:mag_n], linewidth=2)
    lines!(mag_z_ax, mag_time_raw, wahba_static_igrf_body[3,:], linewidth=2)
    lines!(mag_z_ax, mag_time_raw, wahba_dynamic_igrf_body[3,:], linewidth=2)
    
    lines!(mag_scale_ax, mag_time_raw, [norm(mag0_val_raw[:,k]) for k in 1:mag_n], linewidth=2)
    lines!(mag_scale_ax, mag_time_raw, [norm(mag1_val_raw[:,k]) for k in 1:mag_n], linewidth=2)
    lines!(mag_scale_ax, mag_time_raw, [norm(field_igrf_teme[k]) for k in 1:mag_n], linewidth=2)
    lines!(mag_scale_ax, mag_time_raw, [norm(wahba_static_igrf_body[:,k]) for k in 1:mag_n], linewidth=2)
    
    # vspan!(mag_x_ax, mag_time_raw[skiprange[1]], mag_time_raw[skiprange[2]], color=:black, alpha=0.8)

    hlines!(mag_diff_ax, 2.8, linestyle=:dash, color=:red)

    Legend(layout_mag[1:3, 2], [mag0_leg, mag1_leg, igrf_leg, wahba_static_leg, wahba_dynamic_leg], ["MAG0", "MAG1", "IGRF", "Static\nWahba-aligned IGRF", "Dynamic\nWahba-aligned IGRF"])
    # plot!(mag_x_ax, mag_time_raw, mag0_val_raw[1,:], mag1_val_raw[1,:])
    CairoMakie.save("/Users/thanasi/Documents/IMPAX/ADCS/B-field/reference/results/mag_compare.pdf", fig_mag_raw)

    fig_mag_teme = Figure(size=(800, 800), backgroundcolor=:transparent) # create a new figure
    layout_mag = GridLayout(fig_mag_teme[1, 1])
    
    mag_x_ax = Axis(layout_mag[1,1], title="Magnetometer with onboard attitude solution", ylabel="X [µT]")
    mag_y_ax = Axis(layout_mag[2,1], ylabel=:"Y [µT]")
    mag_z_ax = Axis(layout_mag[3,1], ylabel=:"Z [µT]")
    mag_diff_ax = Axis(layout_mag[4,1], ylabel="Error [deg]")
    
    rowsize!(layout_mag, 4, Relative(1/2))

    mag0_ang_err = [acos(dot(field_meas0_teme[k], field_igrf_teme[k]) / norm(field_meas0_teme[k]) / norm(field_igrf_teme[k])) for k in 1:mag_n] .* 180/pi
    mag1_ang_err = [acos(dot(field_meas1_teme[k], field_igrf_teme[k]) / norm(field_meas1_teme[k]) / norm(field_igrf_teme[k])) for k in 1:mag_n] .* 180/pi

    igrf_lab = lines!(mag_x_ax, mag_time_raw, [field_igrf_teme[k][1] for k in 1:mag_n], linewidth=2)
    mag0_lab = lines!(mag_x_ax, mag_time_raw, [field_meas0_teme[k][1] for k in 1:mag_n], linewidth=2)
    mag1_lab = lines!(mag_x_ax, mag_time_raw, [field_meas1_teme[k][1] for k in 1:mag_n], linewidth=2)
    lines!(mag_y_ax, mag_time_raw, [field_igrf_teme[k][2] for k in 1:mag_n], linewidth=2)
    lines!(mag_y_ax, mag_time_raw, [field_meas0_teme[k][2] for k in 1:mag_n], linewidth=2)
    lines!(mag_y_ax, mag_time_raw, [field_meas1_teme[k][2] for k in 1:mag_n], linewidth=2)
    lines!(mag_z_ax, mag_time_raw, [field_igrf_teme[k][3] for k in 1:mag_n], linewidth=2)
    lines!(mag_z_ax, mag_time_raw, [field_meas0_teme[k][3] for k in 1:mag_n], linewidth=2)
    lines!(mag_z_ax, mag_time_raw, [field_meas1_teme[k][3] for k in 1:mag_n], linewidth=2)

    mag0_diff_lab = lines!(mag_diff_ax, mag_time_raw, mag0_ang_err, linewidth = 2, color=:red)
    mag1_diff_lab = lines!(mag_diff_ax, mag_time_raw, mag1_ang_err, linewidth = 2, color=:blue)
    mag_req_lab = hlines!(mag_diff_ax, 2.8, linestyle=:dash, color=:black)
    
    Legend(layout_mag[1:3,2], [igrf_lab, mag0_lab, mag1_lab], ["IGRF", "MAG0", "MAG1"])
    Legend(layout_mag[4,2], [mag0_diff_lab, mag1_diff_lab, mag_req_lab], ["MAG0 to IGRF\nangular error", "MAG1 to IGRF\nangular error", "Requirement"])

    CairoMakie.save("/Users/thanasi/Documents/IMPAX/ADCS/B-field/reference/results/mag_proj_compare.pdf", fig_mag_teme)

    # --------------- Quaternion ---------------------------------------------------------------

    fig_att = Figure(size=(800, 800)) # create a new figure

    layout_att = GridLayout(fig_att[1, 1])
    
    att_q3_ax = Axis(layout_att[1,1], title="Raw quaternion", ylabel=L"q_3")
    att_q012_ax = Axis(layout_att[2,1], ylabel=L"\mathbf{q}")
    plt_q0 = plot!(att_q012_ax, att_time_raw, quat_raw[1,:])
    plt_q1 = plot!(att_q012_ax, att_time_raw, quat_raw[2,:])
    plt_q2 = plot!(att_q012_ax, att_time_raw, quat_raw[3,:])
    plt_q3 = plot!(att_q3_ax, att_time_raw, quat_raw[4,:])
    Legend(layout_att[1:2, 2], [plt_q0, plt_q1, plt_q2, plt_q3], ["q0", "q1", "q2", "q3"])

    att_ypr_ax = Axis(layout_att[3:4,1], title="Yaw/Pitch/Roll", ylabel="YPR [deg]")
    plt_y_raw = lines!(att_ypr_ax, att_time_raw, ypr_raw[1,:], linewidth=1, color=:red, alpha=0.5)
    plt_p_raw = lines!(att_ypr_ax, att_time_raw, ypr_raw[2,:], linewidth=1, color=:green, alpha=0.5)
    plt_r_raw = lines!(att_ypr_ax, att_time_raw, ypr_raw[3,:], linewidth=1, color=:blue, alpha=0.5)
    plt_y_conv = lines!(att_ypr_ax, att_time_raw, ypr_conv[1,:], linewidth=1, linestyle=:dot, color=:red, alpha=0.5)
    plt_p_conv = lines!(att_ypr_ax, att_time_raw, ypr_conv[2,:], linewidth=1, linestyle=:dot, color=:green, alpha=0.5)
    plt_r_conv = lines!(att_ypr_ax, att_time_raw, ypr_conv[3,:], linewidth=1, linestyle=:dot, color=:blue, alpha=0.5)
    plt_y_est = lines!(att_ypr_ax, att_time_raw, ypr_est[1,:], linewidth=1, linestyle=:dash, color=:red, alpha=0.5)
    plt_p_est = lines!(att_ypr_ax, att_time_raw, ypr_est[2,:], linewidth=1, linestyle=:dash, color=:green, alpha=0.5)
    plt_r_est = lines!(att_ypr_ax, att_time_raw, ypr_est[3,:], linewidth=1, linestyle=:dash, color=:blue, alpha=0.5)

    Legend(layout_att[3:4, 2], [plt_y_raw, plt_p_raw, plt_r_raw, plt_y_conv, plt_p_conv, plt_r_conv, plt_y_est, plt_p_est, plt_r_est], ["raw yaw", "raw pitch", "raw roll", "calc yaw", "calc pitch", "calc roll", "est yaw", "est pitch", "est roll"])
    CairoMakie.save("/Users/thanasi/Documents/IMPAX/ADCS/B-field/reference/results/att_compare.pdf", fig_att)

    # --------------- Orbit --------------------------------------------------------------------

    GLMakie.activate!(title="OpenSOAP")
    fig_orbit = Figure(size=(800, 800)) # create a new figure
    layout_orbit = GridLayout(fig_orbit[1, 1])
    display(fig_orbit)
    
    texture = load_earth_texture_to_ecef(joinpath("assets","map_diffuse.png"))
    al = AmbientLight(RGBf(0.7,0.7,0.7))
    ax_orbit = LScene(layout_orbit[1:2,1], show_axis=false, scenekw=(lights=[al], clear=true))

    plot_earth_static!(ax_orbit, SatelliteToolboxBase.date_to_jd(pos_time_raw[1])*24*3600, eops, texture)

    # ax_orbit_radius = Axis(layout_orbit[3,1], xlabel="Time", ylabel="Orbital radius difference [m]")

    ax_dcm = LScene(layout_orbit[3,1], show_axis=false, scenekw=(lights=[al], clear=true, aspect=:equal))

    # xlabel!(ax_orbit_radius, "Time")
    # ylabel!(ax_orbit_radius, "Orbit radius change")
    scatter!(ax_orbit, pos_interp_teme, markersize=4, color=:blue)
    lines!(ax_orbit, pos_interp_teme)

    arrows!(ax_orbit, pos_interp_teme, field_igrf_teme, linewidth=r_E*width_scale, lengthscale=r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), color=:blue, alpha=0.5)
    arrows!(ax_orbit, pos_interp_teme, field_meas0_teme, linewidth=r_E*width_scale, lengthscale=r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), color=:magenta, alpha=0.5)
    arrows!(ax_orbit, pos_interp_teme, field_meas1_teme, linewidth=r_E*width_scale, lengthscale=r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), color=:purple, alpha=0.5)

    # arrows!(ax_orbit, pos_interp_teme, Point3f[mag1_val_raw[:,k] for k in 1:att_n], linewidth=r_E*width_scale, lengthscale=r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), color=:red, alpha=0.25)
    # arrows!(ax_orbit, pos_interp_teme, Point3f[mag0_val_raw[:,k] for k in 1:att_n], linewidth=r_E*width_scale, lengthscale=r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), color=:magenta, alpha=0.25)
    # arrows!(ax_orbit, pos_interp_teme, Point3f[C_dBI[:,:,k]'*mag1_val_raw[:,k] for k in 1:att_n], linewidth=r_E*width_scale, lengthscale=r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), color=:red, alpha=0.25)
    # arrows!(ax_orbit, pos_interp_teme, Point3f[C_dBI[:,:,k]'*mag0_val_raw[:,k] for k in 1:att_n], linewidth=r_E*width_scale, lengthscale=r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), color=:magenta, alpha=0.25)
    
    # arrows!(ax_orbit, pos_interp_teme, Point3f[C_dBI[1,:,k] for k in 1:att_n], linewidth=r_E*width_scale, lengthscale=1e2*r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), alpha=1, color=:red)
    # arrows!(ax_orbit, pos_interp_teme, Point3f[C_dBI[2,:,k] for k in 1:att_n], linewidth=r_E*width_scale, lengthscale=1e2*r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), alpha=1, color=:green)
    # arrows!(ax_orbit, pos_interp_teme, Point3f[C_dBI[3,:,k] for k in 1:att_n], linewidth=r_E*width_scale, lengthscale=1e2*r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), alpha=1, color=:blue)
    
    # scatter!(ax_orbit_radius, mag_time_raw, pos_radius .- Statistics.mean(pos_radius), color=pos_col)
    
    # plotting velocity:
    # arrows!(ax_orbit, pos_interp_teme, vel_interp_teme, linewidth=r_E*width_scale, lengthscale=0.1*r_E*length_scale, arrowsize=Vec3f(1.5*r_E*width_scale, 1.5*r_E*width_scale, 2*r_E*width_scale), alpha=1, color=:orange)

    # scatter!(ax_orbit_radius, pos_interp_time, )

    # plotting DCM (inertial to body, via orbit frame):
    lines!(ax_dcm, Point3f[(C_OI[:,:,k]'*C_OB[:,:,k])[:, 1] for k in 1:mag_n], color=:red)
    lines!(ax_dcm, Point3f[(C_OI[:,:,k]'*C_OB[:,:,k])[:, 2] for k in 1:mag_n], color=:green)
    lines!(ax_dcm, Point3f[(C_OI[:,:,k]'*C_OB[:,:,k])[:, 3] for k in 1:mag_n], color=:blue)

    lines!(ax_dcm, Point3f[zeros(3), C_OI[:,:,1]'*C_OB[:,:,1][:, 1]], linestyle=:dash, color=:red)
    lines!(ax_dcm, Point3f[zeros(3), C_OI[:,:,end]'*C_OB[:,:,end][:, 1]], linestyle=:solid, color=:red)
    lines!(ax_dcm, Point3f[zeros(3), C_OI[:,:,1]'*C_OB[:,:,1][:, 2]], linestyle=:dash, color=:green)
    lines!(ax_dcm, Point3f[zeros(3), C_OI[:,:,end]'*C_OB[:,:,end][:, 2]], linestyle=:solid, color=:green)
    lines!(ax_dcm, Point3f[zeros(3), C_OI[:,:,1]'*C_OB[:,:,1][:, 3]], linestyle=:dash, color=:blue)
    lines!(ax_dcm, Point3f[zeros(3), C_OI[:,:,end]'*C_OB[:,:,end][:, 3]], linestyle=:solid, color=:blue)

    return (att_df, pos_df)
end
