using OpenSOAP
using LinearAlgebra
import Makie, GLMakie, CairoMakie
# import Makie: date_to_number, convert_single_argument, wong_colors, categorical_colors
# import Makie: DimConversions, NoDimConversion, DateTimeConversion

using GeometryBasics
import SatelliteToolboxBase, SatelliteToolboxTransformations

import CSV, DataFrames, Dates


import Statistics

import Base.setindex!
using Printf

# mode_colors = Dict(
#     safe::Modes => [110, 110, 110] ./ 255,
#     idle::Modes => [42, 133, 255] ./ 255,
#     pointing::Modes => [42, 42, 42] ./ 255,
#     charging::Modes => [225, 220, 90] ./ 255,
#     science::Modes => [206, 155, 255] ./ 255,
#     downlink::Modes => [157, 226, 107] ./ 255,
# )

Makie.convert_single_argument(dates::AbstractVector{Dates.AbstractTime}) = Makie.date_to_number.(dates)
Makie.convert_single_argument(date::Dates.DateTime) = Makie.date_to_number(date)

function Base.setindex!(conversions::Makie.DimConversions, value, i::Int)
  isnothing(value) && return nothing # ignore no conversions
  conversions[i] === value && return nothing # ignore same conversion
  # only set new conversion if there is none yet
  # if isnothing(conversions[i]) || conversions[i] === NoDimConversion()
    conversions.conversions[i][] = value
    return nothing
  # end
end

function read_in(path::AbstractString)
    return CSV.read(path, DataFrames.DataFrame)
end

function random_walk_time(times::Vector{T}, drift_rate::T) where T<:Real
    res = zeros(length(times))
    for k in eachindex(times)
        if k == 1
            res[k] = times[k]
            continue 
        end
        dt = times[k] - times[k - 1]
        # draw random walk from a uniform distribution, centered on zero, spanning
        res[k] = times[k - 1] + dt*drift_rate*(2*rand() - 1)
    end
    return res
end

function corrupt_pos_forward(k::Int, scale::Real, rhist::Matrix{<:Real}, vhist::Matrix{<:Real})
    nt = length(rhist[1,:])
    rthis = rhist[:,k]
    vthis = vhist[:,k]
    orbitsteps = 2840
    
    @warn "Assuming "*string(orbitsteps)*" steps per orbit!\n" maxlog=1
    
    kstop = min(nt, k+(orbitsteps ÷ 2))
    rpush = rthis + scale*vthis/norm(vthis)
    
    bestnorm = Inf
    bestk = k
    for kf in k:kstop
        newnorm = norm(rhist[:,kf] - rpush)
        if newnorm < bestnorm
            bestk = kf
            bestnorm = newnorm
        end
    end
    
    return (rhist[:,bestk], vhist[:,bestk])
end

# use sim outputs from D2S2/Riccardo/Endurosat
function from_d2s2(readpath::AbstractString, savepath::AbstractString)
    df = read_in(readpath)
    eops = SatelliteToolboxTransformations.fetch_iers_eop()
    earth_texture_source = joinpath("assets", "map_diffuse.png")
    earth_texture = load_earth_texture_to_ecef(earth_texture_source)
    mag = MagneticTarget(
        "B-field",
        [Pair(35*pi/180,90*pi/180), Pair(-90*pi/180,-35*pi/180)],
        eops
    )
    
    # number of timesteps
    nt = DataFrames.nrow(df)
    
    # add a UTC DateTime column
    utcs = zeros(Dates.DateTime, nt)
    jds = zeros(Float64, nt)
    for k in 1:nt
        utcs[k] = Dates.DateTime(df[k,:Date], Dates.dateformat"d/m/yyyy") + Dates.Millisecond(round(1000*df[k, :met]))
        jds[k] = SatelliteToolboxBase.date_to_jd(utcs[k])
    end
    df[!, "utc"] = utcs
    df[!, "jd"] = jds
    duration = Dates.value(Dates.Second(df[end, :utc] - df[1, :utc]))
    
    # position error scale
    perr_rate = 3e3/24/3600 # 3 km per day, in units of m/s
    perr_inflate = 1000
    
    # compute position estimation error
    perr = zeros(Float64, nt)
    r_I = zeros(Float64, 3, nt)
    v_I = zeros(Float64, 3, nt)
    rest_I = zeros(Float64, 3, nt)
    vest_I = zeros(Float64, 3, nt)
    perr_extra = zeros(Float64, nt)
    con_sun_track = zeros(UInt8, nt)
    con_bfield_track = zeros(UInt8, nt)
    con_xyz_wheel = zeros(UInt8, nt)
    lat_mask = zeros(UInt8, nt)
    utcs_lat_mask_lo = Dates.DateTime[]
    utcs_lat_mask_hi = Dates.DateTime[]
    color_hist = [GLMakie.RGBA(0,0,0) for _ in 1:nt]
    for k in 1:nt
        perr[k] = norm([df[k, :true_pos_eci_x] - df[k, :estimated_pos_eci_x],
                        df[k, :true_pos_eci_y] - df[k, :estimated_pos_eci_y],
                        df[k, :true_pos_eci_z] - df[k, :estimated_pos_eci_z]])
        r_I[:,k] = [df[k, :true_pos_eci_x], df[k, :true_pos_eci_y], df[k, :true_pos_eci_z]]
        v_I[:,k] = [df[k, :true_vel_eci_x], df[k, :true_vel_eci_y], df[k, :true_vel_eci_z]]
        
        lat_mask[k] = abs(position_lla(r_I[:,k], jds[k]*3600*24, eops)[1]) > 35*pi/180
        if k > 1
            if lat_mask[k] > lat_mask[k-1]
                push!(utcs_lat_mask_lo, utcs[k])
            elseif lat_mask[k] < lat_mask[k-1]
                push!(utcs_lat_mask_hi, utcs[k])
            end
        end
        if df[k, :control_mode] == "ConSunTrack"
            con_sun_track[k] = 1
            color_hist[k] = GLMakie.RGBA(([225, 220, 90] ./ 255)...)
        elseif df[k, :control_mode] == "ConBfieldTrack"
            con_bfield_track[k] = 1
            color_hist[k] = GLMakie.RGBA(([206, 155, 255] ./ 255)...)
        elseif df[k, :control_mode] == "ConXYZwheel"
            con_xyz_wheel[k] = 1
            color_hist[k] = GLMakie.RGBA(([42, 133, 255] ./ 255)...)
        else
            println("unknown mode")
        end
    end
    # correct stop and start arrays so they are the same length
    if length(utcs_lat_mask_lo) > length(utcs_lat_mask_hi)
        push!(utcs_lat_mask_hi, utcs[end])
    elseif length(utcs_lat_mask_lo) < length(utcs_lat_mask_hi)
        pushfirst!(utcs_lat_mask_lo, utcs[1])
    end
    
    # corrupt position and velocity data
    for k in 1:nt
        # rv = corrupt_pos_forward(k, perr_rate*df[k, :met], r_I, v_I)
        # rest_I[:,k] = rv[1]
        # vest_I[:,k] = rv[2]
        
        # corrupt the true position by assuming TLE errors scale only the along-track component linearly, based on mission elapsed time (TLE epoch)
        ovec_I = normalize(cross(r_I[:,k], v_I[:,k]))
        # rest_I[:,k] = r_I[:,k] + perr_rate * df[k, :met] * v_I[:,k] / norm(v_I[:,k])
        rest_I[:,k] = r_I[:,k] + 9e3 * v_I[:,k] / norm(v_I[:,k])
        vest_I[:,k] = cross(ovec_I, rest_I[:,k])*norm(v_I[:,k]) / norm(rest_I[:,k])
        perr_extra[k] = norm([df[k, :true_pos_eci_x] - rest_I[1,k],
                        df[k, :true_pos_eci_y] - rest_I[2,k],
                        df[k, :true_pos_eci_z] - rest_I[3,k]])
    end
    
    # compare D2S2 B-field in body frame to SatelliteToolbox
    # Use df[!, true_mag_body_x] and the given ORC to Body quaternion to get values back in ORC. 
    #   Then construct ORC to ECI from local velocity and local geodetic nadir.
    #   Then compare directly between SatelliteToolbox outputs in ECI (queried by true D2S2 position) to transformed D2S2 B-field.
    
    bfield_body = zeros(Float64, 3, nt)
    
    for k in 1:nt
        bfield_body[:,k] = [df[k, :true_mag_body_x],df[k, :true_mag_body_y],df[k, :true_mag_body_z]]
        # using estimated quaternion here for ORC to Body transform because true values unreported from D2S2.
        this_q_BO = [df[k, :estimated_orbit2body_quat_q1], df[k, :estimated_orbit2body_quat_q2], df[k, :estimated_orbit2body_quat_q3], df[k, :estimated_orbit2body_quat_q4]]
        
        # is it best to use true values here in construction of ORC to ECI transform?
        v_eci = [df[k, :true_vel_eci_x], df[k, :true_vel_eci_y], df[k, :true_vel_eci_z]]
        r_eci = [df[k, :true_pos_eci_x], df[k, :true_pos_eci_y], df[k, :true_pos_eci_z]]
    end
    
    ax_q_est = zeros(Float64, 3, nt)
    ang_q_est = zeros(Float64, nt)
    ax_q_tru = zeros(Float64, 3, nt)
    ang_q_tru = zeros(Float64, nt)
    ax_q_err = zeros(Float64, nt)
    ang_q_err = zeros(Float64, nt)
    ax_e_err = zeros(Float64, nt)
    ang_e_err = zeros(Float64, nt)
    ax_qe_err = zeros(Float64, nt)
    ang_qe_err = zeros(Float64, nt)
    
    C_IO = zeros(Float64, 3, 3, nt)
    C_q_BO = zeros(Float64, 3, 3, nt)
    C_e_BO = zeros(Float64, 3, 3, nt)
    Cest_IO = zeros(Float64, 3, 3, nt)
    Cest_q_BO = zeros(Float64, 3, 3, nt)
    Cest_e_BO = zeros(Float64, 3, 3, nt)
    
    C_BO_est_resid = zeros(Float64, nt)
    for k in 1:nt
        # estimated Body-to-ORC coords from D2S2
        Cest_q_BO[:,:,k] = r_from_quat([df[k, :estimated_orbit2body_quat_q1], df[k, :estimated_orbit2body_quat_q2], df[k, :estimated_orbit2body_quat_q3], df[k, :estimated_orbit2body_quat_q4]])'
        aa_est = axisangle(Cest_q_BO[:,:,k])
        ax_q_est[:, k] = aa_est[1]
        ang_q_est[k] = aa_est[2]
        
        # true Body-to-ORC coords from D2S2
        C_q_BO[:,:,k] = r_from_quat([df[k, :true_orbit2body_quat_q1], df[k, :true_orbit2body_quat_q2], df[k, :true_orbit2body_quat_q3], df[k, :true_orbit2body_quat_q4]])'
        aa_tru = axisangle(C_q_BO[:,:,k])
        ax_q_tru[:, k] = aa_tru[1]
        ang_q_tru[k] = aa_tru[2]
        
        # Body-to-ORC estimate vs. true errors
        ax_q_err[k] = acos(ax_q_tru[:,k]'*ax_q_est[:, k] / norm(ax_q_tru[:,k]) / norm(ax_q_est[:,k]))
        ang_q_err[k] = ang_q_tru[k] .- ang_q_est[k]
        
        C_BO_est_resid[k] = residualSO3(Cest_q_BO[:,:,k]'*C_q_BO[:,:,k])
        
        # same, but with D2S2 Euler angles
        C_e_BO[:,:,k] = r_euler3(pi/180*df[k, :true_attitude_yaw])'*r_euler1(pi/180*df[k, :true_attitude_roll])'*r_euler2(pi/180*df[k, :true_attitude_pitch])'
        # C_e_BO[:,:,k] = (r_euler3(pi/180*df[k, :true_attitude_yaw])*r_euler1(pi/180*df[k, :true_attitude_roll])*r_euler2(pi/180*df[k, :true_attitude_pitch]))'
        
        Cest_e_BO[:,:,k] = r_euler3(pi/180*df[k, :estimated_attitude_yaw])'*r_euler1(pi/180*df[k, :estimated_attitude_roll])'*r_euler2(pi/180*df[k, :estimated_attitude_pitch])'
        # C_e_BO[:,:,k] = (r_euler2(pi/180*df[k, :true_attitude_pitch])'*r_euler1(pi/180*df[k, :true_attitude_roll])'*r_euler3(pi/180*df[k, :true_attitude_yaw])')'
        # Cest_e_BO[:,:,k] = r_euler3(pi/180*df[k, :estimated_attitude_yaw])*r_euler1(pi/180*df[k, :estimated_attitude_roll])*r_euler2(pi/180*df[k, :estimated_attitude_pitch])
        ax_tru, an_tru = axisangle(C_e_BO[:,:,k])
        ax_est, an_est = axisangle(Cest_e_BO[:,:,k])
        # println(ax_tru, ", ", an_tru)
        
        # error in Euler angle-derived estimate vs. true
        ax_e_err[k] = acos(ax_tru'*ax_est / norm(ax_tru) / norm(ax_est))
        ang_e_err[k] = an_tru - an_est
        
        # ORC-to-ECI construction
        #   note: projecting this to SO3 due to non-perpendicularity of r_I and v_I for imperfectly circular orbits.
        C_IO_X = v_I[:,k]/norm(v_I[:,k])
        # C_IO_Z = -r_I[:,k]/norm(r_I[:,k])
        # C_IO_Y = normalize(cross(C_IO_Z, C_IO_X))
        C_IO_Y = normalize(cross(-r_I[:,k] / norm(r_I[:,k]), v_I[:,k] / norm(v_I[:,k])))
        C_IO_Z = normalize(cross(C_IO_X, C_IO_Y))
        C_IO[:,:,k] = projSO3([C_IO_X C_IO_Y C_IO_Z])
        
        # ORC-to-ECI construction from estimates
        Cest_IO_X = vest_I[:,k]/norm(vest_I[:,k])
        # Cest_IO_Z = -rest_I[:,k]/norm(rest_I[:,k])
        # Cest_IO_Y = cross(Cest_IO_Z, Cest_IO_X)
        Cest_IO_Y = normalize(cross(-rest_I[:,k] / norm(rest_I[:,k]), vest_I[:,k] / norm(vest_I[:,k])))
        Cest_IO_Z = normalize(cross(Cest_IO_X, Cest_IO_Y))
        Cest_IO[:,:,k] = projSO3([Cest_IO_X Cest_IO_Y Cest_IO_Z])
        
        
        # Let's assume that the Euler angles actually represent intertial-to-Body transform.
        C_e_BI = C_e_BO[:,:,k]
        # C_q_BI = C_q_BO[:,:,k]*C_IO[:,:,k]'
        C_q_BI = C_IO[:,:,k]*C_q_BO[:,:,k]'
        ax_e_BI, ang_e_BI = axisangle(C_e_BI)
        ax_q_BI, ang_q_BI = axisangle(C_q_BI)
        
        # difference in Euler angle vs. quaternion truth
        ax_qe_err[k] = acos(min(1.0, ax_tru'*ax_q_tru[:,k] / norm(ax_tru) / norm(ax_q_tru[:,k])))
        ang_qe_err[k] = an_tru - ang_q_tru[k]
        # ax_qe_err[k] = acos(ax_e_BI'*ax_q_BI / norm(ax_e_BI) / norm(ax_q_BI))
        # ang_qe_err[k] = ang_e_BI - ang_q_BI
    end
    
    # Bfield
    bfield_tru_tp_I = zeros(Float64, 3, nt)
    bfield_tru_rc_I = zeros(Float64, 3, nt)
    bfield_est_rc_I = zeros(Float64, 3, nt)
    bfield_est_rc_notle_I = zeros(Float64, 3, nt)
    bfield_tru_rc_B = zeros(Float64, 3, nt)
    bfield_obs_I = zeros(Float64, 3, nt)
    
    bfield_err_tp_rc = zeros(Float64, nt)
    bfield_err_est_rc = zeros(Float64, nt)
    bfield_err_est_rc = zeros(Float64, nt)
    bfield_err_est_rc_notle = zeros(Float64, nt)
    bfield_err_est_tp_rc = zeros(Float64, nt)
    bfield_err_obs_tp_rc = zeros(Float64, nt)
    bfield_err_obs_rc = zeros(Float64, nt)
    bfield_err_pnt_rc = zeros(Float64, nt)
    bfield_err_agg_rc = zeros(Float64, nt)
    
    b_ref = [0.0;0.0;1.0]
    for k in 1:nt
        bfield_tru_tp_I[:,k] = position_eci(mag, r_I[:,k], jds[k]*3600*24; max_degree=13)
        
        bfield_tru_rc_B[:,k] = [df[k,:true_mag_body_x], df[k,:true_mag_body_y], df[k,:true_mag_body_z]]
        bfield_tru_rc_I[:,k] = C_IO[:,:,k]*C_q_BO[:,:,k]'*bfield_tru_rc_B[:,k]
        bfield_est_rc_I[:,k] = Cest_IO[:,:,k]*Cest_q_BO[:,:,k]'*bfield_tru_rc_B[:,k]
        bfield_est_rc_notle_I[:,k] = C_IO[:,:,k]*Cest_q_BO[:,:,k]'*bfield_tru_rc_B[:,k]
        bfield_obs_I[:,k] =    Cest_IO[:,:,k]*C_e_BO[:,:,k]'*bfield_tru_rc_B[:,k]
        
        # normalize everything:
        bfield_tru_tp_I[:,k] = bfield_tru_tp_I[:,k] / norm(bfield_tru_tp_I[:,k])
        bfield_tru_rc_I[:,k] = bfield_tru_rc_I[:,k] / norm(bfield_tru_rc_I[:,k])
        bfield_est_rc_I[:,k] = bfield_est_rc_I[:,k] / norm(bfield_est_rc_I[:,k])
        bfield_tru_rc_B[:,k] = bfield_tru_rc_B[:,k] / norm(bfield_tru_rc_B[:,k])
        bfield_est_rc_notle_I[:,k] = bfield_est_rc_notle_I[:,k] / norm(bfield_est_rc_notle_I[:,k])
        bfield_obs_I[:,k] = bfield_obs_I[:,k] / norm(bfield_obs_I[:,k])
        
        # compute errors
        #   note: all this mess with the dot product is just to provide a *signed* error measurement. 
        bfield_err_agg_rc[k]        = df[k, :magnetic_nadir_pointing_error]*pi/180
        bfield_err_tp_rc[k]         = asin(sign(dot(b_ref, cross(bfield_tru_tp_I[:,k], bfield_tru_rc_I[:,k])))*norm(cross(bfield_tru_tp_I[:,k], bfield_tru_rc_I[:,k])) / norm(bfield_tru_tp_I[:,k]) / norm(bfield_tru_rc_I[:,k]))
        bfield_err_est_rc[k]        = asin(sign(dot(b_ref, cross(bfield_tru_rc_I[:,k], bfield_est_rc_I[:,k])))*norm(cross(bfield_tru_rc_I[:,k], bfield_est_rc_I[:,k])) / norm(bfield_tru_rc_I[:,k]) / norm(bfield_est_rc_I[:,k]))
        bfield_err_est_tp_rc[k]     = asin(sign(dot(b_ref, cross(bfield_tru_tp_I[:,k], bfield_est_rc_I[:,k])))*norm(cross(bfield_tru_tp_I[:,k], bfield_est_rc_I[:,k])) / norm(bfield_tru_tp_I[:,k]) / norm(bfield_est_rc_I[:,k]))
        bfield_err_est_rc_notle[k]  = asin(sign(dot(b_ref, cross(bfield_tru_rc_I[:,k], bfield_est_rc_notle_I[:,k])))*norm(cross(bfield_tru_rc_I[:,k], bfield_est_rc_notle_I[:,k])) / norm(bfield_tru_rc_I[:,k]) / norm(bfield_est_rc_notle_I[:,k]))
        bfield_err_obs_tp_rc[k]     = asin(sign(dot(b_ref, cross(bfield_tru_tp_I[:,k], bfield_obs_I[:,k])))*norm(cross(bfield_tru_tp_I[:,k], bfield_obs_I[:,k])) / norm(bfield_tru_tp_I[:,k]) / norm(bfield_obs_I[:,k]))
        bfield_err_obs_rc[k]        = asin(sign(dot(b_ref, cross(bfield_tru_rc_I[:,k], bfield_obs_I[:,k])))*norm(cross(bfield_tru_rc_I[:,k], bfield_obs_I[:,k])) / norm(bfield_tru_rc_I[:,k]) / norm(bfield_obs_I[:,k]) * 0.99999999999)
        bfield_err_pnt_rc[k]        = asin(sign(dot(b_ref, cross(bfield_est_rc_I[:,k], bfield_obs_I[:,k])))*norm(cross(bfield_est_rc_I[:,k], bfield_obs_I[:,k])) / norm(bfield_est_rc_I[:,k]) / norm(bfield_obs_I[:,k]))
    end
    
    # binning errors
    bfield_err_tp_rc_mask = [bfield_err_tp_rc[k]            for k in 1:nt if 1 == lat_mask[k]]
    bfield_err_est_rc_mask = [bfield_err_est_rc[k]          for k in 1:nt if 1 == lat_mask[k]]
    bfield_err_est_tp_rc_mask = [bfield_err_est_tp_rc[k]    for k in 1:nt if 1 == lat_mask[k]]
    bfield_err_obs_tp_rc_mask = [bfield_err_obs_tp_rc[k]    for k in 1:nt if 1 == lat_mask[k]]
    bfield_err_obs_rc_mask = [bfield_err_obs_rc[k]    for k in 1:nt if 1 == lat_mask[k]]
    bfield_err_pnt_rc_mask = [bfield_err_pnt_rc[k]    for k in 1:nt if 1 == lat_mask[k]]
    bfield_err_est_rc_notle_mask = [bfield_err_est_rc_notle[k] for k in 1:nt if 1 == lat_mask[k]]
    bfield_err_agg_rc_mask = [bfield_err_agg_rc[k]     for k in 1:nt if 1 == lat_mask[k]]
    
    # plot
        
    legend_fs = 10.0f0 
    bfield_bg_col = GLMakie.RGBAf([206, 155, 255] ./ 255..., 0.25)
    hstrokewidth = 2
    
    # CairoMakie.activate!()
    # fig_density = CairoMakie.Figure(size=(600,600))
    # ax_density = CairoMakie.Axis(fig_density[1:2,1], xlabel="Knowledge error [deg]", ylabel="Density")
    
    MakieBackend = GLMakie
    
    MakieBackend.activate!()
    fig = nothing
    ax_globe, ax_perr, ax_qerr, ax_berr, ax_all_errh, ax_stack = begin
        if MakieBackend === GLMakie
            fig = MakieBackend.Figure(size=(1500,1000))
            left = fig[1,1] = GLMakie.GridLayout()
            center = fig[1,2] = GLMakie.GridLayout()
            right = fig[1,3] = GLMakie.GridLayout()
            al = GLMakie.AmbientLight(GLMakie.RGBf(0.7,0.7,0.7))
            lax_globe = GLMakie.LScene(left[1,1], show_axis = false, scenekw = (
                lights = [al],
                clear=true,
                shading=GLMakie.FastShading
            ))
            
            lax_perr = GLMakie.Axis(center[1,1], xlabel="UTC", ylabel="Position error [km]")
            lax_qerr = GLMakie.Axis(center[2,1], xlabel="UTC", ylabel="Attitude knowledge error [deg]")
            lax_berr = GLMakie.Axis(center[3,1], xlabel="UTC", ylabel="B-field error [deg]")
            # ax_mode = GLMakie.Axis(right[5,1], xlabel="UTC", ylabel="Mode")
            lax_all_errh = GLMakie.Axis(right[1,1], xlabel="Error [deg]", ylabel="Density")
            lax_stack = GLMakie.Axis(right[2,1], ylabel="Error [deg]", xgridvisible = false, ygridvisible = false)
            # lax_temp = GLMakie.Axis(right[2,1])
            MakieBackend.linkxaxes!(lax_perr, lax_qerr, lax_berr)
            
            plot_earth_static!(lax_globe, df[1,:jd]*24*3600, eops, earth_texture)
            
            GLMakie.lines!(lax_globe, r_I[1,:], r_I[2,:], r_I[3,:], color=color_hist, linewidth=4)
            GLMakie.lines!(lax_globe, rest_I[1,:], rest_I[2,:], rest_I[3,:], color=:red, linestyle=:dash, linewidth=1)
            
            lax_globe, lax_perr, lax_qerr, lax_berr, lax_all_errh, lax_stack
        elseif MakieBackend === CairoMakie
            fig_perr = CairoMakie.Figure(size=(600,600))
            fig_qerr = CairoMakie.Figure(size=(600,600))
            fig_berr = CairoMakie.Figure(size=(600,600))
            fig_all_errh = CairoMakie.Figure(size=(600,600))
            fig_stack = CairoMakie.Figure(size=(600,600))
            
            # take two rows per plot so we can put the legend below, not too big, later.
            lax_perr = CairoMakie.Axis(fig_perr[1:2,1], xlabel="UTC", ylabel="Position error [km]")
            lax_qerr = CairoMakie.Axis(fig_qerr[1:2,1], xlabel="UTC", ylabel="Attitude knowledge error [deg]")
            lax_berr = CairoMakie.Axis(fig_berr[1:2,1], xlabel="UTC", ylabel="B-field error [deg]")
            # ax_mode = CairoMakie.Axis(right[5,1], xlabel="UTC", ylabel="Mode")
            lax_all_errh = CairoMakie.Axis(fig_all_errh[1:2,1], xlabel="Error [deg]", ylabel="Density")
            lax_stack = CairoMakie.Axis(fig_stack[1:2,1], ylabel="Error [deg]", xgridvisible = false, ygridvisible = false)
            
            nothing, lax_perr, lax_qerr, lax_berr, lax_all_errh, lax_stack
        end
    end
    
    # ax_terr = GLMakie.Axis(right[1,1], xlabel="Index", ylabel="UTC")
    
    
    # binslide = GLMakie.Slider(right[2,1], range = 5:1:100, startvalue = 30, update_while_dragging=true, tellheight=false, tellwidth=false)
    # nbins = GLMakie.lift(binslide.value) do n
    #     println(n)
    #     return n
    # end
    
    
    # GLMakie.lines!(ax_terr, 1:nt, utcs, color=:black, linewidth=1)
    
    MakieBackend.lines!(ax_perr, utcs, perr ./ 1e3, color=:black, linewidth=1, label="D2S2 propagated TLE error")
    MakieBackend.lines!(ax_perr, utcs, perr_extra ./ 1e3, color=:black, linewidth=1, label="Crude scaled along-track error")
    MakieBackend.lines!(ax_perr, utcs, perr_inflate*perr ./ 1e3, color=:blue, linewidth=1, linestyle=:dash, label="Error scaled onto 3 km/day")
    MakieBackend.lines!(ax_perr, [utcs[1], utcs[end]], [0, 3*duration] ./ 3600 ./ 24, color=:red, linestyle=:dash, linewidth=1, label="3 km/day bound")
    # GLMakie.ylims!(ax_perr, [0, 0.005])
    
    
    # debug_depth = GLMakie.Slider(center[3,2], range = 0.0:0.001:0.01, startvalue = 0.001, update_while_dragging=true, horizontal=false, tellheight=false)
    # depth_val = GLMakie.lift(debug_depth.value) do val
    #     println(val)
    #     return val
    # end
    
    MakieBackend.lines!(ax_qerr, utcs, ax_q_err*180/pi, linewidth=1, label="Quat est vs. true, ax")
    MakieBackend.lines!(ax_qerr, utcs, ang_q_err*180/pi, linewidth=1, label="Quat est vs. true, ang")
    MakieBackend.lines!(ax_qerr, utcs, C_BO_est_resid, linewidth=1, color=:black, linestyle=:dash, label="SO(3) residual")
    # GLMakie.lines!(ax_qerr, utcs, ax_qe_err*180/pi, linewidth=1, label="Quat vs. Euler, ax")
    # GLMakie.lines!(ax_qerr, utcs, ang_qe_err*180/pi, linewidth=1, label="Quat vs. Euler, ang")
    MakieBackend.lines!(ax_qerr, utcs, ax_e_err*180/pi, linewidth=1, label="Euler est vs. true, ax")
    MakieBackend.lines!(ax_qerr, utcs, ang_e_err*180/pi, linewidth=1, label="Euler est vs. true, ang")
    MakieBackend.ylims!(ax_qerr, [-0.5, 0.5])
    MakieBackend.vspan!(ax_qerr, utcs_lat_mask_lo, utcs_lat_mask_hi, color=bfield_bg_col, depth_shift=0.001)
    
    MakieBackend.lines!(ax_berr, utcs, bfield_err_tp_rc*180/pi, linewidth=1, label="UMN vs. ES true field")
    # MakieBackend.lines!(ax_berr, utcs, bfield_err_est_rc*180/pi, linewidth=1, label="ES estimate vs. true field")
    # MakieBackend.lines!(ax_berr, utcs, bfield_err_est_rc_notle*180/pi, linewidth=1, label="ES estimate vs. true field, perfect TLE")
    MakieBackend.lines!(ax_berr, utcs, bfield_err_est_tp_rc*180/pi, linewidth=1, label="ES estimate vs. UMN true field")
    MakieBackend.lines!(ax_berr, utcs, bfield_err_obs_tp_rc*180/pi, linewidth=1, label="ES pointing vs. UMN true field")
    # MakieBackend.lines!(ax_berr, utcs, bfield_err_obs_rc*180/pi, linewidth=1, label="ES pointing vs. ES field")
    MakieBackend.lines!(ax_berr, utcs, bfield_err_pnt_rc*180/pi, linewidth=1, label="ES pointing vs. ES field estimate")
    
    # MakieBackend.lines!(ax_berr, utcs, bfield_err_agg_rc*180/pi, linewidth=1, label="ES reported error")
    # MakieBackend.ylims!(ax_berr, [-0.5, 0.5])
    MakieBackend.ylims!(ax_berr, [-0.5, 0.5])
    MakieBackend.vspan!(ax_berr, utcs_lat_mask_lo, utcs_lat_mask_hi, color=bfield_bg_col, depth_shift=0.001)
    
    # GLMakie.lines!(ax_mode, utcs, con_bfield_track, linewidth=1, linestyle=:dash, label="ConBfieldTrack")
    # GLMakie.lines!(ax_mode, utcs, con_sun_track,    linewidth=1, linestyle=:dash, label="ConSunTrack")
    # GLMakie.lines!(ax_mode, utcs, con_xyz_wheel,    linewidth=1, linestyle=:dash, label="ConXYZwheel")
    # GLMakie.lines!(ax_mode, utcs, lat_mask,         linewidth=1, linestyle=:dash, label=GLMakie.L"| \mathrm{Lat} | > 35^\circ")
    # GLMakie.vspan!(ax_mode, utcs_lat_mask_lo, utcs_lat_mask_hi, color=bfield_bg_col, depth_shift=0.001)
    # GLMakie.axislegend(ax_mode, position=:lt, labelsize=legend_fs)
    
    # GLMakie.linkxaxes!(ax_perr, ax_qerr, ax_berr, ax_mode)
    
    
    # MakieBackend.density!(
    #     ax_all_errh, 
    #     bfield_err_tp_rc_mask.*180/pi, 
    #     strokewidth=hstrokewidth, strokecolor=:black, alpha=0.66,
    #     label="UMN vs. ES true field, (μ,σ)=("*@sprintf("%.3f, %.3f", Statistics.mean(bfield_err_tp_rc_mask)*180/pi, Statistics.std(bfield_err_tp_rc_mask)*180/pi)*") deg"
    # )
    
    MakieBackend.density!(
        ax_all_errh, 
        bfield_err_est_rc_mask.*180/pi, 
        strokewidth=hstrokewidth, strokecolor=:black, alpha=0.66,
        label="ES estimate vs. ES true field, (μ,σ)=("*@sprintf("%.3f, %.3f", Statistics.mean(bfield_err_est_rc_mask)*180/pi, Statistics.std(bfield_err_est_rc_mask)*180/pi)*") deg"
    )
    MakieBackend.density!(
        ax_all_errh, 
        bfield_err_est_rc_notle_mask.*180/pi, 
        strokewidth=hstrokewidth, strokecolor=:black, alpha=0.66,
        label="ES estimate vs. ES true field, perfect TLE, (μ,σ)=("*@sprintf("%.3f, %.3f", Statistics.mean(bfield_err_est_rc_notle_mask)*180/pi, Statistics.std(bfield_err_est_rc_notle_mask)*180/pi)*") deg"
    )
    
    # MakieBackend.density!(
    #     ax_all_errh, 
    #     bfield_err_est_tp_rc_mask.*180/pi, 
    #     strokewidth=hstrokewidth, strokecolor=:black, alpha=0.66,
    #     label="ES estimate vs. UMN true field, (μ,σ)=("*@sprintf("%.3f, %.3f", Statistics.mean(bfield_err_est_tp_rc_mask)*180/pi, Statistics.std(bfield_err_est_tp_rc_mask)*180/pi)*") deg"
    # )
    # MakieBackend.density!(
    #     ax_all_errh, 
    #     bfield_err_obs_tp_rc_mask.*180/pi, 
    #     strokewidth=hstrokewidth, strokecolor=:black, alpha=0.66,
    #     label="ES pointing vs. UMN true field, (μ,σ)=("*@sprintf("%.3f, %.3f", Statistics.mean(bfield_err_obs_tp_rc_mask)*180/pi, Statistics.std(bfield_err_obs_tp_rc_mask)*180/pi)*") deg"
    # )
    ## MakieBackend.density!(
    ##     ax_all_errh, 
    ##     bfield_err_obs_rc_mask.*180/pi, 
    ##     strokewidth=hstrokewidth, strokecolor=:black, alpha=0.66,
    ##     label="ES pointing vs. ES true field, (μ,σ)=("*@sprintf("%.3f, %.3f", Statistics.mean(bfield_err_obs_rc_mask)*180/pi, Statistics.std(bfield_err_obs_rc_mask)*180/pi)*") deg"
    ## )
    # MakieBackend.density!(
    #     ax_all_errh, 
    #     bfield_err_pnt_rc_mask.*180/pi, 
    #     strokewidth=hstrokewidth, strokecolor=:black, alpha=0.66,
    #     label="ES pointing vs. ES estimate field, (μ,σ)=("*@sprintf("%.3f, %.3f", Statistics.mean(bfield_err_pnt_rc_mask)*180/pi, Statistics.std(bfield_err_pnt_rc_mask)*180/pi)*") deg"
    # )
    
    ## MakieBackend.density!(
    ##     ax_all_errh, 
    ##     bfield_err_agg_rc_mask.*180/pi, 
    ##     strokewidth=hstrokewidth, strokecolor=:black, alpha=0.66,
    ##     label="ES reported error, (μ,σ)=("*@sprintf("%.3f, %.3f", Statistics.mean(bfield_err_agg_rc_mask)*180/pi, Statistics.std(bfield_err_agg_rc_mask)*180/pi)*") deg"
    ## )
    
    # MakieBackend.vlines!(ax_all_errh, [Statistics.mean(bfield_err_tp_rc_mask)*180/pi])
    
    MakieBackend.vlines!(ax_all_errh, [Statistics.mean(bfield_err_est_rc_mask)*180/pi])
    MakieBackend.vlines!(ax_all_errh, [Statistics.mean(bfield_err_est_rc_notle_mask)*180/pi])
    
    # MakieBackend.vlines!(ax_all_errh, [Statistics.mean(bfield_err_est_tp_rc_mask)*180/pi])
    # MakieBackend.vlines!(ax_all_errh, [Statistics.mean(bfield_err_obs_tp_rc_mask)*180/pi])
    # # MakieBackend.vlines!(ax_all_errh, [Statistics.mean(bfield_err_obs_rc_mask)*180/pi])
    # MakieBackend.vlines!(ax_all_errh, [Statistics.mean(bfield_err_pnt_rc_mask)*180/pi])
    ## MakieBackend.vlines!(ax_all_errh, [Statistics.mean(bfield_err_agg_rc_mask)*180/pi])
    MakieBackend.xlims!(ax_all_errh, [-0.3, 0.3])
    
    err_stack = Dict(
        "pointing jitter bias"=>abs(Statistics.mean(bfield_err_pnt_rc_mask)*180/pi),
        "pointing jitter 3σ"=>3*Statistics.std(bfield_err_pnt_rc_mask)*180/pi,
        "knowledge error bias"=>abs(Statistics.mean(bfield_err_est_tp_rc_mask)*180/pi),
        "knowledge error 3σ"=>3*Statistics.std(bfield_err_est_tp_rc_mask)*180/pi,
        "IGRF model worst-case error"=>1.0,
        "FIRE worst-case internal boresight error"=>0.016,
        "AXIS worst-case internal boresight error"=>0.098,
        "FIRE worst-case mounting error"=>0.1,
        "AXIS worst-case mounting error"=>0.1,
        "FSS worst-case mounting error"=>0.117 # this is 7-sigma or something, per Pavel
    )
    names = collect(keys(err_stack))
    errs = [err_stack[k] for k in names]
    sI = sortperm(errs, rev=true)
    names = names[sI]
    errs = errs[sI]
    err_his = cumsum(errs)
    err_cent = err_his .- 0.5*errs
    
    err_colors = Makie.categorical_colors(:tab10, length(errs))
    MakieBackend.barplot!(ax_stack,
        zeros(length(names)), errs,
        stack=zeros(Int64, length(names)),
        width=0.5,
        color=1:length(errs),
        strokecolor=:black,
        strokewidth=1,
        colormap=:tab10,
        # colorrange=(1,length(errs)),
        label="nothing"
    )
    labelxstart = 0.4
    labelystart = 0.5
    labelyscale = 0.2
    for (k, name) in enumerate(names)
        MakieBackend.lines!(ax_stack,
            [0.21, labelxstart-0.05, labelxstart - 0.01], [err_cent[k], labelystart + k*labelyscale, labelystart + k*labelyscale], linewidth=1, color=err_colors[k],
        )
        MakieBackend.text!(ax_stack,
            labelxstart, labelystart + k*labelyscale, text=name, align=(:left, :center), space=:data, color=err_colors[k],
        )
    end
    MakieBackend.text!(ax_stack,
        labelxstart, 2.8, text="2.8º pointing requirement", align=(:left, :top), space=:data, color=:red
    )
    
    MakieBackend.hlines!(ax_stack, [2.8], color=:red, linewidth=2, label="2.8º pointing requirement")
    MakieBackend.xlims!(ax_stack, [-0.25, 1.05])
    MakieBackend.ylims!(ax_stack, [0.0, 3.0])
    MakieBackend.hidexdecorations!(ax_stack)
    
    @printf "TOTAL ERROR:\t%.4f\n" err_his[end]
    @printf "TOTAL MARGIN:\t%.4f\n" 2.8 - err_his[end]
    for err in keys(err_stack)
        @printf "%-40s : %7.4f\n" err err_stack[err]
    end
    # println(err_stack)
    
    if MakieBackend === GLMakie
        MakieBackend.axislegend(ax_perr, position=:lt, labelsize=legend_fs)
        MakieBackend.axislegend(ax_qerr, position=:lt, labelsize=legend_fs)
        MakieBackend.axislegend(ax_berr, position=:lt, labelsize=legend_fs)
        MakieBackend.axislegend(ax_all_errh, position=:rt, labelsize=legend_fs)
        display(ax_perr.parent)
    elseif MakieBackend === CairoMakie    
        # x limits to zoom on detail
        xlims = [Dates.DateTime(2026,03,21,0,0,0), Dates.DateTime(2026,03,21,4,0,0)]
        for (k,p) in enumerate([ax_perr, ax_qerr, ax_berr, ax_all_errh, ax_stack])
            leg = MakieBackend.Legend(p.parent[3,1], p, labelsize=legend_fs, tellwidth=false)
            # set x limits to show detail
            if p === ax_qerr || p === ax_berr
                CairoMakie.xlims!(p, xlims)
            end
            CairoMakie.save(joinpath(savepath, "fig"*string(k)*".pdf"), p.parent; backend=CairoMakie)
        end
    end
    
    
    return df
end