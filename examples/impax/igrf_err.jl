using OpenSOAP
import SatelliteToolboxBase, SatelliteToolboxTransformations
using LinearAlgebra
using GeometryBasics, GLMakie

function main()
    eops = SatelliteToolboxTransformations.fetch_iers_eop()
    m = MagneticTarget("mag", eops)

    t_jd_s_start = 3600*24*SatelliteToolboxBase.date_to_jd(1980, 1, 1, 12, 0, 0)
    t_jd_s_end = 3600*24*SatelliteToolboxBase.date_to_jd(2030, 1, 1, 12, 0, 0)
    
    lim = 20
    nt = 1000
    lats = LinRange(-pi/2+0.1, pi/2-0.1, lim)
    lons = LinRange(-pi+0.1, pi-0.1, 2*lim)
    alt = 500e3
    times = LinRange(t_jd_s_start, t_jd_s_end, nt)

    # make the field:
    B = Array{Float64, 4}(undef, length(lats), length(lons), length(times), 3)
    for (i,lat) in enumerate(lats)
        for (j,lon) in enumerate(lons)
            for (k,time) in enumerate(times)
                p_ecef = Vector(SatelliteToolboxTransformations.geodetic_to_ecef(lat, lon, alt))
                B[i,j,k,:] = position_ecef(m, p_ecef, time)
            end
        end
    end

    # equator:
    equator_ref = LinearAlgebra.normalize(B[lim, 1, 1, :])
    # pole
    pole_ref = LinearAlgebra.normalize(B[end, 1, 1, :])

    equator_err = Vector{Float64}(undef, length(times))
    pole_err = Vector{Float64}(undef, length(times))
    for k in 1:length(times)
        equator_val = LinearAlgebra.normalize(B[lim, 1, k, :])
        equator_err[k] = acos(min(equator_val'*equator_ref, 1))*180/pi

        pole_val = LinearAlgebra.normalize(B[end, 1, k, :])
        pole_err[k] = acos(min(pole_val'*pole_ref, 1))*180/pi
    end

    GLMakie.activate!(title="OpenSOAP")
    
    fig = Figure(size=(800,800))
    display(fig)
    ax = Axis(fig[1,1])

    plot!(ax, (Vector(times) .- t_jd_s_start)/3600/24, equator_err)
    plot!(ax, (Vector(times) .- t_jd_s_start)/3600/24, pole_err)

end