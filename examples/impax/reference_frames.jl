using OpenSOAP
using LinearAlgebra
using FileIO
using GLMakie, GeometryBasics
import SatelliteToolboxBase
using SatelliteToolboxTransformations
import SatelliteToolboxCelestialBodies

function plot_orbit!(ax::Makie.LScene, t_jd_s, iers_eops; do_ned::Bool=true)
    r_E = SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS
    r_P = SatelliteToolboxBase.EARTH_POLAR_RADIUS
    r_O = r_E + 3000e3

    C_IO = r_euler2(pi/4)*r_euler1(-pi/3)*r_euler3(pi/3)
    C_IF = Matrix(r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, iers_eops))

    arrow_length = 0.2 * r_E
    arrow_width = 0.01 * r_E

    nθ = 200
    r_X = zeros(nθ)
    r_Y = zeros(nθ)
    r_Z = zeros(nθ)
    color = [RGBAf(0,0,0,1) for i in 1:nθ]

    for i in 1:nθ
        p = C_IO*r_O*[
            cos(2*pi*i/nθ + pi),
            sin(2*pi*i/nθ + pi),
            0.0
        ]
        r_X[i] = p[1]
        r_Y[i] = p[2]
        r_Z[i] = p[3]
        color[i] = RGBAf(0.0,0.0,0.0,i/nθ)
    end

    term_I = [r_X[end],r_Y[end],r_Z[end]]

    eci_basis = diagm([1.0,1.0,1.0])
    body_axes = r_random()

    ned_N = C_IF[:,3]
    ned_D = -term_I / norm(term_I)
    ned_E = cross(ned_D, ned_N)
    ned_E = ned_E / norm(ned_E)
    ned_N = cross(ned_E, ned_D)
    ned = hcat(ned_N, ned_E, ned_D)

    axes_tails = [Point3f(term_I) for k in 1:3]
    ned_heads = [Vec3f(arrow_length*ned[:,k]) for k in 1:3]
    ned_label_pos = [1.5 * ned_heads[k] + axes_tails[k] for k in 1:3]

    colors = [:red, :green, :blue]

    lines!(
        ax,
        r_X,
        r_Y,
        r_Z,
        color=color
    )
    scatter!(
        ax,
        term_I[1],
        term_I[2],
        term_I[3],
        color=:black
    )

    if do_ned
        arrows!(
            ax,
            axes_tails,
            ned_heads,
            color=colors,
            linewidth=arrow_width,
            arrowsize=Vec3f([arrow_width, arrow_width, arrow_width*2])
        )
        text!(
            ax,
            ned_label_pos,
            text=[L"N_1", L"N_2", L"N_3"],
            align=(:center, :center),
            fontsize=20,
            color=:black
            # markerspace=:data
        )
    end
    return (r_X, r_Y, r_Z, color)
end

function plot_body!(ax::Makie.LScene, t_jd_s, iers_eops, object::String, position)
    r_E = SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS
    model = load(joinpath("assets","IMPAX_mech_clean.obj"))
    
    loc = [position[1][end], position[2][end], position[3][end]]

    scale = 0.5*r_E/1e2
    arrow_length = 0.2 * r_E
    arrow_width = 0.01 * r_E

    eci_basis = diagm([1.0,1.0,1.0])
    axes_tails = [Point3f(loc) for k in 1:3]
    body_axes = r_random()
    body_heads = [Vec3f(arrow_length*body_axes*eci_basis[:,k]) for k in 1:3]
    body_label_pos = [1.5 * body_heads[k] + axes_tails[k] for k in 1:3]
    colors = [:red, :green, :blue]

    m = mesh!(
        ax,
        model,
        color=:grey,
        alpha=0.5,
        space=:data,
    )

    q = Makie.Quaternion(body_axes)

    scale!(m, scale, scale, scale)
    translate!(m, loc...)
    GLMakie.rotate!(m, q)

    arrows!(
        ax,
        axes_tails,
        body_heads,
        color=colors,
        linewidth=arrow_width,
        arrowsize=Vec3f([arrow_width, arrow_width, arrow_width*2])
    )
    text!(
        ax,
        body_label_pos,
        text=[L"B_1", L"B_2", L"B_3"],
        align=(:center, :center),
        fontsize=20,
        color=:black
        # markerspace=:data
    )
end

function plot_globe!(ax::Makie.LScene, t_jd_s, iers_eops, texture_ecef)
    r_E = SatelliteToolboxBase.EARTH_EQUATORIAL_RADIUS
    r_P = SatelliteToolboxBase.EARTH_POLAR_RADIUS

    arrow_length = 1.25 * r_E
    arrow_width = 0.025 * r_E
    
    C_IF = Matrix(r_ecef_to_eci(ITRF(), J2000(), t_jd_s/24/3600, iers_eops))

    # making surface texture:
    # nθ = 3600
    nθ = length(texture_ecef[:,1])
    nφ = 100
    θ = range(0, stop=2π, length=nθ)
    φ = range(0, stop=π, length=nφ)

    local_r = r_E*cos(35*pi/180)
    local_z = r_P*sin(35*pi/180)
    
    mesh_X = zeros(nθ, nφ)
    mesh_Y = zeros(nθ, nφ)
    mesh_Z = zeros(nθ, nφ)
    for i in 1:nθ
        for j in 1:nφ
            # a point in the mesh
            p = C_IF*[r_E*cos(θ[i])*sin(φ[j]); r_E*sin(θ[i])*sin(φ[j]); r_P*cos(φ[j])]
            mesh_X[i,j] = p[1]
            mesh_Y[i,j] = p[2]
            mesh_Z[i,j] = p[3]
        end
    end

    eci_basis = diagm([1.0,1.0,1.0])

    axes_tails = [Point3f([0,0,0]) for k in 1:3]
    eci_heads = [Vec3f(arrow_length*eci_basis[:,k]) for k in 1:3]
    ecef_heads = [Vec3f(arrow_length*C_IF*eci_basis[:,k]) for k in 1:3]

    ecef_label_pos = [1.2 * pt for pt in ecef_heads]
    eci_label_pos = [1.2 * pt for pt in eci_heads]

    colors = [:red, :green, :blue]

    surface!(
        ax,
        mesh_X,
        mesh_Y,
        mesh_Z,
        color=texture_ecef,
        diffuse=0.8
    )
    arrows!(
        ax,
        axes_tails,
        ecef_heads,
        color=colors,
        linewidth=arrow_width,
        arrowsize=Vec3f([arrow_width, arrow_width, arrow_width*2])
    )
    arrows!(
        ax,
        axes_tails,
        eci_heads,
        color=colors,
        linewidth=arrow_width,
        arrowsize=Vec3f([arrow_width, arrow_width, arrow_width*2])
    )

    text!(
        ax,
        ecef_label_pos[1:2],
        text=[L"F_1", L"F_2"],
        align=(:center, :center),
        fontsize=20,
        color=:black
        # markerspace=:data
    )
    text!(
        ax,
        eci_label_pos,
        text=[L"I_1", L"I_2", L"I_3,\ F_3"],
        align=(:center, :center),
        fontsize=20,
        color=:black
        # markerspace=:data
    )

end

function main()
    eops = SatelliteToolboxTransformations.fetch_iers_eop()
    texture = load_earth_texture_to_ecef(joinpath("assets","map_pol2.png"))
    t_jd_s = 2460830.26875 * 24 * 3600 - 3600*10

    GLMakie.activate!(title="OpenSOAP")
    
    fig = Figure(size=(800,800))
    display(fig)

    sun_I = 5e6*Vec3f(SatelliteToolboxCelestialBodies.sun_position_mod(t_jd_s/3600/24))
    dl = DirectionalLight(RGBf(243/255, 241/255, 218/255), sun_I)
    al = AmbientLight(RGBf(0.2, 0.2, 0.2))
    al = AmbientLight(RGBf(0.8, 0.8, 0.8))
    al = AmbientLight(RGBf(243/255, 241/255, 230/255))

    ax = LScene(
        fig[1,1], 
        show_axis=false, 
        scenekw=(
            # lights=[dl, al], 
            lights = [al],
            # backgroundcolor=:black, 
            clear=true
        )
    )

    plot_globe!(ax, t_jd_s, eops, texture)
    orb = plot_orbit!(ax, t_jd_s, eops; do_ned = false)
    plot_body!(ax, t_jd_s, eops, "", orb)
    save(joinpath("cases", "eci_ecef_fig.png"), fig)
end