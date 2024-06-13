# Test for using Makie observables successfully, e.g. to drive a plot over time using a slider.

# Use results of this to make test/plot.jl useful and scrollable.

using GLMakie

function plot_surf!(ax::Makie.LScene)
    nθ = 31
    nζ = 31

    R = 1
    b = 0.5
    
    θ = LinRange(0, 2π, nθ)
    ζ = LinRange(-1, 1, nζ)
    X = zeros(nθ, nζ)
    Y = zeros(nθ, nζ)
    Z = zeros(nθ, nζ)
    for i in 1:nθ
        for j in 1:nζ
            X[i,j] = R*sqrt(1 + ζ[j]^2/b^2)*cos(θ[i])
            Y[i,j] = R*sqrt(1 + ζ[j]^2/b^2)*sin(θ[i])
            Z[i,j] = ζ[j]
        end
    end

    surface!(ax, X, Y, Z)
end

function get_trajectory()
    nθ = 100

    θ = LinRange(0, 2π, nθ)
    ζ = LinRange(-1, 1, nθ)
    R = 2
    b = 2

    x = zeros(1, nθ)
    y = zeros(1, nθ)
    z = zeros(1, nθ)

    for i in 1:nθ
        x[i] = R*cos(θ[i])
        y[i] = R*sin(θ[i])
        z[i] = b*ζ[i]
    end
    return [x;y;z]
end

function plot_anim()
    GLMakie.activate!()
    r = get_trajectory()
    n = length(r[1,:])

    fig = Figure(size=(1400,840))
    ax = LScene(fig[1,1], show_axis=false)
    plot_surf!(ax)
    lines!(
        ax,
        r[1,:],
        r[2,:],
        r[3,:]
    )

    time_sliders = SliderGrid(fig[2, 1], (label = "time", format = "{:d}", range = 1:1:n, startvalue = 1))

    slider_observables = [s.value for s in time_sliders.sliders]
    print(slider_observables)
    point = lift(slider_observables[1]) do t_i
        println(t_i)
        Point3f(r[:,t_i])
    end
    dot = scatter!(ax, point, markersize=10, color=:red)
    # points = Observable(Point3f[r[:,1]])
    # dot = scatter!(ax, points, markersize=10, color=:red)
    # fps = 30
    # for i = 2:n
    #     # new_point = Point3f(r[:,i])
    #     # points[] = push!(points[], new_point)
    #     display(fig)
    #     sleep(1/fps)
    # end

    display(fig)
end