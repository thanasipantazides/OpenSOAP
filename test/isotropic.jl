
import GLMakie
using OpenSOAP

function coverage()
    GLMakie.activate!()
    fig = GLMakie.Figure()
    display(fig)

    ax = GLMakie.LScene(fig[1,1])
    for _ in 1:10^3
        u = projSO3(rand(3,3) .- 0.5)*normalize(rand(3))
        GLMakie.lines!(ax, [0, u[1]], [0, u[2]], [0, u[3]], color=:black, linewidth=0.5)
    end
end