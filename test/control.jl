using OpenSOAP
using LinearAlgebra
using InfiniteOpt, Ipopt
import HSL_jll
import MathOptInterface as MOI
using GLMakie, GeometryBasics


# A working example to slew between an initial and target attitude, with zero'd angular rates at both boundaries.
# Uses direct collocation.
# Dynamic constraints, SO(3)-membership constraints, boundary conditions enforced.
# Uses IPOPT with HSL linear solver (requires free license).
# Uses DCM parameterization.
# Minimizes RMS control energy cost over the maneuver.
# Will embed in OpenSOAP and integrate with broader sim.
function jump_opt_vec()
    C_BI0 = diagm([1.0,1.0,1.0])
    C_BIf = r_random()
    # w0 = [0.01, 0.02, 0.03]
    w0 = [0.0, 0.0, 0.0]
    wf = [0.0, 0.0, 0.0]

    J = diagm([0.05, 0.1, 0.15])
    Jinv = inv(J)

    mlim = [2.0e-3, 2.0e-3, 2.0e-3]

    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "fixed_variable_treatment", "relax_bounds")
    
    set_attribute(model, "hsllib", HSL_jll.libhsl_path)
    set_attribute(model, "linear_solver", "ma57")
    # set_attribute(model, "linear_solver", "ma97")
    # set_attribute(model, "linear_solver", "ma27")
    # set_attribute(model, "ma27_meminc_factor", 4)
    # set_optimizer_attribute(model, "mumps_mem_percent", 4000)
    
    n = 1000
    Tf = 60
    dt = Tf/n

    manstep = Int64(ceil(n - n/4))

    # splat control bounds
    mlimmat = hcat([mlim for i in 1:n]...)

    psinf = 1e12
    ident3 = diagm([1.0, 1.0, 1.0])
    
    # splat boundary conditions
    Climlo = fill(-1.0, (3,3,n))
    Climhi = fill(1.0, (3,3,n))
    wlimmatlo = fill(-psinf, (3, n))
    wlimmathi = fill(psinf, (3, n))
    Climlo[:,:,1] = C_BI0
    Climhi[:,:,1] = C_BI0
    Climlo[:,:,n] = C_BIf
    Climhi[:,:,n] = C_BIf
    wlimmatlo[:, 1] = w0
    wlimmathi[:, 1] = w0
    # wlimmatlo[:, n] = wf
    # wlimmathi[:, n] = wf
    wlimmatlo[:, manstep:n] .= wf
    wlimmathi[:, manstep:n] .= wf

    # bounded control:
    @variable(model, u[i = 1:3, k = 1:n])

    # angular rates:
    # @variable(model, w_B[i = 1:3, k = 1:n], start = hcat([w0 for i in 1:n]...))
    @variable(model, w_B[i = 1:3, k = 1:n])

    # attitude:
    @variable(model, C_BI[i=1:3, j=1:3, k=1:n])

    # control constraint:
    @constraint(model, u .<=  mlimmat)
    @constraint(model, u .>= -mlimmat)

    # angular rate boundary conditions
    @constraint(model, w_B .<= wlimmathi)
    @constraint(model, w_B .>= wlimmatlo)
    # attitude boundary conditions
    @constraint(model, C_BI .<= Climhi)
    @constraint(model, C_BI .>= Climlo)

    so3_tol = 1e-4
    angrate_tol = 1e-6
    # SO(3) constraints
    @constraints(model, begin
        [k=2:n-1], ident3 .- so3_tol .<= C_BI[:,:,k]'*C_BI[:,:,k] .<= ident3 .+ so3_tol
    end);

    # explicit determinant for 3x3 matrix
    op_det_func(X) = X[1,1]*X[2,2]*X[3,3] + X[1,2]*X[2,3]*X[3,1] + X[1,3]*X[2,1]*X[3,2] - X[1,3]*X[2,2]*X[3,1] - X[1,2]*X[2,1]*X[3,3] - X[1,1]*X[2,3]*X[3,2]

    @operator(model, op_det, 1, op_det_func)
    @constraints(model, begin
        [k=1:n], 1.0 - so3_tol <= op_det_func(C_BI[:,:,k]) <= 1.0 + so3_tol
    end);

    # dynamic constraints
    ddt_w(x::Matrix, k::Int) = (x[:,k] .- x[:,k-1]) ./ dt
    ddt_C(X::Array, w::Array, k::Int) = (X[:,:,k] .- X[:,:,k-1]) ./ dt
    # ddt_mat(X::Array, w::Array, k::Int) = exp(-(
    #     [ 0 -w[3,k-1] w[2,k-1]; 
    #     w[3,k-1] 0 -w[1,k-1]; 
    #     -w[2,k-1] w[1,k-1] 0]) .* dt)
    # ddt_mat(X::Array, w::Vector, k::Int) = exp(-(
    #     [ 0 -w[3] w[2]; 
    #     w[3] 0 -w[1]; 
    #     -w[2] w[1] 0]) .* dt)

    @constraint(
        model,
        [k=2:n], 
        ddt_w(w_B, k) .<= Jinv*(u[:,k-1] - cross(w_B[:,k-1], J*w_B[:,k-1])) .+ angrate_tol
    );
    @constraint(
        model,
        [k=2:n], 
        ddt_w(w_B, k) .>= Jinv*(u[:,k-1] - cross(w_B[:,k-1], J*w_B[:,k-1])) .- angrate_tol
    );
    @constraint(
        model,
        [k=2:n], 
        ddt_C(C_BI, w_B[:,k], k) .<= -[ 0 -w_B[3,k-1] w_B[2,k-1]; 
            w_B[3,k-1] 0 -w_B[1,k-1]; 
            -w_B[2,k-1] w_B[1,k-1] 0]*C_BI[:,:,k-1] .+ so3_tol
    );
    @constraint(
        model,
        [k=2:n],
        ddt_C(C_BI, w_B[:,k], k) .>= -[ 0 -w_B[3,k-1] w_B[2,k-1]; 
        w_B[3,k-1] 0 -w_B[1,k-1]; 
        -w_B[2,k-1] w_B[1,k-1] 0]*C_BI[:,:,k-1] .- so3_tol 
    );

    @objective(model, Min, sum(u[1,:].^2 + u[2,:].^2 + u[3,:].^2))
    # @objective(model, Min, sum(
    #     w_B[1,manstep:n].^2 + 
    #     w_B[2,manstep:n].^2 + 
    #     w_B[3,manstep:n].^2) + sum(
    #     u[1,:].^2 + 
    #     u[2,:].^2 + 
    #     u[3,:].^2
    #     )
    # );

    optimize!(model)  # Solve for the control and state
    
    GLMakie.activate!()
    t = 1:n
    fig = Figure(size=(500,500))
    axu = Axis(
        fig[2,3],
        xlabel="i", 
        ylabel="u"
    )
    axw = Axis(
        fig[3,3],
        xlabel="i", 
        ylabel="w"
    )
    axc = Axis(
        fig[1,3],
        xlabel="i",
    )
    
    axm = LScene(
        fig[1:3,1:2]
    )

    dethist = zeros(n)
    trhist = zeros(n)
    for i in t
        dethist[i] = det(value.(C_BI[:,:,i]))
        trhist[i] = sum((value.(C_BI[:,:,i]))' * value.(C_BI[:,:,i])) .- 2
    end

    set_theme!(theme_light())
    lines!(axu, t, value.(u[1,:]))
    lines!(axu, t, value.(u[2,:]))
    lines!(axu, t, value.(u[3,:]))
    scatter!(axw, t, value.(w_B[1,:]))
    scatter!(axw, t, value.(w_B[2,:]))
    scatter!(axw, t, value.(w_B[3,:]))
    scatter!(axc, t, dethist, label="determinant")
    scatter!(axc, t, trhist, label="orthogonality")

    basis = diagm([1.0,1.0,1.0])
    cols = [:red, :green, :blue]

    println("C_BI Z check max: ", max(value.(C_BI[3,3,:])...))
    println("C_BI det min/max: ", min(dethist...), " / ", max(dethist...))
    println("C_BI orth min/max: ", min(trhist...), " / ", max(trhist...))

    for b in 1:3
        for i in 1:10:n
            L = value.(C_BI[:,:,i]*basis[:,b])
            lines!(
                axm, 
                [0, L[1]],
                [0, L[2]],
                [0, L[3]],
                color=cols[b]
            )
        end
        lines!(
            axm, 
            [0, (value.(C_BI[:,:,n])*basis[:,b])[1]],
            [0, (value.(C_BI[:,:,n])*basis[:,b])[2]],
            [0, (value.(C_BI[:,:,n])*basis[:,b])[3]],
            linestyle=:dash,
            color=:black
        )
    end

    axislegend(axc)
    display(fig) 
end