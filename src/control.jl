using LinearAlgebra
using JuMP, Ipopt
import HSL_jll



function point_between(start_state::State{<:Real}, stop_state::State{<:Real}, dt::Real, params, verbose::Bool=false)
    model = Model(Ipopt.Optimizer)
    set_optimizer_attribute(model, "fixed_variable_treatment", "relax_bounds")

    set_attribute(model, "hsllib", HSL_jll.libhsl_path)
    set_attribute(model, "linear_solver", "ma57")

    if !verbose
        set_silent(model)
    end

    J = convert(Matrix{Float64}, params.mission.spacecraft.mass.inertia)
    Jinv = inv(J)

    # THIS IS FUDGED! Better way is to inscribe an ellipsoid in the agilitoid. Best way is just use the agilitoid.
    mlim = params.mission.spacecraft.attitude.wheels.torque[1:3]

    # stop and start attitude/angular velocity
    C_BI0 = start_state.attitude
    C_BIf = stop_state.attitude
    w0 = start_state.angular_velocity
    wf = stop_state.angular_velocity

    # time stepping
    steps_per_second = 10
    n = Int(ceil(dt * steps_per_second))
    tstep = dt / n

    # splat control bounds
    mlimmat = hcat([mlim for i in 1:n]...)

    # optimization tolerances/limits
    psinf = 1e12
    so3_tol = 1e-4
    angrate_tol = 1e-6
    ident3 = diagm([1.0, 1.0, 1.0])

    # splat boundary conditions
    Climlo = fill(-1.0, (3, 3, n))
    Climhi = fill(1.0, (3, 3, n))
    wlimmatlo = fill(-psinf, (3, n))
    wlimmathi = fill(psinf, (3, n))
    Climlo[:, :, 1] = C_BI0
    Climhi[:, :, 1] = C_BI0
    Climlo[:, :, n] = C_BIf
    Climhi[:, :, n] = C_BIf
    wlimmatlo[:, 1] = w0
    wlimmathi[:, 1] = w0
    wlimmatlo[:, n] .= wf
    wlimmathi[:, n] .= wf

    # find reasonable initial guesses for decision variables
    C_guess = rotinterp(C_BI0, C_BIf, n)

    # bounded control:
    @variable(model, u[i=1:3, k=1:n])

    # angular rates:
    @variable(model, w_B[i=1:3, k=1:n])

    # attitude:
    @variable(model, C_BI[i=1:3, j=1:3, k=1:n])

    # Figure out if this syntax is doing what I expect:
    for i in 1:3
        for j in 1:3
            for k in 1:n
                set_start_value(C_BI[i, j, k], C_guess[i, j, k])
            end
        end
    end

    # control constraint:
    @constraint(model, u .<= mlimmat)
    @constraint(model, u .>= -mlimmat)

    # angular rate boundary conditions
    @constraint(model, w_B .<= wlimmathi)
    @constraint(model, w_B .>= wlimmatlo)
    # attitude boundary conditions
    @constraint(model, C_BI .<= Climhi)
    @constraint(model, C_BI .>= Climlo)


    # SO(3) constraints
    @constraints(model, begin
        [k = 2:n-1], ident3 .- so3_tol .<= C_BI[:, :, k]' * C_BI[:, :, k] .<= ident3 .+ so3_tol
    end)

    # explicit determinant for 3x3 matrix
    op_det_func(X) = X[1, 1] * X[2, 2] * X[3, 3] + X[1, 2] * X[2, 3] * X[3, 1] + X[1, 3] * X[2, 1] * X[3, 2] - X[1, 3] * X[2, 2] * X[3, 1] - X[1, 2] * X[2, 1] * X[3, 3] - X[1, 1] * X[2, 3] * X[3, 2]

    @operator(model, op_det, 1, op_det_func)
    @constraints(model, begin
        [k = 1:n], 1.0 - so3_tol <= op_det_func(C_BI[:, :, k]) <= 1.0 + so3_tol
    end)

    # dynamic constraints
    ddt_w(x::Matrix, k::Int) = (x[:, k] .- x[:, k-1]) ./ tstep
    ddt_C(X::Array, w::Array, k::Int) = (X[:, :, k] .- X[:, :, k-1]) ./ tstep

    @constraint(
        model,
        [k = 2:n],
        ddt_w(w_B, k) .<= Jinv * (u[:, k-1] - cross(w_B[:, k-1], J * w_B[:, k-1])) .+ angrate_tol
    )
    @constraint(
        model,
        [k = 2:n],
        ddt_w(w_B, k) .>= Jinv * (u[:, k-1] - cross(w_B[:, k-1], J * w_B[:, k-1])) .- angrate_tol
    )
    @constraint(
        model,
        [k = 2:n],
        ddt_C(C_BI, w_B[:, k], k) .<= -[0 -w_B[3, k-1] w_B[2, k-1];
            w_B[3, k-1] 0 -w_B[1, k-1];
            -w_B[2, k-1] w_B[1, k-1] 0] * C_BI[:, :, k-1] .+ so3_tol
    )
    @constraint(
        model,
        [k = 2:n],
        ddt_C(C_BI, w_B[:, k], k) .>= -[0 -w_B[3, k-1] w_B[2, k-1];
            w_B[3, k-1] 0 -w_B[1, k-1];
            -w_B[2, k-1] w_B[1, k-1] 0] * C_BI[:, :, k-1] .- so3_tol
    )

    @objective(model, Min, sum(u[1, :] .^ 2 + u[2, :] .^ 2 + u[3, :] .^ 2))

    optimize!(model)  # Solve for the control and state

    result_state = Vector{State}(undef, n)

    for k in 1:n
        result_state[k] = State{Float64}(zeros(3), zeros(3), value.(w_B[:, k]), value.(C_BI[:, :, k]), 0.0, 0.0, 0.0)
    end

    return result_state
end
