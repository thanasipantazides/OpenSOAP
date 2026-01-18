using OpenSOAP

function check_slew()

    C_AI = diagm([1.0,1.0,1.0])
    C_AI = r_euler3(pi/2)
    p_A = [1.0;0.0;0.0]
    r_A = [0.0;1.0;0.0]

    C_BA = r_min_arc(C_AI*p_A, C_AI*r_A)
    @show C_AI
    @show C_BA
    @show C_BA*C_AI

    @assert all((C_BA*C_AI - r_euler3(pi)).^2 .< 1e-3)
end