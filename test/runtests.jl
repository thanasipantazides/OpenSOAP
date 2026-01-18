using OpenSOAP, Test, LinearAlgebra

@testset "rotation" begin
    n = 1000
    dets = zeros(n)
    traces = zeros(n)
    for i = 1:n
        a = rand(3)
        b = rand(3)
        C = r_min_arc(a/norm(a), b/norm(b))
        dets[i] = det(C)
        traces[i] = sum(diag(C'*C)) - 3
    end
    @test all(dets .≈ 1.0)
    @test all(abs.(traces) .< 1e-9)
end
@testset "chained rotations" begin
    C_AI = diagm([1.0,1.0,1.0])
    C_AI = r_euler3(pi/2)
    p_A = [1.0;0.0;0.0]
    r_A = [0.0;1.0;0.0]

    C_BA = r_min_arc(C_AI*p_A, C_AI*r_A)

    @test all((C_BA*C_AI - r_euler3(pi)).^2 .< 1e-3)
    @test all(C_BA*C_AI*p_A - r_A .< 1e-3)
end
