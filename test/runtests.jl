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
