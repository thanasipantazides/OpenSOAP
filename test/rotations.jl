using Test
using OpenSOAP
using GeometryBasics
using LinearAlgebra

@testset "rotation tools" begin

    @testset "residuals" begin
        a = rand(3)
        b = rand(3)
        R = rand(3,3)
        
        tol = 1e-12

        @test residualso3(OpenSOAP.cross(a)) < 1e-12
        @test residualso3(R) > tol

        @test residualSO3(R) > tol
        @test residualSO3(I(3)) < tol
        @test residualSO3(projSO3(R)) < tol
    end
    
    @testset "minimum arc slew" begin
        X = Vec3d(1.0,0.0,0.0)
        Y = Vec3d(0.0,1.0,0.0)
        Z = Vec3d(0.0,0.0,1.0)
        a = rand(3)
        b = rand(3)

        tol = 1e-3
    
        @test norm(r_min_arc(X,Y)*X .- Y) < tol
        @test norm(r_min_arc(Y,Z)*Y .- Z) < tol
        @test norm(r_min_arc(Z,X)*Z .- X) < tol
    
        @test residualSO3(r_min_arc(a,b)) < 1e-12
        @test norm(r_min_arc(a,b)*a .- b) < tol
        @test norm(r_min_arc(b,a)*b .- a) < tol
    
        @test sum(r_min_arc(a,a) .- I(3)) < 1e-12
    end

end