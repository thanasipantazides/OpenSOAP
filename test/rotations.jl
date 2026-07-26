using Test
using OpenSOAP
using GeometryBasics
using LinearAlgebra

import GLMakie

@testset "rotation tools" begin
    tol = 1e-12
    a = rand(3)
    b = rand(3)
    R = rand(3,3)

    @testset "cross" begin
        @test norm(OpenSOAP.uncross(OpenSOAP.cross(a)) .- a) < tol
    end
    
    @testset "residuals" begin
        @test residualso3(OpenSOAP.cross(a)) < tol
        @test residualso3(R) > tol

        @test residualSO3(R) > tol
        @test residualSO3(I(3)) < tol
        @test residualSO3(projSO3(R)) < tol
    end
    
    @testset "Lie group exp and log maps" begin
        @test residualSO3(exp(OpenSOAP.cross(a))*projSO3(R)) < tol
        @test residualso3(log(projSO3(R))) < tol
    end
    
    @testset "minimum arc slew" begin
        X = Vec3d(1.0,0.0,0.0)
        Y = Vec3d(0.0,1.0,0.0)
        Z = Vec3d(0.0,0.0,1.0)
        a = normalize(rand(3))
        b = normalize(rand(3))

        # tol = 1e-9
    
        @test norm(r_min_arc(X,Y)*X .- Y) < tol
        @test norm(r_min_arc(Y,Z)*Y .- Z) < tol
        @test norm(r_min_arc(Z,X)*Z .- X) < tol

        for _ in 1:10^3
            u = projSO3(rand(3,3) .- 0.5)*normalize(rand(3))
            v = projSO3(rand(3,3) .- 0.5)*normalize(rand(3))
            
            @test residualSO3(r_min_arc(u,v)) < tol
            @test norm(r_min_arc(u,v)*u .- v) < tol
            @test norm(r_min_arc(v,u)*v .- u) < tol
        
            @test sum(r_min_arc(u,u) .- I(3)) < tol
            @test sum(r_min_arc(u,v)*r_min_arc(v,u) .- I(3)) < tol
        end
    end
end
