using Test

@testset "OpenSOAP" begin
    include("rotations.jl")
    include("bounded.jl")
    include("loading.jl")
    include("network.jl")
end
