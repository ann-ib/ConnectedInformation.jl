using Test

include("../src/entropy.jl")

@testset "ConnectedInformation" begin
    include("ConnectedInformation_test.jl")
    include("entropy_test.jl")
end