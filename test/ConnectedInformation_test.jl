using ConnectedInformation
using Test

using .Entropy

@testset "ConnectedInformation.jl" begin
    @testset "estimate_connected_information" begin
        @testset "parity" begin
            parity = zeros(Float64, 2, 2, 2, 2)
            parity[1, 1, 1, 1] = 1/8
            parity[1, 1, 2, 2] = 1/8
            parity[1, 2, 1, 2] = 1/8
            parity[1, 2, 2, 1] = 1/8
            parity[2, 1, 1, 2] = 1/8
            parity[2, 1, 2, 1] = 1/8
            parity[2, 2, 1, 1] = 1/8
            parity[2, 2, 2, 2] = 1/8
            entropies = caluclate_all_entropies(parity)
            
            distr_cards = [2, 2, 2, 2]
            connected_information = estimate_connected_information(4, distr_cards, entropies)
            @test isnan(connected_information[1])
            @test connected_information[2:end] ≈ [0.0, 0.0, 1.0]
        end

        @testset "random values" begin
            distr_cards = [3, 2, 4]
            rnd_values = Array{Float64,3}(undef, distr_cards...)
            rnd_values[:,:,1] = [0.0732422  0.0258789;
                                 0.124512   0.0473633;
                                 0.0322266  0.000976562]
            rnd_values[:,:,2] = [0.0102539  0.00732422;
                                 0.0317383  0.120605;
                                 0.0791016  0.0595703;]
            rnd_values[:,:,3] = [0.00390625  0.0249023;
                                 0.0292969   0.00927734;
                                 0.0356445   0.0893555]
            rnd_values[:,:,4] = [0.0209961  0.0439453;
                                 0.0839844  0.0219727;
                                 0.0136719  0.0102539]
            entropies = caluclate_all_entropies(rnd_values)
            connected_information = estimate_connected_information(3, distr_cards, entropies)
            @test isnan(connected_information[1])
            # TODO: calculate real maximum value
            @test connected_information[2:end] ≈ [0.20138705009070534, 0.17931841065010445]
        end
    end

    @testset "estimate_max_enthropy" begin
        @testset "parity" begin
            distr_cards = [2, 2, 2, 2]

            one_dim_entropies = Dict{Vector{Int64}, Float64}([1] => 1, [2] => 1, [3] => 1, [4] => 1)
            @test estimate_max_enthropy(1, distr_cards, one_dim_entropies) ≈ 4.0
            
            parity = zeros(Float64, 2, 2, 2, 2)
            parity[1, 1, 1, 1] = 1/8
            parity[1, 1, 2, 2] = 1/8
            parity[1, 2, 1, 2] = 1/8
            parity[1, 2, 2, 1] = 1/8
            parity[2, 1, 1, 2] = 1/8
            parity[2, 1, 2, 1] = 1/8
            parity[2, 2, 1, 1] = 1/8
            parity[2, 2, 2, 2] = 1/8
            parity_entropies = caluclate_all_entropies(parity)
            @test estimate_max_enthropy(4, distr_cards, parity_entropies) == 3.0
        end

        @testset "random values" begin
            
        end
    end
end