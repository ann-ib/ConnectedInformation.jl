using ConnectedInformation
using Test

@testset "ConnectedInformation.jl" begin
    @testset "connected information approximation" begin
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
        @show(connected_information)
    end

    @testset "max entropy approximation" begin
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

    @testset "calculation of all entropies" begin
        parity = zeros(Float64, 2, 2, 2, 2)
        parity[1, 1, 1, 1] = 1/8        
        parity[1, 1, 2, 2] = 1/8
        parity[1, 2, 1, 2] = 1/8
        parity[1, 2, 2, 1] = 1/8
        parity[2, 1, 1, 2] = 1/8
        parity[2, 1, 2, 1] = 1/8
        parity[2, 2, 1, 1] = 1/8
        parity[2, 2, 2, 2] = 1/8

        expected_entropies = Dict([] => 0,
                                  [1] => 1.0,
                                  [2] => 1.0,
                                  [3] => 1.0,
                                  [4] => 1.0,
                                  [1, 2] => 2.0,
                                  [1, 3] => 2.0,
                                  [1, 4] => 2.0,
                                  [2, 3] => 2.0,
                                  [2, 4] => 2.0,
                                  [3, 4] => 2.0,
                                  [1, 2, 3] => 3.0,
                                  [1, 2, 4] => 3.0,
                                  [1, 3, 4] => 3.0,
                                  [2, 3, 4] => 3.0,
                                  [1, 2, 3, 4] => 3.0)

        @test caluclate_all_entropies(parity) == expected_entropies
    end

    @testset "entropy calculation" begin
        eight_size_vector = [0.14, 0.07, 0.04, 0.13, 0.15, 0.08, 0.17, 0.22]

        empty = Float64[]
        one_distr = [0.25, 0.3, 0.05, 0.4]
        two_distrs = reshape(eight_size_vector, 2, 4)
        three_distrs = reshape(eight_size_vector, 2, 2, 2)
        test_cases = [empty, one_distr, two_distrs, three_distrs]

        exp_empty = 0.0
        exp_one_distr = 1.76595732 
        exp_two_distrs = 2.851277264
        exp_three_distrs = 2.851277264
        exp_results = [exp_empty, exp_one_distr, exp_two_distrs, exp_three_distrs]

        @testset for (probability_table, exp_result) in zip(test_cases, exp_results)
            @test calculate_entropy(probability_table) ≈ exp_result
        end
    end
end
