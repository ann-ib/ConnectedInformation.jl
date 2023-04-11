using Test

using .Entropy: caluclate_all_entropies, calculate_entropy

@testset "entropy.jl" begin
    @testset "caluclate_all_entropies" begin
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

    @testset "calculate_entropy" begin
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