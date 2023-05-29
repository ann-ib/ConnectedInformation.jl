using MAT
using Serialization

include("../src/ConnectedInformation.jl")
using .ConnectedInformation

include("../src/entropy.jl")
using .Entropy

include("../src/helpers.jl")
using .Helpers


upper_max_entropies_results = Vector{Vector{Float64}}()
lower_con_inf_results = Vector{Vector{Float64}}()
for i in 1:1000
    prob_table = generate_prob_table(2, 2, 2, 2)
    entropies = calculate_all_entropies(prob_table)
    estimated_max_entropies = estimate_max_entropies(4, [2, 2, 2, 2], entropies)
    push!(upper_max_entropies_results, estimated_max_entropies)
    connected_information = calculate_connected_information(4, estimated_max_entropies)
    push!(lower_con_inf_results, connected_information)
end

filename = "random_1000"

matwrite(
    "resources/random_distributions_results/$(filename).mat",
    Dict(
        "upper_max_entropies_results" => upper_max_entropies_results,
        "lower_con_inf_results" => lower_con_inf_results))

serialize(
    "resources/random_distributions_results/$(filename).ser",
    Dict(
        "upper_max_entropies_results" => upper_max_entropies_results,
        "lower_con_inf_results" => lower_con_inf_results))
