using MAT
using Serialization

include("../src/entropy.jl")
using .Entropy


include("../src/helpers.jl")
using .Helpers

DISTRIBUTIONS_NUM = 1000

prob_tables = Vector{Array{Float64}}()
entropies = Vector{Dict{Vector{Int64}, Float64}}()
for i in 1:DISTRIBUTIONS_NUM
    prob_table = generate_prob_table(2, 2, 2, 2)
    prob_table_entropies = calculate_all_entropies(prob_table)
    push!(prob_tables, prob_table)
    push!(entropies, prob_table_entropies)
end

filename = "random_nsb_1000"

matwrite(
    "resources/comparison_to_nsb/$(filename).mat",
    Dict("prob_tables" => prob_tables))

serialize(
    "resources/comparison_to_nsb/$(filename).ser",
    Dict(
        "prob_tables" => prob_tables,
        "entropies" => entropies))