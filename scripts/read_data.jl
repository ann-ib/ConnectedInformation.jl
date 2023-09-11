using MAT
using Serialization

include("../src/ConnectedInformation.jl")
include("../src/entropy.jl")
include("../src/helpers.jl")

import .ConnectedInformation: estimate_max_entropies, estimate_connected_information
import .Entropy: calculate_all_entropies
import .Helpers: discretize_data

SURR = true  # compute for pure data or surrogate 
bins_n = 2  # diskretization number

# Read
vars = matread("resources/connectivities_estimates_withsurr_twoAALnets.mat")
DMNlist = vars["DMNlist"]
variables_n = length(DMNlist)
dmn_data = vars["dmn_results"]["data"]  # size N x 10
dmn_surr = vars["dmn_results"]["surr"]  # size N x 10

if SURR
    cur_data = dmn_surr
else
    cur_data = dmn_data
end

dmn_diskretized = discretize_data(cur_data, variables_n, bins_n)

# Calculate probability table 
data_size = length(cur_data[:, 1])
distr_cards = [bins_n for i=1:variables_n]
prob_table = zeros(Float64, distr_cards...)
el_count = zeros(Int64, distr_cards...)
for row_i=1:data_size
    row = dmn_diskretized[row_i, :]
    el_count[row...] += 1
end

for i=1:length(el_count)
    prob_table[i] = el_count[i]/data_size
end

# Calculate entropies
@show "Calculating entropies"
entropies = calculate_all_entropies(prob_table)

# Claculate approximations
@show "Calculating max entropies"
max_entropies_upper = estimate_max_entropies(10, distr_cards, entropies)
max_entropies_lower = estimate_max_entropies(10, distr_cards, entropies; lower_bound = true)
@show "Calculating connected information"
connected_information_upper = estimate_connected_information(10, distr_cards, entropies)
connected_information_lower = estimate_connected_information(10, distr_cards, entropies; lower_bound = true)

# Wrtie
# filename = "dmn_diskretized_$(if SURR; "surr"; else; "data"; end)_$(K)"

# matwrite(
#     "resources/$(filename).mat",
#     Dict(
#         "diskretized" => dmn_diskretized,
#         "prob_table" => prob_table,
#         "max_entropies_upper" => max_entropies_upper,
#         "max_entropies_lower" => max_entropies_lower,
#         "connected_information_upper" => connected_information_upper,
#         "connected_information_lower" => connected_information_lower))

# serialize(
#     "resources/$(filename).ser",
#     Dict(
#         "diskretized" => dmn_diskretized,
#         "prob_table" => prob_table,
#         "entropies" => entropies,
#         "max_entropies_upper" => max_entropies_upper,
#         "max_entropies_lower" => max_entropies_lower,
#         "connected_information_upper" => connected_information_upper,
#         "connected_information_lower" => connected_information_lower))