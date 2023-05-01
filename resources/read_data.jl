using MAT
using Serialization

using ConnectedInformation

include("../src/entropy.jl")
using .Entropy: calculate_all_entropies

SURR = true  # compute for pure data or surrogate 
K = 2  # diskretization number

# Read
vars = matread("resources/connectivities_estimates_withsurr_twoAALnets.mat")
DMNlist = vars["DMNlist"]
dmn_data = vars["dmn_results"]["data"]  # size N x 10
dmn_surr = vars["dmn_results"]["surr"]  # size N x 10

if SURR
    cur_data = dmn_surr
else
    cur_data = dmn_data
end

# Calculate bins' boundaries
bins_n = K
variables_n = length(DMNlist)
upper_boundaries = zeros(Float64, bins_n, variables_n)
lower_boundaries = zeros(Float64, bins_n, variables_n)
data_size = length(cur_data[:, 1])
elements_per_bin = data_size / bins_n
for i=1:variables_n
    column = cur_data[:,i]
    sorted_column = sort(column)
    (left_i, right_i) = (0, 0)
    for j=1:bins_n
        left_i = right_i + 1
        right_i = floor(Int, j * elements_per_bin)
        if j == bins_n
            right_i = data_size
        end
        lower_boundaries[j, i] = sorted_column[left_i]
        upper_boundaries[j, i] = sorted_column[right_i]
    end
end

# Diskretize data
dmn_diskretized = zeros(Int64, size(cur_data))
for i=1:length(DMNlist)
    for j=1:data_size
        for t=1:bins_n
            if lower_boundaries[t, i] <= cur_data[j, i] <= upper_boundaries[t, i]
                dmn_diskretized[j, i] = t
                break
            end
        end
    end
end

# Calculate probability table 
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
filename = "dmn_diskretized_$(if SURR; "surr"; else; "data"; end)_$(K)"

matwrite(
    "resources/$(filename).mat",
    Dict(
        "diskretized" => dmn_diskretized,
        "prob_table" => prob_table,
        "max_entropies_upper" => max_entropies_upper,
        "max_entropies_lower" => max_entropies_lower,
        "connected_information_upper" => connected_information_upper,
        "connected_information_lower" => connected_information_lower))

serialize(
    "resources/$(filename).ser",
    Dict(
        "diskretized" => dmn_diskretized,
        "prob_table" => prob_table,
        "entropies" => entropies,
        "max_entropies_upper" => max_entropies_upper,
        "max_entropies_lower" => max_entropies_lower,
        "connected_information_upper" => connected_information_upper,
        "connected_information_lower" => connected_information_lower))