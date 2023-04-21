using MAT

include("../src/entropy.jl")
using .Entropy: calculate_all_entropies

# Read
vars = matread("resources/connectivities_estimates_withsurr_twoAALnets.mat")
DMNlist = vars["DMNlist"]
dmn_surr = vars["dmn_results"]["surr"]  # size N x 10

# Calculate bins' boundaries
bins_n = 4
variables_n = length(DMNlist)
upper_boundaries = zeros(Float64, bins_n, variables_n)
lower_boundaries = zeros(Float64, bins_n, variables_n)
data_size = length(dmn_surr[:, 1])
elements_per_bin = data_size / bins_n
for i=1:variables_n
    column = dmn_surr[:,i]
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
        # @show left_i, right_i
    end
end

# Diskretize data
dmn_diskretized_surr = zeros(Int64, size(dmn_surr))
for i=1:length(DMNlist)
    for j=1:data_size
        for t=1:bins_n
            if lower_boundaries[t, i] <= dmn_surr[j, i] <= upper_boundaries[t, i]
                dmn_diskretized_surr[j, i] = t
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
    row = dmn_diskretized_surr[row_i, :]
    el_count[row...] += 1
end

for i=1:length(el_count)
    prob_table[i] = el_count[i]/data_size
end

@show "Calculating entropies"

# Calculate entropies
entropies = calculate_all_entropies(prob_table)

# Claculate Connected Information
connected_information = estimate_connected_information(10, distr_cards, entropies)

# Wrtie
matwrite(
    "resources/dmn_diskretized_surr_4.mat", 
    Dict(
        "dmn_diskretized_surr" => dmn_diskretized_surr,
        "dmn_prob_table" => prob_table,
        "dmn_connected_information" => connected_information))
