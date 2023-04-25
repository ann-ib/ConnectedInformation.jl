using MAT

using ConnectedInformation

include("../src/entropy.jl")
using .Entropy: calculate_all_entropies

# Read
vars = matread("resources/connectivities_estimates_withsurr_twoAALnets.mat")
DMNlist = vars["DMNlist"]
dmn_data = vars["dmn_results"]["data"]  # size N x 10
dmn_surr = vars["dmn_results"]["surr"]  # size N x 10

# Calculate bins' boundaries
bins_n = 4
variables_n = length(DMNlist)
upper_boundaries = zeros(Float64, bins_n, variables_n)
lower_boundaries = zeros(Float64, bins_n, variables_n)
data_size = length(dmn_surr[:, 1])
# data_size = length(dmn_data[:, 1])
elements_per_bin = data_size / bins_n
for i=1:variables_n
    column = dmn_surr[:,i]
    # column = dmn_data[:,i]
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
# dmn_diskretized_data = zeros(Int64, size(dmn_data))
for i=1:length(DMNlist)
    for j=1:data_size
        for t=1:bins_n
            if lower_boundaries[t, i] <= dmn_surr[j, i] <= upper_boundaries[t, i]
            # if lower_boundaries[t, i] <= dmn_data[j, i] <= upper_boundaries[t, i]
                dmn_diskretized_surr[j, i] = t
                # dmn_diskretized_data[j, i] = t
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
    # row = dmn_diskretized_data[row_i, :]
    el_count[row...] += 1
end

for i=1:length(el_count)
    prob_table[i] = el_count[i]/data_size
end

# Calculate entropies
@show "Calculating entropies"
entropies = calculate_all_entropies(prob_table)
@show entropies 
# For 4 bins: 
#   Entropy for data: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] => 13.999890148403557
#   Entropy for surr: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] => 14.032503344201544
# For 6 bins:
#   Entropy for data: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] => 14.891790717456699
#   Entropy for surr: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] => 14.920331924191213


# Claculate Connected Information
@show "Calculating connected information"
connected_information = estimate_connected_information(10, distr_cards, entropies)
@show connected_information

# Connected information for data (4 bins)
# connected_information = [
#     NaN, 2.5681297562362424, 0.03634687464144193, 0.5756762091823333, 
#     0.01875476469014359, 0.3010863924518361, 0.0, 0.5546088181557671, 
#     0.556051386135417, 1.3894556501032618]

# Connected information for surrogate data (4 bins)
# connected_information = [
#     NaN, 2.5395243767783313, 0.026449042291297076, 0.580634812329663, 
#     0.014735425704273553, 0.29555065008044323, 0.0, 0.5443055604757809, 
#     0.5725793609165137, 1.3937174272221533]

# Connected information for data (6 bins)
# connected_information = [
#     NaN, 3.0150545111866762, 0.036558343057553344, 0.6967680665352631,
#     0.19541830869415833, 0.6017043277840166, 0.6146646815975956, 
#     1.8280439960335464, 1.5504137887368046, 2.4192082661292478]

# Connected information for surrogate data (6 bins)
# connected_information = [
#     NaN, 2.97281254152216, 0.03663785939826525, 0.6973971733843101,
#     0.17451342333183817, 0.5965618182873058, 0.6300059278948389,
#     1.8235096124054841, 1.5835037303086565, 2.4143509964874887]

# # Wrtie
# matwrite(
#     "resources/dmn_diskretized_surr_6.mat", 
#     # "resources/dmn_diskretized_data_6.mat",
#     Dict(
#         "dmn_diskretized_surr" => dmn_diskretized_surr,
#         # "dmn_diskretized_data" => dmn_diskretized_data,
#         "dmn_prob_table" => prob_table,
#         "dmn_connected_information" => connected_information))
