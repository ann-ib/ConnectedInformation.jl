using MAT
using Serialization

include("../src/ConnectedInformation.jl")
include("../src/helpers.jl")
include("../src/nsb_helpers.jl")

using Revise

import .ConnectedInformation: estimate_max_entropies, estimate_connected_information
import .Helpers: discretize_data
import .NSBHelpers: calc_nsb_entropies_single, calc_nsb_entropies_bulk

# Read data
vars = matread("resources/connectivities_estimates_withsurr_twoAALnets.mat")
DMNlist = vars["DMNlist"]
dmn_data = vars["dmn_results"]["data"]  # size N x 10
dmn_surr = vars["dmn_results"]["surr"]  # size N x 10

data_size = size(dmn_data, 1)
variables_n = length(DMNlist)
bins_n = 2
data = dmn_data

discretized_data = discretize_data(data, variables_n, bins_n, false)

nsb_entropies = calc_nsb_entropies_bulk(discretized_data, bins_n, variables_n)

# Claculate approximations
distr_cards = [bins_n for i=1:variables_n]
@show "Calculating max entropies"
max_entropies = estimate_max_entropies(
    variables_n, distr_cards, nsb_entropies)
@show "Calculating connected information"
connected_information = estimate_connected_information(
    variables_n, distr_cards, nsb_entropies)


# Wrtie
filename = "dmn_diskretized_nsb_bulk_$(bins_n)"

matwrite(
    "resources/$(filename).mat",
    Dict(
        "discretized_data" => discretized_data,
        "max_entropies" => max_entropies,
        "connected_information" => connected_information))

serialize(
    "resources/$(filename).ser",
    Dict(
        "discretized_data" => discretized_data,
        "nsb_entropies" => nsb_entropies,
        "max_entropies" => max_entropies,
        "connected_information" => connected_information))
