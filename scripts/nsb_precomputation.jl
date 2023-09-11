using MAT

include("../src/ConnectedInformation.jl")
include("../src/helpers.jl")
include("../src/nsb_helpers.jl")

using Revise

import .ConnectedInformation
import .Helpers: discretize_data
import .NSBHelpers: calc_nsb_entropies

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

entropies = calc_nsb_entropies(discretized_data, bins_n, variables_n)

print(entropies)
