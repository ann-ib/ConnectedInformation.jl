module Entropy 

using Combinatorics

export caluclate_all_entropies, calculate_entropy

"""
    caluclate_all_entropies(probability_table::Array{Float64})

Return all possible joint entropies for distributions characterized by 
`probability_table`.

# Example

TBD
"""
function caluclate_all_entropies(probability_table::Array{Float64})::Dict{Vector{Int64}, Float64}
    entropies = Dict()
    distributions_n = ndims(probability_table)
    for choosen_distributions in powerset(collect(1:distributions_n))  
        cur_prob_table = _reduce_prob_table(probability_table, choosen_distributions)
        entropies[choosen_distributions] = calculate_entropy(cur_prob_table)
    end
    return entropies
end


"""
    calculate_entropy(probability_table::Array{Float64})

Calculates joint entropy of random variables for given `probability_table` 
(marginal distribution probabilities).

# Example

TBD
"""
function calculate_entropy(probability_table::Array{Float64})::Float64
    entropy = 0
    for p in probability_table
        if p == 0 
            continue
        end
        entropy -= p * log2(p)
    end
    return entropy
end 


"""
    _reduce_prob_table(probability_table::Array{Float64}, choosen_distributions::Vector{Int64})

Return probability table derived from given `probability_table` 
for `choosen_distributions`.

# Example

TBD
"""
function _reduce_prob_table(probability_table::Array{Float64}, choosen_distributions::Vector{Int64})::Array{Float64}
    distributions_n = ndims(probability_table)
    if distributions_n == size(choosen_distributions)[1]
        return probability_table
    end
    distributions_to_drop = Tuple(setdiff((1:distributions_n), choosen_distributions))
    cur_prob_table = sum(probability_table, dims=distributions_to_drop)
    cur_prob_table = dropdims(cur_prob_table, dims=distributions_to_drop)
    return cur_prob_table
end

end
