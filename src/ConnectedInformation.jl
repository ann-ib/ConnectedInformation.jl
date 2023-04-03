module ConnectedInformation

using JuMP
using MosekTools
using LinearAlgebra
using Combinatorics

export estimate_connected_information, estimate_max_enthropy, calculate_entropy, caluclate_all_entropies


"""
    estimate_connected_information(to_order::Int64, 
                                   distrs_card::Vector{Int64}, 
                                   entropy_constraints::Dict{Vector{Int64}, Float64})

Return vector with connected information approximation of order from 2 
up to `to_order` for given entropies `entropy_constraints` characterizing
distributions with given `distrs_card` cardinalities.
``

# Example

TBD
"""
function estimate_connected_information(to_order::Int64, 
                                        distrs_card::Vector{Int64}, 
                                        entropy_constraints::Dict{Vector{Int64}, Float64})::Vector{Float64}
    # TBD: make to_order argument optional, calculate it by default from marginals dict

    con_inf = Vector{Float64}(undef, to_order)
    entropies = Vector{Float64}(undef, to_order)
    con_inf[1] = NaN

    # Calculate maximum entropy approximation consistent with marginals of order 1 
    cur_constraints = filter(p -> (length(first(p))) <= 1, entropy_constraints)
    entropies[1] = estimate_max_enthropy(1, distrs_card, cur_constraints)

    # Calculate maximum entropy and connected informatoin approximations up to |to_order| order
    for i in 2:to_order
        cur_constraints = filter(p -> (length(first(p))) <= i, entropy_constraints)
        entropies[i] = estimate_max_enthropy(i, distrs_card, cur_constraints)
        # Approximation is not precise, entropy with more constraints could have smaller value
        if entropies[i] > entropies[i - 1]
            con_inf[i] = 0
        else 
            con_inf[i] = entropies[i - 1] - entropies[i]
        end
    end

    return con_inf
end


"""
    estimate_max_enthropy(k::Int64, distrs_card::Vector{<:Int64}, 
                          entropy_constraints::Dict{Vector{Int64}, Float64})

Return approximation of maximum entrophy for given entropies 
`entropy_constraints` characterizing distributions with given 
`distrs_card` cardinalities.

# Example

TBD
"""
function estimate_max_enthropy(k::Int64, distrs_card::Vector{Int64}, 
                               entropy_constraints::Dict{Vector{Int64}, Float64})::Float64
    # ℎ(∅) = 0
    entropy_constraints[[]] = 0 
    
    distributions_n = length(distrs_card)
    subset_to_index = Dict()
    subset_to_index[[]] = 1

    model = Model(optimizer_with_attributes(Mosek.Optimizer, "QUIET" => true))

    # non-negativity constraints
    # ℎ(𝐴) ≥ 0, ∀𝐴 ∈ 𝒫(𝑁)
    @variable(model, h[1:(2^distributions_n + 1)] >= 0)  

    # ∀𝐴 ∈ 𝒫(𝑁)
    for subset_A in powerset(collect(1:distributions_n), 1)  
        subset_to_index[subset_A] = length(subset_to_index) + 1 
        subset_A_index = subset_to_index[subset_A]
    
        # lower bound constraints (part of ℎ ∈ Γₙ)
        # ℎ(𝐴) ≤ ℎ(𝐵), ∀𝐴,𝐵 ⊆ 𝑁 : 𝐴 ⊂ 𝐵
        subset_B_size = length(subset_A) - 1
        for subset_B in powerset(subset_A, subset_B_size, subset_B_size)  # subset of |subset_A| with one less element
            subset_B_index = subset_to_index[subset_B]
            @constraint(model, h[subset_B_index] <= h[subset_A_index])
        end

        # upper bound constraints (part of ℎ ∈ Γₙ)
        # ℎ(𝐴) ≤ ℎ(𝐵) + ℎ(𝐴/𝐵), ∀𝐵 ⊆ 𝐴
        subset_union = subset_A
        for subset_C in powerset(subset_union, ceil(Int64, length(subset_union) / 2), length(subset_union))  # subsets with at least half of |subset_A| elements
            subset_D = setdiff(subset_union, subset_C)
            subset_intersect = intersect(subset_C, subset_D)
            

            subset_C_index = subset_to_index[subset_C]
            subset_D_index = subset_to_index[subset_D]
            subset_union_index = subset_to_index[subset_union]
            subset_intersect_index = subset_to_index[subset_intersect]
            
            @constraint(model, h[subset_intersect_index] + h[subset_union_index] <= h[subset_C_index] + h[subset_D_index])
        end

        # given entropy constraints
        # ℎ(𝐴) = ℎₚ(𝐴), ∀𝐴 ∈ 𝒫ₖ(𝑁)
        if length(subset_A) <= k
            if haskey(entropy_constraints, subset_A)
                @constraint(model, h[subset_A_index] == entropy_constraints[subset_A])
            end
        end

        # upper bound constraints 
        # ℎ(𝐴) ≤ 𝑙𝑜𝑔₂∣𝒳ₐ∣, ∀𝐴 ∈ 𝒫(𝑁)∖𝒫ₖ(𝑁)
        if length(subset_A) > k
            cardinality = _calculate_cardinality(distrs_card, subset_A)
            @constraint(model, h[subset_A_index] <= log(2, cardinality))
        end
    end
    
    
    # Maximize 𝐻(𝑋₁, 𝑋₂, …, 𝑋ₙ)
    @objective(model, Max, h[subset_to_index[collect(1:distributions_n)]])  
    optimize!(model);
    return objective_value(model)
end


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


function _calculate_cardinality(distrs_card::Vector{<:Int64}, choosen_distributions::Vector{Int64})
    if length(choosen_distributions) == 0
        return 0
    end
    cardinality = 1
    for d in choosen_distributions
        cardinality *= distrs_card[d]
    end
    return cardinality
end

end
