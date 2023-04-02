module PrimalProblem

using JuMP
using Ipopt
using LinearAlgebra
using Combinatorics

export solve_primal_prioblem

"""
For 3 variables distribution, 3 possible value each   
"""
function solve_primal_prioblem()
    # Example
    p = fill(1/32, 3, 3, 3)
    p[3, 3, 1] = 1/16
    p[3, 1, 3] = 1/16
    p[2, 3, 1] = 1/16
    p[3, 1, 2] = 1/16
    p[1, 3, 2] = 1/16

    model = Model(optimizer_with_attributes(Ipopt.Optimizer, "QUIET" => true))

    register(model, :calculate_entropy1, 1, calculate_entropy; autodiff = true)
    register(model, :calculate_entropy2, 2, calculate_entropy; autodiff = true)

    @variable(model, q[1:27] >= 0)  # linear inequality constraints q(X) >= 0
    @constraint(model, sum(q) == 1)  # linear equality constraint ∑q(X) == 1
    for choosen_distributions in powerset([1,2,3], 1, 2)
        p_entropy = calculate_entropy(p, choosen_distributions)
        @NLconstraint(model, 
            calculate_entropy2(reshape(q, (3, 3, 3)), choosen_distributions) == p_entropy)
    end

    @NLobjective(model, Max, calculate_entropy1(q[1:27]))  
    optimize!(model);
    return objective_value(model)
end

function _reduce_prob_table(probability_table::Array{<:Real}, choosen_distributions::Vector{Int64})::Array{<:Real}
    distributions_n = ndims(probability_table)
    if distributions_n == size(choosen_distributions)[1]
        return probability_table
    end
    distributions_to_drop = Tuple(setdiff((1:distributions_n), choosen_distributions))
    cur_prob_table = sum(probability_table, dims=distributions_to_drop)
    cur_prob_table = dropdims(cur_prob_table, dims=distributions_to_drop)
    return cur_prob_table
end

function calculate_entropy(probability_table::Array{<:Real}, choosen_distributions::Vector{Int64})::Real
    cur_prob_table = _reduce_prob_table(probability_table, choosen_distributions)
    entropy = 0
    for p in cur_prob_table
        if p == 0 
            continue
        end
        entropy -= p * log2(p)
    end
    return entropy
end 

function calculate_entropy(probability_table::Array{<:Real})::Real
    entropy = 0
    for p in probability_table
        if p == 0 
            continue
        end
        entropy -= p * log2(p)
    end
    return entropy
end 

end