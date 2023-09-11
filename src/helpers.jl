module Helpers

export generate_prob_table, discretize_data

function generate_prob_table(dims...)
    rnd_n::Int64 = 1
    for dim in dims
        rnd_n *= dim
    end
    tmp_probs_v = zeros(Float64, rnd_n + 1)
    tmp_probs_v[rnd_n + 1] = 1
    tmp_probs_v[2:rnd_n] = rand(Float64, (rnd_n - 1))
    sort!(tmp_probs_v)
    probabilities = [tmp_probs_v[i] - tmp_probs_v[i - 1] for i=2:(rnd_n + 1)]
    return reshape(probabilities, dims)
end

function discretize_data(data, variables_n, bins_n = 2, serialize = false) 

    # Calculate bins' boundaries
    upper_boundaries = zeros(Float64, bins_n, variables_n)
    lower_boundaries = zeros(Float64, bins_n, variables_n)
    data_size = length(data[:, 1])
    elements_per_bin = data_size / bins_n
    for i=1:variables_n
        column = data[:,i]
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
    diskretized = zeros(Int16, size(data))
    for i=1:variables_n
        for j=1:data_size
            for t=1:bins_n
                if lower_boundaries[t, i] <= data[j, i] <= upper_boundaries[t, i]
                    diskretized[j, i] = t
                    break
                end
            end
        end
    end

    # Serialize
    if serialize
        matwrite("discretized_data.mat", Dict("diskretized" => diskretized))
    end
    return diskretized
end

end