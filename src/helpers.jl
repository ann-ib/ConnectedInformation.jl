module Helpers

export generate_prob_table

function generate_prob_table(dims...)
    rnd_n = 1
    for dim in dims
        rnd_n *= dim
    end
    tmp_probs_v = zeros(rnd_n + 1)
    tmp_probs_v[rnd_n + 1] = 1
    tmp_probs_v[2:rnd_n] = rand(Float16, (rnd_n - 1))
    sort!(tmp_probs_v)
    probabilities = [tmp_probs_v[i] - tmp_probs_v[i - 1] for i=2:(rnd_n + 1)]
    return reshape(probabilities, dims)
end

end