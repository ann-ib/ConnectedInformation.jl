module NSBHelpers

using Combinatorics
using MAT

export calc_nsb_entropies_single, calc_nsb_entropies_bulk

nsb_path = abspath("src/nsb_octave")

function _call_nsb_octave(data)::Float64
    print("Call with $data.\n")
    K = length(data)
    nx = data[data .> 0]
    kx = ones(Int16, length(nx));
    p=Pipe()
    run(pipeline(
        `octave --eval "addpath('$nsb_path'); 
                        S_nsb = find_nsb_entropy([$(join(kx, ' '))], [$(join(nx, ' '))], $K, 0.1, 1); 
                        disp(S_nsb)"`,
        stdout=p,
        stderr=devnull))
    close(p.in)
    S_nsb = parse(Float64, readchomp(pipeline(p, `tail -n 1`)))
    return S_nsb
end

function _serialize_data(disc_data_table, choosen_variables)
    if choosen_variables == []
        return []
    end
    variables_n = ndims(disc_data_table)
    if variables_n == length(choosen_variables)
        return disc_data_table
    end
    vars_to_drop = Tuple(setdiff((1:variables_n), choosen_variables))
    cur_table = sum(disc_data_table, dims=vars_to_drop)
    cur_table = dropdims(cur_table, dims=vars_to_drop)
    return vec(cur_table)
end

function calc_nsb_entropies_single(disc_data, bins_n, variables_n)::Dict{Vector{Int64}, Float64}
    entropies = Dict()
    entropies[[]] = 0
    data_size = size(disc_data, 1)

    # Calculate discretized data table 
    distr_cards = [bins_n for i=1:variables_n]
    disc_data_table = zeros(Int64, distr_cards...)
    for row_i=1:data_size
        row = disc_data[row_i, :]
        disc_data_table[row...] += 1
    end

    # Go through all combinations (2 ^ variables_n)
    mem = Dict()
    for choosen_variables in powerset(collect(1:variables_n), 1)  
        nsb_serialized_data = _serialize_data(disc_data_table, choosen_variables)
        if !haskey(mem, nsb_serialized_data)
            entropy = _call_nsb_octave(nsb_serialized_data)
            bin_entropy = entropy * log2(ℯ)
            mem[nsb_serialized_data] = bin_entropy
        end
        entropies[choosen_variables] = mem[nsb_serialized_data]
    end
    return entropies
end

function _call_nsb_octave_bulk(disc_data, bins_n, variables_n)::Dict{Vector{Int64}, Float64}
    in_filename = "tmp/nsb_discretized_data.mat"
    out_filename = "tmp/nsb_all_entropies.mat"

    matwrite(
        in_filename,
        Dict(
            "disc_data" => disc_data))

    print("Call bulk.\n")
    run(pipeline(
        `octave --eval "warning('off','all');
                        addpath('$nsb_path'); 
                        entropies = find_all_nsb_entropies('$in_filename', '$out_filename', $bins_n, $variables_n);"`,
        stdout=devnull))

    nsb_entropies = matread(out_filename)["entropies"]
    entropies = Dict()
    entropies[[]] = 0
    for (key, value) in nsb_entropies
        variables = map(x -> parse(Int64, x), split(key))
        entropies[variables] = value * log2(ℯ)
    end

    return entropies
end


function calc_nsb_entropies_bulk(disc_data, bins_n, variables_n)::Dict{Vector{Int64}, Float64}
    entropies = _call_nsb_octave_bulk(disc_data, bins_n, variables_n)
    return entropies
end
end