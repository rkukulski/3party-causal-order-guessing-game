using Convex, SCS
using LinearAlgebra
using QuantumInformation
using Combinatorics
using SparseArrays
const MOI = Convex.MOI


# Function returning a matrix for permuting elements of 
# vector defined on dims according to permutation
function permutesystems_matrix(permutation, dims)
    M = sparse(zeros( prod(dims), prod(dims)))

    for i = 1:prod(dims)
        temp = zeros(Int, length(dims))
        j = i - 1
        for idx = length(dims):-1:1
            current = mod(j, dims[idx])
            temp[idx] = current
            j = Int((j - current) / dims[idx])
        end
        j = temp[permutation[1]]
        for idx = 2:length(dims)
            j = j * dims[permutation[idx]] + temp[permutation[idx]]
        end
        j = j+1
        M[j, i] = 1
    end
    return M
end

# Function permuting subsystems of vector x defined on systems dims according to permutation
function permute_for_vectors(x, permutation, dims)
    M = permutesystems_matrix(permutation, dims)
    return M*x
end

# Function calculating partial trace over system for vector x defined on systems dims
function partialtrace_for_vectors(x, system, dims)
    y = permute_for_vectors(x, [2, 1, 3], [prod(dims[1:(system-1)]), dims[system], prod(dims[(system+1):length(dims)])])
    chunk = div(prod(dims), dims[system])
    return sum(y[(chunk*i + 1):(chunk*(i+1))] for i=0:(dims[system]-1))
end

# Program for computing the optimal probability of causal order discrimination using
# classical non-signaling strategies.
function optimal_non_signaling_classical_strategy_sdp()
    # Instead of π ∈ S_3 we use variable j ∈ {1,...,6}

    # List of variables D_j, which can be represented by their diagonal elements
    listD = [Variable(256) for j=1:6]

    # Positivity constraint
    constraints = [element >= 0 for element in listD[1]]
    for j in 2:6
        constraints += [element >= 0 for element in listD[j]]
    end

    # System order for Vectors D_j reads: 
    # S_P, S_I^A, S_O^A, S_I^B, S_O^B, S_I^C, S_O^C, S_F
    sumD = sum(listD)
    
    # Linear constraints
    # 1.
    constraints += [sumD == kron(partialtrace_for_vectors(sumD, 2, [128, 2]), [1/2, 1/2])]
    
    # System order for reduced sumD reads: 
    # S_P, S_I^A, S_O^A, S_I^B, S_O^B, S_I^C, S_O^C
    sumD = partialtrace_for_vectors(sumD, 2, [128, 2])
    # 2.
    constraints += [partialtrace_for_vectors(sumD, 2, [4, 2, 16]) == 
    permute_for_vectors(kron(partialtrace_for_vectors(sumD, 2, [2, 4, 16]), [1/2, 1/2]), [1, 3, 2], [2, 16, 2])]

    # # 3.
    constraints += [partialtrace_for_vectors(sumD, 2, [16, 2, 4]) == 
    permute_for_vectors(kron(partialtrace_for_vectors(sumD, 2, [8, 4, 4]), [1/2, 1/2]), [1, 3, 2], [8, 4, 2])]

    # # 4.
    constraints += [partialtrace_for_vectors(sumD, 2, [64, 2]) == 
    kron(partialtrace_for_vectors(sumD, 2, [32, 4]), [1/2, 1/2])]

    # 5.
    constraints += [sum(sumD) == 16]

    # Objective function f
    f = 0
    list_of_permutation_process_matrices = [perm for perm in permutations(1:3)]
    vec_id = [1, 0, 0, 1]
    temp_W_pi = kron(vec_id, vec_id, vec_id, vec_id)
    for j=1:6
        # reminder of the order: 
        # S_P, S_I^A, S_O^A, S_I^B, S_O^B, S_I^C, S_O^C, S_F
        permj = list_of_permutation_process_matrices[j]
        permutation = [1, 2, 3, 4, 5, 6, 7, 8]
        permutation[2*permj[1]] = 2
        permutation[2*permj[2]] = permutation[2*permj[1]+1]+1
        permutation[2*permj[3]] = permutation[2*permj[2]+1]+1
        permutation[8] = permutation[2*permj[3]+1]+1
        Wj = permute_for_vectors(temp_W_pi, permutation, fill(2, 8))
        f += sum(Wj .* listD[j])
    end
    f = f/6
    
    # SDP problem
    problem = maximize(f, constraints)
    solve!(
        problem,
        MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => 1e-8, "eps_rel" => 1e-8);
        silent_solver = true
    )
    
    if string(problem.status) == "OPTIMAL"
        return problem.optval
    else
        return 0
    end
end

optimal_non_signaling_classical_strategy_sdp()

