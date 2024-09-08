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
            j = Int((j - current) // dims[idx])
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

# Program for computing the optimal quantum initial state for which LOSE strategy
# gives perfect discrimination of order guessing.  
function optimal_lose_strategy()
    # List of system permutation matrices M_π, where the order of systems is
    # S, A_I, B_I, C_I
    L = [
        permutesystems_matrix([4, 1, 2, 3], [2, 2, 2, 2]),
        permutesystems_matrix([3, 1, 4, 2], [2, 2, 2, 2]),
        permutesystems_matrix([4, 3, 1, 2], [2, 2, 2, 2]),
        permutesystems_matrix([2, 4, 1, 3], [2, 2, 2, 2]),
        permutesystems_matrix([3, 4, 2, 1], [2, 2, 2, 2]),
        permutesystems_matrix([2, 3, 4, 1], [2, 2, 2, 2])
    ]

    # List of M_π' * M_π
    LL = []
    for i=1:6
        for j=(i+1):6
            append!(LL, [L[j]'*L[i]])
        end
    end

    #Definition of quantum state
    sigma = ComplexVariable(2^4, 2^4)
    constraints = [sigma in :SDP]
    # If the state sigma exists then tr = 1 otherwise tr = 0
    constraints += [real(tr(sigma)) <= 1]

    for i=1:15
        constraints += [tr(sigma * LL[i]) == 0]
    end

    #Objective function
    f = real(tr(sigma))
    
    problem = maximize(f, constraints)
    solve!(
        problem,
        MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => 1e-8, "eps_rel" => 1e-8);
        silent_solver = true
    )
    return problem.optval, sigma.value
end

solution = optimal_lose_strategy()
# trace of sigma is equal 1, so there is a solution for perfect discrimination
print("Trace of the numerical state = $(solution[1]) \n")
# Numerical approximation of sigma
sigma_numerical = solution[2]
# Analitical version of sigma
tau0 = proj([1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
tau1 = proj([0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0])
tau2 = proj([0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0])
tau3 = proj([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0])
tau4 = proj([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])
sigma = 1//12*I(16) -( 1//15*(tau0 + tau4) + 1//60*(tau1 + tau3) + 1//90*tau2)
# Norm of the difference sigma - sigma_numerical
print("Distance between analitical and numerical solution = $(norm(sigma - sigma_numerical)) \n")

# List of system permutation matrices M_π, where the order of systems is
# S, A_I, B_I, C_I
L = [
    permutesystems_matrix([4, 1, 2, 3], [2, 2, 2, 2]),
    permutesystems_matrix([3, 1, 4, 2], [2, 2, 2, 2]),
    permutesystems_matrix([4, 3, 1, 2], [2, 2, 2, 2]),
    permutesystems_matrix([2, 4, 1, 3], [2, 2, 2, 2]),
    permutesystems_matrix([3, 4, 2, 1], [2, 2, 2, 2]),
    permutesystems_matrix([2, 3, 4, 1], [2, 2, 2, 2])
]
# List of M_π' * M_π
LL = []
for i=1:6
    for j=(i+1):6
        append!(LL, [L[j]'*L[i]])
    end
end
# Check if tr(M_π' * M_π) = 4
print("Number of matrices that has trace equal 4 is $(sum([tr(M) == 4 for M in LL])) out of $(length(LL))")