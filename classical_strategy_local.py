import itertools

# Program for computing the optimal LOSR strategy. The output is a tuple containing 
# the maximal number of distinct tuples for optimal strategy divided by six, giving us
# the probability of guessing order and the form of optimal strategy for each party. 


def optimal_losr_strategy():
    # List of extremal operations given as a tuple, e.g. (a^m(0), a^m(1))
    avaliable_operations =  [(0, 0), (0, 1), (1, 0), (1, 1)]

    max_number_of_distinct_tuples = 0
    optimal_strategy = []
    
    for strategy in itertools.product(avaliable_operations, repeat=3):
        # List of outputs for a given strategy interatted over permutations S_3
        outputs = set()

        for pi in itertools.permutations([0, 1, 2]):
            # Position in output_pi: Alice = 0, Bob = 1, Charlie = 2, System = 3
            output_pi = [0, 0, 0, 0]
            current_state = 0
            
            # Composition pi(a^m, b^m, c^m)(0)
            for i in range(3):
                output_pi[pi[i]] = current_state
                current_state = strategy[pi[i]][current_state]
            output_pi[3] = current_state

            outputs.add(tuple(output_pi))

        number_of_distinct_tuples = len(outputs)
        if number_of_distinct_tuples > max_number_of_distinct_tuples:
            max_number_of_distinct_tuples = number_of_distinct_tuples
            optimal_strategy = strategy
    return max_number_of_distinct_tuples/6, optimal_strategy

solved = optimal_losr_strategy()
print(f"The probability of order guessing: {solved[0]}")
print(f"Optimal strategy: a^m={solved[1][0]}, b^m={solved[1][1]}, c^m={solved[1][2]}")

