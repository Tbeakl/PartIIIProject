# Fragment-template power-analysis attacks against microcontroller implementations of the 32-bit stream cipher ChaCha

This is a repository containing the source code for performing a fragment-template attack against different implementations of ChaCha.

## C source code

This contains the implementations of ChaCha which are attacked when run on a microcontroller.  They are designed to be run with a ChipWhisperer and are based around the simple serial base provided.  When run for the actual attacks they were compiled with the gcc compiler option `-Os`.

All of the implementations accept the same commands:

- `k` - Sets the key of ChaCha
- `n` - Sets the nonce of ChaCha
- `c` - Sets the counter of ChaCha
- `p` - Sets the plaintext to encrypt
- `x` - Reset the cipher (clear the plaintext, ciphertext and state)
- `i` - Output the input state
- `m` - Output the plaintext
- `e` - Perform the encryption
- `f` - Output the final state
- `o` - Output the ciphertext

The different source files refer to different implementations which were attacked

- `chacha-8bit-on-32bit.c` - An implementation of ChaCha which only uses 8 bit values
- `chacha-32bit.c` - An implementation of ChaCha which keeps most of the state in registers and uses 32-bit values, taken from [Lightweight Cryptography Primitives](https://rweather.github.io/lightweight-crypto/index.html)
- `chacha-32bit-volatile.c` - A modified version of the normal 32-bit implementation which has its state marked as volatile so every operation needs to access SRAM.

## Python source code

This consists of two main parts `initial_belief_propagation` which contains some initial experimentation with belief propagation on trees.  `ChipWhisperer` contains the code required for controlling the ChipWhisperer, oscilloscope and clock generator to make the power recordings.  `make_full_set_of_recordings.py` can be used to generate a full set of recordings for performing an attack.

## Julia source code

This makes up the bulk of the code which performs the processing steps for the attacks including in simulation.

### Belief propagation

This contains the basic code for belief propagation:

- `node.jl` - Contains the definition of the structures of types which make up the factor graph and nodes.
- `messages.jl` - This contains the basic code for performing belief propagation around a factor graph the main functions of interest are:
    - `variable_to_factor_messages` - takes an input variable and outputs updated messages to its connected factors, the messages are damped by the damping value.
    - `factor_to_variable_messages` - takes an input factor and outputs updated messages to its connected variables, which are damped by the damping value.  There are several implementations which are chosen between by multiple dispatch at runtime depending on the type of factor allowing for specialised calculations.
    - `marginal` - calculates the marginal distribution of a given variable.
- `dynamic_message_scheduling.jl` - performs a greedy dynamic approach for trying to reduce entropy in the factor graph

### ChaCha factor graph

This contains the code for creating the specific factor graph of ChaCha and adding leakages to it alongside being able to visualise the flow of information and quality of values throughout the graph.

- `chacha_factor_graph.jl` - This creates the factor graph for an execution of ChaCha in particular the `chacha_factor_graph!` function.
- `add_leakage_to_graph.jl` - Adds leakage probabilities to the graph when provided with a function for creating probabilities.
- `heatmap_visualisation.jl` - Allows the creation of a matrix representing the names of all the variables and provides functionality for creating a similar similarly sized matrix with either the current entropy of the corresponding variables or the change in entropy.
- `heatmap_visualisation_intermediate_values.jl` - Allows the creation of the same heatmap as the entropy but when the values relate to intermediate values of the execution such as the logarithmic guessing entropy, first order success rate or number of interesting cycles.

### Encryption

This contains the functions for performing ChaCha, outputting the intermediate values from the execution, converting hamming weights to probabilities, and for performing key enumeration (or estimating the amount of key enumeration required).

- `chacha.jl` - Contains the code for a normal implementation of ChaCha in Julia.
- `leakage_functions.jl` - Contains the code for leaking out the intermediate values and the associated functions for generating leakages.
- `hamming_weight_probability_calculations.jl` - Contains code for converting a leaked hamming weight into a probability distribution over values.
- `key_enumeration.jl` - Contains the code for enumerating through the key candidates in likelihood order.
- `rank_estimation.jl` - Contains the code for estimating the number of key candidates which need to be enumerated to find the correct solution.

### Evaluation

This contains the code for generating the figures which are included inside of the report

### Plotting traces

The files plot out raw power traces.

### Profiling

The code here is concerned with creating and evaluating the quality of the templates.

The basic steps when attempting to make new traces are:

- `make_intermediate_values_traces.jl` - This creates matrices of the intermediate values of the power traces which are used in the training stage.
- `make_mean_trace.jl` - This creates a mean trace which all future traces can be aligned against so that slight differences in trigger time do not affect the results much.
- `test_all_other_traces_highly_correlated.jl` - This checks that all of the traces are correlated to the mean trace so are of correct executions.
- `plot_correlations_to_find_clock_signal.jl` - This plots various statistical relationships of the traces to the intermediate values to find how the clock signal falls in the alignment.
- `plot_correlations_to_find_relevant_clock_cycles.jl` - This plots and outputs the correlation with a linear model on a per clock cycle basis.
- `make_bitmaps_of_relevant_cycles.jl` - This creates bitmaps of where the correlation with a linear model is above a certain threshold.
- `make_downsampled_matrix.jl` - This creates a matrix of downsampled traces for faster loading and computations in the future.
- `make_templates.jl` - This creates the templates from the downsampled traces.

The other files provide the qualities of different sets of templates as well as providing statistics about their performance.

### Tests

These contain various different tests for checking different parts of the program, such as that ChaCha has been correctly implemented, the factor graph correctly implements ChaCha, that belief propagation works around various known structures and that the success rates of simulated templates are close to their actual success rates.

### Attacks

This contains various different versions where SASCA is performed in different scenarios when given different types of leakage (whether simulated or real).  They all have a similar structure to each other where the `*_experimentation.jl` file performs SASCA in a specific instance with the variables such as `number_of_bits` controlling how the algorithm executes, the other files associated with each are either for providing functions for inputting the probability distributions into the factor graph or performing a set of attacks against lots of traces.