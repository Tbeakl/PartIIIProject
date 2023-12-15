from node import Variable, Factor
from factorGraph import FG
from labeledArray import LabeledArray
from typing import Dict, List, Tuple
from messages import Messages
import numpy as np
import math

def make_add_including_carry_prob_array(num_of_bits: int):
    output = np.zeros((2, 1 << num_of_bits, 1 << num_of_bits, 1 << (num_of_bits + 1)))
    for carry_in in range(2):
        for input_a in range(1 << num_of_bits):
            for input_b in range(1 << num_of_bits):
                for output_val in range(1 << (num_of_bits + 1)):
                    output[carry_in, input_a, input_b, output_val] = (
                        input_a + input_b + carry_in
                    ) == output_val
    return output.astype("float64")


def make_top_bit_prob_array(num_of_bits: int):
    output = np.zeros(((1 << (num_of_bits + 1)), 2))
    for val_in in range(1 << (num_of_bits + 1)):
        for output_val in range(2):
            if (output_val == 1 and val_in >= 1 << num_of_bits) or (
                output_val == 0 and val_in < (1 << num_of_bits)
            ):
                output[val_in, output_val] = 1.0
    return output.astype("float64")


def make_bottom_bits_prob_array(num_of_bits: int):
    output = np.zeros(((1 << (num_of_bits + 1)), 1 << num_of_bits))
    for val_in in range(1 << (num_of_bits + 1)):
        for output_val in range(1 << num_of_bits):
            output[val_in, output_val] = val_in % (1 << num_of_bits) == output_val
    return output.astype("float64")


def make_32_bit_adder(number_of_bits_per_cluster: int) -> Tuple[Dict[str, Variable], Dict[str, Factor], Dict[str, LabeledArray]]:
    number_of_clusters: int = math.ceil(32 // number_of_bits_per_cluster)
    variables: Dict[str, Variable] = {}
    factors: Dict[str, Factor] = {}
    distributions: Dict[str, LabeledArray] = {}

    full_add_dist = make_add_including_carry_prob_array(number_of_bits_per_cluster)
    full_add_to_output = make_bottom_bits_prob_array(number_of_bits_per_cluster)
    full_add_carry = make_top_bit_prob_array(number_of_bits_per_cluster)

    # Create the base variables and factors for the addition
    for i in range(number_of_clusters):
        variables[f"carry_{i}"] = Variable(f"carry_{i}")
        variables[f"input_a_{i}"] = Variable(f"input_a_{i}")
        variables[f"input_b_{i}"] = Variable(f"input_b_{i}")
        variables[f"output_temp_{i}"] = Variable(f"output_temp_{i}")
        variables[f"output_{i}"] = Variable(f"output_{i}")

        factors[f"f_add_{i}"] = Factor(f"f_add_{i}")
        factors[f"f_add_output_{i}"] = Factor(f"f_add_output_{i}")
        factors[f"f_add_carry_{i}"] = Factor(f"f_add_carry_{i}")

    variables[f"carry_{number_of_clusters}"] = Variable(f"carry_{number_of_clusters}")

    # Create the factors and their distributions for the internals of the factor graph
    for i in range(number_of_clusters):
        variables[f"carry_{i}"].add_neighbour(factors[f"f_add_{i}"])
        variables[f"input_a_{i}"].add_neighbour(factors[f"f_add_{i}"])
        variables[f"input_b_{i}"].add_neighbour(factors[f"f_add_{i}"])
        variables[f"output_temp_{i}"].add_neighbour(factors[f"f_add_{i}"])

        factors[f"f_add_{i}"].add_neighbour(variables[f"carry_{i}"])
        factors[f"f_add_{i}"].add_neighbour(variables[f"input_a_{i}"])
        factors[f"f_add_{i}"].add_neighbour(variables[f"input_b_{i}"])
        factors[f"f_add_{i}"].add_neighbour(variables[f"output_temp_{i}"])

        variables[f"output_temp_{i}"].add_neighbour(factors[f"f_add_output_{i}"])
        variables[f"output_{i}"].add_neighbour(factors[f"f_add_output_{i}"])

        factors[f"f_add_output_{i}"].add_neighbour(variables[f"output_temp_{i}"])
        factors[f"f_add_output_{i}"].add_neighbour(variables[f"output_{i}"])

        variables[f"output_temp_{i}"].add_neighbour(factors[f"f_add_carry_{i}"])
        variables[f"carry_{i+1}"].add_neighbour(factors[f"f_add_carry_{i}"])

        factors[f"f_add_carry_{i}"].add_neighbour(variables[f"output_temp_{i}"])
        factors[f"f_add_carry_{i}"].add_neighbour(variables[f"carry_{i+1}"])

        distributions[f"f_add_{i}"] = LabeledArray(full_add_dist, [f"carry_{i}", f"input_a_{i}", f"input_b_{i}", f"output_temp_{i}"])
        distributions[f"f_add_output_{i}"] = LabeledArray(full_add_to_output, [f"output_temp_{i}", f"output_{i}"])
        distributions[f"f_add_carry_{i}"] = LabeledArray(full_add_carry, [f"output_temp_{i}", f"carry_{i+1}"])

    return variables, factors, distributions

def add_base_probabilities(number_of_bits_per_cluster: int, variables: Dict[str, Variable], factors: Dict[str, Factor], distributions: Dict[str, LabeledArray]):
    number_of_clusters: int = math.ceil(32 // number_of_bits_per_cluster)

    p_carry_in_array = np.zeros((2, 1))
    p_carry_in_array[0] = 1.
    factors["f_carry_0"] = Factor("f_carry_0")
    variables["carry_0"].add_neighbour(factors["f_carry_0"])
    factors["f_carry_0"].add_neighbour(variables["carry_0"])
    distributions["f_carry_0"] = LabeledArray(p_carry_in_array, ["carry_0"])

    # Put in a distribution for input A, just setting the first, second and fourth to 1
    for i in range(number_of_clusters):
        factors[f"f_input_a_{i}"] = Factor(f"f_input_a_{i}")
        variables[f"input_a_{i}"].add_neighbour(factors[f"f_input_a_{i}"])
        factors[f"f_input_a_{i}"].add_neighbour(variables[f"input_a_{i}"])
        p_input_array = np.zeros((1 << number_of_bits_per_cluster, 1))
        if i in [0,1,3]:
            p_input_array[1] = 1.
        else:
            p_input_array[0] = 1.
        distributions[f"f_input_a_{i}"] = LabeledArray(p_input_array, [f"input_a_{i}"])

    # Put in distribution for input B, just set to all 0
    for i in range(number_of_clusters):
        factors[f"f_input_b_{i}"] = Factor(f"f_input_b_{i}")
        variables[f"input_b_{i}"].add_neighbour(factors[f"f_input_b_{i}"])
        factors[f"f_input_b_{i}"].add_neighbour(variables[f"input_b_{i}"])
        p_input_array = np.zeros((1 << number_of_bits_per_cluster, 1))
        if i in []:
            p_input_array[1] = 1.
        else:
            p_input_array[0] = 1.
        distributions[f"f_input_b_{i}"] = LabeledArray(p_input_array, [f"input_b_{i}"])

NUMBER_OF_BITS = 2

variables, factors, distributions = make_32_bit_adder(NUMBER_OF_BITS)
add_base_probabilities(NUMBER_OF_BITS, variables=variables, factors=factors, distributions=distributions)

fg = FG(factors=factors, variables=variables)

fg.set_data(
    distributions
)

m = Messages()

for i in range(math.ceil(32 // NUMBER_OF_BITS)):
    print(f"Output {i}: {m.marginal(fg.variable_from_name(f'output_{i}'))}")