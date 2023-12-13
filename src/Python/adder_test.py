from node import Variable, Factor
from factorGraph import FG
from labeledArray import LabeledArray
from typing import Dict, List
from messages import Messages
import numpy as np


def make_add_mod_prob_array(num_of_bits: int):
    output = np.zeros((2, 1 << num_of_bits, 1 << num_of_bits, 1 << num_of_bits))
    for carry_in in range(2):
        for input_a in range(1 << num_of_bits):
            for input_b in range(1 << num_of_bits):
                for output_val in range(1 << num_of_bits):
                    output[carry_in, input_a, input_b, output_val] = (
                        input_a + input_b + carry_in
                    ) % (1 << num_of_bits) == output_val
    return output.astype("float64")


def make_carry_out_mod_prob_array(num_of_bits: int):
    output = np.zeros((2, 1 << num_of_bits, 1 << num_of_bits, 2))
    MAX_VALUE = (1 << num_of_bits) - 1
    for carry_in in range(2):
        for input_a in range(1 << num_of_bits):
            for input_b in range(1 << num_of_bits):
                output[carry_in, input_a, input_b, 1] = (
                    input_a + input_b + carry_in
                ) > MAX_VALUE
    return output.astype("float64")


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


variables: Dict[str, Variable] = {
    "carry_in": Variable("carry_in"),
    "input_a": Variable("input_a"),
    "input_b": Variable("input_b"),
    "output_temp": Variable("output_temp"),
    "output": Variable("output"),
    "carry_out": Variable("carry_out"),
}

factors: Dict[str, Factor] = {
    "carry_in_dist": Factor("carry_in_dist"),
    "input_a_dist": Factor("input_a_dist"),
    "input_b_dist": Factor("input_b_dist"),
    "f_add": Factor("f_add"),
    "f_add_carry_out": Factor("f_add_carry_out"),
    "f_add_output": Factor("f_add_output"),
}

NUMBER_OF_BITS = 4

variables["carry_in"].add_neighbour(factors["carry_in_dist"])
variables["carry_in"].add_neighbour(factors["f_add"])

variables["input_a"].add_neighbour(factors["input_a_dist"])
variables["input_a"].add_neighbour(factors["f_add"])

variables["input_b"].add_neighbour(factors["input_b_dist"])
variables["input_b"].add_neighbour(factors["f_add"])

variables["output_temp"].add_neighbour(factors["f_add"])
variables["output_temp"].add_neighbour(factors["f_add_carry_out"])
variables["output_temp"].add_neighbour(factors["f_add_output"])

variables["output"].add_neighbour(factors["f_add_output"])
variables["carry_out"].add_neighbour(factors["f_add_carry_out"])

factors["carry_in_dist"].add_neighbour(variables["carry_in"])
factors["input_a_dist"].add_neighbour(variables["input_a"])
factors["input_b_dist"].add_neighbour(variables["input_b"])

factors["f_add"].add_neighbour(variables["carry_in"])
factors["f_add"].add_neighbour(variables["input_a"])
factors["f_add"].add_neighbour(variables["input_b"])
factors["f_add"].add_neighbour(variables["output_temp"])

factors["f_add_carry_out"].add_neighbour(variables["output_temp"])
factors["f_add_carry_out"].add_neighbour(variables["carry_out"])
factors["f_add_output"].add_neighbour(variables["output_temp"])
factors["f_add_output"].add_neighbour(variables["output"])


p_carry_in_array = np.zeros((2, 1))
p_input_a_array = np.zeros((1 << NUMBER_OF_BITS, 1))
p_input_b_array = np.zeros((1 << NUMBER_OF_BITS, 1))

p_carry_in_array[0] = 1.0
p_input_a_array[1] = 1.0
p_input_b_array[1] = 1.0

p_carry_in = LabeledArray(p_carry_in_array, ["carry_in"])
p_input_a = LabeledArray(p_input_a_array, ["input_a"])
p_input_b = LabeledArray(p_input_b_array, ["input_b"])

p_add = LabeledArray(
    make_add_including_carry_prob_array(NUMBER_OF_BITS),
    ["carry_in", "input_a", "input_b", "output_temp"],
)
p_add_output = LabeledArray(
    make_bottom_bits_prob_array(NUMBER_OF_BITS), ["output_temp", "output"]
)
p_add_carry = LabeledArray(
    make_top_bit_prob_array(NUMBER_OF_BITS), ["output_temp", "carry_out"]
)

fg = FG(factors=factors, variables=variables)

fg.set_data(
    {
        "carry_in_dist": p_carry_in,
        "input_a_dist": p_input_a,
        "input_b_dist": p_input_b,
        "f_add": p_add,
        "f_add_output": p_add_output,
        "f_add_carry_out": p_add_carry,
    }
)

m = Messages()
print(f"Output: {m.marginal(fg.variable_from_name('output'))}")
print(f"Carry out: {m.marginal(fg.variable_from_name('carry_out'))}")
