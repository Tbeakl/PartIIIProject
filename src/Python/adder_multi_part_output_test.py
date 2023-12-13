from node import Variable, Factor
from factorGraph import FG
from labeledArray import LabeledArray
from typing import Dict, List
from messages import Messages
import numpy as np

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
    "input_a_0": Variable("input_a_0"),
    "input_b_0": Variable("input_b_0"),
    "output_temp_0": Variable("output_temp_0"),
    "output_0": Variable("output_0"),
    "carry_out_0": Variable("carry_out_0"),

    "input_a_1": Variable("input_a_1"),
    "input_b_1": Variable("input_b_1"),
    "output_temp_1": Variable("output_temp_1"),
    "output_1": Variable("output_1"),
    "carry_out_1": Variable("carry_out_1"),
}

factors: Dict[str, Factor] = {
    "carry_in_dist": Factor("carry_in_dist"),
    "input_a_0_dist": Factor("input_a_0_dist"),
    "output_0_dist": Factor("output_0_dist"),
    "f_add_0": Factor("f_add_0"),
    "f_add_carry_out_0": Factor("f_add_carry_out_0"),
    "f_add_output_0": Factor("f_add_output_0"),
    "input_a_1_dist": Factor("input_a_1_dist"),
    "output_1_dist": Factor("output_1_dist"),
    "f_add_1": Factor("f_add_1"),
    "f_add_carry_out_1": Factor("f_add_carry_out_1"),
    "f_add_output_1": Factor("f_add_output_1"),
}

NUMBER_OF_BITS = 4

variables["carry_in"].add_neighbour(factors["carry_in_dist"])
variables["carry_in"].add_neighbour(factors["f_add_0"])

variables["input_a_0"].add_neighbour(factors["input_a_0_dist"])
variables["input_a_0"].add_neighbour(factors["f_add_0"])

variables["input_b_0"].add_neighbour(factors["f_add_0"])

variables["output_temp_0"].add_neighbour(factors["f_add_0"])
variables["output_temp_0"].add_neighbour(factors["f_add_carry_out_0"])
variables["output_temp_0"].add_neighbour(factors["f_add_output_0"])

variables["output_0"].add_neighbour(factors["f_add_output_0"])
variables["output_0"].add_neighbour(factors["output_0_dist"])
variables["carry_out_0"].add_neighbour(factors["f_add_carry_out_0"])

factors["carry_in_dist"].add_neighbour(variables["carry_in"])
factors["input_a_0_dist"].add_neighbour(variables["input_a_0"])
factors["output_0_dist"].add_neighbour(variables["output_0"])

factors["f_add_0"].add_neighbour(variables["carry_in"])
factors["f_add_0"].add_neighbour(variables["input_a_0"])
factors["f_add_0"].add_neighbour(variables["input_b_0"])
factors["f_add_0"].add_neighbour(variables["output_temp_0"])

factors["f_add_carry_out_0"].add_neighbour(variables["output_temp_0"])
factors["f_add_carry_out_0"].add_neighbour(variables["carry_out_0"])
factors["f_add_output_0"].add_neighbour(variables["output_temp_0"])
factors["f_add_output_0"].add_neighbour(variables["output_0"])


variables["carry_out_0"].add_neighbour(factors["f_add_1"])

variables["input_a_1"].add_neighbour(factors["input_a_1_dist"])
variables["input_a_1"].add_neighbour(factors["f_add_1"])

variables["input_b_1"].add_neighbour(factors["f_add_1"])

variables["output_temp_1"].add_neighbour(factors["f_add_1"])
variables["output_temp_1"].add_neighbour(factors["f_add_carry_out_1"])
variables["output_temp_1"].add_neighbour(factors["f_add_output_1"])

variables["output_1"].add_neighbour(factors["f_add_output_1"])
variables["output_1"].add_neighbour(factors["output_1_dist"])
variables["carry_out_1"].add_neighbour(factors["f_add_carry_out_1"])

factors["input_a_1_dist"].add_neighbour(variables["input_a_1"])
factors["output_1_dist"].add_neighbour(variables["output_1"])

factors["f_add_1"].add_neighbour(variables["carry_out_0"])
factors["f_add_1"].add_neighbour(variables["input_a_1"])
factors["f_add_1"].add_neighbour(variables["input_b_1"])
factors["f_add_1"].add_neighbour(variables["output_temp_1"])

factors["f_add_carry_out_1"].add_neighbour(variables["output_temp_1"])
factors["f_add_carry_out_1"].add_neighbour(variables["carry_out_1"])
factors["f_add_output_1"].add_neighbour(variables["output_temp_1"])
factors["f_add_output_1"].add_neighbour(variables["output_1"])


p_carry_in_array = np.zeros((2, 1))
p_input_a_0_array = np.zeros((1 << NUMBER_OF_BITS, 1))
p_output_0_array = np.zeros((1 << NUMBER_OF_BITS, 1))
p_input_a_1_array = np.zeros((1 << NUMBER_OF_BITS, 1))
p_output_1_array = np.zeros((1 << NUMBER_OF_BITS, 1))


p_carry_in_array[0] = 1.0
p_input_a_0_array[5] = 1.0
p_input_a_1_array[0] = 1.0

p_output_0_array[12] = 1.0
p_output_1_array[1] = 1.0


p_carry_in = LabeledArray(p_carry_in_array, ["carry_in"])
p_input_a_0 = LabeledArray(p_input_a_0_array, ["input_a_0"])
p_output_0 = LabeledArray(p_output_0_array, ["output_0"])
p_input_a_1 = LabeledArray(p_input_a_1_array, ["input_a_1"])
p_output_1 = LabeledArray(p_output_1_array, ["output_1"])

p_add_0 = LabeledArray(
    make_add_including_carry_prob_array(NUMBER_OF_BITS),
    ["carry_in", "input_a_0", "input_b_0", "output_temp_0"],
)
p_add_output_0 = LabeledArray(
    make_bottom_bits_prob_array(NUMBER_OF_BITS), ["output_temp_0", "output_0"]
)
p_add_carry_0 = LabeledArray(
    make_top_bit_prob_array(NUMBER_OF_BITS), ["output_temp_0", "carry_out_0"]
)

p_add_1 = LabeledArray(
    make_add_including_carry_prob_array(NUMBER_OF_BITS),
    ["carry_out_0", "input_a_1", "input_b_1", "output_temp_1"],
)
p_add_output_1 = LabeledArray(
    make_bottom_bits_prob_array(NUMBER_OF_BITS), ["output_temp_1", "output_1"]
)
p_add_carry_1 = LabeledArray(
    make_top_bit_prob_array(NUMBER_OF_BITS), ["output_temp_1", "carry_out_1"]
)

fg = FG(factors=factors, variables=variables)

fg.set_data(
    {
        "carry_in_dist": p_carry_in,
        "input_a_0_dist": p_input_a_0,
        "output_0_dist": p_output_0,
        "f_add_0": p_add_0,
        "f_add_output_0": p_add_output_0,
        "f_add_carry_out_0": p_add_carry_0,
        "input_a_1_dist": p_input_a_1,
        "output_1_dist": p_output_1,
        "f_add_1": p_add_1,
        "f_add_output_1": p_add_output_1,
        "f_add_carry_out_1": p_add_carry_1,
    }
)

m = Messages()
print(f"Input B 0 : {m.marginal(fg.variable_from_name('input_b_0'))}")
print(f"Carry out 0: {m.marginal(fg.variable_from_name('carry_out_0'))}")

print(f"Input B 1 : {m.marginal(fg.variable_from_name('input_b_1'))}")
print(f"Carry out 1: {m.marginal(fg.variable_from_name('carry_out_1'))}")
