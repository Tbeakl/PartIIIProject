from node import Variable, Factor
from factorGraph import FG
from labeledArray import LabeledArray
from typing import Dict, List
from messages import Messages
import numpy as np

variables: Dict[str, Variable] = {
    "a": Variable("a"),
    "b": Variable("b"),
    "c": Variable("c"),
}

factors: Dict[str, Factor] = {
    "f1": Factor("f1"),
    "f2": Factor("f2"),
    "f_xor": Factor("f_xor"),
}

NUMBER_OF_BITS = 4

def make_xor_prob_array(num_of_bits: int):
    output = np.zeros((1 << num_of_bits, 1 << num_of_bits, 1 << num_of_bits))
    for input_a in range(1 << num_of_bits):
        for input_b in range(1 << num_of_bits):
            for output_val in range(1 << num_of_bits):
                output[input_a, input_b, output_val] = ((input_a ^ input_b) == output_val)
    return output.astype('float64')

variables["a"].add_neighbour(factors["f1"])
variables["a"].add_neighbour(factors["f_xor"])

variables["b"].add_neighbour(factors["f2"])
variables["b"].add_neighbour(factors["f_xor"])

variables["c"].add_neighbour(factors["f_xor"])

factors["f1"].add_neighbour(variables["a"])

factors["f2"].add_neighbour(variables["b"])

factors["f_xor"].add_neighbour(variables["a"])
factors["f_xor"].add_neighbour(variables["b"])
factors["f_xor"].add_neighbour(variables["c"])

p_a_array = np.zeros((1 << NUMBER_OF_BITS, 1))
p_b_array = np.zeros((1 << NUMBER_OF_BITS, 1))

p_a_array[3] = 1.
p_b_array[10] = 1.

p_a = LabeledArray(p_a_array, ["a"])
p_b = LabeledArray(p_b_array, ["b"])

p_xor = LabeledArray(make_xor_prob_array(NUMBER_OF_BITS), ["a", "b","c"])

fg = FG(factors=factors, variables=variables)

fg.set_data({"f1": p_a,
             "f2": p_b,
             "f_xor": p_xor})

m = Messages()
print(m.marginal(fg.variable_from_name('c')))