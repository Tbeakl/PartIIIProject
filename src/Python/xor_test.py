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

p_a = LabeledArray(np.array([[1.], [0.]]), ["a"])
p_b = LabeledArray(np.array([[1.], [0.]]), ["b"])
# p_c = LabeledArray(np.array([[0.9], [0.1]]), ["c"])

p_xor = LabeledArray(np.array([[[1.,0.],[0.,1.]],[[0.,1.],[1.,0.]]]), ["a", "b","c"])

fg = FG(factors=factors, variables=variables)

fg.set_data({"f1": p_a,
             "f2": p_b,
             "f_xor": p_xor})

m = Messages()
print(m.marginal(fg.variable_from_name('c')))