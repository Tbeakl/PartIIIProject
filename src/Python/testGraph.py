from node import Variable, Factor
from factorGraph import FG
from labeledArray import LabeledArray
from typing import Dict, List
from messages import Messages
import numpy as np

variables: Dict[str, Variable] = {
    "h1": Variable("h1"),
    "h2": Variable("h2"),
    "v1": Variable("v1"),
    "v2": Variable("v2"),
}

factors: List[Factor] = [
    Factor("p(h1)"),
    Factor("p(h2|h1)"),
    Factor("p(v1|h1)"),
    Factor("p(v2|h2)"),
]

variables["h1"].add_neighbour(factors[0])
variables["h1"].add_neighbour(factors[1])
variables["h1"].add_neighbour(factors[2])

variables["h2"].add_neighbour(factors[1])
variables["h2"].add_neighbour(factors[3])

variables["v1"].add_neighbour(factors[2])

variables["v2"].add_neighbour(factors[3])

factors[0].add_neighbour(variables["h1"])

factors[1].add_neighbour(variables["h2"])
factors[1].add_neighbour(variables["h1"])

factors[2].add_neighbour(variables["v1"])
factors[2].add_neighbour(variables["h1"])

factors[3].add_neighbour(variables["v2"])
factors[3].add_neighbour(variables["h2"])

p_h1 = LabeledArray(np.array([[0.2], [0.8]]), ["h1"])
p_h2_given_h1 = LabeledArray(np.array([[0.5, 0.2], [0.5, 0.8]]), ["h2", "h1"])
p_v1_given_h1 = LabeledArray(np.array([[0.6, 0.1], [0.4, 0.9]]), ["v1", "h1"])
p_v2_given_h2 = LabeledArray(p_v1_given_h1.array, ["v2", "h2"])

fg = FG(factors=factors, variables=variables)

fg.set_data({"p(h1)": p_h1,
             "p(h2|h1)": p_h2_given_h1,
             "p(v1|h1)": p_v1_given_h1,
             "p(v2|h2)": p_v2_given_h2})

m = Messages()
print(m.marginal(fg.variable_from_name('v2')))