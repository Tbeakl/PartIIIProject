from labeledArray import LabeledArray
from typing import Dict, List
from node import Factor, Variable

class FG:
    def __init__(self, factors: Dict[str, Factor], variables: Dict[str, Variable]) -> None:
        self._factors: Dict[str, Factor] = factors
        self._variables: Dict[str, Variable] = variables

    def set_data(self, data: Dict[str, LabeledArray]):
        var_dims = {}
        for factor_name in self._factors:
            factor = self._factors[factor_name]
            factor_data = data[factor.name]

            if set(factor_data.axes_labels) != set(v.name for v in factor.neighbours):
                missing_axes = set(v.name for v in factor.neighbours) - set(factor_data.axes_labels)
                raise ValueError(f"data[{factor.name}] is missing axes: {missing_axes}")
            
            for var_name, dim in zip(factor_data.axes_labels, factor_data.array.shape):
                if var_name not in var_dims:
                    var_dims[var_name] = dim
                
                if var_dims[var_name] != dim:
                    raise ValueError(f"data[{factor.name}] axes is wrong size, {dim}. Expected {var_dims[var_name]}")
                
            factor.data = factor_data

    def variable_from_name(self, var_name: str) -> Variable:
        return self._variables[var_name]