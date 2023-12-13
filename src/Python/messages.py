from node import Variable, Factor
from labeledArray import LabeledArray, name_to_axis_mapping, other_axes_from_labeled_axes
import numpy as np


class Messages:
    def __init__(self) -> None:
        self.messages = {}

    def _tile_to_shape_along_axis(self, arr, target_shape, target_axis):
        raw_axes = list(range(len(target_shape)))
        tile_dimensions = [target_shape[a] for a in raw_axes if a != target_axis]
        if len(arr.shape) == 0:
            tile_dimensions += [target_shape[target_axis]]
        elif len(arr.shape) == 1:
            assert arr.shape[0] == target_shape[target_axis]
            tile_dimensions += [1]
        else:
            raise NotImplementedError()
        tiled = np.tile(arr, tile_dimensions)

        shifted_axes = raw_axes[:target_axis] + [raw_axes[-1]] + raw_axes[target_axis:-1]
        transposed = np.transpose(tiled, shifted_axes)
        return transposed
    
    def _tile_to_other_dist_along_axis_name(self, tiling_labeled_array, target_array):
        target_axis_label = tiling_labeled_array.axes_labels[0]
        return LabeledArray(
            self._tile_to_shape_along_axis(tiling_labeled_array.array, target_array.array.shape, name_to_axis_mapping(target_array)[target_axis_label]),
            axes_labels=target_array.axes_labels
        )

    def _variable_to_factor_messages(self, variable: Variable, factor: Factor):
        incoming_messages = [
            self.factor_to_variable_message(neighbour_factor, variable)
            for neighbour_factor in variable.neighbours
            if neighbour_factor.name != factor.name
        ]

        return np.prod(incoming_messages, axis=0)

    def _factor_to_variable_messages(self, factor: Factor, variable: Variable):
        # Compute the product
        factor_dist = np.copy(factor.data.array)
        for neighbour_variable in factor.neighbours:
            if neighbour_variable.name == variable.name:
                continue
            incoming_message = self.variable_to_factor_messages(neighbour_variable, factor)
            factor_dist *= self._tile_to_other_dist_along_axis_name(
                LabeledArray(incoming_message, [neighbour_variable.name]),
                factor.data
            ).array

        other_axes = other_axes_from_labeled_axes(factor.data, variable.name)
        return np.squeeze(np.sum(factor_dist, axis=other_axes))

    def marginal(self, variable: Variable):
        unnorm_p = np.prod([
            self.factor_to_variable_message(neighbour_factor, variable)
            for neighbour_factor in variable.neighbours
        ], axis=0)
        
        return unnorm_p/np.sum(unnorm_p)

    def variable_to_factor_messages(self, variable: Variable, factor: Factor):
        message_name = (variable.name, factor.name)
        if message_name not in self.messages:
            self.messages[message_name] = self._variable_to_factor_messages(variable, factor)
        return self.messages[message_name]

    def factor_to_variable_message(self, factor: Factor, variable: Variable):
        message_name = (factor.name, variable.name)
        if message_name not in self.messages:
            self.messages[message_name] = self._factor_to_variable_messages(factor, variable)
        return self.messages[message_name]
