from collections import namedtuple

LabeledArray = namedtuple('LabeledArray', [
    'array',
    'axes_labels',
])


def name_to_axis_mapping(labeled_array):
    return {
        name: axis
        for axis, name in enumerate(labeled_array.axes_labels)
    }

def other_axes_from_labeled_axes(labeled_array, axis_label):
    # returns the indexes of the axes that are not axis label
    return tuple(
        axis
        for axis, name in enumerate(labeled_array.axes_labels)
        if name != axis_label
    )