import numpy as np

def _tile_to_shape_along_axis(arr, target_shape, target_axis):
    # This should be relatively simple I think
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
    # print(tiled)
    # print(tile_dimensions)
    # print(tiled.shape)
    shifted_axes = raw_axes[:target_axis] + [raw_axes[-1]] + raw_axes[target_axis:-1]
    transposed = np.transpose(tiled, shifted_axes)
    return transposed

array = np.array([1,2,3,4])
target_shape = (4,3,2,1)
target_axis = 0
# _tile_to_shape_along_axis(array, target_shape, target_axis)
print(_tile_to_shape_along_axis(array, target_shape, target_axis))
print(_tile_to_shape_along_axis(array, target_shape, target_axis).shape)