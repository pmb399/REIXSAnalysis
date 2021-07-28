import numpy as np
from scipy.interpolate import interp2d

def grid_data(mcpRSXS_scales, mcpRSXS_scatters, img, x_low=None, x_high=None, y_low=None, y_high=None, new_x=None, new_y=None):

    mcpRSXS_scale = mcpRSXS_scales[img]
    mcpRSXS_scatter = mcpRSXS_scatters[img]

    if x_low != None:
        if x_low < mcpRSXS_scale[0, :].min():
            xmin = mcpRSXS_scale[0, :].min()
        else:
            xmin = x_low
    else:
        xmin = mcpRSXS_scale[0, :].min()

    if x_high != None:
        if x_high < mcpRSXS_scale[0, :].max():
            xmax = x_high
        else:
            xmax = mcpRSXS_scale[0, :].max()
    else:
        xmax = mcpRSXS_scale[0, :].max()

    if y_low != None:
        if y_low < mcpRSXS_scale[1, :].min():
            ymin = mcpRSXS_scale[1, :].min()
        else:
            ymin = y_low
    else:
        ymin = mcpRSXS_scale[1, :].min()

    if y_high != None:
        if y_high < mcpRSXS_scale[1, :].max():
            ymax = y_high
        else:
            ymax = mcpRSXS_scale[1, :].max()
    else:
        ymax = mcpRSXS_scale[1, :].max()

    xpoints = int(
        np.abs(np.ceil((xmax-xmin)/np.abs(np.diff(mcpRSXS_scale[0, :])).min())))
    ypoints = int(
        np.abs(np.ceil((ymax-ymin)/np.abs(np.diff(mcpRSXS_scale[1, :])).min())))

    f = interp2d(mcpRSXS_scale[0], mcpRSXS_scale[1], mcpRSXS_scatter)

    if type(new_x) == type(None):
        new_x = np.linspace(xmin, xmax, xpoints)
    if type(new_y) == type(None):
        new_y = np.linspace(ymin, ymax, ypoints)

    new_z = f(new_x, new_y)
    new_z = np.transpose(new_z)

    return xmin, xmax, ymin, ymax, new_x, new_y, new_z


def img_to_sca(mcpRSXS_scales, mcpRSXS_scatters, img, x_low=None, x_high=None, y_low=None, y_high=None, axis=0, new_x=None, new_y=None):

    xmin, xmax, ymin, ymax, new_x, new_y, new_z = grid_data(
        mcpRSXS_scales, mcpRSXS_scatters, img, x_low, x_high, y_low, y_high, new_x, new_y)

    newz_sum = np.sum(new_z, axis=axis)

    if axis == 0:
        caxis = new_y
    elif axis == 1:
        caxis = new_x
    else:
        raise TypeError("Invalid axis")

    return caxis, newz_sum


def grid_stack(mcpRSXS_scales, mcpRSXS_scatters):

    z_stack = dict()
    x_min = dict()
    x_max = dict()
    y_min = dict()
    y_max = dict()

    for i in range(len(mcpRSXS_scatters)):
        xmin, xmax, ymin, ymax, new_x, new_y, new_z = grid_data(
            mcpRSXS_scales, mcpRSXS_scatters, i)

        z_stack[i] = new_z
        x_min[i] = xmin
        x_max[i] = xmax
        y_min[i] = ymin
        y_max[i] = ymax

    return z_stack, x_min, x_max, y_min, y_max


def check_array_ordering(mcpRSXS_scales, dictkey):
    if mcpRSXS_scales[dictkey][0][0] < mcpRSXS_scales[dictkey][0][1]:
        sorter_0 = 'asc'
    else:
        sorter_0 = 'desc'

    if mcpRSXS_scales[dictkey][1][0] < mcpRSXS_scales[dictkey][1][1]:
        sorter_1 = 'asc'
    else:
        sorter_1 = 'desc'

    return sorter_0, sorter_1


def get_interpolation_boundaries(mcpRSXS_scales):

    sorter_0_vals_min = []
    sorter_0_vals_max = []
    sorter_1_vals_min = []
    sorter_1_vals_max = []
    min_diff_0_vals = []
    min_diff_1_vals = []

    sorter_0, sorter_1 = check_array_ordering(mcpRSXS_scales, 0)

    for k, v in mcpRSXS_scales.items():
        sorter_0_vals_min.append(v[0, 0])
        sorter_1_vals_min.append(v[1, 0])

        sorter_0_vals_max.append(v[0, -1])
        sorter_1_vals_max.append(v[1, -1])

        min_diff_0_vals.append(np.abs(np.diff(v[0])).min())
        min_diff_1_vals.append(np.abs(np.diff(v[1])).min())

        min_diff_0 = min(min_diff_0_vals)
        min_diff_1 = min(min_diff_1_vals)

        if sorter_0 == 'asc':
            sorter_0_min = max(sorter_0_vals_min)
            sorter_0_max = min(sorter_0_vals_max)
        else:
            sorter_0_min = min(sorter_0_vals_min)
            sorter_0_max = max(sorter_0_vals_max)

        if sorter_1 == 'asc':
            sorter_1_min = max(sorter_1_vals_min)
            sorter_1_max = min(sorter_1_vals_max)

        else:
            sorter_1_min = min(sorter_1_vals_min)
            sorter_1_max = max(sorter_1_vals_max)

    return sorter_0_min, sorter_0_max, sorter_1_min, sorter_1_max, min_diff_0, min_diff_1


def stack_to_mca(mcpRSXS_scales, mcpRSXS_scatters, x_low=None, x_high=None, y_low=None, y_high=None, axis=0):

    sorter_0_min, sorter_0_max, sorter_1_min, sorter_1_max, min_diff_0, min_diff_1 = get_interpolation_boundaries(
        mcpRSXS_scales)

    if x_low != None:
        if x_low < sorter_0_min:
            xmin = sorter_0_min
        else:
            xmin = x_low
    else:
        xmin = sorter_0_min

    if x_high != None:
        if x_high < sorter_0_max:
            xmax = x_high
        else:
            xmax = sorter_0_max
    else:
        xmax = sorter_0_max

    if y_low != None:
        if y_low < sorter_1_min:
            ymin = sorter_1_min
        else:
            ymin = y_low
    else:
        ymin = sorter_1_min

    if y_high != None:
        if y_high < sorter_1_max:
            ymax = y_high
        else:
            ymax = sorter_1_max
    else:
        ymax = sorter_1_max

    xpoints = int(np.abs(np.ceil((xmax-xmin)/min_diff_0)))
    ypoints = int(np.abs(np.ceil((ymax-ymin)/min_diff_1)))

    new_x = np.linspace(xmin, xmax, xpoints)
    new_y = np.linspace(ymin, ymax, ypoints)

    mca = []
    for i in range(len(mcpRSXS_scatters)):
        caxis, newz_sum = img_to_sca(
            mcpRSXS_scales, mcpRSXS_scatters, i, xmin, xmax, ymin, ymax, axis, new_x, new_y)
        mca.append(newz_sum)

    if axis == 0:
        new_x = np.arange(0, len(mcpRSXS_scales))
        image = np.transpose(np.stack(mca))

    elif axis == 1:
        new_y = np.arange(0, len(mcpRSXS_scatters))
        image = np.stack(mca)
    else:
        raise TypeError("Invalid axis.")

    return new_x, new_y, image