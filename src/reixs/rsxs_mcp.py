from .simplemath import apply_offset, take_derivative1d
from .ReadData import REIXS
import numpy as np


def loadRSXS1dROIscans(basedir, file, images, axis, *args, x=[None, None], y=[None, None], norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None, deriv=None):

    data = dict()
    REIXSobj = REIXS(basedir, file)

    if axis not in [0, 1]:
        raise UserWarning("Invalid axis!")

    if len(args) == 0:
        raise UserWarning("You did not specify a scan")

    for arg in args:
        data[arg] = REIXSobj.Scan(arg)
        data[arg].imagesca = dict()
        data[arg].caxis = dict()

        if axis == 0:
            data[arg].caxis_label = data[arg].mcpRSXS_axes[1]
        else:
            data[arg].caxis_label = data[arg].mcpRSXS_axes[0]

        for image in images:
            data[arg].caxis[image], data[arg].imagesca[image] = data[arg].RSXS_1dROI(image,
                                                                                     x_low=x[0], x_high=x[1], y_low=y[0], y_high=y[1], axis=axis)
            data[arg].imagesca[image] = apply_offset(
                data[arg].imagesca[image], yoffset, ycoffset)
            data[arg].caxis[image] = apply_offset(
                data[arg].caxis[image], xoffset, xcoffset)

            if norm == True:
                data[arg].imagesca[image] = data[arg].imagesca[image] / \
                    np.max(data[arg].imagesca[image])

            if deriv != None:
                data[arg].imagesca[image] = take_derivative1d(
                    data[arg].imagesca[image], data[arg].caxis[image], deriv)
                if norm == True:
                    data[arg].imagesca[image] = data[arg].imagesca[image] / \
                        data[arg].imagesca[image].max()

    return data


def loadRSXS2dROIscans(basedir, file, axis, *args, x=[None, None], y=[None, None], norm=False, xoffset=None, xcoffset=None, yoffset=None, ycoffset=None):

    data = dict()
    REIXSobj = REIXS(basedir, file)
    for arg in args:
        data[arg] = REIXSobj.Scan(arg)
        data[arg].imagesca = dict()
        data[arg].caxis_labels = data[arg].mcpRSXS_axes

        if axis == 0:
            data[arg].caxis_labels[0] = "Index"
        elif axis == 1:
            data[arg].caxis_labels[-1] = "Index"
        else:
            raise UserWarning("Invalid axis.")

        new_x, new_y, image = data[arg].RSXS_2dROI(
            x_low=x[0], x_high=x[1], y_low=y[0], y_high=y[0], axis=axis)
        data[arg].mcp_x = apply_offset(new_x, xoffset, xcoffset)
        data[arg].mcp_y = apply_offset(new_y, yoffset, ycoffset)
        data[arg].imagemca = image

        if norm == True:
            data[arg].imagemca = data[arg].imagemca / \
                np.max(data[arg].imagemca)

    return data


def loadRSXSImageStack(basedir, file, *args):
    data = dict()
    REIXSobj = REIXS(basedir, file)
    for arg in args:
        data[arg] = REIXSobj.Scan(arg)
        data[arg].caxis_labels = data[arg].mcpRSXS_axes

        z_stack, x_min, x_max, y_min, y_max = data[arg].RSXS_Images(
        )

        data[arg].imagemca = z_stack
        data[arg].x_min = x_min
        data[arg].x_max = x_max
        data[arg].y_min = y_min
        data[arg].y_max = y_max

    return data
