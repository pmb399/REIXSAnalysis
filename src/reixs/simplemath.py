import numpy as np


def apply_offset(stream, offset=None, coffset=None):
    if offset != None:
        offsetarray = np.array(offset)
        coeff = np.polyfit(
            offsetarray[:, 0], offsetarray[:, 1], deg=len(offsetarray)-1)

        # Make sure that constant shift (only one tuple provided is handled as offset)
        if len(coeff) == 1:
            shift = offsetarray[0, 1] - offsetarray[0, 0]
            stream = stream+shift
        else:
            stream = np.polyval(coeff, stream)

    else:
        pass

    if coffset != None:
        return stream+coffset
    else:
        return stream


def take_derivative1d(x, y, deg):
    if deg == 1:
        return np.gradient(y, x)
    elif deg == 2:
        dy = np.gradient(y, x)
        return np.gradient(dy, x)
    else:
        raise TypeError("No other derivatives implemented.")


def take_derivative2d(z, x, y, deg):
    if deg == 1:
        return np.gradient(z, x, y)
    elif deg == 2:
        dz = np.gradient(z, x, y)
        return np.gradient(dz, x, y)
    else:
        raise TypeError("No other derivatives implemented.")
