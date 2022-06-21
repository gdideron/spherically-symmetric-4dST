import numpy as np
from sdfprocess import load_sdfs
from operations import deriv_1_o4


def char_determinant(lambda2, ffile, gfile, output_file="char-speeds.sdf"):
    data, coord, times, nan_time = load_sdfs([ffile,gfile])
    try:
       np.assert_array_equal(coord[0],coord[1])
    except AssertionError as error:
        raise error
    print(coord)
    coord = coord[0]
    print(coord)
    f = data[0]
    g = data[1]
    f2 = f*f 
    g2 = g*g 
    fx = deriv_1_o4(f,coord)
    gx = deriv_1_o4(g,coord)
    return 1 + lambda2*(f2+g2+2*fx+4*f/coord+(f*g+gx)**2*lambda2)
