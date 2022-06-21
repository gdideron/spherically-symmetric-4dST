import numpy as np
import matplotlib.pyplot as plt
from lib.weast.sdf import writeArrayToSDF
from sdfprocess import load_sdfs
from operations import deriv_1_o4

def char_determinant(lambda_1, ffile, gfile,
        output_file="char-speeds.sdf",output=True):
    data, coord, times = load_sdfs([ffile,gfile])
    try:
       np.testing.assert_array_equal(coord[0],coord[1])
    except AssertionError as error:
        raise error
    x = coord[0]
    xinv = 1./x
    xinv2 = xinv*xinv
    B = 1 + 2*xinv
    C = 1 - 2*xinv
    D = 1 - 12*xinv2
    f = data[0]
    g = data[1]
    f2 = f*f 
    g2 = g*g 
    fg = f*g
    fx = deriv_1_o4(f,x)
    gx = deriv_1_o4(g,x)
    g_tt = -B + lambda_1*( -3*B*B*g2 + 12*xinv*B*fg + D*f2)
    g_tx = 2*xinv + lambda_1*(6*xinv*B*g2 + 2*D*fg - 6*xinv*C*f2)
    g_xx = C + lambda_1*(D*g2 - 12*xinv*C*fg - 3*C*C*f2)
    det = np.zeros(f.shape)
    det = -g_tx*g_tx + g_tt*g_xx
    if output:
        out = open(output_file,"wb")
        writeArrayToSDF( det, np.array([0,1]), out, 
                         dataName="Characteristic determinant",
                         coordName="r", coordData=x, timeArray=times[0])
    else:
        return det 

def char_speed(lambda_1, direction, ffile, gfile, output_file="char-speed.sdf"):
    data, coord, times = load_sdfs([ffile,gfile])
    try:
       np.testing.assert_array_equal(coord[0],coord[1])
    except AssertionError as error:
        raise error
    x = coord[0]
    xinv = 1./x
    xinv2 = xinv*xinv
    B = 1 + 2*xinv
    C = 1 - 2*xinv
    D = 1 - 12*xinv2
    f = data[0]
    g = data[1]
    f2 = f*f 
    g2 = g*g 
    fg = f*g
    fx = deriv_1_o4(f,x)
    gx = deriv_1_o4(g,x)
    g_tt = -B + lambda_1*( -3*B*B*g2 + 12*xinv*B*fg + D*f2)
    g_tx = 2*xinv + lambda_1*(6*xinv*B*g2 + 2*D*fg - 6*xinv*C*f2)
    g_xx = C + lambda_1*(D*g2 - 12*xinv*C*fg - 3*C*C*f2)
    speed = (-g_tx + direction*(g_tx*g_tx - g_tt*g_xx))/g_tt 
    out = open(output_file,"wb")
    writeArrayToSDF( speed, np.array([0,1]), out, 
                     dataName="char_speed{}".format(direction),
                     coordName="r", coordData=x, timeArray=times[0])
