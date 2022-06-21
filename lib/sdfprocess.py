"""
2021
GD Guillaume Dideron
Functions for analysing sdf files
"""
import os
from os.path import dirname, realpath
import sys
import numpy as np
from matplotlib import pyplot as plt
from .sdf import SDFData,writeArrayToSDF
from .operations import time_order,L2_norm,get_grid_size


pixel_density = 300
def load_sdfs(filenames):
    files = []
    data = []
    coord = []
    times = [] 
    for i in range(len(filenames)):
        filestream = open(filenames[i],"rb")
        sdf = SDFData(filestream)
        data.append(sdf.getData())
        coord.append(sdf.getCoord())
        times.append(sdf.timeArray)
        filestream.close()
    return data,coord,times
    

def convergence_order( filenames,check_time=True,
                       output_file = "convergence-order.sdf",
                       plot = False, name = None ):


#    if len(times) == 0:
#        times = sdf.timeArray
#        time_steps_number = len(sdf.timeArray)
#    else:
#        if len(sdf.timeArray) < time_steps_number:
#            time_steps_number = len(sdf.timeArray)
#            times = times[:time_steps_number]
#            nan_time = times[-1]
#        np.testing.assert_almost_equal(times,sdf.timeArray)
    try:
        assert len(filenames) == 3,"Must have 3 step sizes, has " +\
        "{}.".format(len(filenames))
    except AttributeError:
        raise TypeError
    data,coord,times = load_sdfs(filenames)
    min_time_step_number = min([len(row) for row in times])
    time = np.array(times[0][:min_time_step_number])
    if check_time:
        for time_array in times:
            np.testing.assert_array_equal(time_array[:min_time_step_number],time)
    data = [datum[1:min_time_step_number] for datum in data]
    out = open(output_file,"wb")
    order = time_order(*data,2,*coord) 
    bbox = np.array([time[1], time[-1]])
    writeArrayToSDF( order, bbox, out, dataName="order",
                     coordName="t", coordData=time[1::] )
    
    if plot:
       plt.figure(dpi=pixel_density)
       plt.ylabel("Convergence Order")
       plt.xlabel("Time")
       plt.plot(times[1:],order,'k.')
       if name != None:
           plt.title(name)
       plt.savefig("convergence-order.png",format="png")
       plt.clf()


def scri_norm(filenames,output_file="log_scri_norm.sdf",dataname="log_norm"):
    data,coord,time = load_sdfs(filenames)
    data = data[0][0]
    coord = coord[0]
    time = time[0]
    norm = np.log(abs(data-data[-1]))
    out = open(output_file,"wb")
    bbox = np.array([coord[0],coord[-1]])
    norm = np.array(norm)
    writeArrayToSDF( norm, bbox, out, dataName=dataname,
                     coordName="t", coordData=coord)

def norms(filenames,output_file="norms.sdf",dataname="norms",plot=False,name=None):
    norms = []
    data,coord,times = load_sdfs(filenames)
    if plot:
        plt.figure(dpi=pixel_density)
        plt.title(name)
        plt.xlabel("x")
        plt.ylabel("Norm")
    for i,datum in enumerate(data):
        grid_size = get_grid_size(coord[i])
        norms.append(L2_norm(datum,grid_size))
        if plot:
            plt.plot(times,norms[i],'.',label="h={:e}".format(grid_size))
    if plot:
        plt.legend(loc="best")
        plt.savefig("norms.png",format="png")
        plt.clf()
    out = open(output_file,"wb")
    for norm,time in zip(norms,times):
        bbox = np.array([time[0],time[-1]])
        writeArrayToSDF( norm, bbox, out, dataName=dataname,
                         coordName="t", coordData=time )


def plot_sdf(filename,output,coordinate=None,ylabel=None,title=None):
    data,coord,times = load_sdfs([filename])
    plt.figure(dpi=pixel_density) 
    plt.title(title)
    plt.xlabel(coordinate)
    plt.ylabel(ylabel)
    plt.plot(coord[0],data[0][0],'k.')
    plt.savefig(output,format="png")
