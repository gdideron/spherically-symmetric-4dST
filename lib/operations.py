"""
2021
GD Guillaume Dideron
Ensure that the data for different step sizes is aligned in both time and
space.
"""
import numpy as np


class Error(Exception):
    pass

class NotNumpyArray(Error):
    def __init__(self):
        message = "Data and coordinates must be given in Numpy arrays."
        super().__init__(message=message)

class ShapeError(Error):
    def __init__(self):
        message = "Data for different grid sizes have incompatible shapes."
        super().__init__(message=message)

def data_check( data_1,data_2,data_3,grid_factor=2,
                coord_1=None,coord_2=None,coord_3=None):
    """Checks data type and array shape for all step sizes"""
    error_shape = "Data for all step sizes must be of same shape"
    data = [data_1,data_2,data_3]
    coord = [coord_1,coord_2,coord_3]
    coord_loaded = False 
    if any(c is not None for c in coord):
        if all(c is not None for c in coord):
            coord_loaded = True
        else:
            print( "Coordinates must be given for all " + 
                   "or none of the data: Assuming none given.")
    for i in range(len(data)):
        try:
            assert isinstance(data[i],np.ndarray) 
        except AssertionError:
            raise NotNumpyArray
        try:
            assert data[i].shape[0] == data_1.shape[0]
        except AssertionError:
            raise ShapeError

        # Checks if data have compatible sizes, 
        # but sample grid points could still be misaligned
        factor = int(grid_factor**i)
        assert data[i][:,::factor].shape[1] == data_1.shape[1]
        # Can check alignment if coordinates given
        if coord_loaded:
            assert coord[::factor] == coord


def time_order( data_1,data_2,data_3,grid_factor,
                coord_1=None,coord_2=None,coord_3=None ):
    """Returns time array of global convergence order"""
    data_check( data_1,data_2,data_3,
                coord_1=None,coord_2=None,coord_3=None )
    grid_size_1 = get_grid_size(coord_2)
    grid_size_2 = get_grid_size(coord_1)
    difference_1 = L2_norm(data_3[:,::grid_factor] - data_2,grid_size_1)
    difference_2 = L2_norm(data_2[:,::grid_factor] - data_1,grid_size_2)
    multiplier = -1.0/np.log(grid_factor)
    order_time_array = multiplier*np.log(difference_1/difference_2) 
    return order_time_array


def get_grid_size(coord):
    try :
        assert isinstance(coord,np.ndarray)
    except AssertionError:
        print(coord)
        raise NotNumpyArray
    grid_size = coord[1] - coord[0]
    #assert np.all((coord[1:] - coord[:-1]) == grid_size), "Grid size not fixed"
    return grid_size


def L2_norm(array,step_size):
    """For now, only for 1d and fixed step size"""
    return np.linalg.norm(array,axis=1)*np.sqrt(step_size)


def deriv_1_o4(array,coord):
    """
    First order derivative to fourth order accuracy for constant size grid.
    Uses sided operator at boundary
    """
    dx = get_grid_size(coord)
    back_2 = np.roll(array,2,axis=1)
    back_1 = np.roll(array,1,axis=1)
    front_1 = np.roll(array,-1,axis=1)
    front_2 = np.roll(array,-2,axis=1)
    deriv = 1.0/12*(back_2 - 8*back_1 + 8*front_1 - front_2)/dx
    deriv[:,0] = 1.0/12*( - 25*array[:,0] + 48*array[:,1] - 36*array[:,2] 
                        + 16*array[:,3] - 3*array[:,4] )/dx
    deriv[:,1] = 1.0/12*( - 3*array[:,0] - 10*array[:,1] + 18*array[:,2] 
                        - 6*array[:,3] + array[:,4] )/dx
    deriv[:,-2] = 1.0/12*( 3*array[:,-1] + 10*array[:,-2] - 18*array[:,-3] 
                        + 6*array[:,-4] - array[:,-5] )/dx
    deriv[:,-1] = 1.0/12*( 25*array[:,-1] - 48*array[:,-2] + 36*array[:,-3] 
                        - 16*array[:,-4] + 3*array[:,-5] )/dx
    return deriv


