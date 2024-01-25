import os
import sys              # sys.stdout.flush
from ctypes import *
import numpy as np

current_dir = os.getcwd()
_halo = CDLL(current_dir + "/libhalo.so")
_halo.int_rtbp.argtypes = (c_double, POINTER(c_double), c_double, c_double,
        c_double)

def int_rtbp(t, x, tol, hmin, hmax):
    global _halo
    array_type = c_double * 6
    x = (array_type)(*x)
    result = _halo.int_rtbp(c_double(t), x, c_double(tol), 
            c_double(hmin), c_double(hmax))
    return np.array(x)

# Set numpy to use this lambda function for every float it prints out
np.set_printoptions(formatter={'float': lambda x: "{0:0.16e}".format(x)})

# Josep passed me the i.c. as x,y,z,xd,yd,zd
xJ = np.array([-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00])

# My vectorfield expects the order x,px,y,py,z,pz. 
# Recall that xd=px+y; yd=py-x; zd=pz;
# Hence px=xd-y; py=yd+x; pz=zd;
x = np.array([xJ[0], xJ[3]-xJ[1], xJ[1], xJ[4]+xJ[0], xJ[2], xJ[5]])

t = 0.3059226605957322E+01  # integration time
tol = 1.e-14        # tolerance to integration error
n = 128             # number of points in the orbit
h = t/n             # required step size for n points
hmin = h            # min step size
hmax = h            # max step size

print("The initial condition is: ", x)
print("The step size is: ", h, flush=True)
y = int_rtbp(t, x, tol, hmin, hmax)
print("The final condition is: ", y)
print("Error = Initial condition - Final condition = ", np.subtract(x, y))
