# Check that the stable direction w computed in matlab (file eig_diferencial.m)
# is correct by checking that it is indeed (aproximately) invariant

import os
import sys              # sys.stdout.flush
from ctypes import *
import numpy as np
import math             # acos

current_dir = os.getcwd()
_halo = CDLL(current_dir + "/libhalo.so")
_halo.int_rtbp.argtypes = (c_double, POINTER(c_double), c_double, c_double,
        c_double, c_int)

def int_rtbp(t, x, tol, hmin, hmax, bPrint):
    global _halo
    array_type = c_double * 6
    x = (array_type)(*x)
    result = _halo.int_rtbp(c_double(t), x, c_double(tol), 
            c_double(hmin), c_double(hmax), c_int(bPrint))
    return np.array(x)

# Set numpy to use this lambda function for every float it prints out
np.set_printoptions(formatter={'float': lambda x: "{0:0.16e}".format(x)})

# Set i.c. to L_1 + w*\epsilon, where w = stable direction.
L1 = np.array([-0.99003, 0, 0, 0, 0, 0])
w = np.array([0.323988273471917, 0.173290594261497, 0, -0.820113803552770, -0.438651704448293, 0])
w = w/np.linalg.norm(w);    # normalize w to size 1
ic = L1 + w*0.001

# My vectorfield expects the order x,px,y,py,z,pz. 
# Recall that xd=px+y; yd=py-x; zd=pz;
# Hence px=xd-y; py=yd+x; pz=zd;
x = np.array([ic[0], ic[3]-ic[1], ic[1], ic[4]+ic[0], ic[2], ic[5]])

t = 1.e-3           # integration time
tol = 1.e-14        # tolerance to integration error
h = t/10            # step size
hmin = h/10         # min step size
hmax = t            # max step size
bPrint = 0          # set print flag to True

print("The initial condition is: ", ic)
print("The step size is: ", h, flush=True)
y = int_rtbp(t, x, tol, hmin, hmax, bPrint)

# Final condition in pos-vel (x,y,z,dx,dy,dz)
fc = np.array([y[0], y[2], y[4], y[1]+y[2], y[3]-y[0], y[5]])

v = np.subtract(ic, fc)
v = v/np.linalg.norm(v)

# Angle between v and w
alpha = math.acos(np.dot(w, v)/(np.linalg.norm(w)*np.linalg.norm(v)))

print("The final condition is: ", fc)
#print("v = Initial condition - Final condition = ", v)
print("w = ", w)
print("v = ", v)
print("Angle between v and w = ", alpha)
