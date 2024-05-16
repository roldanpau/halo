## \file shadowing_NN.py
#  \brief Construct a trajectory that ``shadows'' an LPO orbit.
#
#  Construct a trajectory that ``shadows'' an LPO orbit for 20 years by making
#  tiny adjustments in velocity every CORREC_TIME (usually 90 days).
#  
#  We extend the time interval within the LPO region by querying the Neural
#  Network about which dv correction to use (function correction_NN). The NN
#  has been trained in neuralnet.ipynb.
#  
#  USAGE:  ./shadowing
#  

import os
import sys              # sys.stdout.flush
from ctypes import *
import numpy as np

from cv_module import * 	# posmom_to_posvel, posvel_to_posmom

from correction_module import * 	# correction_NN

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

## \brief Period of nominal halo orbit. 
#
# The period is close to \f$\pi\f$, i.e. approximately 180 days. 
#
T = 0.3059226605957322E+01

twentyYrs = 20*2*np.pi  ##< 20 yrs (in normalized units)
#twentyYrs = 2*np.pi  ##< 20 yrs (in normalized units)

GOLDEN_FRACT = 0.381966    ##< Fraction of a circle occupied by the golden angle

CORREC_TIME = 3*GOLDEN_FRACT*T	##< Perform one correction dv every CORREC_TIME
SHADOW_TIME = 3*T	##< Correction stays in LPO region longer than SHADOW_TIME

def print_array(arr):
    """
    prints a 1-D numpy array in a nicer format
    """
    for elem in arr:
        print("{}".format(elem), end=" ")
    print("")

def apply_correction_st(q90, dv):

    # w = stable direction
    # Masde says: Take only first 3 (position) components
    w = np.array([0, 0, 0, 0.323988273471917, 0.173290594261497, 0])

    w_mod = np.linalg.norm(w, 2)       # Modulus of stable direction

    # Set maneuver to dv*w, where w is the unitary stable dir
    w = w/w_mod
    Dv = dv*w

    # Apply maneuver to q90
    q90_new = q90 + Dv

    return q90_new

# I.c. close to the nominal halo orbit, given as x,y,z,xd,yd,zd
X1 = np.array([-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00])

#print("The initial condition is: ", X1)

time = 0.0 ##< Total time of extended orbit (in LPO region)
while(time < twentyYrs):

    # My vectorfield expects the order x,px,y,py,z,pz
    q = posvel_to_posmom(X1)

    # Integrate orbit for CORREC_TIME
    q = int_rtbp(CORREC_TIME, q, 1.e-15, 1.e-5, 1.e-1, 0)

    q90 = posmom_to_posvel(q)
    print("q_90: ", q90)
    #print_array(q90)

    # Find correction maneuver, according to the neural network
    dv = correction_NN(q90)
    dv = dv[0,0]
    print("dv: ", dv)

    # Apply correction maneuver, and iterate again.
    X1 = apply_correction_st(q90, dv)
    print("X1: ", X1)

    time = time + CORREC_TIME

#print("The final condition is: ", X1)
