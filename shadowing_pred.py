## \file shadowing_pred.py
#  \brief Construct a trajectory that ``shadows'' an LPO orbit.
#
#  Construct a trajectory that ``shadows'' an LPO orbit for 20 years by making
#  tiny adjustments in velocity every CORREC_TIME (around 180 days).
#  
#  We extend the time interval within the LPO region by querying the predictive
#  model (either a polynomial regression model or a neural network) about which
#  dv correction to use (function correction_NN/correction_regression). The
#  model has been trained in supervised.ipynb or neuralnet.ipynb.
#  
#  USAGE:  
#   (Remember to `conda activate halo' environment first)
#   python shadowing_pred.py
#   python shadowing_pred.py >shadowing_pred.res 2>shadowing_pred.err

import os
import sys                      # sys.stdout.flush
from ctypes import *
import numpy as np
from cv_module import * 	    # posmom_to_posvel, posvel_to_posmom
from correction_module import * # correction_regression, correction_NN

## dimension of RTBP system
DIM = 6     

## \brief Period of nominal halo orbit. 
#
# The period is close to \f$\pi\f$, i.e. approximately 180 days. 
#
T = 0.3059226605957322E+01

## 20 yrs (in normalized units)
twentyYrs = 20*2*np.pi  

# GOLDEN_FRACT =20.381966    ##< Fraction of a circle occupied by the golden angle

## Perform one correction dv every CORREC_TIME
CORREC_TIME = T	
## Correction stays in LPO region longer than SHADOW_TIME
SHADOW_TIME = 3*T	

# Interface to C functions 

current_dir = os.getcwd()
_halo = CDLL(current_dir + "/libhalo.so")

_halo.int_rtbp.argtypes = (c_double, POINTER(c_double), c_double, c_double,
        c_double, c_int)

## Interface to corresponding C function
def int_rtbp(t, x, tol, hmin, hmax, bPrint):
    global _halo
    array_type = c_double * DIM
    x = (array_type)(*x)
    result = _halo.int_rtbp(c_double(t), x, c_double(tol), 
            c_double(hmin), c_double(hmax), c_int(bPrint))
    return np.array(x)

_halo.correction_opt.restype = c_double
_halo.correction_opt.argtypes = (POINTER(c_double), c_double, c_double,
        c_int, POINTER(c_double), POINTER(c_double))

## Interface to corresponding C function
def correction_opt(q_Masde, CORREC_TIME, SHADOW_TIME, corr, q90, q90_new):
    global _halo
    array_type = c_double * DIM
    q_Masde = (array_type)(*q_Masde)
    q90 = (array_type)(*q90)
    q90_new = (array_type)(*q90_new)
    dv = _halo.correction_opt(q_Masde, c_double(CORREC_TIME),
            c_double(SHADOW_TIME), c_int(corr), q90, q90_new)
    return dv, np.array(q90), np.array(q90_new)

# Set numpy to use this lambda function for every float it prints out
np.set_printoptions(formatter={'float': lambda x: "{0:0.16e}".format(x)})

def print_array(arr):
    """
    prints a 1-D numpy array in a nicer format
    """
    for elem in arr:
        print("{}".format(elem), end=" ")
    print("")

## Apply correction dv to q90
#
# This function applies the maneuver dv to the initial condition \f$
# q_{90}=(x,y,z,v_x,v_y,v_z) \f$ in the STable direction \f$ E^s \f$. The IC
# after applying the maneuver, \f$ q90 + \Delta v\f$, is returned in q90_new.
#
# @param[in]      q90         I.c.  (pos-vel coordinates)
# @param[in]      dv          Magnitude of the maneuver (may be negative)

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

## I.c. close to the nominal halo orbit, given as x,y,z,xd,yd,zd
X1 = np.array([-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00])

#print("The initial condition is: ", X1)

## Total time of extended orbit (in LPO region)
time = 0.0 

# The first 3 corrections are not printed, since we are initially on the
# halo and they are not significant.
while(time < 3*CORREC_TIME):

    # Find correction maneuver according to optimal method
    q90 = np.empty([DIM,])      # Make space for q90
    q90_new = np.empty([DIM,])
    [dv, q90, q90_new] = correction_opt(X1, CORREC_TIME, SHADOW_TIME, 1, q90,
            q90_new)

    # Apply correction maneuver, and iterate again.
    X1 = apply_correction_st(q90, dv)
    print_array(X1)
    time = time + CORREC_TIME

while(time < twentyYrs):

    # No need to integrate orbit for CORREC_TIME, since this is done inside
    # correction_opt

    # Find correction maneuver according to optimal method
    q90 = np.empty([DIM,])      # Make space for q90
    q90_new = np.empty([DIM,])
    [dv, q90, q90_new] = correction_opt(X1, CORREC_TIME, SHADOW_TIME, 1, q90,
            q90_new)

    # Find correction maneuver according to the predictor (either regression or
    # NN)
    dv_pred = correction_regression(q90)
    dv_pred = dv_pred[0,0]

    print("Years: ", time/(2*np.pi), "\tdv_opt: ", dv, "\tdv_pred: ", dv_pred)

    # Apply correction maneuver, and iterate again.
    X1 = apply_correction_st(q90, dv_pred)
    #print_array(X1)
    time = time + CORREC_TIME

#print("The final condition is: ", X1)
