## \file control_module.py
#  \brief Functions to control a (perturbed) trajectory around a nominal LPO orbit.
#
#  Assume that a trajectory close to a nominal LPO orbit is perturbed due e.g.
#  to integration errors (of size 10^{-5} -- 10^{-4}). Apply control to this
#  trajectory by making small adjustments in velocity every CORREC_TIME (around
#  180 days), and maintain it close to the nominal orbit for longer than 20
#  years.
#
#  We extend the time interval within the LPO region by querying the predictive
#  model (either a polynomial regression model or a neural network) about which
#  dv correction to use (function correction_NN/correction_regression). The
#  model has been trained in supervised.ipynb or neuralnet.ipynb.
#

import os
import sys                      # sys.stdout.flush
from ctypes import *
import numpy as np
from cv_module import * 	    # posmom_to_posvel, posvel_to_posmom
from correction_module import * # correction_regression, correction_NN
from numpy.linalg import norm

## dimension of RTBP system
DIM = 6     

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

_halo.int_correction_opt.restype = c_double
_halo.int_correction_opt.argtypes = (POINTER(c_double), c_double, c_double,
        c_int, POINTER(c_double), POINTER(c_double))

## Interface to corresponding C function
def int_correction_opt(q_Masde, CORREC_TIME, SHADOW_TIME, corr, q90, q90_new):
    global _halo
    array_type = c_double * DIM
    q_Masde = (array_type)(*q_Masde)
    q90 = (array_type)(*q90)
    q90_new = (array_type)(*q90_new)
    dv = _halo.int_correction_opt(q_Masde, c_double(CORREC_TIME),
            c_double(SHADOW_TIME), c_int(corr), q90, q90_new)
    return dv, np.array(q90), np.array(q90_new)

_halo.correction_opt.restype = c_double
_halo.correction_opt.argtypes = (POINTER(c_double), c_double, c_int,
        POINTER(c_double))

## Interface to corresponding C function
def correction_opt(q, SHADOW_TIME, corr, q_new):
    global _halo
    array_type = c_double * DIM
    q = (array_type)(*q)
    q_new = (array_type)(*q_new)
    dv = _halo.correction_opt(q, c_double(SHADOW_TIME), c_int(corr), q_new)
    return dv, np.array(q_new)

# Set numpy to use this lambda function for every float it prints out
np.set_printoptions(formatter={'float': lambda x: "{0:0.16e}".format(x)})

def print_array(arr):
    """
    prints a 1-D numpy array in a nicer format
    """
    for elem in arr:
        print("{:.3e}".format(elem), end=" ")
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


# Control a (perturbed) trajectory around a nominal LPO orbit.
#
#  Assume that a trajectory close to a nominal LPO orbit is perturbed due e.g.
#  to integration errors (of size 10^{-5} -- 10^{-4}). Apply control to this
#  trajectory by making small adjustments in velocity every CORREC_TIME (around
#  180 days), and maintain it close to the nominal orbit for longer than 20
#  years.
#
#  We extend the time interval within the LPO region by querying the predictive
#  model (either a polynomial regression model or a neural network) about which
#  dv correction to use (function correction_NN/correction_regression). The
#  model has been trained in supervised.ipynb or neuralnet.ipynb.
#
# If procedure is unable to keep the trajectory close to LPO orbit for the
# stipulated period of yrs, funcion returns bOK = False. Otherwise, bOK = True.
#
# @param[in]      X1    I.c.  (pos-vel coordinates)
# @param[in]      yrs   Length of shadowing (in years)
# @param[out]     bOK   Return flag (True or False)    

def control(X1, yrs):
    ## \brief Period of nominal halo orbit. 
    #
    # The period is close to \f$\pi\f$, i.e. approximately 180 days. 
    #
    T = 0.3059226605957322E+01

    ## 20 yrs (in normalized units)
    twentyYrs = yrs*2*np.pi  

    # GOLDEN_FRACT =20.381966    ##< Fraction of a circle occupied by the golden angle

    ## Perform one correction dv every CORREC_TIME
    CORREC_TIME = T	
    ## Correction stays in LPO region longer than SHADOW_TIME
    SHADOW_TIME = 2*T	

    #print("The initial condition is: ", X1)

    ## Total time of extended orbit (in LPO region)
    time = 0.0 

    q = np.empty([DIM,])      # Make space for q
    q90 = np.empty([DIM,])      # Make space for q90
    q90_new = np.empty([DIM,])  # Make space for q90_new

    print(X1)

    # Perturb X1 randomly
    #R = 1.e-4   # Size of perturbation
    #X = X1 + R*(np.random.rand(1, DIM)-0.5).flatten()

    # Find correction maneuver according to predictor
    #dv = correction_regression(X)
    #if(dv==0):
    #    return False

    # Apply correction maneuver.
    # This is not strictly necessary for correction_opt, where dv is
    # actually applied and modified q90 is returned as q90_new.
    # However, it IS necessary for correction_regression or correction_NN.
    #X = apply_correction_st(X, dv)

    X = X1

    while(time < twentyYrs):

        # Perturb X randomly
        #X = X + R*(np.random.rand(1, DIM)-0.5)

        # Integrate orbit for CORREC_TIME
        q = posvel_to_posmom(X)
        q = int_rtbp(CORREC_TIME, q, 1.e-15, 1.e-5, CORREC_TIME/10, 0)
        q90 = posmom_to_posvel(q)

        time = time + CORREC_TIME

        # Find correction maneuver according to the predictor (either regression or
        # NN)
        dv = correction_regression(q90)
        if(dv==0):
            return False

        # Apply correction maneuver, and iterate again.
        # This is not strictly necessary for correction_opt, where dv is
        # actually applied and modified q90 is returned as q90_new.
        # However, it IS necessary for correction_regression or correction_NN.
        X = apply_correction_st(q90, dv)
        
        # The first 3 corrections are not printed, since we are initially on the
        # halo and they are not significant.
        #if(time > 3*CORREC_TIME)
        print("Time: ", '%0.2f' % (time/(2*np.pi)), "yrs")
        #print("state before maneuver: ", end=" ")
        #print_array(q90) 
        print("distance to LPO: ", '%.0e' % (norm(q90-X1)))
        print("dv predicted: ", '%.4e' % dv)
        print("state after maneuver: ", end=" ")
        print_array(X) 
        print("")

    #print("The final condition is: ", X)
    return True
