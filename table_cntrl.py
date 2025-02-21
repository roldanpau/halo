## \file table_ctrl.py
# \brief Generate table of optimal control maneuvers for supervised learning
#
# Generate a table of labelled examples of the form
# x y z dx dy dz dv
# associating to each state (x,y,z,dx,dy,dz) the correction dv (magnitude
# applied in the direction of current velocity).
# 
# We generate these examples by taking a bunch of initial conditions inside a
# small cube of radius R around a state X1 on the nominal halo, and computing
# the optimal control maneuver dv needed to keep the trajectory close to the
# halo.
#
#  USAGE:  
#   (Remember to `conda activate halo' environment first)
#   python table_ctrl.py
#   python table_ctrl.py > maneuvers_ctrl.csv

import os
import sys                      # sys.stdout.flush, sys.exit
from ctypes import *
import numpy as np
from shadowing_module import correction_opt

## dimension of RTBP system
DIM = 6

## Radius of the sampling cube
R = 1.e-4

## I.c. on the nominal halo orbit (at t=0), given as x,y,z,xd,yd,zd
#X1 = np.array([-0.9916647163367744E+00,  0.0000000000000000E+00,
#    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
#    0.0000000000000000E+00])

## I.c. at QUARTER PERIOD (T/4) of the nominal halo orbit
## (x,px,y,py,z,pz) = -9.905881e-01 -1.456353e-03 4.517090e-03 -9.909699e-01 -2.167137e-04 -2.007546e-03
## In (x,y,z, dx=px+y, dy=py-x, dz) coords.
#X1 = np.array({-9.905881e-01, 4.517090e-03, -2.167137e-04,
#	-1.456353e-03+4.517090e-03, -9.909699e-01-(-9.905881e-01),
#	-2.007546e-03])

## I.c. at T/8 of the nominal halo orbit
## (x,px,y,py,z,pz) = -9.914575e-01 -2.013749e-03 3.327175e-03 -9.850071e-01 5.350106e-04 -1.703465e-03
## In (x,y,z, dx=px+y, dy=py-x, dz) coords.
X1 = np.array([-9.914575e-01, 3.327175e-03, 5.350106e-04,
    -2.013749e-03+3.327175e-03, -9.850071e-01-(-9.914575e-01),
    -1.703465e-03])

## Amount of points inside the cube
amount = 1000

# Sample small cube of radius R around X1 randomly
cube = R*(np.random.rand(amount, DIM)-0.5)+X1

## \brief Period of nominal halo orbit. 
#
# The period is close to \f$\pi\f$, i.e. approximately 180 days. 
#
T = 0.3059226605957322E+01

## Correction stays in LPO region longer than SHADOW_TIME
SHADOW_TIME = 2*T	

q_new = np.empty([DIM,])  # Make space for q_new

print("x,y,z,dx,dy,dz,dv_st")
for idx, x in enumerate(cube):  # For each point inside the cube
    # Find correction maneuver according to optimal method
    [dv, q_new] = correction_opt(x, SHADOW_TIME, 1, q_new)
    if(dv==0):
        sys.exit("correction_opt returned an error")

    print(*x, dv, sep=',')
