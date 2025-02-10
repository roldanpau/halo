## \file controllable.py
# \brief Check if all initial conditions inside a cube are controllable.
#
# Take a bunch of initial conditions inside a small cube of radius R around a
# state X1 on the nominal halo, and check if they are all controllable, in the
# sense that we can construct a trajectory that stays close to the halo for 20
# years by making tiny adjustments in velocity periodically.
#  
#  USAGE:  
#   (Remember to `conda activate halo' environment first)
#   (IF NECESSARY, remember to generate the correct regression predictor by
#   running the supervised.ipynb notebook first)
#   python controllable.py
#   python controllable.py >controllable.res 2>controllable.err

import os
import sys                      # sys.stdout.flush
from ctypes import *
import numpy as np
from shadowing_module import *  # shadowing

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
amount = 100

# Sample small cube of radius R around X1 randomly
cube = R*(np.random.rand(amount, DIM)-0.5)+X1

for idx, x in enumerate(cube):  # For each point inside the cube
    print(idx,x)
    bOK = shadowing(x,10)  # Try shadowing for 10 years
    if(bOK):
        print("Shadowing of i.c.", idx, "for", 10, "years was successful!")
    else:
        print("Shadowing of i.c.", idx, "for", 10, "years was UNSUCCESSFUL!")
        break
