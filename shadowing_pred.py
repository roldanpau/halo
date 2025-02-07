## \file shadowing_pred.py
#  \brief Construct a trajectory that ``shadows'' an LPO orbit.
#
#  Construct a trajectory that ``shadows'' an LPO orbit for 20 years by making
#  tiny adjustments in velocity every CORREC_TIME (around 180 days).
#  
#  USAGE:  
#   (Remember to `conda activate halo' environment first)
#   (Remember to generate the correct regression predictor by running the
#   supervised.ipynb notebook first)
#   python shadowing_pred.py
#   python shadowing_pred.py >shadowing_pred.res 2>shadowing_pred.err

import os
import sys                      # sys.stdout.flush
from ctypes import *
import numpy as np
from shadowing_module import *  # shadowing

## I.c. on the nominal halo orbit (at t=0), given as x,y,z,xd,yd,zd
X1 = np.array([-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00])

## Do NOT apply corrections inmediatly, first integrate for CORREC_TIME
bNow = 0
bOK = shadowing(X1,20,bNow)
if(bOK):
    print("Shadowing of X1 for", 20, "years was successful")
