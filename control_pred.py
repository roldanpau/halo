## \file control_pred.py
#  \brief Control a (perturbed) trajectory around a nominal LPO orbit.
#
#  Assume that a trajectory close to a nominal LPO orbit is perturbed due e.g.
#  to integration errors (of size 10^{-5} -- 10^{-4}). Apply control to this
#  trajectory by making small adjustments in velocity every CORREC_TIME (around
#  180 days), and maintain it close to the nominal orbit for longer than 20
#  years.
#  
#  USAGE:  
#   (Remember to `conda activate halo' environment first)
#   (Remember to generate the correct regression predictor by running the
#   supervised_cntrl.ipynb notebook first)
#   python control_pred.py
#   python control_pred.py >control_pred.res 2>control_pred.err

import os
import sys                      # sys.stdout.flush
from ctypes import *
import numpy as np
from control_module import *  # control

## I.c. on the nominal halo orbit (at t=0), given as x,y,z,xd,yd,zd
X1 = np.array([-0.9916647163367744E+00,  0.0000000000000000E+00,
    0.8983543483564242E-03, -0.0000000000000000E+00,  0.9931014021976879E-02,
    0.0000000000000000E+00])

bOK = control(X1,20)
if(bOK):
    print("Shadowing of X1 for", 20, "years was successful")
