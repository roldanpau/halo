## \file shadowing_pred.py
#  \brief Construct a trajectory that ``shadows'' an LPO orbit.
#
#  Construct a trajectory that ``shadows'' an LPO orbit for 20 years by making
#  tiny adjustments in velocity every CORREC_TIME (around 180 days).
#
#  This program is OBSOLETE, since shadowing() does not use
#  correction_NN/correction_supervised anymore, it uses correction_opt!!!
#  The reason we don't do shadowing using a prediction model is that the
#  prediction model was derived from a table of control maneuvers and has a
#  RMSE of 10^{-6}, which is not good enough to predict shadowing maneuvers
#  (with stdev 10^{-11}).
#  
#  USAGE:  
#   (Remember to `conda activate halo' environment first)
#   (Remember to generate the correct regression predictor by running the
#   supervised_cntrl.ipynb notebook first)
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

bOK = shadowing(X1,20)
if(bOK):
    print("Shadowing of X1 for", 20, "years was successful")
