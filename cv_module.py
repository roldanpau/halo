## \file cv_module.py
#
# \brief General changes of variables.

import numpy as np

## dimension of RTBP system
DIM = 6

## Transform coordinates from position-momenta to position-velocities.
#
#
# @param[in]  vx 	position-momenta (x,px,y,py,z,pz).
# @param[out] vy 	position-velocities (x,y,z,dx,dy,dz).

def posmom_to_posvel(vx):
   vy = np.empty([DIM,])
   x = vx[0]; px = vx[1]
   y = vx[2]; py = vx[3]
   z = vx[4]; pz = vx[5]
   dx = px+y
   dy = py-x
   dz = pz
   vy[0]=x; vy[3]=dx
   vy[1]=y; vy[4]=dy
   vy[2]=z; vy[5]=dz
   return vy

## Transform coordinates from position-velocities to position-momenta.
#
#
# @param[in]  vx 	position-velocities (x,y,z,dx,dy,dz).
# @param[out] vy 	position-momenta (x,px,y,py,z,pz).

def posvel_to_posmom(vx):
   vy = np.empty([DIM,])
   x = vx[0]; dx = vx[3]
   y = vx[1]; dy = vx[4]
   z = vx[2]; dz = vx[5]
   px = dx-y
   py = dy+x
   pz = dz
   vy[0]=x; vy[1]=px
   vy[2]=y; vy[3]=py
   vy[4]=z; vy[5]=pz
   return vy
