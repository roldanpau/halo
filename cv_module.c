/* cv

   General changes of variables.
*/

/*
   Transform coordinates from position-momenta to position-velocities.

   params
   vx (input)	position-momenta (x,px,y,py,z,pz).
   vy (output)	position-velocities (x,y,z,dx,dy,dz).
*/
void posmom_to_posvel(const double vx[6], double vy[6])
{
   double x = vx[0]; double px = vx[1];
   double y = vx[2]; double py = vx[3];
   double z = vx[4]; double pz = vx[5];
   double dx = px+y;
   double dy = py-x;
   double dz = pz;
   vy[0]=x; vy[3]=dx;
   vy[1]=y; vy[4]=dy;
   vy[2]=z; vy[5]=dz;
}

void posvel_to_posmom(const double vx[6], double vy[6])
{
   double x = vx[0]; double dx = vx[3];
   double y = vx[1]; double dy = vx[4];
   double z = vx[2]; double dz = vx[5];
   double px = dx-y;
   double py = dy+x;
   double pz = dz;
   vy[0]=x; vy[1]=px;
   vy[2]=y; vy[3]=py;
   vy[4]=z; vy[5]=pz;
}
