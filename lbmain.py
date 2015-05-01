import numpy as np

import lbfunctionsD as lbfD


#define shape grid and pressure gradient
xdim,ydim = 30,7
pres_grad = 1
timesteps = 20

#define initial velocity density set
ones = np.ones((ydim,xdim))
ones[0] = 0
ones[-1] = 0
ones = np.swapaxes(ones,0,1)
ones[0] = 0
ones[-1] = 0
first_fill = np.swapaxes(ones,0,1)
old_vel_dens_9 = np.array([first_fill,first_fill,first_fill,first_fill,first_fill,first_fill,first_fill,first_fill,first_fill])

# Actual algorithm
for tmstp in range(timesteps):
  new_vel_dens_9 = lbfD.update_pos_inc_boundary(old_vel_dens_9)
  vel_vec = lbfD.find_vel_vec_at_pos_inc_presgrad(new_vel_dens_9,xdim,ydim,pres_grad)
  vel_dens_9_eq = lbfD.calc_Ni_eq(vel_vec,new_vel_dens_9,rho,m,c)
  old_vel_dens_9 = lbfD.relax_dens(tau,vel_dens_9_eq,new_vel_dens_9)




'''
zeros = np.zeros((xdim,ydim))
zeros[0] = 1
first_fill = np.swapaxes(zeros,0,1)
first_fill[0] = 0
first_fill[-1] = 0
ffswp = np.swapaxes(first_fill,0,1)
ffswproll = np.roll(ffswp,1,axis=0)
first_fill = np.swapaxes(ffswproll,1,0)
 '''
