import numpy as np

import lbfunctionsD as lbfD
import lbfunctionsw


#define shape grid and pressure gradient
xdim,ydim = 10,4
pres_grad = 0.01
timesteps = 20
tau = 0.001
rho,m,c = 1,1,1

#define initial velocity density set
first_fill = np.ones((ydim,xdim))
first_fill[0,:] = 0
first_fill[-1,:] = 0
'''
ones[0] = 0
ones[-1] = 0
ones = np.swapaxes(ones,0,1)
ones[0] = 0
ones[-1] = 0
first_fill = np.swapaxes(ones,0,1)
#'''
old_n_i = np.array([first_fill,first_fill,first_fill,first_fill,first_fill,first_fill,first_fill,first_fill,first_fill])

# Actual algorithm
for tmstp in range(timesteps):
  print tmstp
  new_n_i = lbfunctionsw.move_densities(old_n_i)#lbfD.update_pos_inc_boundary(old_n_i)
  vel_vec = lbfD.find_vel_vec_at_pos_inc_presgrad(new_n_i,xdim,ydim,pres_grad)
  n_i_eq = lbfD.calc_Ni_eq(vel_vec,new_n_i,rho,m,c)
  old_n_i = lbfD.relax_dens(tau,n_i_eq,new_n_i)




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
