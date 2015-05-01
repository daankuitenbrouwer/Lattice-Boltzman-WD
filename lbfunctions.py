import numpy as np
import matplotlib.pyplot as plt

def image(matrix,name):
  plt.figure(1)
  plt.subplot(111)
  plt.imshow(matrix,cmap="Greys",interpolation="nearest")
  plt.title(name)
  plt.colorbar()
  plt.show()
  return

def update_pos_inc_boundary(old_n_i):#vel_dens_matrix must have 2 more sites in both x and y direction to account for boundary
  image(old_n_i[0],'old')
  new_n_i = np.empty_like(old_n_i)
  new_n_i[0] = old_n_i[0]
  #image(new_n_i[0],'0')
  # dir 1
  new_n_i[1] = np.roll(old_n_i[1],1,axis=1) 
  # dir 2
  new_n_i[2] = np.roll(np.roll(old_n_i[2],1,axis=1),-1,axis=0)
  #image(new_n_i[2],'2before')
  new_n_i[6][0] = new_n_i[2][0]
  new_n_i[2][0] = 0
  
  # dir 3
  new_n_i[3] = np.roll(old_n_i[3],-1,axis=0)
  new_n_i[7][0] = new_n_i[3][0]
  new_n_i[3][0] = 0
  
  # dir 4
  new_n_i[4] = np.roll(np.roll(old_n_i[4],-1,axis=0),-1,axis=1)
  new_n_i[8][0] = new_n_i[4][0]
  new_n_i[4][0] = 0
  
  # dir 5
  new_n_i[5] = np.roll(old_n_i[5],-1,axis=1)
  
  # dir 6
  new_n_i[6] += np.roll(np.roll(old_n_i[6],-1,axis=1),1,axis=0)
  #image(new_n_i[6],'6')
  new_n_i[2][-1] = new_n_i[6][-1]
  new_n_i[6][-1] = 0
  
  # dir 7
  new_n_i[7] += np.roll(old_n_i[7],1,axis=0)
  new_n_i[3][-1] = new_n_i[7][-1]
  new_n_i[7][-1] = 0
  
  # dir 8
  new_n_i[8] += np.roll(np.roll(old_n_i[8],1,axis=0),1,axis=1)
  new_n_i[4][-1] = new_n_i[8][-1]
  new_n_i[8][-1] = 0
  '''
  image(new_n_i[1],'1')
  image(new_n_i[2],'2')
  image(new_n_i[3],'3')
  image(new_n_i[4],'4')
  image(new_n_i[5],'5')
  image(new_n_i[6],'6')
  image(new_n_i[7],'7')
  image(new_n_i[8],'8')
   #'''
  return new_n_i
  
def find_vel_vec_at_pos_inc_presgrad(new_n_i,xdim,ydim,pres_grad):
  vel_vec = np.array([np.zeros((ydim,xdim)),np.zeros((ydim,xdim))])
  inv_Ni = 1./(new_n_i[0] + new_n_i[1] + new_n_i[2] + new_n_i[3] + new_n_i[4] + new_n_i[5] + new_n_i[6] + new_n_i[7] + new_n_i[8])
  vel_vec[0] = pres_grad + (new_n_i[1] + new_n_i[2] - new_n_i[4] - new_n_i[5] - new_n_i[6] + new_n_i[8])*inv_Ni
  vel_vec[1] = (new_n_i[2] + new_n_i[3] + new_n_i[4] - new_n_i[6] - new_n_i[7] - new_n_i[8])*inv_Ni
  image(vel_vec[0],'x')
  image(vel_vec[1],'y')
  return vel_vec
  
def calc_Ni_eq(vel_vec,new_n_i,rho,m,c):
  n_i_eq = np.empty_like(new_n_i)
  wi_rho_invm = (rho/m)*[4./9.,1./9.,1./36.]
  4th_term = (vel_vec[0]**2 + vel_vec[1]**2)*(3./(2.*c**2))
  n_i_eq[0] = wi_rho_invm[0] + 4th_term
  for ind,i in enumerate[2,-4,-6,8]:
    n_i_eq[abs(i)] = wi_rho_invm[1]*(1 + (3./c**2.)*(i/abs(i))*vel_vec[0] + (9./(2.*c**4.))*((-1)**ind)*vel_vec[1]*vel_vec[0] - (3./2.*c**2)*vel_vec[0]**2)
  for i in [1,-5]:
    n_i_eq[abs(i)] = wi_rho_invm[2]*(1 + (3./c**2)*(i/abs(i))*vel_vec[0] - (3./2.*c**2)*vel_vec[0]**2)
  n_i_eq[3] = wi_rho_invm[2]*(1 - (3./2.*c**2)*vel_vec[0]**2)
  n_i_eq[7] = wi_rho_invm[2]*(1 - (3./2.*c**2)*vel_vec[0]**2)
  return n_i_eq
  
def relax_dens(tau,n_i_eq,new_n_i):
  final_n_i = np.empty_like(new_n_i)
  for i in range(9):
    final_n_i[i] = (1 - (1./tau))*new_n_i[i] + n_i_eq[i]/tau # how do I sum per layer in the 3d matrix?
  return final_n_i
  
  
  
  
  '''
  #taking all horizontal into account
  tot_dens_matrix += vel_dens_matrix[0] + np.roll(vel_dens_matrix[1],1,axis=1) + np.roll(vel_dens_matrix[5],-1,axis=1) #no bc needed
  
  image(tot_dens_matrix,'dens_matrix')
  dens_matrix = np.swapaxes(dens_matrix,0,1)
  dens_matrix[-1] = 0
  dens_matrix = np.swapaxes(dens_matrix,1,0)
  dens_matrix[0] = 0
  dens_matrix[-1] = 0
  image(dens_matrix,'after bc')
  #direction 2
  dens_matrix = dens_matrix + 
  image(dens_matrix,'dens_matrix')
  dens_matrix[-1] = 0
  dens_matrix[1] += np.roll(dens_matrix[0],-1)
  dens_matrix[0] = 0
  image(dens_matrix,'after2')
  #direction 3
  dens_matrix = dens_matrix + np.roll(vel_dens_matrix[3],1,axis=0)
  dens_matrix[-1] = 0
  dens_matrix[1] += dens_matrix[0]
  dens_matrix[0] = 0
  image(dens_matrix,'after3')
  #direction 4
  dens_matrix += np.roll(np.roll(vel_dens_matrix[4],1,axis=0),-1,axis=1)
  dens_matrix = np.swapaxes(dens_matrix,0,1)
  dens_matrix[-1] = 0
  dens_matrix = np.swapaxes(dens_matrix,1,0)
  dens_matrix[-1] = 0
  dens_matrix[1] += np.roll(dens_matrix[0],1)
  dens_matrix[0] = 0
  image(dens_matrix,'after4')
  # to do: direction 678
  
  
  #print dens_matrix
  
  
  return dens_matrix
  
  #need to find out how to roll to the left.
   '''
