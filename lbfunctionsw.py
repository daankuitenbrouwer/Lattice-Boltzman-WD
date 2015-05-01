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

#Initializes density array, width is number of gridpoints inside
#Orthogonal to pipewall, length number of points parallel to wall.
def init_dens(width,length):
    dens = np.ones((width,length,9))
    dens[0,:,:]=0
    dens[-1,:,:]=0
    return dens

def init_mask(width,length):
    mask = np.zeros((width,length),dtype="bool")
    mask[0,:]=1
    mask[-1,:]=1
    return mask

def init_test(width,length):
    dens = np.ones((width,length,9))
    dens[0,:,:]=0
    dens[-1,:,:]=0
    dens[:,0,:]=0
    dens[:,-1,:]=0
    return dens

#vel_dens_matrix must have 2 more sites 
#in both x and y direction to account for boundary
def move_densities(old_dens):
  temp_dens = np.zeros(np.shape(old_dens))
  temp_dens[:,:,0] = old_dens[:,:,0]
  # let appropriate densities move in 1st direction
  temp_dens[:,:,1] = np.roll(old_dens[:,:,1],1,axis=1) 
  # let densities move in 2nd direction, 3rd direction etc. 
  temp_dens[:,:,2] = np.roll(np.roll(old_dens[:,:,2],1,axis=1),-1,axis=0)
  temp_dens[:,:,3] = np.roll(old_dens[:,:,3],-1,axis=0)
  temp_dens[:,:,4] = np.roll(np.roll(old_dens[:,:,4],-1,axis=0),-1,axis=1)
  temp_dens[:,:,5] = np.roll(old_dens[:,:,5],-1,axis=1)
  temp_dens[:,:,6] = np.roll(np.roll(old_dens[:,:,6],-1,axis=1),1,axis=0)
  temp_dens[:,:,7] = np.roll(old_dens[:,:,7],1,axis=0)
  temp_dens[:,:,8] = np.roll(np.roll(old_dens[:,:,8],1,axis=0),1,axis=1)
  return temp_dens
  
#fix boundary conditions by switching appropriate
#densities   
def fix_bc(dens,mask):    
    dens[:,:,1][mask],dens[:,:,5][mask] = dens[:,:,5][mask],dens[:,:,1][mask]
    dens[:,:,2][mask],dens[:,:,6][mask] = dens[:,:,6][mask],dens[:,:,2][mask]
    dens[:,:,3][mask],dens[:,:,7][mask] = dens[:,:,7][mask],dens[:,:,3][mask]
    dens[:,:,4][mask],dens[:,:,8][mask] = dens[:,:,8][mask],dens[:,:,4][mask]
    return dens
    
#u is the average velocity  
def calc_u(dens,pres_grad):
  total_dens = np.sum(dens,axis=2) + 0.00000001
  
  return total_dens
  
def calc_Ni_eq(vel_vec,new_n_i,rho,m,c):
  n_i_eq = np.empty_like(new_n_i)
  wi_rho_invm = (rho/m)*[4./9.,1./9.,1./36.]
  n_i_eq[0] = wi_rho_invm[0]
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
