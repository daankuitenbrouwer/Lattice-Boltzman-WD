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
  
def fix_bc(dens,mask):    
    dens[:,:,1][mask],dens[:,:,5][mask] = dens[:,:,5][mask],dens[:,:,1][mask]
    dens[:,:,2][mask],dens[:,:,6][mask] = dens[:,:,6][mask],dens[:,:,2][mask]
    dens[:,:,3][mask],dens[:,:,7][mask] = dens[:,:,7][mask],dens[:,:,3][mask]
    dens[:,:,4][mask],dens[:,:,8][mask] = dens[:,:,8][mask],dens[:,:,4][mask]
    return dens
  
def find_vel_vec_at_pos_inc_presgrad(dens,xdim,ydim,pres_grad):
  print 'here'
  vel_vec = np.array([np.zeros((ydim,xdim)),np.zeros((ydim,xdim))])
  inv_Ni = 1./(np.sum(dens,axis=2) + 0.0000001)
  dir_matrix = np.transpose(np.array([[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]]))
  print np.shape(dir_matrix),np.shape(dens),'dir,dens'
  vel_vec = np.dot(dens,dir_matrix)
  
  print vel_vec,'velvec',np.shape(vel_vec),'shpvelvec',np.shape(inv_Ni),'shpaeinve'
  return vel_vec
  
def calc_Ni_eq(vel_vec,new_dens,rho,m,c):
  dens_eq = np.empty_like(new_dens)
  wi_rho_invm = (rho/m)*[1./9.,1./36.,4./9.]
  fourth_term = (vel_vec[0]**2 + vel_vec[1]**2)*(3./(2.*c**2))
  dens_eq[0] = wi_rho_invm[0] + fourth_term
  dir_matrix = [[0,0],[1,0],[1,1],[0,1],[-1,1],[-1,0],[-1,-1],[0,-1],[1,-1]]
  for i in range(8):
    dens_eq[i+1] = wi_rho_invm[i%2]*(1 + (3./c**2.)*(dir_matrix[i+1][0]*vel_vec[0] + dir_matrix[i+1][1]*vel_vec[1]) + (9./(2.*c**4.))*(  (vel_vec[0]**2)*(dir_matrix[i+1][0]**2) + 2*dir_matrix[i+1][0]*dir_matrix[i+1][1]*vel_vec[0]*vel_vec[1] + (vel_vec[1]**2)*(dir_matrix[i+1][1]**2)) - fourth_term  )
  return dens_eq
  
def relax_dens(tau,dens_eq,new_dens):
  final_dens = np.empty_like(new_dens)
  for i in range(9):
    final_dens[i] = (1 - (1./tau))*new_dens[i] + dens_eq[i]/tau # how do I sum per layer in the 3d matrix?
  return final_dens
  
  
  
  
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
