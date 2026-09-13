
# coding: utf-8

# In[1]:


from htmd.ui import *


# In[4]:


# coordinates
m = Molecule('x.pdb')
m.read('output.coor')
m.filter('not resid 264')
m.write('x_out.coor', type='coor')


# In[5]:


# velocities
m_orig = Molecule('x.pdb')
vel = np.fromfile('output.vel', dtype=np.float64).reshape(-1, 3)
keep = m_orig.atomselect('not resid 264')
kept_indices = np.where(keep)[0]
vel_filtered = vel[kept_indices]
vel_filtered.astype(np.float64).tofile('x_out.vel')
