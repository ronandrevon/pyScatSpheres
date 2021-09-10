from utils import *                                 ;imp.reload(dsp)
# from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import utils as ut ;imp.reload(ut)
plt.close('all')


# df_name = 'data/df_qdotSpheresN2.pkl'
# kas0   = np.array([0.5,1,2,3,4,5,7,10,15,20,30])
# kdkas0 = np.array([2.1,3,5,10,25,20])
# kps0   = np.array([1.2,1.1,1.01,1.001,1.0001])
# kas,kdkas,kps = np.meshgrid(kas0,kdkas0,kps0)
# kas,kdkas,kps =kas.flatten(),kdkas.flatten(),kps.flatten()
# Ns=np.array([2]*kas.size,dtype=int)

ka,kp,kdka = 7,1.001,3
df_name = 'data/df_qdot_ka%d_kp%d_Ns.pkl' %(ka,abs(np.round(np.log10(kp-1))))
Ns   = np.array([2,3,4,6,8,10,15,20,30,45,50,70,100])

# ka,kp,kdka,Ns = 15,np.array([1.01,1.001,1.0001]),3,100
# df_name = 'data/df_qdot_ka%d_N_%d_np.pkl' %(ka,Ns)#bs(np.round(np.log10(kp-1))))


kas,kps,kdkas,Ns = ut.get_param_list([ka,kp,kdka,Ns])
df = ut.solve_set(df_name,kas,kps,kdkas,Ns,new=1,n=10)
