from utils import*;imp.reload(dsp)
from pyScatSpheres import qdot_sphere_array as qsa;imp.reload(qsa)
from pyScatSpheres import qdot_sphere_array_arb as qsa_arb;imp.reload(qsa_arb)
from pyScatSpheres import spherical_utils as spu
import scipy.special as spe

ka,n_p,nmax = 1,1.5,3
args = {'ka':ka,'kd_z':0,'kd_y':0,'kp':n_p,'nmax':nmax}


# qdot1_a0 = qsa_arb.QdotSphereArray(alpha=0 ,**args)
qdot1_a1 = qsa_arb.QdotSphereArray(alpha=45,opt2=1,**args)
# qdot1_a1_opt0 = qsa_arb.QdotSphereArray(alpha=45,opt2=1,**args)

phi_a = np.pi/2
ls = np.hstack([[l]*(2*l+1)for l in range(nmax+1)])
ms = np.hstack([list(np.arange(l+1))+list(np.flipud(np.arange(-l,0))) for l in range(nmax+1)])

# analytical zeros incidence
jl0,jl1   = np.array([spu.jn( l,np.array([ka,n_p*ka]))  for l in ls]).T
djl0,djl1 = np.array([spu.jnp(l,np.array([ka,n_p*ka]))  for l in ls]).T
hl0  = np.array([spu.hn1(l,ka)  for l in ls])
dhl0 = np.array([spu.hn1p(l,ka) for l in ls])
al_0 = (hl0*djl0-    dhl0*jl0)/(n_p*djl1*hl0-jl1*dhl0)
bl_0 = (jl1*djl0-n_p*djl1*jl0)/(n_p*djl1*hl0-jl1*dhl0)

Ylma = np.conj([spe.sph_harm(m,l,phi_a,qdot1_a1.alpha) for l,m in zip(ls,ms)])
clm = 4*np.pi*1j**ls*Ylma
# Ylm0 = np.conj([spe.sph_harm(0,l,phi_a,0) for l in ls])
# Ylma/Ylm0

ap1_0 = abs(qdot1_a1.ap-clm*al_0)
print(ap1_0.sum())#,ap1_0.imag.sum())
