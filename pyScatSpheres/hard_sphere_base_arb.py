import scipy.special as spe
import numpy as np
from . import displayStandards as dsp 
from . import spherical_utils as spu 

class HardSphereArrayBaseArb():
	def __init__(self,ka=1,kd_z=2*np.ones((1,1)),kd_y=2*np.ones((1,1)),kp=np.inf,k=1,nmax=7,Cp=None,solve=True,copt=1):
        """
        - n,m,ka,kp,kd : nb spheres,wave number inside spheres, normalized radius, normalized distances
        - k : incident wave number
        - nmax : max included order expansion
        """
        self.k,self.nmax = k,nmax+1
        self.ka,self.kp,self.kd_z,self.kd_y= ka,kp,kd_z,kd_y
        self.kd_z*self.kd_y=N
        n,m=kd_z.shape[0],kd_y.shape[0]
        self.kd_pz=np.zeros(n)
        self.kd_pz[0]=kd_z[0]
        for i in range(1,n):
            self.kd_pz[i]=kd_z[i]+self.kd_pz[i-1] 
        self.d_pz  = self.kd_pz/k
        self.kd_py=np.zeros(m)
        self.kd_py[0]=kd_y[0]
        for i in range(1,m):
            self.kd_py[i]=kd_y[i]+self.kd_py[i-1] 
        self.d_py  = self.kd_py/k
        self.kd_p=np.sqrt(kd_pz**2+kd_py**2)
        self.d_p=np.sqrt(d_pz**2+d_py**2)

        if isinstance(Cp,np.ndarray):self.set_Cp(Cp)
        if solve:self.solve(copt=copt)

    def _params(self):
        ss = '$N=%d$, $n_{ref}=%g$, $ka=%.1f$ ,kd_x=[' %(self.kd_z*self.kd_y,self.kp, self.ka)
        for i in self.kd_z :
            pp='%.1f,' %i 
            ss+=pp
        ss+='], kd_y=['
        for i in self.kd_y :
            pp='%.1f,' %i 
            ss+=pp
        ss+=']'
        return ss

