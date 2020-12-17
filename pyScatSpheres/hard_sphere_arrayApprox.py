# import importlib as imp
import scipy.special as spe
import numpy as np
from . import displayStandards as dsp  #;imp.reload(dsp)
from . import spherical_utils as spu  #;imp.reload(spu)
from . import hard_sphere_base as hsb #;imp.reload(hsb)

Pl = lambda l,theta: np.sqrt(4*np.pi/(2*l+1))*spe.sph_harm(0,l,0,theta)

class HardSphereArray(hsb.HardSphereArrayBase):
    def solve(self,copt=1,v=0):
        N,ka,kd,nmax = self.N,self.ka,self.kd,self.nmax
        kdp = self.kd_p

        ll = np.arange(nmax)
        bl = 1J**ll*(2*ll+1)
        ul = -np.array([spu.jn(l,ka)/spu.hn1(l,ka) for l in range(nmax)])
        cpl0 = np.sqrt(4*np.pi/(2*ll+1))*bl*ul
        Yl = np.array([spe.sph_harm(0,l,0,[0,np.pi]) for l in range(nmax)])

        f = np.zeros(2,dtype=complex)#np.exp(1j*kdp)
        for l in range(nmax):
            f += (-1J)**(l+1)*cpl0[l]*Yl[l]

        self.L,self.T = np.zeros((N),dtype=complex),np.zeros((N,N),dtype=complex)
        self.Cp = np.zeros((N),dtype=complex)
        if copt:
            for p in range(self.N):
                fp = f*np.exp(1J*kdp[p])  #;print('p=%d, kdp=%.1f, f_p=' %(p,kdp[p]),fp)

                # print(colors.yellow+'p=%d'%p+colors.black)
                q = np.setdiff1d(np.arange(N),p)
                kdpq = np.abs(kdp[p]-kdp[q])    #; print('f1:',2*f/kdpl)
                theta_qp = np.array(q>p,dtype=int) #0,0,.,pi,pi
                theta_pq = np.array(q<p,dtype=int) #pi,pi,.,0,0
                gp1 = np.sum(fp[theta_pq])/kdpq
                gp2 = np.sum(fp[theta_qp])/kdpq
                # incident field scattering term
                self.L[p] = gp2
                # rear spheres scattering terms
                r  = np.arange(p)
                self.T[p,r] = np.exp(-1J*kdp[r])*gp1  #;print('q:',q)
                # forward spheres scattering terms
                r  = np.arange(p+1,N)
                self.T[p,r] = np.exp(+1J*kdp[r])*gp2  #;print('q:',q)

            if v:print("...solving...")
            self.Cp = np.linalg.solve((N-1)*np.identity(N)-self.T,self.L)
        return self.Cp

    def _cpl(self):
        k,ka,kdp,nmax,N = self.k,self.ka,self.kd_p,self.nmax,self.N
        ll = np.arange(nmax)
        bl = 1J**ll*(2*ll+1)
        ul = -np.array([spu.jn(l,ka)/spu.hn1(l,ka) for l in range(nmax)])
        cpl0 = np.sqrt(4*np.pi/(2*ll+1))*bl*ul
        self.cpl = np.zeros((N*nmax),dtype=complex)
        for l in range(nmax):
            for p in range(self.N):
                q0,q1 = np.arange(p-1),np.arange(p+1,self.N)
                self.cpl[p*nmax+l]=cpl0[l]*np.exp(1J*kdp[p])*(
                    1 + np.sum(self.Cp[q0]*np.exp(-1J*kdp[q0])) +
                    (-1)**l*np.sum(self.Cp[q1]*np.exp(+1J*kdp[q1]))
                    )
        return self.cpl


    def get_ff(self,theta=np.linspace(0,np.pi,361)):
        '''far field scattering amplitude'''
        k,ka,kdp,nmax,N = self.k,self.ka,self.kd_p,self.nmax,self.N
        ct = np.cos(theta)

        ll = np.arange(nmax)
        bl = 1J**ll*(2*ll+1)
        ul = -np.array([spu.jn(l,ka)/spu.hn1(l,ka) for l in range(nmax)])
        cpl0 = np.sqrt(4*np.pi/(2*ll+1))*bl*ul

        f = np.zeros(theta.shape,dtype=complex)
        for l in range(nmax):
            Yl = spe.sph_harm(0,l,0,theta)
            for p in range(self.N):
                q0,q1 = np.arange(p-1),np.arange(p+1,self.N)
                f+=(-1J)**(l+1)*cpl0[l]*np.exp(1J*kdp[p])*(
                    Yl +
                    np.sum(self.Cp[q0]*np.exp(-1J*kdp[q0]))*Yl +
                    np.sum(self.Cp[q1]*np.exp(+1J*kdp[q1]))*np.flipud(Yl)
                    )*np.exp(-1J*kdp[p]*ct)
        return f
