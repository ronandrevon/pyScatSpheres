# import importlib as imp
import scipy.special as spe
import numpy as np, pandas as pd
from . import displayStandards as dsp #;imp.reload(dsp)
from . import glob_colors as colors
from . import spherical_utils as spu  #;imp.reload(spu)
from . import hard_sphere_base as hsb #;imp.reload(hsb)

class QdotSphereArray(hsb.HardSphereArrayBase):
    def solve(self,nmax=-1,copt=1,opt2=0,v=0):
        ''' Finds the unknown coefficients
        - nmax : max inlcuded order
        - copt : solve coupled problem
        - opt2 : get
        '''
        if nmax>0:self.nmax=nmax+1
        N,n_p,ka,kd,kdp,nmax = self.N,self.kp,self.ka,self.kd,self.kd_p,self.nmax
        fz_q = spu.hn1

        #preprocess,
        ll = np.arange(nmax)
        cl = 1J**ll*(2*ll+1)*np.sqrt(4*np.pi/(2*ll+1))
        jl0,jl1   = np.array([spu.jn( l,np.array([ka,n_p*ka]))  for l in ll]).T
        djl0,djl1 = np.array([spu.jnp(l,np.array([ka,n_p*ka]))  for l in ll]).T
        hl0  = np.array([spu.hn1(l,ka)  for l in ll])
        dhl0 = np.array([spu.hn1p(l,ka) for l in ll])
        if opt2:
            al = cl*(hl0*djl0-    dhl0*jl0)/(n_p*djl1*hl0-jl1*dhl0)
            bl = cl*(jl1*djl0-n_p*djl1*jl0)/(n_p*djl1*hl0-jl1*dhl0)
        else:
            al = cl*jl0
            bl = cl*djl0

        #### incident plane wave response
        L = np.hstack([np.tile(al,N),np.tile(bl,N)])
        for p in range(N):
            idpl = p*nmax+np.arange(nmax) #;print('p=%d'%p)#;print('idp:' , idpl)
            L[idpl]        *= np.exp(1J*kdp[p])
            L[idpl+N*nmax] *= np.exp(1J*kdp[p])
        if opt2:
            self.ap,self.bp = L[:N*nmax],L[N*nmax:]
            return self.ap,self.bp

        #### assembling
        T = np.zeros((2*N*nmax,2*N*nmax),dtype=complex)
        P = np.zeros((2*N*nmax,2*N*nmax),dtype=complex)
        if v:print("...assembling coupling...")
        for p in range(N):
            print(colors.yellow+'p=%d' %p+colors.black)
            idp1,idp2 = slice(p*nmax,(p+1)*nmax),slice((N+p)*nmax,(N+p+1)*nmax)
            P[idp1,idp1] =  np.diag(jl1)
            P[idp1,idp2] = -np.diag(hl0)
            P[idp2,idp1] =  np.diag(djl1*n_p)
            P[idp2,idp2] = -np.diag(dhl0)
            if copt:
                for q in range(p):
                    kdpq = kdp[p]-kdp[q]
                    idqb = slice((N+q)*nmax,(N+q+1)*nmax)
                    aln0 = spu.a_ln(nmax-1,nmax-1, fz_q,kdpq,np.pi)
                    T[idp1,idqb] =  jl0[:,None] * aln0
                    T[idp2,idqb] = djl0[:,None] * aln0
                for q in range(p+1,N):
                    idqb = slice((N+q)*nmax,(N+q+1)*nmax)
                    kdpq = kdp[q]-kdp[p]
                    aln1 = spu.a_ln(nmax-1,nmax-1, fz_q,kdpq,0)
                    T[idp1,idqb] =  jl0[:,None]* aln1
                    T[idp2,idqb] = djl0[:,None]* aln1
        #### Solving
        if v:print("...solving...")
        cp = np.linalg.solve(P-copt*T,L)
        self.ap,self.bp = cp[:N*nmax],cp[N*nmax:]

        return self.ap,self.bp

    def set_Cp(self,cp):
        # print(self.N,cp)
        self.ap,self.bp = cp[:self.N*self.nmax],cp[self.N*self.nmax:]


    def compute_f(self,r,theta,phi,ftype='t',Gopt=0,idp=None):
        ''' computes scattering amplitude f
        - r,theta,phi : np.ndarray each - coordinates
        - ftype : str - 't'(total), 's'(scattered), 'i'(incident),
        - Gopt : compute gradient
        - idp : index of sphere to show (None or -1 => all)
        return :
        - np.real(E_ftype)
        '''
        k,n_p,d_p,nmax,N = self.k,self.kp,self.d_p,self.nmax,self.N
        x,y,z = spu.sphere2cart(r,theta,phi)
        idp = self._check_idp(idp) #;print(idp)

        #incident wave
        if ftype in 'ita':
            fi = np.zeros(r.shape,dtype=complex)
            gi = np.zeros(r.shape,dtype=complex)
            # plane wave at sphere idp : use local spherical decomposition
            if isinstance(idp,list):
                # print('incident field at p sphere')
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[idp[0]])
                for l in range(nmax):
                    fi += 1J**l*(2*l+1)*spu.jn(l,k*r_p)*spu.Pl(l,theta_p)
                fi *= np.exp(1J*k*self.d_p[idp])
            # Otherwise actual plane wave propagating along z
            else:
                # print('full incident field')
                fi = np.exp(1J*k*z)
                gi = 1J*k*np.cos(theta)*np.exp(1J*k*z)
            #remove incident field inside pth spheres
            for p in range(N):
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
                idx_i = r_p<self.ka
                fi[idx_i] = 0
                gi[idx_i] = 0

        #scattered fields
        if ftype in 'sta':
            fs = np.zeros(r.shape,dtype=complex)
            gs = np.zeros(r.shape,dtype=complex)
            #outgoing field
            for p in idp:
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
                idx_o = r_p>=self.ka     #;print(idx_o.shape)
                for l in range(nmax):
                    Yl = spe.sph_harm(0,l,0,theta_p)
                    fs[idx_o] += self.bp[p*nmax+l]*spu.hn1(l, k*r_p[idx_o])*Yl[idx_o]
                    gs[idx_o] += self.bp[p*nmax+l]*spu.hn1p(l,k*r_p[idx_o])*Yl[idx_o]
            #remove scattered field inside spheres
            for p in range(N):
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
                idx_i = r_p<self.ka
                fs[idx_i] = 0
                gs[idx_i] = 0
            #inside field
            for p in idp:
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
                idx_i = r_p<self.ka     #;print(idx_i.shape)
                for l in range(nmax):
                    Yl = spe.sph_harm(0,l,0,theta_p)
                    fs[idx_i] += self.ap[p*nmax+l]*spu.jn( l,     n_p*k*r_p[idx_i])*Yl[idx_i]
                    gs[idx_i] += self.ap[p*nmax+l]*n_p*spu.jnp( l,n_p*k*r_p[idx_i])*Yl[idx_i]

        #total field
        if Gopt:
            # print('     return Gradient');
            if   ftype=='t' :return gs+gi
            elif ftype=='s' :return gs
            elif ftype=='i' :return gi
            elif ftype=='a' :return gi,gs
        else:
            if   ftype=='t' :return fs+fi
            elif ftype=='s' :return fs
            elif ftype=='i' :return fi
            elif ftype=='a' :return fi,fs

    def get_ff(self,theta):
        k,d_p,nmax,N = self.k,self.d_p,self.nmax,self.N
        ct = np.cos(theta)

        f = np.zeros(theta.shape,dtype=complex)
        for l in range(nmax):
            Yl = spe.sph_harm(0,l,0,theta)
            for p in np.arange(N):
                f += (-1J)**(l+1)*Yl*self.bp[p*nmax+l]*np.exp(-1J*self.kd_p[p]*ct)
        return f
