import importlib as imp
import scipy.special as spe
import numpy as np, pandas as pd
from . import displayStandards as dsp #;imp.reload(dsp)
from . import glob_colors as colors
from . import spherical_utils as spu  #;imp.reload(spu)
from . import hard_sphere_base as hsb ;imp.reload(hsb)

class QdotSphereArray(hsb.HardSphereArrayBase):
    def solve(self,nmax=-1,copt=1,opt2=0,optT=0,v=1):
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
            self.vl=bl/cl
            self.ul=al/cl
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
        self.Nu = 2*N*nmax
        self.Nap = N*nmax
        T = np.zeros((2*N*nmax,2*N*nmax),dtype=complex)
        P = np.zeros((2*N*nmax,2*N*nmax),dtype=complex)
        if v:print(colors.blue+"...assembling coupling..."+colors.black)
        for p in range(N):
            if v>1:print(colors.yellow+'p=%d' %p+colors.black)
            idp1,idp2 = slice(p*nmax,(p+1)*nmax),slice((N+p)*nmax,(N+p+1)*nmax)
            P[idp1,idp1] =  np.diag(jl1)
            P[idp1,idp2] = -np.diag(hl0)
            P[idp2,idp1] =  np.diag(djl1*n_p)
            P[idp2,idp2] = -np.diag(dhl0)
            if copt==1:
                q0,q1 = np.arange(p),np.arange(p+1,N)
                q = np.hstack([q0,q1])                               #;print(q)
                kdpq = abs(kdp[p]-kdp[q])                            #;print(kdpq)
                theta_pq = np.hstack([[np.pi]*q0.size,[0]*q1.size])  #;print(theta_pq)
                aln0 = spu.get_aln(nmax-1,nmax-1, fz_q,kdpq,theta_pq)#;print(aln0)
                idqb = np.hstack([np.arange((N+q0)*nmax,(N+q0+1)*nmax) for q0 in q])
                idp1,idp2 = np.arange(p*nmax,(p+1)*nmax),np.arange((N+p)*nmax,(N+p+1)*nmax)
                # print(aln0.shape,idp1.shape,idqb.shape)# print(idp1,idqb)
                # print(np.ix_(idp1,idqb))
                T[np.ix_(idp1,idqb)] =  jl0[:,None] * aln0
                T[np.ix_(idp2,idqb)] = djl0[:,None] * aln0
            elif copt==2:
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
        if optT:
            # print('optT')
            j,i = np.meshgrid(np.arange(self.Nu),np.arange(self.Nu))
            T[i+self.Nap<j] = 0
            T[(i>=self.Nap) & (i+self.Nap<j+self.Nap)] = 0
        #### Solving
        self.T=T
        self.L=L
        self.P=P
        if v:print(colors.blue+"...solving %dx%d system..."%(self.Nu,self.Nu)+colors.black)
        self.P1=np.linalg.inv(P)
        cp = np.linalg.solve(P-(copt>0)*T,L)
        self.ap,self.bp = cp[:N*nmax],cp[N*nmax:]
        return self.ap,self.bp

    def get_cp0(self):
        cp = self.P1.dot(self.L)
        self.ap,self.bp = cp[:self.N*self.nmax],cp[self.N*self.nmax:]
        return self.ap,self.bp
    def get_cp2(self):
        T=self.T.copy()
        j,i = np.meshgrid(np.arange(self.Nu),np.arange(self.Nu))
        T[i+self.Nap<j] = 0
        T[(i>=self.Nap) & (i+self.Nap<j+self.Nap)] = 0
        bp0 = self.P1.dot(self.L)
        cp = (np.identity(self.Nu)+self.P1.dot(T)).dot(bp0)
        self.ap,self.bp = cp[:self.N*self.nmax],cp[self.N*self.nmax:]
        return self.ap,self.bp

    def get_cpn(self,n=2):
        T=self.T.copy()
        j,i = np.meshgrid(np.arange(self.Nu),np.arange(self.Nu))
        bp0 = self.P1.dot(self.L)
        T[i+self.Nap<j] = 0
        T[(i>=self.Nap) & (i+self.Nap<j+self.Nap)] = 0
        T = self.P1.dot(T)
        ap,bp = [],[]
        Tn = np.identity(self.Nu)
        for i in range(n):
            cp = Tn.dot(bp0)#;print(cp.shape)
            ap += [cp[:self.N*self.nmax]]
            bp += [cp[self.N*self.nmax:]]
            Tn = Tn.dot(T)
        return ap,bp

    def get_cpa(self):
        j,i = np.meshgrid(np.arange(self.Nu),np.arange(self.Nu))
        T=self.T.copy()
        T[i+self.Nap<j] = 0
        T[(i>=self.Nap) & (i+self.Nap<j+self.Nap)] = 0
        cp = np.linalg.solve(self.P-T,self.L)
        self.ap,self.bp = cp[:self.N*self.nmax],cp[self.N*self.nmax:]
        return self.ap,self.bp

    def set_Cp(self,cp):
        # print(self.N,cp)
        self.ap,self.bp = cp[:self.N*self.nmax],cp[self.N*self.nmax:]
        self.cp=self.bp
    def get_ffn(self,theta,bp):
        k,d_p,nmax,N = self.k,self.d_p,self.nmax,self.N
        ct = np.cos(theta)

        f = np.zeros(theta.shape,dtype=complex)
        for l in range(nmax):
            Yl = spe.sph_harm(0,l,0,theta)
            for p in np.arange(N):
                f += (-1J)**(l+1)*Yl*bp[p*nmax+l]*np.exp(-1J*self.kd_p[p]*ct)
        return f

    def compute_f(self,r,theta,phi,ftype='t',Gopt='',idp=None):
        ''' computes scattering amplitude f
        - r,theta,phi : np.ndarray each - coordinates
        - ftype : str - 't'(total), 's'(scattered), 'i'(incident),
        - Gopt : compute gradient(G), Flux(F)
        - idp : index of sphere to show (None or -1 => all)
        return :
            - Field
        '''
        print('...computing near field...')
        k,n_p,d_p,nmax,N = self.k,self.kp,self.d_p,self.nmax,self.N
        x,y,z = spu.sphere2cart(r,theta,phi)
        self._check_idp(idp) #;print(idp)

        #incident wave
        if ftype in 'ita':
            fi = np.zeros(r.shape,dtype=complex)
            gi = np.zeros(r.shape,dtype=complex)
            # plane wave at sphere idp : use local spherical decomposition
            if isinstance(idp,int):
                # print('incident field at p sphere')
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[idp])
                for l in range(nmax):
                    fi += 1J**l*(2*l+1)*spu.jn(l,k*r_p)*spu.Pl(l,theta_p)
                    gi += 1J**l*(2*l+1)*spu.jnp(l,k*r_p)*spu.Pl(l,theta_p)*k
                fi *= np.exp(1J*k*self.d_p[idp])
                gi *= np.exp(1J*k*self.d_p[idp])
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
        # idp = range(N)
        if ftype in 'sta':
            fs = np.zeros(r.shape,dtype=complex)
            gs = np.zeros(r.shape,dtype=complex)
            #outgoing field
            for p in range(N):
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
                idx_o = r_p>=self.ka     #;print(idx_o.shape)
                for l in range(nmax):
                    Yl = spe.sph_harm(0,l,0,theta_p)
                    fs[idx_o] += self.bp[p*nmax+l]*spu.hn1(l, k*r_p[idx_o])*Yl[idx_o]
                    gs[idx_o] += self.bp[p*nmax+l]*spu.hn1p(l,k*r_p[idx_o])*Yl[idx_o]
                    #gs SHOULD USE THE TRANSLATION TO GET RADIAL DERIVATIVE r_p
            #remove scattered field inside spheres
            for p in range(N):
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
                idx_i = r_p<self.ka
                fs[idx_i] = 0
                gs[idx_i] = 0
            #inside field
            for p in range(N):
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
                idx_i = r_p<self.ka     #;print(idx_i.shape)
                for l in range(nmax):
                    Yl = spe.sph_harm(0,l,0,theta_p)
                    fs[idx_i] += self.ap[p*nmax+l]*spu.jn( l,     n_p*k*r_p[idx_i])*Yl[idx_i]
                    gs[idx_i] += self.ap[p*nmax+l]*n_p*spu.jnp( l,n_p*k*r_p[idx_i])*Yl[idx_i]

        #total field
        if Gopt=='G':
            # print('Radial derivative');
            if   ftype=='t' :return gs+gi
            elif ftype=='s' :return gs
            elif ftype=='i' :return gi
            elif ftype=='a' :return gi,gs
        elif Gopt=='F':
            # print('Flux');
            if   ftype=='t' :return np.conj(fs+fi)*(gs+gi)
            elif ftype=='s' :return np.conj(fs)*gs
            elif ftype=='i' :return np.conj(fi)*gi
            elif ftype=='a' :return np.conj(fi)*gi,np.conj(fs)*gs
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


    def show_T(self):
        dsp.stddisp(im=[abs(self.T)],pOpt='im')
