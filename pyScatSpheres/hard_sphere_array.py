import importlib as imp
import scipy.special as spe
import numpy as np, pandas as pd
# import utils.displayStandards as dsp ;imp.reload(dsp)
import utils.glob_colors as colors
from . import spherical_utils as spu ;imp.reload(spu)
from . import hard_sphere_base as hsb ;imp.reload(hsb)

Pl = lambda l,theta: np.sqrt(4*np.pi/(2*l+1))*spe.sph_harm(0,l,0,theta)


class HardSphereArray(hsb.HardSphereArrayBase):
    def solve(self,nmax=-1,copt=1,v=1):
        if nmax>0:self.nmax=nmax+1
        N,ka,kd,kdp,nmax = self.N,self.ka,self.kd,self.kd_p,self.nmax
        fz_q = spu.hn1

        #preprocess
        ll = np.arange(nmax)
        bl = 1J**ll*(2*ll+1)*np.sqrt(4*np.pi/(2*ll+1))
        ul = -np.array([spu.jn(l,ka)/spu.hn1(l,ka) for l in range(nmax)])

        #### incident plane wave response
        L = np.tile(bl*ul,N)
        for p in range(N):
            idpl = p*nmax+np.arange(nmax) #;print('p=%d'%p)#;print('idp:' , idpl)
            L[idpl] *= np.exp(1J*kdp[p])

        #### coupling
        T = np.zeros((N*nmax,N*nmax),dtype=complex)
        if copt:
            ul = ul[:,None]
            if v:print("...assembling coupling...")
            for p in range(N):
                idpl = p*nmax+np.arange(nmax) #;print('p=%d'%p)#;print('idp:' , idpl)
                for q in range(p):
                    kdpq = kdp[p]-kdp[q]
                    aln0 = a_ln(nmax-1,nmax-1, fz_q,kdpq,np.pi)
                    T[p*nmax:(p+1)*nmax,q*nmax:(q+1)*nmax] = ul*aln0
                for q in range(p+1,N):
                    kdpq = kdp[q]-kdp[p]
                    aln1 = a_ln(nmax-1,nmax-1, fz_q,kdpq,0*np.pi)
                    T[p*nmax:(p+1)*nmax,q*nmax:(q+1)*nmax] = ul*aln1

        if v:print("...solving...")
        self.Cp = np.linalg.solve(np.identity(N*nmax)-copt*T,L)
        return self.Cp

    def set_Cp(self,Cp):self.Cp=Cp

    def compute_f(self,r,theta,phi,ftype='t',Gopt=0,idp=None):
        ''' computes scattering amplitude f
        - r,theta,phi : np.ndarray each - coordinates
        - ftype : str - 't'(total), 's'(scattered), 'i'(incident)
        - idp : index of sphere to show (None or -1 =>all)
        return :
        - np.real(E_ftype)
        '''
        k,d_p,nmax,N = self.k,self.d_p,self.nmax,self.N
        x,y,z = spu.sphere2cart(r,theta,phi)
        idp = self._check_idp(idp) #;print(idp)

        #incident wave
        if ftype in 'ita':
            fi = np.zeros(r.shape,dtype=complex)
            # plane wave at sphere idp : use local spherical decomposition
            if isinstance(idp,list):
                # print('incident field at p sphere')
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[idp[0]])
                for l in range(nmax):
                    fi += 1J**l*(2*l+1)*spu.jn(l,k*r_p)*Pl(l,theta_p)
                fi *= np.exp(1J*k*self.d_p[idp])
            # Otherwise actual plane wave propagating along z
            else:
                # print('full incident field')
                fi = np.exp(1J*k*z)

        #scattered fields
        if ftype in 'sta':
            fs = np.zeros(r.shape,dtype=complex)
            for p in idp:
                r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
                for l in range(nmax):
                    fs += self.Cp[p*nmax+l]*spu.hn1(l,k*r_p)*spe.sph_harm(0,l,0,theta_p)

        #remove values inside pth spheres
        for p in range(N):
            r_p,theta_p,phi_p = spu.cart2sphere(x,y,z-self.d_p[p])
            idx_p = r_p<self.ka/k
            if ftype in 'sta' : fs[idx_p] = 0
            if ftype in 'ita' : fi[idx_p] = 0

        #total field
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
                f += (-1J)**(l+1)*Yl*self.Cp[p*nmax+l]*np.exp(-1J*self.kd_p[p]*ct)
        return f

def a_ln(lmax,nmax, fz_q,r_d,theta_d):
    '''Scalar translational addition theorem coefficients
    - Nmax : order of expansion
    - fz_q : type of expansion (bessel or hankel)
    - r_d,theta_d : translation vector
    '''
    l = np.arange(lmax+nmax+2)
    Yq  = np.sqrt((2*l+1)/(4*np.pi))*spe.lpn(lmax+nmax+1,np.cos(theta_d))[0]
    # Yq = np.array([spe.sph_harm(0,q,0,theta_d) for q in range(lmax+nmax)])
    zq = np.array([fz_q(q,r_d) for q in range(lmax+nmax+2)])

    aln = np.zeros((lmax+1,nmax+1),dtype=complex)
    for l in range(lmax+1):
        # print('l=%d' %l)
        for n in range(nmax+1):
            q = np.arange(abs(l-n),n+l+1)   #;cput = time.time()
            Glnq = np.array([np.float(spu.gaunt(l,n,q, 0,0,0)) for q in q])
            # Glnq = spu.w3j(l,n,q, 0,0,0)
            aln[l,n] = 4*np.pi*np.sum((-1J)**(l-n-q)*zq[q]*Yq[q]*Glnq) #;print('n=%d' %n, q,zq[q],Yq[q],Glnq)
    return aln


def sweep_ka(df_name,N=2,nmax=7,kas=[0.5,1,2,4],kds=None,nkds=50,kdr=(2,10),v=1):
    ''' solve for a series of values kas,kds and save the Cp
    - df_name : name of data to create or update
    - kas : array of ka to solve
    - kds : array of kd to solve
    - kdr,nkds : kds=np.linspace(kdr[0]*ka,kdr[1]*ka,nkds) (only if kds=None)
    '''
    kas = np.array(kas)
    nkas = kas.size
    kds_var=not (isinstance(kds,np.ndarray) or isinstance(kds,list))
    cols = ['N','ka','kd','nmax','Cp','sigma']
    df = pd.DataFrame([],columns=cols)
    for i,ka in enumerate(kas):
        if v:print('...%d : solving ka=%.1f...' %(i,ka))
        if kds_var: kds = np.linspace(kdr[0]*ka,kdr[1]*ka,nkds)
        for j,kd in enumerate(kds):
            if kd>2*ka:
                s = HardSphereArray(N,ka,kd,nmax,solve=True) #s;print(s.Cp)
                di = dict(zip(cols,[s.N,s.ka,s.kd,s.nmax,s.Cp,s.get_s()])) #;print(di)
                df = df.append(di,ignore_index=True)

    try :
        df0 = pd.read_pickle(df_name)
        df.update(df0)
    except FileNotFoundError:
        print('creating new file')
    df.to_pickle(df_name)
    print(colors.green+'DataFrame saved : \n'+colors.yellow+df_name+colors.black)
    return df
