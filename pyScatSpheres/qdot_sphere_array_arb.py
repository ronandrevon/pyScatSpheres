import importlib as imp
import scipy.special as spe
import numpy as np, pandas as pd
from . import displayStandards as dsp #;imp.reload(dsp)
from . import glob_colors as colors
from . import spherical_utils as spu       ;imp.reload(spu)
from . import hard_sphere_base_arb as hsba ;imp.reload(hsba)
from . import harmo_sphe as hs             ;imp.reload(hs)
import time

class QdotSphereArray(hsba.HardSphereArrayBaseArb):
    def solve(self,nmax=-1,copt=1,opt2=1,v=1):
        print(colors.red,"Arbitrary",colors.black)
        if nmax>0:self.nmax=nmax+1
        if self.N==1:copt=0
        if opt2:self.solve2(copt)
        else:self.solve1()

    def solve2(self,copt=0,v=1):
        ''' Finds the unknown coefficients alternative assembling '''
        print(colors.red,"alternative assembling",colors.black)
        k=1
        N,n_p,ka = self.N,self.kp,self.ka,
        # kd_z,kd_y = self.kd_z,self.kd_y
        kdpz,kdpy,nmax = self.kd_pz,self.kd_py,self.nmax
        ka = ka[0]
        #### ul,vl : nmax x 1 (if all ka,n_p identical)
        ll = np.arange(nmax)
        jl0,jl1   = np.array([spu.jn( l,np.array([ka,n_p*ka]))  for l in ll]).T
        djl0,djl1 = np.array([spu.jnp(l,np.array([ka,n_p*ka]))  for l in ll]).T
        hl0  = np.array([spu.hn1(l,ka)  for l in ll])
        dhl0 = np.array([spu.hn1p(l,ka) for l in ll])
        # print(hl0.shape,djl0.shape)
        self.ul = (hl0*djl0-    dhl0*jl0)/(n_p*djl1*hl0-jl1*dhl0)
        self.vl = (jl1*djl0-n_p*djl1*jl0)/(n_p*djl1*hl0-jl1*dhl0)

        ###----------------------------------------------------------------------------------------
        #   Construction du second membre L
        ###----------------------------------------------------------------------------------------
        ls = np.hstack([[l]*(2*l+1)for l in range(nmax)])
        ms = np.hstack([list(np.arange(l+1))+list(np.flipud(np.arange(-l,0))) for l in range(nmax)])
        phi_a = np.pi/2
        Ylm_a = [spe.sph_harm(m,l,phi_a,self.alpha) for l,m in zip(ls,ms)]
        self.clm = 4*np.pi*1J**ls*np.conj(Ylm_a)

        #### size : nmax**2*N
        self.upl = np.tile(self.ul[ls],N)
        self.vpl = np.tile(self.vl[ls],N)
        self.kdp,self.theta_p,self.phi_p=spu.cart2sphere(0*kdpz,kdpy,kdpz)
        self.xi_p = np.sin(self.theta_p)*np.sin(self.phi_p)*np.sin(self.alpha)+np.cos(self.theta_p)*np.cos(self.alpha)
        self.ep = np.exp(1J*self.kdp*self.xi_p)
        # print(self.clm.shape,ep.shape)
        self.L0 = np.outer(self.ep,self.clm).flatten()
        self.L = np.hstack([self.upl*self.L0,self.vpl*self.L0])

        ###----------------------------------------------------------------------------------------
        #   Construction de T
        ###----------------------------------------------------------------------------------------
        self.Nap = N*nmax**2
        if copt:
            Tu=np.zeros((self.Nap,)*2,dtype=complex)
            Tv=np.zeros((self.Nap,)*2,dtype=complex)
            i = 0
            for p in range(N):
                print(colors.blue+'p=%d' %p+colors.black)
                q = np.hstack([np.arange(p),np.arange(p+1,N)])
                idx = np.hstack([q0*nmax**2+np.arange(nmax**2) for q0 in q])
                kd_pq    = np.sqrt((kdpz[q]-kdpz[p])**2 + (kdpy[q]-kdpy[p])**2) #;print(kd_pq)
                theta_pq = np.arccos((kdpz[q]-kdpz[p])/kd_pq)                   #;print(np.rad2deg(theta_pq))
                phi_pq   = -np.arctan2((kdpy[q]-kdpy[p]),0)                     #;print(np.rad2deg(phi_pq))
                # phi_pq   = np.pi/2*np.ones(kd_pq.shape)
                for l,m in zip(ls,ms):
                    # print(colors.yellow+'l=%d,m=%d' %(l,m)+colors.black)
                    #a_lmmunu : N-1 x nmax**2
                    a_lmnumu = spu.get_almnumu(l,m,kd_pq/k,theta_pq,phi_pq, Nmax=nmax)
                    # a_numulm = spu.get_anumulm(l,m,kd_pq/k,theta_pq,phi_pq, Nmax=nmax)
                    # a_lmnumu = np.ones((N-1, nmax**2))
                    a_lmnumu = a_lmnumu.flatten()
                    Tu[i,idx] = a_lmnumu*self.ul[l]
                    Tv[i,idx] = a_lmnumu*self.vl[l]
                    i+=1
            T = np.vstack([Tu,Tv])
            self.T = np.hstack([0*T,T])

            ###----------------------------------------------------------------------------------------
            #   Solving
            ###----------------------------------------------------------------------------------------
            if v:print(colors.blue+"...solving..."+colors.black)
            I = np.identity(2*N*nmax**2)
            cp = np.linalg.solve(I-self.T,self.L)
        else:
            cp = self.L

        self.ap,self.bp = cp[:self.Nap],cp[self.Nap:]
        return self.ap,self.bp

    def compute_f(self,r,theta,phi,ftype='t',Gopt='',idp=None):
        ''' computes scattering amplitude f
        - r,theta,phi : np.ndarray each - coordinates
        - ftype : str - 't'(total), 's'(scattered), 'i'(incident),
        - Gopt : compute gradient(G), Flux(F)
        - idp : index of sphere to show (None or -1 => all)
        return :
            - Field
        '''
        k,n_p,d_pz,d_py,nmax,N = self.k,self.kp,self.d_pz,self.d_py,self.nmax,self.N
        x,y,z = spu.sphere2cart(r,theta,phi)
        idp = self._check_idp(idp) #;print(idp)
        # print(y.min(),y.max())

        #### incident wave -------------------------------------------------------------------------
        if ftype in 'ita':
            fi = np.zeros(r.shape,dtype=complex)
            gi = np.zeros(r.shape,dtype=complex)
            # plane wave at sphere idp : use local spherical decomposition
            fi = np.exp(1J*k*(z*np.cos(self.alpha)+y*np.sin(self.alpha)))
            # gi = 1J*k*np.cos(theta)*fi
            #remove incident field inside pth spheres
            for p in idp:
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[p],z-self.d_pz[p])
                idx_i = r_p<self.ka[p]
                fi[idx_i] = 0
                # gi[idx_i] = 0

        #### scattered fields --------------------------------------------------------------------
        ls = np.hstack([[l]*(2*l+1)for l in range(self.nmax)])
        ms = np.hstack([list(np.arange(l+1))+list(np.flipud(np.arange(-l,0))) for l in range(self.nmax)])
        if ftype in 'sta':
            fs = np.zeros(r.shape,dtype=complex)
            # gs = np.zeros(r.shape,dtype=complex)
            #### outgoing field
            for p in idp:
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[p],z-self.d_pz[p])
                idx_o = r_p>=self.ka[p]     #;print(idx_o.shape)
                # print(idx_o)
                # A=hs.harmonique(nmax-1,theta_p,'tab',phi_p)
                i = p*self.nmax**2
                for l,m in zip(ls,ms):
                    hl = spu.hn1(l, k*r_p[idx_o])
                    Ylm=spe.sph_harm(m,l,phi_p[idx_o],theta_p[idx_o])
                    fs[idx_o] += self.bp[i]*hl*Ylm
                    i+=1
            #remove scattered field inside the q spheres
            for p in idp:# range(N):
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[p],z-self.d_pz[p])
                idx_i = r_p<self.ka[p]
                fs[idx_i] = 0
                # gs[idx_i] = 0

            #### inside field
            # A=hs.harmonique(nmax-1,theta_p,'tab',phi_p)
            for p in idp:#range(N):
                i = p*self.nmax**2
                # print(colors.blue+'p=%d' %p+colors.black)
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[p],z-self.d_pz[p])
                idx_i = r_p<self.ka[p]     #;print(idx_i.shape)
                for l,m in zip(ls,ms):
                    # print(colors.yellow+'i=%d,l=%d,m=%d' %(i,l,m)+colors.black)
                    jl = spu.jn( l, n_p*k*r_p[idx_i])
                    Ylm=spe.sph_harm(m,l,phi_p[idx_i],theta_p[idx_i])
                    # fs[idx_i] += self.ap[p*nmax**2+i]*jl*Ylm
                    fs[idx_i] += self.ap[i]*jl*Ylm
                    i+=1


        ### total field -----------------------------------------------------------------------------------------
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



    def solve1(self,v=1):
        ''' Finds the unknown coefficients'''
        N,n_p,ka,kd_z,kd_y,kdpz,kdpy,nmax = self.N,self.kp,self.ka,self.kd_z,self.kd_y,self.kd_pz,self.kd_py,self.nmax

        print(colors.red,"Non linéaire",colors.black)

        alpha=self.alpha

        ###----------------------------------------------------------------------------------------
        #   Construction du second membre L
        ###----------------------------------------------------------------------------------------
        #Initialisation
        # Ylm_a = [spu.sph_harm(l,m,np.pi/2,alpha) for m in idm]
        # Ylm_a = hs.harmonique(nmax-1,self.alpha,'tab',np.pi/2)
        L=np.zeros(2*(N*(nmax)**2),dtype=complex) #cf calcul de ma feuilleeeeee trop dur les maths ahhaha grrr
        self.ap0 = np.zeros((N*(nmax)**2),dtype=complex)
        for p in range(N) : #pas faire boucle sur sur boule mais pour l'instant on va faire ca sinon embrouiller moi tete
            m=0
            r_p,theta_p,phi_p = spu.cart2sphere(0,self.kd_py[p],self.kd_pz[p])
            for l in range(nmax) :

                    #Calcul coeff Cl
                    # cl = 1J**l*(2*l+1)*np.sqrt(4*np.pi/(2*l+1))#c depend que de l

                    #Calcul des fct de bessel
                    jl0=spu.jn(l,ka[p]) #j'espere que c'est bon ecrit comme ca ? peut etre c'est juste ka???????
                    djl0=spu.jnp(l,ka[p])
                    exp1=np.exp(1J*self.kdp[p]*self.xi_p[p]) #ca ca dépend que de N

                    idm = np.arange(2*l+1)
                    ims = np.hstack([np.arange(l+1),np.flipud(np.arange(-l,0))])
                    Ylm_a = [spe.sph_harm(m,l,phi_a,self.alpha) for m in ims]
                    clm = 4*np.pi*1J**l*np.conj(Ylm_a)#[l,idm]#c depend que de l
                    L[m+idm+p*nmax**2            ]=clm*exp1*jl0
                    L[m+idm+N*(nmax)**2+p*nmax**2]=clm*exp1*djl0
                    self.ap0[m+idm+p*nmax**2 ]    =clm*exp1
                    #Calcul des exp
                    zeta=np.sin(theta_p)*np.sin(phi_p)*np.sin(alpha)+np.cos(theta_p)*np.cos(alpha)
                    exp1=np.exp(1J*r_p*self.k*zeta) #ca ca dépend que de N


                    #Assemblage matrice
                    # !!L[m+p*nmax**2:m+2*l+1+p*nmax**2]=cl*exp1*jl0
                    # !!L[m+N*(nmax)**2+p*nmax**2:m+2*l+1+N*(nmax)**2+p*nmax**2]=cl*exp1*djl0 #hihi j'espere que les indices sont bons sinon je pleure
                    # L[m+p*nmax**2]=cl*exp1*jl0
                    # L[m+N*(nmax)**2+p*nmax**2]=cl*exp1*djl0
                    m=m+2*l+1

        ###----------------------------------------------------------------------------------------
        #   Construction de P
        ###----------------------------------------------------------------------------------------
        #Initialisation
        P = np.zeros((2*(N*(nmax)**2),2*(N*(nmax)**2)),dtype=complex)
        #print(P.shape)
        if v:print(colors.blue+"...assembling coupling..."+colors.black)
        m=0
        for p in range(N):
            if v>1:print(colors.yellow+'p=%d' %p+colors.black)
            #calcul des fonctions bessel hankel
            ll = np.arange(nmax)
            jl0,jl1   = np.array([spu.jn( l,np.array([ka[p],n_p*ka[p]]))  for l in ll]).T
            djl0,djl1 = np.array([spu.jnp(l,np.array([ka[p],n_p*ka[p]]))  for l in ll]).T
            hl0  = np.array([spu.hn1(l,ka[p])  for l in ll])
            dhl0 = np.array([spu.hn1p(l,ka[p]) for l in ll])
            for l in range(nmax):
                #remplissage de la matrice
                P[m:m+2*l+1,m:m+2*l+1]=np.diag([jl1[l]]*(2*l+1)) #en haut a gauche EST CE QUE C'EST INCIDE P OU L ????????
                P[m:m+2*l+1,m+N*(nmax)**2:m+2*l+1+N*(nmax)**2]=-np.diag([hl0[l]]*(2*l+1))#en haut a droite
                P[m+N*(nmax)**2:m+2*l+1+N*(nmax)**2,m:m+2*l+1]=np.diag([djl1[l]*n_p]*(2*l+1)) #en bas a gauche
                P[m+N*(nmax)**2:m+2*l+1+N*(nmax)**2,m+N*(nmax)**2:m+2*l+1+N*(nmax)**2]=-np.diag([dhl0[l]]*(2*l+1)) #en bas a droite

                m=m+2*l+1

        ###----------------------------------------------------------------------------------------
        #   Construction de T
        ###----------------------------------------------------------------------------------------

        ### Construction de T1 et T2 ----------------------------------------------------------
        T1=np.zeros(((N*(nmax)**2),(N*(nmax)**2)),dtype=complex)
        T2=np.zeros(((N*(nmax)**2),(N*(nmax)**2)),dtype=complex)
        fz_q = spu.hn1

        #p_indic=0
        #l_indic=0
        #m_indic=0
        ll = np.arange(nmax)
        for p in range(N):
            jl0  = np.array([spu.jn( l,ka[p])  for l in ll])
            djl0 = np.array([spu.jnp(l,ka[p])  for l in ll])
            # jl0,jl1   = np.array([spu.jn( l,np.array([ka[p],n_p*ka[p]]))  for l in ll]).T
            # djl0,djl1 = np.array([spu.jnp(l,np.array([ka[p],n_p*ka[p]]))  for l in ll]).T
            # hl0  = np.array([spu.hn1(l,ka[p])  for l in ll])
            # dhl0 = np.array([spu.hn1p(l,ka[p]) for l in ll])
            for l in range (nmax):
                for m in range (2*l+1):
                    vect=np.zeros((N*nmax**2),dtype=complex)
                    #c'est pour le coeff a la ligne lm de l'atome p
                    # vect=np.zeros((1,N*nmax**2),dtype=complex)
                    #coeff de la ligne qu'on regarde
                    # for q in list(range(p))+list(range(p+1,N)):
                    for q in range(N):
                        for nu in range(nmax):
                            for mu in range(2*nu+1):
                                """
                                print("\n")
                                print(colors.red,"indice m : ",colors.black)
                                print(p*nmax**2+(l)**2+m)
                                print(colors.green,"indice mu :",colors.black)
                                print(q*nmax**2+(nu)**2+mu)
                                """
                                if (q!=p):
                                    #print(q*nmax**2+(nu)**2+mu)
                                    x=kdpz[p]-kdpz[q]
                                    y=kdpy[p]-kdpy[q]
                                    kdpq = np.sqrt(x**2+y**2)
                                    theta=np.arccos((kdpz[q]-kdpz[p])/kdpq)
                                    phi=np.arctan2(-y,0)
                                    if np.isnan(theta):print(kdpq,kdpz[q],kdpz[p],theta)
                                    #theta=np.arctan2(y,x)
                                    #else:
                                        #theta=np.arctan2(y,x)
                                    #if (kdpz[p]-kdpz[q]<0):
                                        #theta+=np.pi/2
                                    if (q!=p):
                                        """
                                        print("boule ref :")
                                        print(p)
                                        print("boule actuelle :")
                                        print(q)
                                        print("theta :")
                                        print(theta)
                                        """
                                        m1=m
                                        mu1=mu
                                        if (m>l):
                                            m1=l-m
                                        if (mu>nu):
                                            mu1=nu-mu
                                        
                                        a=spu.a_lmnp(l,m1,nu,mu1,fz_q,kdpq,theta,phi)
                                        """
                                        print("m1:")
                                        print(m1)
                                        print("mu1")
                                        print(mu1)
                                        print("a :")
                                        print(a)
                                        """

                                        vect[0,q*nmax**2+(nu)**2+mu]=a

                                    #print(q*nmax**2+(nu)**2+mu)
                        #p*nmax**2 c'est pour passer d'un p a l'autre
                        #2*l+1 c'est pour passer d'un l a l'autre

                        T1[p*nmax**2+(l)**2+m,:]=vect*jl0[l]
                        #if((mu==0)&(m==0)):print(T1[p*nmax**2+(l)**2+m,:])
                        T2[p*nmax**2+(l)**2+m,:]=vect*djl0[l]




        ### Assemblage de T ---------------------------------------------------------------
        T=np.zeros((2*(N*(nmax)**2),2*(N*(nmax)**2)),dtype=complex)

        T[0:N*nmax**2,N*nmax**2:2*N*nmax**2]=T1
        T[N*nmax**2:2*N*nmax**2,N*nmax**2:2*N*nmax**2]=T2

        ###----------------------------------------------------------------------------------------
        #   Résolution
        ###----------------------------------------------------------------------------------------
        #print(colors.red,"pas linéaire\n",colors.black)

        #dsp.stddisp(im=[np.log10(abs(T))], pOpt='im',cmap='Spectral',imOpt='c')
        #dsp.plt.show()
        idx=np.array([0,1,4,5,8,9])
        idx1 = np.hstack([idx,idx+12])

        #V=(P-T)[np.ix_(idx1,idx1)]
        #print(V)


        if v:print(colors.blue+"...solving..."+colors.black)
        cp = np.linalg.solve(P-T,L)
        self.ap,self.bp = cp[:N*nmax**2],cp[N*nmax**2:]
        # self.V=V
        self.L=L
        self.P=P
        self.T=T

        return self.ap,self.bp


    def compute_f2(self,r,theta,phi,ftype='t',Gopt='',idp=None):
        ''' computes scattering amplitude f
        - r,theta,phi : np.ndarray each - coordinates
        - ftype : str - 't'(total), 's'(scattered), 'i'(incident),
        - Gopt : compute gradient(G), Flux(F)
        - idp : index of sphere to show (None or -1 => all)
        return :
            - Field
        '''
        k,n_p,d_pz,d_py,nmax,N = self.k,self.kp,self.d_pz,self.d_py,self.nmax,self.N
        x,y,z = spu.sphere2cart(r,theta,phi)
        idp = self._check_idp(idp) #;print(idp)
        # print(y.min(),y.max())

        ### incident wave -------------------------------------------------------------------------
        if ftype in 'ita':
            fi = np.zeros(r.shape,dtype=complex)
            gi = np.zeros(r.shape,dtype=complex)
            # plane wave at sphere idp : use local spherical decomposition
            if 0:#isinstance(idp,int):
                # print('incident field at p sphere')
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[idp],z-self.d_pz[idp])
                for l in range(nmax):
                    fi += 1J**l*(2*l+1)*spu.jn(l,k*r_p)*spu.Pl(l,theta_p)
                    gi += 1J**l*(2*l+1)*spu.jnp(l,k*r_p)*spu.Pl(l,theta_p)*k
                fi *= np.exp(1J*k*self.d_pz[idp])
                gi *= np.exp(1J*k*self.d_pz[idp])
            # Otherwise actual plane wave propagating along z
            else:
                # print('full incident field')
                # fi = np.exp(1J*k*z)
                # ku = np.array([np.cos(self.alpha),0,np.sin(self.alpha)])
                fi = np.exp(1J*k*(z*np.cos(self.alpha)+y*np.sin(self.alpha)))
                gi = 1J*k*np.cos(theta)*fi
            #remove incident field inside pth spheres
            for p in idp:
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[p],z-self.d_pz[p])
                idx_i = r_p<self.ka[p]
                # print(fi[idx_i])
                fi[idx_i] = 0
                gi[idx_i] = 0

        ### scattered fields --------------------------------------------------------------------
        # idp = range(N)
        if ftype in 'sta':
            fs = np.zeros(r.shape,dtype=complex)
            gs = np.zeros(r.shape,dtype=complex)
            #print(coor)
            #outgoing field

            for p in idp:
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[p],z-self.d_pz[p])
                idx_o = r_p>=self.ka[p]     #;print(idx_o.shape)
                # print(idx_o)

                # A=hs.harmonique(nmax-1,theta_p,'tab',phi_p)
                """
                for l in range(nmax):
                    Yl=A[:,0][l]

                    fs[idx_o] += self.bp[p*nmax**2+l]*spu.hn1(l, k*r_p[idx_o])*Yl[idx_o]
                    gs[idx_o] += self.bp[p*nmax**2+l]*spu.hn1p(l,k*r_p[idx_o])*Yl[idx_o]

                    #t1=time.time()
                    #A1=spu.hn1p(l,k*r_p[idx_o])
                    #tf=time.time()-t1
                    #print(tf)
                """
                for l in range(nmax):
                    hl = spu.hn1(l, k*r_p[idx_o])
                    for m in range(2*l+1):
                        #q*nmax**2+(nu)**2+mu
                        #Yl=A[l,m]

                        m1=m
                        if (m>l):
                            m1=l-m
                        Ylm=spe.sph_harm(m1,l,phi_p[idx_o],theta_p[idx_o])
                        #print(Yl.shape)
                        #print(idx_o.shape)
                        # Yl=A[:,m1][l]
                        #print(Yl.shape)
                        #print(idx_o.shape)

                        fs[idx_o] += self.bp[p*nmax**2+l**2+m]*hl*Ylm
                        # fs[idx_o] += self.bp[p*nmax**2+l**2+m]*spu.hn1(l, k*r_p[idx_o])*Ylm[idx_o]
                        # gs[idx_o] += self.bp[p*nmax**2+l**2+m]*spu.hn1p(l,k*r_p[idx_o])*Ylm[idx_o]

                    #gs SHOULD USE THE TRANSLATION TO GET RADIAL DERIVATIVE r_p
            #remove scattered field inside spheres
            for p in idp:# range(N):
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[p],z-self.d_pz[p])
                idx_i = r_p<self.ka[p]
                fs[idx_i] = 0
                gs[idx_i] = 0

            #inside field
            for p in idp:#range(N):
                r_p,theta_p,phi_p = spu.cart2sphere(x,y-self.d_py[p],z-self.d_pz[p])
                idx_i = r_p<self.ka[p]     #;print(idx_i.shape)

                # A=hs.harmonique(nmax-1,theta_p,'tab',phi_p)
                for l in range(nmax):
                    jl = spu.jn( l, n_p*k*r_p[idx_i])
                    for m in range(2*l+1):
                        m1=m
                        if (m>l):
                            m1=l-m
                        #Yl=hs.harmonique(l,theta_p,m1)
                        #print(time.time()-ti2)
                        #Yl=hs.harmonique(l,theta_p,0)
                        Ylm=spe.sph_harm(m1,l,phi_p[idx_i],theta_p[idx_i])
                        # Yl=A[:,m1][l]
                        # print('id=%d,l=%d,m=%d' %(p*nmax**2+l**2+m,l,m1))
                        fs[idx_i] += self.ap[p*nmax**2+l**2+m]*jl*Ylm
                        #fs[idx_i] += self.ap[p*nmax**2+l**2+m]*spu.jn( l, n_p*k*r_p[idx_i])*Ylm[idx_i]
                        #gs[idx_i] += self.ap[p*nmax**2+l**2+m]*n_p*spu.jnp( l,n_p*k*r_p[idx_i])*Ylm[idx_i]

                        #t1=time.time()
                        #A1=spu.jn( l,     n_p*k*r_p[idx_i])
                        #tf=time.time()-t1
                        #print(tf)


        ### total field -----------------------------------------------------------------------------------------
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
