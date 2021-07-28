import importlib as imp
import scipy.special as spe
import numpy as np, pandas as pd
from . import displayStandards as dsp #;imp.reload(dsp)
from . import glob_colors as colors
from . import spherical_utils as spu  #;imp.reload(spu)
from . import hard_sphere_base_arb as hsba ;imp.reload(hsba)
from . import harmo_sphe as hs
import time
class QdotSphereArray(hsba.HardSphereArrayBaseArb):
    def solve(self,nmax=-1,copt=1,opt2=0,v=1):
        ''' Finds the unknown coefficients
        - nmax : max inlcuded order
        - copt : solve coupled problem
        - opt2 : get
        '''
        if nmax>0:self.nmax=nmax+1
        N,n_p,ka,kd_z,kd_y,kdpz,kdpy,nmax = self.N,self.kp,self.ka,self.kd_z,self.kd_y,self.kd_pz,self.kd_py,self.nmax

        print(colors.red,"Non linéaire",colors.black)

        ###----------------------------------------------------------------------------------------
        #   Construction du second membre L
        ###----------------------------------------------------------------------------------------
        #Initialisation 
        L=np.zeros(2*(N*(nmax)**2),dtype=complex) #cf calcul de ma feuilleeeeee trop dur les maths ahhaha grrr
        for p in range(N) : #pas faire boucle sur sur boule mais pour l'instant on va faire ca sinon embrouiller moi tete
            m=0
            for l in range(nmax) :

                    #Calcul coeff Cl
                    cl = 1J**l*(2*l+1)*np.sqrt(4*np.pi/(2*l+1))#c depend que de l 

                    #Calcul des fct de bessel
                    jl0=spu.jn(l,ka) #j'espere que c'est bon ecrit comme ca ? peut etre c'est juste ka???????
                    djl0=spu.jnp(l,ka)

                    #Calcul des exp
                    exp1=np.exp(1J*kdpz[p]) #ca ca dépend que de N

                    #Assemblage matrice
                    #L[m+p*nmax**2:m+2*l+1+p*nmax**2]=cl*exp1*jl0
                    L[m+p*nmax**2]=cl*exp1*jl0
                    #L[m+N*(nmax)**2+p*nmax**2:m+2*l+1+N*(nmax)**2+p*nmax**2]=cl*exp1*djl0 #hihi j'espere que les indices sont bons sinon je pleure
                    L[m+N*(nmax)**2+p*nmax**2]=cl*exp1*djl0
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
            jl0,jl1   = np.array([spu.jn( l,np.array([ka,n_p*ka]))  for l in ll]).T
            djl0,djl1 = np.array([spu.jnp(l,np.array([ka,n_p*ka]))  for l in ll]).T
            hl0  = np.array([spu.hn1(l,ka)  for l in ll])
            dhl0 = np.array([spu.hn1p(l,ka) for l in ll])
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

        ### Vecteur colonne des j_nu et dj_nu ---------------------------------------------
        """
        jnu=np.zeros(((nmax)**2,1))
        djnu=np.zeros(((nmax)**2,1))
        m=0
        for l in range(nmax) :
            idp1=slice(m,m+2*l+1)
            jnu[m:m+2*l+1]=jl0[l]*np.ones((2*l+1,1))
            djnu[idp1]=djl0[l]*np.ones((2*l+1,1))
            m=m+2*l+1
        """
 

        ### Construction de T1 et T2 ----------------------------------------------------------
        T1=np.zeros(((N*(nmax)**2),(N*(nmax)**2)),dtype=complex)
        T2=np.zeros(((N*(nmax)**2),(N*(nmax)**2)),dtype=complex)
        fz_q = spu.hn1

        #p_indic=0
        #l_indic=0
        #m_indic=0
        for p in range(N):
            for l in range (nmax):
                    for m in range (2*l+1):
                        #c'est pour le coeff a la ligne lm de l'atome p
                        vect=np.zeros((1,N*nmax**2),dtype=complex)
                        #coeff de la ligne qu'on regarde
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
                                    #print(q*nmax**2+(nu)**2+mu)
                                    #coor=np.array([abs(kdpz[p]-kdpz[q]),abs(kdpy[p]-kdpy[q])])
                                    x=kdpz[p]-kdpz[q]
                                    y=kdpy[p]-kdpy[q]
                                    #theta=np.arccos(coor[0]/kdpq)
                                    kdpq = np.sqrt(x**2+y**2)
                                    #if (kdpz[q]-kdpz[p]>=0):
                                    theta=np.arccos((kdpz[q]-kdpz[p])/kdpq)
                                        #theta=np.arctan2(y,x)
                                    #else:
                                        #theta=np.arctan2(y,x)
                                    #if (kdpz[p]-kdpz[q]<0):
                                        #theta+=np.pi/2
                                    if (q!=p):
                                        
                                        print("boule ref :")
                                        print(p)
                                        print("boule actuelle :")
                                        print(q)
                                        print("theta :")
                                        print(theta)
                                        
                                        m1=m
                                        mu1=mu
                                        if (m>l):
                                            m1=l-m
                                        if (mu>nu):
                                            mu1=nu-mu

                                        a=spu.a_lmnp(l,m1,nu,mu1,fz_q,kdpq,theta,0)
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

        """
        T1[2:4,6:8]=np.zeros((2,2))
        T1[2:4,10:12]=np.zeros((2,2))
        T1[6:8,2:4]=np.zeros((2,2))
        T1[6:8,10:12]=np.zeros((2,2))
        T1[10:12,2:4]=np.zeros((2,2))
        T1[10:12,6:8]=np.zeros((2,2))

        T2[2:4,6:8]=np.zeros((2,2))
        T2[2:4,10:12]=np.zeros((2,2))
        T2[6:8,2:4]=np.zeros((2,2))
        T2[6:8,10:12]=np.zeros((2,2))
        T2[10:12,2:4]=np.zeros((2,2))
        T2[10:12,6:8]=np.zeros((2,2))
        """
        #print(jnu)
        #print(T1.shape)
        T[0:N*nmax**2,N*nmax**2:2*N*nmax**2]=T1
        T[N*nmax**2:2*N*nmax**2,N*nmax**2:2*N*nmax**2]=T2
        #la je mets les endroits ou c'est 0
        #for i in range (12):
            #for j in range(12):
                #T[p*nmax**2:(p+1)*nmax**2,p*nmax**2:(p+1)*nmax**2]=np.zeros((nmax**2,nmax**2))
        
        #print(T[2:4,2:4])

        ###----------------------------------------------------------------------------------------
        #   Résolution
        ###----------------------------------------------------------------------------------------
        #print(colors.red,"pas linéaire\n",colors.black)
        
        dsp.stddisp(im=[np.log10(abs(P))], pOpt='im',cmap='Spectral')
        #dsp.plt.show()
        idx=np.array([0,1,4,5,8,9])
        idx1 = np.hstack([idx,idx+12])

        V=(P-T)[np.ix_(idx1,idx1)]
        print(V)


        if v:print(colors.blue+"...solving..."+colors.black)
        cp = np.linalg.solve(P-T,L)
        self.ap,self.bp = cp[:N*nmax**2],cp[N*nmax**2:]
        self.V=V

        return self.ap,self.bp 

### Notes à moi meme : 
#- pb pour mu=nu 
#- pb theta