import importlib as imp
import scipy.special as spe
import numpy as np, pandas as pd
from . import displayStandards as dsp #;imp.reload(dsp)
from . import glob_colors as colors
from . import spherical_utils as spu  #;imp.reload(spu)
from . import hard_sphere_base_arb as hsba ;imp.reload(hsba)
import harmo_sphe as hs
import time
class QdotSphereArray(hsba.HardSphereArrayBaseArb):
    def solve(self,nmax=-1,copt=1,opt2=0,v=1):
        ''' Finds the unknown coefficients
        - nmax : max inlcuded order
        - copt : solve coupled problem
        - opt2 : get
        '''
        if nmax>0:self.nmax=nmax+1
        N,n_p,ka,kd_z,kd_y,kdpz,kdpy,nmax,kdp = self.kd_z*self.kd_y,self.kp,self.ka,self.kd_z,self.kd_y,self.kd_pz,self.kd_py,self.nmax,self.kd_p

        ###----------------------------------------------------------------------------------------
        #   Construction du second membre L
        ###----------------------------------------------------------------------------------------
        #Initialisation 
        L=np.zeros(2*(N*nmax*(nmax+1)**2)) #cf calcul de ma feuilleeeeee trop dur les maths ahhaha grrr
        for p in range(N) : #pas faire boucle sur sur boule mais pour l'instant on va faire ca sinon embrouiller moi tete
            m=0
            for l in range(nmax) :

                    #Calcul coeff Cl
                    cl = 1J**l*(2*l+1)*np.sqrt(4*np.pi/(2*l+1))#c depend que de l 

                    #Calcul des fct de bessel
                    jl0=spu.jn(l,ka[p]) #j'espere que c'est bon ecrit comme ca ? peut etre c'est juste ka???????
                    djl0=spu.jnp(l,ka[p])

                    #Calcul des exp
                    exp1=np.exp(1J*kdp[p]) #ca ca dépend que de N

                    #Assemblage matrice
                    L[m:m+2*l+1]=cl*exp1*jl0
                    L[m+N*nmax*(nmax+1)**2+1:m+2*l+nmax*(nmax+1)**2+1]=cl*exp1*djl0 #hihi j'espere que les indices sont bons sinon je pleure
                    m=m+2*l+1

        ###----------------------------------------------------------------------------------------
        #   Construction de P 
        ###----------------------------------------------------------------------------------------
        #Initialisation
        P = np.zeros((2*(N*nmax*(nmax+1)**2),2*(N*nmax*(nmax+1)**2)),dtype=complex)
        if v:print(colors.blue+"...assembling coupling..."+colors.black)
        for p in range(N):
            if v>1:print(colors.yellow+'p=%d' %p+colors.black)
            m=0
            for l in range(nmax):
                #calcul des fonctions bessel hankel 
                jl0,jl1   = np.array([spu.jn( l,np.array([ka,n_p*ka]))  for l in ll]).T
                djl0,djl1 = np.array([spu.jnp(l,np.array([ka,n_p*ka]))  for l in ll]).T
                hl0  = np.array([spu.hn1(l,ka)  for l in ll])
                dhl0 = np.array([spu.hn1p(l,ka) for l in ll])

                #remplissage de la matrice
                P[m:m+2*l+1,m:m+2*l+1]=np.diag([jl1[p]]*(2*l+1)) #en haut a gauche EST CE QUE C'EST INCIDE P OU L ???????? 
                P[m:m+2*l+1,m+N*nmax*(nmax+1)**2+1:m+2*l+nmax*(nmax+1)**2+1]=-np.diag([hl0[p]]*(2*l+1))#en haut a droite
                P[m+N*nmax*(nmax+1)**2+1:m+2*l+nmax*(nmax+1)**2+1,m:m+2*l+1]=np.diag([djl1[p]]*(2*l+1)) #en bas a gauche
                P[m+N*nmax*(nmax+1)**2+1:m+2*l+nmax*(nmax+1)**2+1,m+N*nmax*(nmax+1)**2+1:m+2*l+nmax*(nmax+1)**2+1]=-np.diag([dhl0[p]]*(2*l+1)) #en bas a droite

                m=m+2*l+1

        ###----------------------------------------------------------------------------------------
        #   Construction de T
        ###----------------------------------------------------------------------------------------
        T = np.zeros((2*(N*nmax*(nmax+1)**2),2*(N*nmax*(nmax+1)**2)),dtype=complex)
        T1=np.zeros(N*nmax*(nmax+1)**2)
        T2=np.zeros(N*nmax*(nmax+1)**2)

        #vecteur ligne 1
        Un=np.array([1]*(nmax*(nmax+1)**2))

        #vecteur colonne des j_nu et dj_nu
        jnu=np.zeros(nmax*(nmax+1)**2)
        djnu=np.zeros(nmax*(nmax+1)**2)
        m=0
        for l in range(nmax) :
            jnu[m:m+2*l+1]=np.array([jl0[p]]*(2*l+1))#je crois que ya confusion sur l'indice ?D:::  je sais paaaaas
            djnu[m:m+2*l+1]=np.array([djl0[p]]*(2*l+1))
            m=m+2*l+1

        #Construction de A
        fz_q = spu.hn1
        
        for p in range(N):
            A=np.zeros((nmax*(nmax+1)**2,nmax*(nmax+1)**2))
            for l in range(nmax) :
                for m in range (2*l+1):#ca fait beaucoup de boucles oupsssss
                    #somme sur q!=p des matrices carrées
                    As=np.zeros((nmax*(nmax+1)**2,nmax*(nmax+1)**2))
                    for lk in range(nmax) :
                        for mk in range(2*lk+1):
                            for q in range(p):
                                coor=np.array([abs(kdpz[p]-kdpz[q]),abs(kdpy[p]-kdpy[q])])
                                kdpq = np.linalg.norm(coor)
                                theta=np.arccos(coor[0]/kdpq)#coor polaire ?
                                As1[lk,mk]=spu.a_lnmp(l,m,lk,mk,fz_q,kdpq,theta,0)#il manque l,m??
                            for q in range(p+1,N):
                                coor=np.array([abs(kdpz[p]-kdpz[q]),abs(kdpy[p]-kdpy[q])])
                                kdpq = np.linalg.norm(coor)
                                theta=np.arccos(coor[0]/kdpq)#coor polaire ?
                                As2[lk,mk] = spu.a_lnmp(l,m,lk,mk,fz_q,kdpq,theta,0)#il manque l,m??
                    A+=As1+As2

            T1[p]=np.dot(Un,A)
            T1[p]=np.dot(T1[p],jnu)
            T2[p]=np.dot(T2[p],djnu)

        #assemblage dans la matrice 
        idp1=slice(0,N*nmax*(nmax+1)**2+1)
        idp2=slice(N*nmax*(nmax+1)**2+1,2*(N*nmax*(nmax+1)**2)+1)
        T[idp1,idp2]=np.diag(T1)
        T[idp2,idp2]=np.diag(T2)

        ###----------------------------------------------------------------------------------------
        #   Résolution
        ###----------------------------------------------------------------------------------------
        if v:print(colors.blue+"...solving..."+colors.black)
        cp = np.linalg.solve(P-copt*T,L)
        print(P.shape)
        self.ap,self.bp = cp[:N*nmax],cp[N*nmax:]

        return self.ap,self.bp
