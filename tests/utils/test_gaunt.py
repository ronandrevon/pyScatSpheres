from utils import *
import scipy,copy,time
from pyScatSpheres import spherical_utils as spu   ;imp.reload(spu)
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)


def speed_copt():
    args = {'kp':1.25,'kd':20,'ka':2,'nmax':10,'N':20,'v':1}

    t0 = time.time()
    qdot2_1 = qsa.QdotSphereArray(copt=1,**args)
    print('vectorial : ',time.time()-t0)

    t0 = time.time()
    qdot2_2 = qsa.QdotSphereArray(copt=2,**args)
    print('serial : ',time.time()-t0)

    print(abs(qdot2_1.bp - qdot2_2.bp).sum())
    print(abs(qdot2_1.ap - qdot2_2.ap).sum())


def speed_gaunt(lmax=10,N=3):
    ds = 1+3*np.arange(N)
    ts = np.hstack([[0]*3,[np.pi]*(N-3)])

    t0 = time.time()
    aln0 = spu.get_aln(lmax-1,lmax-1,spu.hn1,ds,ts)
    print('vectorial : ',time.time()-t0)

    t0 = time.time()
    aln1 = np.zeros(aln0.shape,dtype=complex)
    for i,(d,t) in enumerate(zip(ds,ts)):
        aln = spu.a_ln(lmax-1,lmax-1,spu.hn1,d,t)
        aln1[:,i*lmax+np.arange(lmax)] = aln
    print('serial : ' ,time.time()-t0)

    print(abs(aln0-aln1).sum())
    # print(aln0.shape,aln1.shape)
# speed_gaunt(lmax=10,N=10)
speed_copt()
