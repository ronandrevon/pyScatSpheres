from utils import *                                 ;imp.reload(dsp)
import scipy,copy,os,scipy.special as spe
from scipy.signal import find_peaks
from scipy.integrate import trapz,quad
from scipy.interpolate import interp1d
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import utils as ut ;imp.reload(ut)
plt.close('all')
name='figures/qdotSphereArray2approx_kakp_'
df_name = 'data/df_qdotSpheres_2.pkl'
opt = 'sp'
opts = '' # s(resolve) e(set_error)

#### load dataframe with precomputed solutions
df = pd.read_pickle(df_name)
if 'e' in opts:df = ut.set_errors(df_name)


#### gather the errors from df
df['kdka'] = df.kd/df.ka
df = df.loc[df.kdka<100]
kas0=np.array(df.ka.unique(),dtype=int)
kps0=np.array(df.kp.unique(),dtype=float)
kds0=np.array(df.kdka.unique(),dtype=float)

kas0=kas0[kas0>1]
kd0=3
dfkd = df.loc[(df.kdka==kd0) & (df.N==2)]
errs_pa = np.zeros((kas0.size,kps0.size))
errs_p0 = np.zeros((kas0.size,kps0.size))
for i,ka in enumerate(kas0):
    dfkap = dfkd.loc[dfkd.ka==ka]
    errs_pa[i,:] = np.log10(np.array(dfkap.err_a.values,dtype=float))
    errs_p0[i,:] = np.log10(np.array(dfkap.err_0.values,dtype=float))

lkps0 = np.arange(kps0.size);fkp = interp1d(kps0,lkps0)
lkas0 = np.arange(kas0.size);fka = interp1d(kas0,lkas0)

#### Carbone pts
E = 200              #keV
lam = cst.keV2lam(E) #A
r,V = np.loadtxt('data/C_V.txt').T
fV  = interp1d(r,V/1e3)
kasC = np.array(kas0[kas0>=6],dtype=float)
r0 = kasC/(2*np.pi/lam)
kpsC = np.sqrt(1+fV(r0)/E)
plts = [fkp(kpsC),fka(kasC),'bo-','C@200keV']

argsM={'pOpt':'im','cmap':'RdYlGn_r','labs':['$n_p$','$ka$'],
    'plots':plts,'ms':10,'caxis':[-5,0],'fonts':'2',
    'xyTicks':[lkps0,lkas0],'xyTickLabs':[np.array(kps0,dtype=str),np.array(kas0,dtype=str)],
    'opt':opt}
dsp.stddisp(im=[lkps0,lkas0,errs_pa],#title='error forward   kd/ka=%.4f' %kd0,
    name=name+'forward.png',**argsM)
dsp.stddisp(im=[lkps0,lkas0,errs_p0],#title='error uncoupled kd/ka=%.4f' %kd0,
    name=name+'uncoupled.png',**argsM)


# idx = abs(ff)/ffmax>1e-3
# err = abs((ff[idx]-ffu[idx])/ff[idx]).mean()
# print(err)

# if err<0.05:qdot2.show_f()

# plts = [[kas0,kderrs,'b-o','$err=%1.E$' %err]]
# txts = [[ka,kdka,'%.2E' %err,'b'] for ka,kdka,err in zip(kas0,kderrs,errs)]
# dsp.stddisp(plts,texts=txts,labs=['$ka$','$kdka$'],fonts={'text':15})

# s2s,s0s = [],[]
# for i,N in enumerate(Ns):
#     s2 = qsa.QdotSphereArray(N=N,ka=ka,kp=kp,kd=kd*ka,nmax=nmax);s2.solve(v=1)
#     s0 = copy.copy(s2); s0.solve(copt=0)
#     fig,ax = s0.show_ff(npts=2000,fopts='m',opt='',leg=0)
#     s2.show_ff(npts=2000,ax=ax,fopts='m',legElt={'uncoupled':'k--'},name=name+'N%d_' %i,opt='p')
#     s2s+=[s2]
#     s0s+=[s0]
