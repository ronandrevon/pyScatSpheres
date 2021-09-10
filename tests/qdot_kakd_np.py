from utils import *                                 ;imp.reload(dsp)
import scipy,copy,os,scipy.special as spe
from scipy.signal import find_peaks
from scipy.integrate import trapz,quad
from scipy.interpolate import interp1d
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from pyScatSpheres import utils as ut ;imp.reload(ut)
plt.close('all')
name='figures/qdotSphereArray2approx_kakp_'
df_name = 'data/df_qdotSpheresN2.pkl'
opt = 'p'

opts = 'e' # s(solve) e(set_error)
optsA='n' #d(ka,kd) n(ka,kp) e(err vs param)

if 's' in opts:
    # kas0   = np.array([0.5,1,2,3,4,5,7,10,15,20,30])
    # kdkas0 = np.array([2.1,3,5,10,25,20])
    # kps0   = np.array([1.2,1.1,1.01,1.001,1.0001])
    # kas,kdkas,kps = np.meshgrid(kas0,kdkas0,kps0)
    # kas,kdkas,kps =kas.flatten(),kdkas.flatten(),kps.flatten()
    # Ns=np.array([2]*kas.size,dtype=int)

    kas0   = np.array([5,7,10,15])
    kps0   = np.array([1.1,1.01,1.001,1.0001])
    Ns0    = np.array([2,3,4,6,8,10,15,20,30]) #np.stack([np.arange(2,11,1),np.arange(15,26,5)])
    kas,Ns,kps = np.meshgrid(kas0,Ns0,kps0)
    kas,Ns,kps = kas.flatten(),Ns.flatten(),kps.flatten()
    kdkas = np.array([3]*kas.size)

    df = ut.solve_set(df_name,kas,kps,kdkas,Ns)

df = pd.read_pickle(df_name)
if 'e' in opts:df = ut.set_errors(df_name)

# thetas = np.deg2rad(np.linspace(0,180,2001))
# ff = qdot.get_ff(thetas)
# # idp = np.hstack([ p*qdot.nmax+np.arange(4,qdot.nmax) for p in range(qdot.N)])
# idp = abs(qdot.bp)<1e-4
# qdot.bp[idp]=0#; qdot.show_ff()
# ff0 = qdot.get_ff(thetas)
# dsp.stddisp([[thetas,np.log10(abs(ff0-ff)),'b-']])
# qdot.test_convergence(nmaxs=[5,7,10])

df['kdka'] = df.kd/df.ka
df = df.loc[df.kdka<100]
kas0=np.array(df.ka.unique(),dtype=float)
kps0=np.array(df.kp.unique(),dtype=float)
kds0=np.array(df.kdka.unique(),dtype=float)
argsM = {'pOpt':'im','cmap':'RdYlGn_r'}


ikp = 3
kp = kps0[ikp]
# kderrs,errs = np.zeros(kas0.shape),np.zeros(kas0.shape)
dfkp = df.loc[df.kp==kp]
nkas = kas0.size
errs_a = np.zeros((nkas,kds0.size))
errs_0 = np.zeros((nkas,kds0.size))
for i,ka in enumerate(kas0):
    dfka = dfkp.loc[dfkp.ka==ka]
    errs_a[i,:] = np.log10(np.array(dfka.err_a.values,dtype=float))
    errs_0[i,:] = np.log10(np.array(dfka.err_0.values,dtype=float))


xylims=['y',min(errs_0.min(),errs_a.min()),max(errs_0.max(),errs_a.max())]
# cs = dsp.getCs('Spectral',nkas)
# args = {'labs':['$kd/ka$',r'$log_{10}(err)$'],'xylims':xylims,'lw':2}
# plts0 = [[kds0,errs_0[i,:],[cs[i],'-o'],'$ka=%.1f$'%ka] for (i,ka) in enumerate(kas0) ]
# pltsa = [[kds0,errs_a[i,:],[cs[i],'-o'],'$ka=%.1f$'%ka] for (i,ka) in enumerate(kas0) ]
# dsp.stddisp(plts0,title='error uncoupled, $k_p=%.4f$' %kp       ,**args)
# dsp.stddisp(pltsa,title='error forward coupling, $k_p=%.4f$' %kp,**args)

# kds,kas =np.meshgrid(kas0,kds0)
#### maps ka,kdka
lkds0 = np.arange(kps0.size)#;fkp = interp1d(kps0,lkps0)
lkas0 = np.arange(kas0.size)#;fka = interp1d(kas0,lkas0)
argsM.update({'labs':['$kd/ka$','$ka$'],'caxis':[-8,0],
    'xyTicks':[lkds0,lkas0],'xyTickLabs':[np.array(kds0,dtype=str),np.array(kas0,dtype=str)],
    'opt':'ps'})
dsp.stddisp(im=[lkds0,lkas0,errs_a],title='error forward   kp=%.4f' %kp,
    name=name+'kakd_forward%d.png' %ikp,**argsM)#caxis=[errs_a.min(),errs_a.max()]
dsp.stddisp(im=[lkds0,lkas0,errs_0],title='error uncoupled kp=%.4f' %kp,
    name=name+'kakd_uncoupled%d.png' %ikp,**argsM)#caxis=[errs_0.min(),errs_0.max()]
