from utils import *                                ;imp.reload(dsp)
from pyScatSpheres import qdot_sphere_array as qsa ;imp.reload(qsa)
from scipy.integrate import trapz,quad
from scipy.interpolate import interp1d
import scipy.special as spe
plt.close('all')
name='figures/TmatrixSingleS.eps'

E = 200
k0 = 2*np.pi/cst.keV2lam(E)
sigAng = lambda sig,kas:np.log10(sig*np.pi*(kas/k0)**2)

##### precomputed np curves
df_name = '../tests/data/qdotSphereArray1.pkl'
df = pd.read_pickle(df_name)
kas,kps = df.ka.unique(),df.kp.unique()[:-1]
# ras = ka/k0
plts,cs = [],dsp.getCs('jet',kps.size)
for ik,kp in enumerate(kps):
    dfn = df[df.kp==kp]
    sig_ka = np.array([4*np.sum(np.abs(r.bp)**2)/r.ka**2 for i,r in dfn.iterrows()]) #;print(sig_th)
    plts += [[kas,sigAng(sig_ka,kas),cs[ik],'$n_p=%.4f$' %kp]]


#### Carbone shells points @200keV
r,V = np.loadtxt('../tests/data/C_V.txt').T
kpC  = interp1d(r*k0,np.sqrt(1+V/1e3/E))
kas0 = np.hstack([[5.5],np.arange(10,101,10)])
kps0 = kpC(kas0)
qdotC = [qsa.QdotSphereArray(N=1,ka=ka1,kp=kp1,kd=0,copt=0,solve=1,nmax=int(ka1)+3) for ka1,kp1 in zip(kas0,kps0)]
sigC = np.array([qdot.get_s(npts=3600,norm=1) for qdot in qdotC])

plts+= [[kas0,sigAng(sigC,kas0),'b-o','$C@200keV$']]

dsp.stddisp(plts,labs=[r"$ka$",r"$\log_{10}(\sigma/\AA^2)$"],
    lw=2,fonts='2',xylims=[0,100,-8,1],#logOpt='y',
    name=name,opt='ps')
