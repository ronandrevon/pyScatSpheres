import importlib as imp
from utils import *
from pyScatSpheres import hard_sphere_arrayApprox as hsaA   ;imp.reload(hsaA)
from pyScatSpheres import hard_sphere_array as hsa          ;imp.reload(hsa)
plt.close('all')
path  ='../docs/figures/'
name=path+'hardSphereArray2approx_'


def sphere_array2Approx(N=2,nmax=7,ka=2.0,kd=8,fopts='ue', **kwargs):
    '''fopts : 'u' (uncoupled), e'(exact)'''
    theta  = np.linspace(0,np.pi,361)

    s2a = hsaA.HardSphereArray(N,ka,kd,nmax)
    cpa = s2a.solve(copt=1)#;print(cpa)
    fa  = s2a.get_ff(theta)
    theta_d = np.rad2deg(theta)
    plts= [[theta_d,np.abs(fa) ,'r-' ,'$approx_{coupled}$']]

    #exact
    if 'e' in fopts:
        s2  = hsa.HardSphereArray(N,ka,kd,nmax)
        cp  = s2.solve(copt=1)#;print(cp)
        fe  = s2.get_ff(theta)
        plts+= [[theta_d,np.abs(fe),'b-','$exact$']]
    #uncoupled
    if 'u' in fopts:
        cp  = s2a.solve(copt=0)
        fa0 = s2a.get_ff(theta)
        plts+= [[theta_d,np.abs(fa0) ,'b--' ,'$uncoupled$']]




    dsp.stddisp(plts,labs=[r"$\theta(^{\circ})$",r"$|f(\theta)|$"],lw=2,
        title='Scattering amplitudes for $ka=%.1f$, $kd=%.1f$' %(ka,kd),xylims=['x',0,180],
        **kwargs)

    return s2a


# s2  = hsa.HardSphereArray(N=1,ka=10,kd=0,nmax=15,solve=True);s2.show_ff()
s2a=sphere_array2Approx(N=2,ka=10.0,kd=25,nmax=15,fopts='ue',name=name+'fka0.svg',opt='p')
# s2a=sphere_array2Approx(N=2,ka=5.0,kd=40,nmax=40,fopts='ue',name=name+'fka1.svg',opt='p')
# ss = sweep_ka(N=2,nmax=20 ,kas=[2.0],nkds=10,kdr=(10,11),Nmax=5,opts='cA',name=name,opt='p')
