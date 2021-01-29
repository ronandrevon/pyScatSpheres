import importlib as imp
from optics.scattering import spherePEC_arrayApprox as peca ;imp.reload(peca)
from optics.scattering import spherePEC as sPEC             ;imp.reload(sPEC)
from optics.scattering import spherical_utils as spu        ;imp.reload(spu)
from utils import *
plt.close('all')
path='../../docs/figures/'
impath='../../docs/resources/hamid1990/'
name=path+'PECspheresApprox_'
opts= 'C'

k = 1
Nmax = 5
lam = 2*np.pi/k

def PEC2(ka,kd,Nmax):
    f = np.zeros((2),dtype=complex)
    for n in range(1,Nmax+1):
        bn = spu.jn(n,ka)/spu.hn1(n,ka)
        cn = (ka*spu.jnp(n,ka)+spu.jn(n,ka))/(ka*spu.hn1p(n,ka)+spu.hn1(n,ka))
        f[0] +=        -(2*n+1)*(bn+cn)/2
        f[1] += -(-1)**n*(2*n+1)*(bn-cn)/2
    f *= 1J/kd
    ed = np.exp(1J*kd)      #;print('f2:',f);print('f2*ed:',f*ed)
    T = ed*np.array([[0,f[1]],[f[1],0]],dtype=complex)
    L = np.array([f[0],ed*f[1]],dtype=complex)
    Cp = np.linalg.solve(np.identity(2) - T,L)
    return Cp

if '0' in opts:
    a_p,d = 2, 5
    rcs0 = sPEC.get_monoRCS(ka=k*a_p,Nmax=Nmax)
    s2    = peca.PECarray(a_p, np.array([0,d]), lam,Nmax)
    rcs1  = s2.compute_monoRCS()
    print('|Cp|:',abs(s2.Cp));print('phi(Cp):',np.rad2deg(np.angle(s2.Cp)))
    ka,kd = a_p,d
    Cp = PEC2(ka,kd,Nmax)
    # print(s2.Cp);print(Cp)
    # s2.Cp = Cp
    # rcs2  = s2.compute_monoRCS()
    # print('rcs0:',rcs0)
    # print('rcs1:',rcs1)
    # print('rcs2:',rcs2)

if 'C' in opts:
    ka,kds = 1.5,np.linspace(1,10,101)
    Cp = np.array([PEC2(ka,kd,Nmax) for kd in kds])
    plts=[]
    plts += [[Cp[:,0].real,Cp[:,0].imag,'b-o','$C_1$'],[Cp[0,0].real,Cp[0,0].imag,'bs','']]
    plts += [[Cp[:,1].real,Cp[:,1].imag,'r-o','$C_2$'],[Cp[0,1].real,Cp[0,1].imag,'rs','']]
    dsp.stddisp(plts,labs=['$Re$','$Im$'],lw=2)



#single mono_rcs
if '1' in opts:
    ka  = np.linspace(0.1,2,30)
    s2  = [peca.PECarray(ap, np.array([0,2*ap]), lam,Nmax) for ap in ka]
    rcs  = np.array([s.compute_monoRCS(idp=[0],Copt=0) for s in s2])
    ka0  = np.linspace(ka[0],ka[-1],100)
    rcs0 = sPEC.get_monoRCS(ka0,Nmax)
    plts = [[ka,rcs,'bo'],[ka0,rcs0,'r']]
    dsp.stddisp(plts,lw=2,labs=[r'$ka$',r'$\sigma/\pi a_0^2$'],
        opt='p')#,xylims=[1,11,0,rcs.max()])
#fig2b
if '2' in opts:
    N = 5
    a_p = 0.5
    d_ps = np.linspace(1,11,500)
    rcs0 = sPEC.get_monoRCS(a_p,Nmax)

    print('...solving...')
    s2 = [peca.PECarray(a_p, dp*np.arange(N), lam,Nmax) for dp in d_ps]
    print('...rcs...')
    rcs = np.array([s.compute_monoRCS() for s in s2])

    plts=[[d_ps,rcs,'b'],[d_ps[[0,-1]],rcs0*N**2*np.array([1,1]),'k--']]
    # dsp.stddisp(plts,lw=2,labs=[r'$kd$',r'$\sigma/\pi a_0^2$'],
    #     opt='p',xylims=[1,d_ps[-1],0,19],xyTicks=[2,3])

    im = impath+'fig2b0.png'
    fig = dsp.image_bg(im,plots=plts,lw=2,labs=[r'$kd$',r'$\sigma/\pi a_0^2$'],
        xylims=[1,11,0,19],xyTicks=[2,3],
        opt='ps',name=name+'fig2b.png')
#fig3d
if '3' in opts:
    N,ka,kd = 5,0.5,4
    s5 = peca.PECarray(ka, kd*np.arange(N), lam,Nmax)
    fig=s5.biRCS_show(npts=100,is_3d=False,phi=[0,90],bargs={'Copt':0},im=impath+'fig3d0.png',
        xylims=[0,180,0,12],xyTicks=[30,2],
        opt='ps',name=name+'fig3d.png')
