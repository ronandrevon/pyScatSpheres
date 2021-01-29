import importlib as imp
from optics.scattering import qdot_sphere_array as qsa ;imp.reload(qsa)
from optics.scattering import spherical_utils as spu   ;imp.reload(spu)
from optics.scattering import spherical_utils as sphu  ;imp.reload(sphu)
from utils import *
imp.reload(dsp)
from scipy.integrate import trapz,quad
import scipy.special as spe
import scipy.fft as fft

plt.close('all')
path='../../docs/figures/'
name=path+'qdotSphereSingle_'

opts = 'bQF' #b(Born) #w(weak phase)  Q(qa)F(f)M(multi-shell)
k0 = 1
# kas = np.array([2,5,10])
kas = np.array([10,30,60])
eps = 0.1
# eps = 0.01
cm = 'Spectral'

# theoretical exact scattering
kp = np.sqrt(1+eps)
# theta = np.linspace(1e-3,np.pi,361)
# theta_deg = np.rad2deg(theta)
theta_deg = np.linspace(1e-3,30,361)
theta = np.deg2rad(theta_deg)
s1 = [qsa.QdotSphereArray(N=1,ka=ka,kp=kp,kd=0,nmax=80,solve=1,copt=0) for ka in kas]
f1 = np.array([s.get_ff(theta) for s in s1])


ct,st = np.cos(theta),np.sin(theta)
cs = dsp.getCs(cm,kas.size)
ka_max = kas.max()
q = 2*np.sin(theta/2)

# Th Born approximation Radial fourier transform of the potential
if 'b' in opts:
    if 'Q' in opts:
        qa = np.linspace(0,2*ka_max,300)
        qa_x,qa_y = np.meshgrid(qa,qa)
        qas = np.sqrt(qa_x**2+qa_y**2)
        f = (np.sin(qas)/(qas)-np.cos(qas))/qas**2

        plts = [[ka*st,ka*(1-ct),cs[i],'ka=%.1f' %ka] for i,ka in enumerate(kas)]
        dsp.stddisp(plts,im=[qa_x,qa_y,np.abs(f)],lw=2,labs=['$aq_x$','$aq_z$'],title='$|f(qa)|$',
            xylims=[0,2*ka_max,0,2*ka_max],cmap='viridis',
            name=name+'fqa.png',opt='ps')

    ### plot scattering amplitudes
    if 'F' in opts:
        f = np.array([ eps*ka**3*(np.sin(ka*q)/(ka*q)-np.cos(ka*q))/(ka*q)**2 for i,ka in enumerate(kas)])
        f  = np.abs(f); f /= f.max()
        f1 = np.abs(f1);f1/= f1.max()

        plts=[]
        plts += [[theta_deg, f[i],[cs[i],'o'  ],'$ka=%.1f$' %ka] for i,ka in enumerate(kas)]
        plts += [[theta_deg,f1[i],[cs[i],'--' ],''             ] for i,ka in enumerate(kas)]
        legElt = {'exact':'k--','born':'ko'}
        inset={'axpos':[0.3,0.3,0.4,0.6],'xylims':[65,180,-0.005,0.025],'ms':3,'ec':'k'}
        # inset={'axpos':[0.27,0.33,0.45,0.55],'xylims':[0,60,-0.01,0.15],'ms':3,'ec':'k'}
        dsp.stddisp(plts,lw=2,labs=[r'$\theta(^{\circ})$',r'$|f(\theta)|$'],inset=inset,iOpt='GtX',
            xylims=['x',0,180],legElt=legElt,title = r'$\epsilon=%.2f$, $n_{ref}=%.3f$' %(eps,kp),
            name=name+'fka1.svg',opt='ps')

    if 'M' in opts:
        E,lam = 200,0.025
        # r0 = np.hstack([[0.025,0.05] ,np.arange(1,6)*0.1+0.05])
        # V0 = np.hstack([[15,2.20],[0.8,0.3,0.2,0.15,0.05]])
        r0 = np.arange(50)*0.01+0.01
        V0 = 0.016/r0**1.5
        kas = r0*2*np.pi/lam
        eps = V0/E

        csf = dsp.getCs('jet',kas.size)
        img = 'data/potential_C.png'
        rs = np.hstack([0,r0])
        patches = [dsp.Rectangle((rs[i],0),rs[i+1]-rs[i],V,fill=1,color=csf[i],ec='b',alpha=0.5) for i,V in enumerate(V0)]
        fig=dsp.image_bg(img,patches=patches,xylims=[0,0.5,0,20],xyTicks=[0.05,2],
            name=name+'shells_pot1.png',opt='ps')
        # dsp.fig_pos(fig)

        eps = np.hstack([eps,0])
        fi = [ (eps[i]-eps[i+1])*ka**3*(np.sin(ka*q)/(ka*q)-np.cos(ka*q))/(ka*q)**2 for i,ka in enumerate(kas)]
        f = np.array(fi).sum(axis=0)

        # epss = np.flipud(np.cumsum(np.flipud(eps[1:])))
        # f = fi[-1]
        # for i in range(len(fi)-1):
        #     f += (eps[i]-epss[i])/eps[i]*fi[i]
        plts = [[theta_deg,np.abs(f),'k','total' ]]
        plts += [[theta_deg,fi[i],[csf[i],'--'],''] for i,ka in enumerate(kas) ]
        # plts += [[theta_deg,fi[i],[csf[i],'--'],'$ka=%d$' %ka ] for i,ka in enumerate(kas) ]
        dsp.stddisp(plts,lw=2,labs=[r'$\theta(^{\circ})$',r'$|f(\theta)|$'],#inset=inset,iOpt='GtX',
            xylims=['x',0,30],#legElt=legElt,title = r'$\epsilon=%.2f$, $n_{ref}=%.3f$' %(eps,kp),
            name=name+'shells_fka1.svg',opt='ps')

sigV0 = k0*eps/2
#weak phase
if 'w' in opts:
    # projected potential and tranmsission via fourier transform
    # if 'F' in opts:
        # N = int(npts/2)
        # s = qmax*np.linspace(-1,1,npts)
        # qx,qy = np.meshgrid(s,s)
        # qs = np.sqrt(qx**2+qy**2)
        # Vq = np.abs(ka/qs**2*(np.sin(ka*qs)/(ka*qs)-np.cos(ka*qs)))
        # # plts = [st,(1-ct),'k--']
        # dsp.stddisp(im=[qx,qy,Vq],labs=['$q_x$','$q_y$'],lw=2)
        # Vz = scipy.fft.ifft2(Vq)
        # # Vz = scipy.fft.fftshift(Vz)
        # # tranmsission function
        # u = np.arange(-N,N)/(2*qmax)
        # x,y = np.meshgrid(u,u)
        # T = np.exp(1J*sig*Vz)
        # # dsp.stddisp(im=[x,y,np.real(T)],labs=['$x$','$y$'],lw=2)
        # # dsp.stddisp(im=[x,y,np.imag(T)],labs=['$x$','$y$'],lw=2)
        # dsp.stddisp(im=[x,y,np.abs(  scipy.fft.fftshift(T))],labs=['$x$','$y$'],lw=2,title='mag')
        # dsp.stddisp(im=[x,y,np.angle(scipy.fft.fftshift(T))],labs=['$x$','$y$'],lw=2,title='phase')
        # #diffraction pattern
        # F = scipy.fft.fft2(T)
        # F[0] = 0
        # # F = scipy.fft.fftshift(F)
        # dsp.stddisp(im=[qx,qy,np.abs(F)],labs=['$q_x$','$q_y$'],lw=2)


    # radial fourier transform of the transmission function
    if 'F' in opts:
        for ia,ka in enumerate(kas):
            fvr = lambda rho,q :np.real((1-np.exp(2*1J*sigV0*np.sqrt(ka**2-rho**2)))*spe.jv(0,q*rho)*rho)
            fvi = lambda rho,q :np.imag((1-np.exp(2*1J*sigV0*np.sqrt(ka**2-rho**2)))*spe.jv(0,q*rho)*rho)
            fqr = np.array([quad(fvr,0,ka,(q0,))[0] for q0 in q])
            fqi = np.array([quad(fvi,0,ka,(q0,))[0] for q0 in q])
            fq = np.abs(fqr+1J*fqi); fq_max=fq.max()

            # plts = [[theta_deg, fq/fq_max,'k-','eikonal']]
            # plts+= [[theta_deg,fqr/fq_max,'b-.']]
            # plts+= [[theta_deg,fqi/fq_max,'r-.']]
            f0a = np.abs(eps*ka**3*(np.sin(ka*q)/(ka*q)-np.cos(ka*q))/(ka*q)**2); f0a /= f0a.max()
            f1a = np.abs(f1[ia]);f1a/= f1a.max()
            f2a = fq/fq_max
            plts=[]
            plts+= [[theta_deg,f1a,['k','-'],'exact']]
            plts+= [[theta_deg,f2a,['r','--'],'MS']]
            plts+= [[theta_deg,f0a,['b','-.'],'Born']]
            # legElt = {'exact':'k--','born':'ko'}
            dsp.stddisp(plts,lw=2,labs=[r'$\theta(^{\circ})$',r'$|f(\theta)|$'],#inset=inset,iOpt='GtX',
                xylims=['x',0,theta_deg.max()],title = r'$ka=%.1f$, $\epsilon=%.2f$, $n_{ref}=%.3f$' %(ka,eps,kp),#legElt=legElt,
                name=name+'feka%d.svg' %ia,opt='ps')


    ### num DFT to get Fourier transform of projected potential
    # if 'K' in opts:
    # ka = 1
    # x0 = 10
    # xi = np.linspace(-x0,x0,500)
    # x,y,z = np.meshgrid(xi,xi,xi)
    # r = np.sqrt(x**2+y**2+z**2)
    # f =  np.zeros(x.shape)
    # f[r<=1] = ka
    # F = scipy.fft.fftshift(scipy.fft.fftn(f))
    #
    # f0 = np.abs(f[:,:,100]);plt.figure();plt.imshow(f0,extent=[-x0,x0,-x0,x0]);plt.show()
    # F0 = np.abs(F[:,:,100]);plt.figure();plt.imshow(F0)
    # plt.show()
    #
    # q = scipy.fft.fftshift(scipy.fft.fftfreq(xi.size,1/(200)))
    # f0 = np.abs(F[:,100,100])
    # qs = np.linspace(1e-2,q.max(),1000)
    # f = np.abs(1/qs**2*(np.sin(ka*qs)/(ka*qs)-np.cos(ka*qs))) #; f[q==0]=0
    # plts = [[q,f0/f0.max(),'bo']]
    # plts+= [[qs,f/f.max(),'r-']]
    # dsp.stddisp(plts,lw=2)
