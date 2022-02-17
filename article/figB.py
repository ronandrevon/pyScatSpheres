'''Computes the distance of a point to a parabola,
if the parabola is almost flat, the excitation error is very close to
the approximation used in Cowley1957'''
from utils import*;imp.reload(dsp)
plt.close('all')


def multiple_scattering(h=5,h0=0,h1=2,**kwargs):
    n,N=7,3
    K,dz = 10,2
    hs,ls = np.meshgrid(np.arange(n),np.arange(N)*dz)
    hs,ls = hs.flatten(),ls.flatten()
    hr,lr = [h0,h1,h-h1-h0],np.arange(N)*dz
    hn = np.cumsum(hr)

    #### ewald paraboloids
    scat=([ls,hs,10,'k','o'],)
    scat+=([lr,hr,50,'r','o'],)
    # t = np.linspace(0,np.deg2rad(30),100)
    # ct,st = K*(np.cos(t)-1),K*np.sin(t)
    u = np.linspace(0,n-0.5,100)
    f = lambda u:-1/(2*K)*u**2
    w = f(u)
    plts = [[w+i,u,'b-',''] for i in np.arange(N)*dz]
    fig,ax = dsp.stddisp(plts,scat=scat,lw=2,opt='')

    #### texts for  h_r,zeta_n
    wm,fsT=0.1,40
    txts_z = [[lr[i] +0*wm,hn[i] +wm,r'$\zeta_%d$' %(i+1),'b'] for i in range(N-1)]
    txts_z+= [[lr[-1]+0*wm,hn[-1]+wm,r'$\zeta_h=\zeta$','b']]
    dsp.stddisp(ax=ax,texts=txts_z,targs={'fontsize':fsT,'ha':'right','va':'bottom'},opt='')
    txts_h = [[lr[i] +wm,hr[i],r'$h_%d$' %(i+1),'r'] for i in range(N-1)]
    txts_h+= [[lr[-1]+wm,hr[-1],r'$h_N$','r'] ]
    dsp.stddisp(ax=ax,texts=txts_h,targs={'fontsize':fsT,'ha':'left','va':'top','linespacing':2},opt='')
    #### arrows,zetas
    eps=0.925
    hns = np.hstack([0,hn])
    arrows = [[(i-1)*dz,hns[i],eps*2,eps*hr[i],'r'] for i in range(N)]
    plts = [[[lr[i],lr[i]+f(hn[i])],[hn[i]]*2,'b--',''] for i in np.arange(N)]
    plts+= [[ [i,i],[0,n],'k--',''] for i in np.arange(N)*dz]
    txts = [[i*dz,n+0.1,r'$z_%d$' %(i+1),'k'] for i in np.arange(N)]
    ax.plot(lr,hn,'bo',ms=15,mew=2,mfc='none')
    dsp.stddisp(plts,ax=ax,texts=txts,arrows=arrows,lw=1.5,
        pOpt='XpG',axPos=[0,0,1,1],xylims=[-2.5,2*N-1,-1,n+1],
        xyTicks=2,fonts={'text':fsT},#xyTickLabs=[[],[]],
        **kwargs)


def ff_kin(N=30,dp=80 ):
    '''dp in units of lambda'''
    lam =1;k0=1/lam

    theta_d = np.linspace(1e-3,10,10000)
    theta = np.deg2rad(theta_d)
    q = 2*k0*np.sin(theta/2)/lam
    beta = dp*q*np.sin(theta/2)
    beta0 = 1
    # beta = dp*(1-np.cos(theta))
    f = np.sin(np.pi*N*beta)/np.sin(np.pi*beta)
    f2 = np.sin(np.pi*N*beta)/(np.pi*beta)
    f3 = np.sin(np.pi*N*(beta-beta0))/(np.pi*(beta-beta0))
    plts = [[q,abs(f),'b-',r'$sin(\pi\beta H)/sin(\pi\beta)$']]
    plts+= [[q,abs(f2),'r--',r'$sin(\pi\beta H)/\pi\beta$']]
    plts+= [[q,abs(f3),'g--','']]
    dsp.stddisp(plts,labs=[r'$q$','f'],xylims=[0,q.max(),0,abs(f).max()],lw=2)


def show_zeta(**kwargs):
    from  scipy.optimize import fsolve
    k0 = 5
    u0,w0=1,1.5
    f = lambda ua,u0,w0:(u0-ua)+1/k0*ua*(w0-1/(2*k0)*ua**2)
    ua = fsolve(f,u0,args=(u0,w0))[0]
    # print(f(ua,u0,w0))
    t = np.array([1,1/k0*ua])
    t/=np.linalg.norm(t)
    fe = lambda u:1/(2*k0)*u**2
    fc = lambda u:k0-np.sqrt(k0**2-u**2)

    h = np.array([u0,w0])
    A = np.array([ua,fe(ua)])
    ka = 0.1
    a0,a1 = A - ka*t,A + ka*t

    u = np.linspace(0.8,2.5,1000)
    plts = [
        [u,fe(u),'b','Ewald parabola'],
        [u,fc(u),'b--','Ewald circle'],
        [A[0],A[1],'bo',''],
        [u0,w0,'ro',''],
        [u0,fe(u0),'ro',''],
        [[u0,u0],[fe(u0),w0],'r--',''],#'$\zeta$'],
        # [[u0,A[0]],[A[1],A[1]],'k--',''],
        [[u0,A[0]],[w0,A[1]],'r-' ,''],#'$\hat d$'],
        [[a0[0],a1[0]],[a0[1],a1[1]],'c-',''],#'$\hat t$']
        ]
    wm,wmr = 0.075,0.15# 0.025
    txts = [
        [u0+wm,w0+wm     ,'$(u_n,w_n)$','r'],
        [A[0]+wm,A[1]+wm ,'$P$','b'],
        # [a1[0]+wm,a1[1]+wm,'$\hat t$','c'],
        [(u0+A[0])/2+wmr,(w0+A[1])/2+wm,'$\zeta^{(BW)}$','r'],
        [u0-2*wm,(w0+fe(u0))/2+wm,'$\zeta_n^{(MS)}$','r'],
        ]
    dsp.stddisp(plts,texts=txts,lw=2,ms=10,labs=['$u$','$w$'],pOpt='pe',
        fonts = {'text':40},xylims=[0,3,0,5],axPos=[0.1,0.07 ,0.85,0.85],
        **kwargs)

figpath='figures/'
show_zeta(name=figpath+'parabola.eps',opt='ps')
# ff_kin(N=30,dp=80 )
# multiple_scattering(h=5,h0=0,h1=0,name=figpath+'scat1_0.svg',opt='p')
# multiple_scattering(h=5,h0=5,h1=0,name=figpath+'scat1_1.svg',opt='p')
# multiple_scattering(h=5,h0=0,h1=5,name=figpath+'scat1_2.svg',opt='p')
# multiple_scattering(h=5,h0=0,h1=2,name=figpath+'scat2_0.svg',opt='ps')
# multiple_scattering(h=5,h0=2,h1=0,name=figpath+'scat2_1.svg',opt='ps')
# multiple_scattering(h=5,h0=2,h1=3,name=figpath+'scat2_2.svg',opt='ps')
# multiple_scattering(h=5,h0=5,h1=1,name=figpath+'scat3_0.svg',opt='p')
# multiple_scattering(h=5,h0=2,h1=1,name=figpath+'scat3_1.svg',opt='ps')
# multiple_scattering(h=5,h0=2,h1=1,name=figpath+'scat3_1.eps',opt='ps')
