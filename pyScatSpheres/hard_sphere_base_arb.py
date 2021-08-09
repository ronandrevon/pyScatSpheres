import importlib as imp
import scipy.special as spe
import numpy as np
from . import displayStandards as dsp
from . import glob_colors as colors
from . import spherical_utils as spu ;imp.reload(spu)

class HardSphereArrayBaseArb():
    def __init__(self,ka=1,kd_z=2*np.ones((1,1)),kd_y=2*np.ones((1,1)),
        kp=np.inf,alpha=0,k=1,nmax=7,Cp=None,solve=True,copt=1,opt2=1):
        """
        - n,m,ka,kp,kd : nb spheres,wave number inside spheres, normalized radius, normalized distances
        - k : incident wave number
        - nmax : max included order expansion
        """
        if type(kd_z) in [float,int,np.int64,np.float64]:kd_z,kd_y = [kd_z],[kd_y]
        if type(kd_z) in [list]:kd_z,kd_y = np.array(kd_z),np.array(kd_y)
        self.N = kd_z.size
        if type(ka) in [float,int,np.int64,np.float64]:ka=np.array([ka]*self.N)
        self.k,self.nmax = 1,nmax+1
        self.ka,self.kp,self.kd_z,self.kd_y= ka,kp,kd_z,kd_y
        n,m=kd_z.shape[0],kd_y.shape[0]
        self.kd_pz=np.zeros(n)
        self.kd_pz[0]=kd_z[0]
        for i in range(1,n):
            self.kd_pz[i]=kd_z[i]+self.kd_pz[i-1]
        self.d_pz  = self.kd_pz/k
        self.kd_py=np.zeros(m)
        self.kd_py[0]=kd_y[0]
        for i in range(1,m):
            self.kd_py[i]=kd_y[i]+self.kd_py[i-1]
        self.d_py  = self.kd_py/k
        #self.kd_p=np.sqrt(kd_pz**2+kd_py**2)
        #self.d_p=np.sqrt(d_pz**2+d_py**2)
        self.alpha_deg=alpha
        self.alpha=np.deg2rad(alpha)
        if isinstance(Cp,np.ndarray):self.set_Cp(Cp)
        if solve:self.solve(copt=copt,opt2=opt2)

    def _check_idp(self,idp=None):
        if type(idp) in [list,int] :

            if not all(np.array(idp)>=0) and all(np.array(idp)<self.N):
                raise Exception(colors.red+'sphere index not within 0<=idp<N'+colors.black)
        else:
            idp=range(self.N)
        return idp

    def _params(self):
        ss = '$N=%d$, $n_{ref}=%g$, $ka=[$' %(self.kd_z.size,self.kp)
        for i in self.ka :
            pp='%.1f,' %i
            ss+=pp
        ss+='], kd=['
        for i in self.kd_z :
            pp='%.1f,' %i
            ss+=pp
        ss+='], kd_y=['
        for i in self.kd_y :
            pp='%.1f,' %i
            ss+=pp
        ss+=']'
        return ss

    def show_f(self,cart_opt=True,npts=200,r=None,opts='tT',idp=None,fz=np.abs,
        name='',def_args=True,short_title=False,v=0,**kwargs):
        '''
        opts : i(incident), s(scattered), t(total), P(sphere idp),T(All), G(Gradient), F(flux)
        '''
        print(colors.blue+"...compute field %s..." %opts+colors.black)
        ka,kd_z,kd_y,N = self.ka,self.kd_z,self.kd_y,self.N
#LIGNE A CHANGER !!
        if not isinstance(r,tuple) and not isinstance(r,list) and not isinstance(r,np.ndarray):
            #r=(-4*ka+N*kd_y,4*ka+N*kd_z,-N*kd_y+2*ka,N*kd_z+2*ka) #CHANGER ICIIIIIIIIIII
            # r=(1e-3,4*ka.max()+N*kd_y,-2*ka.max(),N*kd_z+2*ka.max())
            r=(1e-3,4*ka.max()+N*kd_y,-2*ka.max(),N*kd_z+2*ka.max())
        # print(r)
        # r,theta,y,z = spu.polar_mesh2D(cart_opt,npts,r)
        r0=r
        x,y,z,r,theta,phi = spu.spherical_mesh2D('yz',npts,r0)
        # print(r0,x.sum(),y.min(),y.max(),z.min(),z.max())
        k,N,nmax = self.k,self.N,self.nmax
        # ap,dp_z,dp_y = self.ka/k,self.kd_z/k,self.kd_y/k
        ap,dp_z,dp_y = self.ka/k,self.kd_pz/k,self.kd_py/k

        #self._check_idp(idp);
        args = {}
        if def_args:
            args = {'lw':2,'labs':['$y$','$z$'],'imOpt':'c','axPos':'V','fonts':{'title':25},}

        if not ('T'  in opts or 'P' in opts) : opts+='T'
        t = np.linspace(-np.pi,np.pi,100)
        ct,st = np.cos(t), np.sin(t)
#LIGNE A CHANGER !!
        #plts = [ [ap*ct, dp*p+ap*st,'k-',''] for p in range(N)]
        plts = [ [dp_y[p]+ap[p]*ct, dp_z[p]+ap[p]*st,'k-',''] for p in range(N)]
        fstr =  r'$%s \psi(r,\theta)$' %(['',r'\partial_r']['G' in opts])
        Gopt =  ''.join([c for c in opts if c in 'GF'] )
        name+=['f','df']['G' in opts]
        if 'P' in opts:
            p=idp[0]
            if v:print("...compute field at spheres p=%d..." %p)
            fi,fs = self.compute_f(r,theta,phi,ftype='a',Gopt=Gopt,idp=idp)
            if 'i' in opts:
                f = fz(fi)
                fig,ax=dsp.stddisp(plts,im=[y,z,f],
                    title = r"Incident %s at sphere %d, $N_{max}=%d$" %(fstr,p,nmax),
                    name=name+'i_sphere%d.png' %p,**args,**kwargs)
            if 's' in opts:
                f = fz(fs)
                fig,ax=dsp.stddisp(plts,im=[y,z,f],
                    title = r"Scattered %s at sphere %d, $N_{max}=%d$" %(fstr,p,nmax),
                    name=name+'s_sphere%d.png' %p,**args,**kwargs)
            if 't' in opts:
                f = fz(fi+fs)
                fig,ax=dsp.stddisp(plts,im=[y,z,f],
                    title = r"Total %s at sphere %d, $N_{max}=%d$" %(fstr,p,nmax),
                    name=name+'t_sphere%d.png' %p,**args,**kwargs)
        if 'T' in opts:
            if v:print("...compute field from all spheres ...")
            #fi,fs = self.compute_f(r,theta,0,ftype='a',Gopt=Gopt)
            fi = self.compute_f(r,theta,phi,ftype='i',Gopt=Gopt)
            fs = self.compute_f(r,theta,phi,ftype='s',Gopt=Gopt)
            params = self._params()
            if short_title:params,fstr = '',[r'$\Psi$',r'$\partial_r\Psi$']['G' in opts]
            # if ''
            if 'i' in opts:
                f = fz(fi)
                fig,ax=dsp.stddisp(plts,im=[y,z,f],
                    title = r"Incident %s, %s " %(fstr,params),
                    name=name+'i.png',**args,**kwargs)
            if 's' in opts:
                f = fz(fs)
                fig,ax=dsp.stddisp(plts,im=[y,z,f],
                    title = r"Scattered %s, %s" %(fstr,params),
                    name=name+'s.png',**args,**kwargs)
            if 't' in opts:
                f = fz(fi+fs)
                fig,ax=dsp.stddisp(plts,im=[y,z,f],
                    title = r"Total %s, %s" %(fstr,params),
                    name=name+'t.png',**args,**kwargs)
        return fig,ax# return f.max(),f.min()
