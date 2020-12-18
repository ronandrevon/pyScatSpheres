'''Gui for scattering the linear array of hard spheres'''
# import importlib as imp
import os,time,copy,easygui,matplotlib, pandas as pd, numpy as np, matplotlib.pyplot as plt
from matplotlib.widgets import CheckButtons
from . import displayStandards as dsp   #;imp.reload(dsp)
from . import glob_colors as colors     #;imp.reload(colors)
from . import qdot_sphere_array as qsa  #;imp.reload(qsa)
from . import hard_sphere_array as hsa  #;imp.reload(hsa)
from . import gui_config as config      #;imp.reload(config)


wm,hm = 0.05,0.07
props = {'lw':2,'fonts':{'title':15,'lab':13,'leg':15,'tick':10},'opt':''}
ax_props = [
    {'axPos':[0.0+wm,0.5+hm,0.5-2*wm,0.5-1.5*hm],},
    # {'axPos':[0.0+wm,0.5+hm,0.5-2*wm,0.5-1.5*hm],'labs':[r'Normalized radius $ka$',r'Normalized Cross section $\sigma/ka^2$']},
    {'axPos':[0.0+wm,0.0+hm,0.5-2*wm,0.5-1.5*hm],'labs':[r'$\Re$',r'$\Im$'  ]},
    {'axPos':[0.5+wm,0.5+hm,0.5-2*wm,0.5-1.5*hm],},#'xylims':['x',0,180]        },
    {'axPos':[0.5+wm,0.0+hm,0.5-2*wm,0.5-1.5*hm],'labs':['$y$','$z$'        ]},
]
for d in ax_props:d.update(props)
labs0 = {'kdka' :r'Normalized distance $kd/ka$',
         'kd'   :r'Normalized distance $kd$',
         'ka'   :r'Normalized radius $ka$',
         'sigka':r'Normalized Cross section $\sigma/ka^2$',
         'sigma':r'Cross section $\sigma$',}
clab0 = {'kdka':r'$kd/ka$','kd':'$kd$','ka':'$ka$','sigka':r'$\sigma/ka^2$','sigma':r'$\sigma$'}
btns_pos={'open':[0.01,0.99],'axis':[0.035,0.99], 'settings':[0.057,0.99],'help':[0.093,0.99]}

class GUI_handler():
    def __init__(self,df_name,Qdot=1,cols=['kd','sigka','ka']):
        self.fig=plt.figure()
        manager = plt.get_current_fig_manager()
        manager.window.showMaximized()
        self.ax = [self.fig.add_subplot(position=ax_props[i]['axPos']) for i in range(len(ax_props))]
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self)
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self)
        self.cid = self.fig.canvas.mpl_connect('motion_notify_event', self)

        self.bg = self.fig.canvas.copy_from_bbox(self.ax[0].bbox)
        self.init(df_name,Qdot,cols,v=1)
        plt.show(block=True )
        # plt.show(block=False)
        # self.refresh()

    def __call__(self, event):
        #### Mouse events
        if type(event)==matplotlib.backend_bases.MouseEvent:
            #### axes
            if event.inaxes==self.ax[0]:
                xydat  = np.array([event.xdata, event.ydata])             #; print(xydat)
                xydat  = self._normalize_dat(xydat)                       #; print(xydat)
                id0    = np.linalg.norm(self.dats-xydat,axis=1).argmin()  #; print(id0)
                xydat0 = self.df.loc[id0][self.cols0].values              #; print(xydat0)
                dist   = np.linalg.norm(self._normalize_dat(xydat0)-xydat)#; print(dist)
                self.update(button=event.button==1,hover=dist<0.1,id0=id0)#; print('event');print('   xydat:',xydat);print('   xydat0:',xydat0, 'id0',id0);print('   dist:',dist)
            #### buttons
            if event.button==1:
                for i,btn in enumerate(self.btns):
                    if btn.contains(event)[0]:self.btns_cb[i]()
            self._plot_btn(int(self.plot_btn.contains(event)[0]),int(event.button==1))
        #### keys
        elif type(event)==matplotlib.backend_bases.KeyEvent:
            if event.key=='f1': self.help()
            if event.key=='f5': self.refresh()
            if event.key==' ' : self.refresh()
            if event.key=='enter': self.change_settings()
            if event.key in ['up','right','left','down']:
                self.arrowActionsParam(event.key)

    ################################################################################################
    # callbacks
    ################################################################################################
    def _ckbtns(self,label):
        '''checkbox callback'''
        if label=='not\ncoupled':
            self.uncoupled = not self.uncoupled
            if self.uncoupled:
                self.s0 = copy.copy(self.s)
                self.s0.solve(copt=0)
        else:
            l = {'re':0,'im':1,r'$\phi$':2,'$||$':3,'$||^2$':4}[label]
            self.fopts[l] = not self.fopts[l]
        self.ax[2].cla()
        self.plot_f()
        self.refresh()

    def _plot_btn(self,hover,btn):
        plotBtnTxt = [[['b','click to plot nf'],['b','click to plot nf']],
                      [['g','click to plot nf'],['r','... plotting ...']]]
        c,txt = plotBtnTxt[hover][btn]
        self.plot_btn.set_backgroundcolor(c)
        self.plot_btn.set_text(txt.ljust(20))
        self.fig.draw_artist(self.plot_btn)
        self.fig.canvas.update()
        self.fig.canvas.flush_events()
        if hover and btn :
            self.plot_fields(field=1)
            c,txt = plotBtnTxt[0][0]
            self.plot_btn.set_backgroundcolor(c)
            self.plot_btn.set_text(txt.ljust(20))
            self.refresh()

    def arrowActionsParam(self,key):
        #### move up or down across the colored curves
        if key == 'up' or key=='down':
            d = self.df.loc[self.id0] #active data point
            idcb0 = np.abs(d[self.cb0col]-self.cb0vals).argmin()      #;print(idcb0)
            if   key == 'up'  :ik = min(self.cb0vals.size-1,idcb0+1)  #;print(ik)
            elif key == 'down':ik = max(0,idcb0-1)                    #;print(ik)
            df0 = self.df[self.df[self.cb0col]==self.cb0vals[ik]]     #;print(ik)
            id0 = np.abs(df0[self.cols0[0]]-d[self.cols0[0]]).argmin()#;print(id0)
            id0 = df0.iloc[id0].name
        #### move left or right along the curves
        elif key == 'right':
            id0 = min(self.id0+1,self.df.shape[0]-1)
        elif key == 'left':
            id0 = max(self.id0-1,0)
        if not id0==self.id0:
            self.update(button=1,hover=0,id0=id0)
        else:
            print('no data point at that value')

    def open(self):
        # datapath = os.path.realpath(os.path.dirname(__file__)+'/../../tests/scattering/data')+'/*.pkl' #; print(datapath)
        datapath = os.path.realpath(os.path.dirname(self.df_name))+'/*.pkl' #; print(datapath)
        df_name = easygui.fileopenbox(msg='open datafile', title='open', default=datapath,
            filetypes=["*.pkl"], multiple=False)
        if df_name:
            print(df_name)
            Qdot = 'qdot' in df_name
            self.load_data(df_name,Qdot,self.cols0+[self.cb0col])
            self.init_plots()

    def set_ax0(self):
        fieldNames = ['xlab','ylab','clab']
        fieldValues = [self.cols0[0],self.cols0[1],self.cb0col]
        msg = "change axis labels. Available values are :\n %s" %list(labs0.keys())
        fieldValues = config.multenterbox(msg,"ax0",fieldValues,fieldNames)
        if fieldValues:
            xlab,ylab,self.cb0col = fieldValues
            self.cols0 = [xlab,ylab]
            self.set_dats()
            self.ax[0].cla()
            self.h0  = self.ax[0].plot(0,0,'s',color='k' ,markersize=12,animated=True)[0]
            self.h1  = self.ax[0].plot(0,0,'o',color='g' ,markersize=10,animated=True)[0]
            self.update_cb0()
            self.plot_s()
            self.refresh()

    def change_settings(self):
        fieldNames = ['opts3', 'npts','range nY','range nZ', 'cmap(sigma)','cmap(l)','cmap(near field)','cmax ']
        fieldValues = [self.opts3,str(self.npts),str(self.r[0]),'%.2f,%.2f' %(self.r[1],self.r[2]), self.cmS,self.cmF,self.cmL, str(self.cmax)]
        fieldValues = config.multenterbox("Change settings","settings", fieldValues,fieldNames)
        if fieldValues:
            self.opts3,npts,ry,rz,self.cmS,self.cmF,self.cmL,cmax = fieldValues
            self.npts = max(10,min(int(npts),1000))
            self.r = np.array([ry]+rz.split(','),dtype=float)
            # self.cmax = max(0,int(cmax))
            # self.cbs[3].set_caxis([-self.cmax,self.cmax])
            self.refresh()

    def help(self):
        print(colors.green+'Commands : '+colors.black)
        print(colors.yellow+config.help_cmd+colors.black)
        easygui.msgbox(config.help_msg+config.help_cmd,title='help') #text='More info at : \n'+url)

    ##############################################################################################################
    ##### displays
    ##############################################################################################################
    ##ax0
    def plot_s(self):
        # print("...plotting cross sections")
        plts = []
        self.csS = dsp.getCs(self.cmS,self.cb0vals.size+1)
        for i,val in enumerate(self.cb0vals):
            dfa = self.df[self.df[self.cb0col]==val]
            plts += [[dfa[self.cols0[0]], dfa[self.cols0[1]],[self.csS[i],'-o' ],'']] #,'$ka=%.1f$' %ka]]
        xM,yM = self.df[self.cols0].max()
        xm,ym = self.df[self.cols0].min()
        labs = labs0[self.cols0[0]],labs0[self.cols0[1]]
        dsp.stddisp(plts,ax=self.ax[0],labs=labs,xylims=[0.99*xm,1.02*xM,0.99*ym,1.02*yM],**ax_props[0])#,pargs={'picker':1,'pickradius':1})

    ##ax1
    def plot_cp(self):
        # print('...plotting coefficients...')
        nmax = self.s.nmax
        idp = self.ip*nmax
        idl = np.arange(nmax)
        val0 = self.df.loc[self.id0][self.cb0col]
        coeff = ['Cp','bp'][self.Qdot]

        #### all coefficients on same curve as active point
        self.csL = dsp.getCs(self.cmL,nmax)
        df = self.df[self.df[self.cb0col]==val0] #;print(df)
        cp = np.vstack(df[coeff].values)         #;print(cp.shape) cause df.Cp.values=[Cp0,Cp1,...] is a list
        mm = np.abs(cp).max()
        plts  = [[ cp[:,idp+l].real ,cp[:,idp+l].imag,[self.csL[l],'-o' ]] for i,l in enumerate(idl)]
        if self.uncoupled:
            cp0 = self.s0.__getattribute__(coeff)
            plts+=[ [cp0[l].real,cp0[l].imag,[self.csL[l],'s'],''] for l in idl ]

        #### active point coefficients
        cpi = self.df.loc[self.id0][coeff][self.ip*nmax+idl] #;print(cpi)
        plts+=[[cpi.real,cpi.imag,'go','']]
        tle = '$b_{pl}$ for $%s=%g$ and varying $%s$ ' %(self.cb0col,val0,self.cols0[0])
        dsp.stddisp(plts,ax=self.ax[1],xylims=[-mm,mm,-mm,mm],title=tle,**ax_props[1])
    ##ax2
    def plot_f(self):
        # print('...plotting scatering amplitude...')
        fopts = ''.join(np.array(['r','i','a','m','2'])[self.fopts])
        legElt={}
        if self.uncoupled:
            self.s0.show_ff(npts=361,fopts=fopts,ax=self.ax[2],leg=0,title=None,**ax_props[2])
            legElt={'uncoupled':'k--'}
        self.s.show_ff(npts=361,fopts=fopts,ax=self.ax[2],legElt=legElt,**ax_props[2])
    ##ax3
    def plot_fields(self, field=1):
        ka,kd = self.s.ka,self.s.kd
        em = 0.05
        if self.r[0]==-1:
            ry = 2+kd/ka+self.r[1]+self.r[2]
        else:
            ry = 1+self.r[0]
        r = np.array([1e-3,ry*ka,-ka*(1+self.r[1]),kd+ka*(1+self.r[2])])
        xylims = [0,r[1]-em,r[2],r[3]-em]
        if field:
            # print('...plotting near fields...')
            self.s.show_f(cart_opt=True,npts=self.npts,r=r,opts=self.opts3,def_args=False,fz=np.real,
                xylims=xylims,cmap=self.cmF,imOpt='c',caxis=[-self.cmax,self.cmax],cb_pos=[0.955,hm,0.015,0.5-1.5*hm],fig=self.fig,
                cb=self.cbs[3],ax=self.ax[3],**ax_props[3])
        else:
            t = np.linspace(-np.pi/2,np.pi/2,100)
            ct,st = np.cos(t), np.sin(t)
            plts = [ [ka*ct, kd*p+ka*st,'k-',''] for p in range(self.s.N)]
            dsp.stddisp(plts,xylims=xylims,ax=self.ax[3],**ax_props[3])

    ##############################################################################################################
    #### initialization
    ##############################################################################################################
    def init(self,df_name,Qdot,cols,v=0):
        print(colors.green+'Initialization'+colors.black)
        t=time.time()
        self.init_params()          #;print("init_params %.3f" %(time.time()-t));t=time.time()
        self.load_data(df_name,Qdot,cols)#;print("load %.3f" %(time.time()-t));t=time.time()
        self.init_widgets()                ;print("\twidget     :%.3f" %(time.time()-t));t=time.time()
        self.init_colorbars()              ;print("\tcb         :%.3f" %(time.time()-t));t=time.time()
        self.init_plots()           #;print("init_plots %.3f" %(time.time()-t));t=time.time()
        print(colors.green+'Completed'+colors.black)

    def load_data(self,df_name,Qdot,cols=['kdka','sigka','ka']):
        self.df_name=df_name
        self.Qdot=Qdot
        if Qdot:
            self.SphereArray = qsa.QdotSphereArray
        else:
            self.SphereArray = hsa.HardSphereArray
        self.df  = pd.read_pickle(df_name)       #; print(self.df)
        self.df['sigka'] = self.df.sigma/self.df.ka**2
        self.df['kdka'] = self.df.kd/self.df.ka
        self.cols0  = cols[:2]
        self.cb0col = cols[2]
        self.set_dats()
        d = self.df.loc[0]
        self.s = self.SphereArray(N=d.N,ka=d.ka,kp=d.kp,nmax=d.nmax-1,solve=True)
        self.s0 = copy.copy(self.s)
        self.s0.solve(copt=0)

    def set_dats(self):
        dats = self.df[self.cols0].values.T
        self.datsmM  = (dats[0].min(), dats[0].max(),  dats[1].min(), dats[1].max()) #;print(self.datsmM )
        self.dats    = self._normalize_dat(dats).T
        self.cb0vals = self.df[self.cb0col].unique()
    def _normalize_dat(self,dats):
        xm,xM,ym,yM = self.datsmM
        dats[0] = (dats[0] - xm)/(xM-xm)
        dats[1] = (dats[1] - ym)/(yM-ym)
        return dats

    def init_params(self):
        self.opts3 = 'tT'  #options for field to plot
        self.npts = 200
        self.r = (-1,1,1)   #range for extend of field to plot
        self.id0=10        #active case
        self.cmS,self.cmF,self.cmL = 'viridis', 'jet','Spectral'
        self.ip=0          #active sphere coeffs
        self.uncoupled = False
        self.fopts = [False]*3+[True,False]   #riam2
        self.cmax = 1

    def init_widgets(self):
        # interactive buttons,checkbox
        self.btns_cb=[self.open,self.set_ax0,self.change_settings,self.help]
        self.btns=[]
        for k,pos in btns_pos.items():
            self.btns += [self.fig.text(pos[0],pos[1],k,size=12,ha='left',va='top',backgroundcolor=(0.8,0.8,0.8))]
        self.plot_btn = self.fig.text(0.99,0.5,'click to plot nf'.ljust(20)+'!',size=12,ha='right',va='bottom',backgroundcolor='b',alpha=0.9)
        self.ckbtns = CheckButtons(self.fig.add_axes([0.955,0.5+hm,0.1,0.5-1.5*hm]), ['not\ncoupled','re','im',r'$\phi$','$||$','$||^2$'],[self.uncoupled]+self.fopts)
        self.ckbtns.on_clicked(self._ckbtns)

    def init_colorbars(self):
        # print('...init colorbars...')
        args = {'fig':self.fig,'is_3d':0,'cb':None,'cs':None}#,'shading':'auto'}
        cm, nmax = self.cmax, self.s.nmax
        nkps,kpM = self.cb0vals.size,self.cb0vals.max()

        self.cbs = np.zeros((4),dtype=object)
        self.cbs[0] = dsp.add_colorbar(ax=self.ax[0],cmap=self.cmS,caxis=[0,nkps],cb_pos=[0.455,0.5+hm,0.015,0.5-1.5*hm], values=np.arange(nkps), boundaries=np.arange(nkps)-0.5,ticks=np.arange(nkps),**args)
        self.cbs[1] = dsp.add_colorbar(ax=self.ax[1],cmap=self.cmL,caxis=[0,nmax],cb_pos=[0.455,hm    ,0.015,0.5-1.5*hm], values=np.arange(nmax), boundaries=np.arange(nmax)-0.5,ticks=np.arange(nmax),**args)
        self.cbs[3] = dsp.add_colorbar(ax=self.ax[3],cmap=self.cmF,caxis=[-cm,cm],cb_pos=[0.955,hm    ,0.015,0.5-1.5*hm], **args)
        self.cbs[1].ax.set_xlabel('$ l$'  ,color='k',fontsize=15)
        self.cbs[3].ax.set_xlabel('$ |f|$',color='k',fontsize=15)
        self.update_cb0()

        #interactive markers
        self.h0  = self.ax[0].plot(0,0,'s',color='k' ,markersize=12,animated=True)[0]
        self.h1  = self.ax[0].plot(0,0,'o',color='g' ,markersize=12,animated=True)[0]
        self.h11 = self.ax[1].plot(0,0,'o',color='g' ,markersize=5)[0]

    def update_cb0(self):
        nkps,kpM = self.cb0vals.size,self.cb0vals.max()
        tlabs = ['%.1f' %v for v in self.cb0vals]
        # tlabs = ['%g' %v for v in self.cb0vals]
        self.cbs[0].cmap = self.cmS
        self.cbs[0].max = nkps
        self.cbs[0].values = np.arange(nkps)
        self.cbs[0].boundaries = np.arange(nkps+1)-0.5
        self.cbs[0].ticks = np.arange(nkps)
        self.cbs[0].ax.set_xlabel('%s' %clab0[self.cb0col],color='k',fontsize=15)
        self.cbs[0].set_ticklabels(tlabs)
        self.cbs[0].update_ticks()

    # def update_cb3(self):
    def init_plots(self):
        self.clear()
        t=time.time()
        self.plot_s()                      ;print("\tplot_s     :%.3f" %(time.time()-t));t=time.time()
        self.plot_cp()                     ;print("\tplot_cp    :%.3f" %(time.time()-t));t=time.time()
        self.plot_fields(field=0)          ;print("\tplot_fields:%.3f" %(time.time()-t));t=time.time()
        self.update(button=1,hover=0,id0=0);print("\tupdate     :%.3f" %(time.time()-t));t=time.time()

    ##############################################################################################################
    ##### Update
    ##############################################################################################################
    def clear(self):
        for ax in self.ax:ax.cla()
        for i,ax in enumerate(self.ax):dsp.stddisp(ax=ax,**ax_props[i])
    def refresh(self):
        self.fig.canvas.draw()
        # self.bg = self.fig.canvas.copy_from_bbox(self.ax[0].bbox)
        self.bg = self.fig.canvas.copy_from_bbox(self.fig.bbox)
        print(colors.green+'figure refreshed'+colors.black)
    def update_ax0(self,i=0):
        self.fig.canvas.restore_region(self.bg)
        self.ax[i].draw_artist(self.h0)
        self.ax[i].draw_artist(self.h1)
        self.fig.canvas.update() #; self.fig.canvas.blit(self.fig.bbox)
        self.fig.canvas.flush_events()

    def update(self,button,hover,id0):
        d = self.df.loc[id0]
        col = self.csS[np.abs(self.cb0vals-d[self.cb0col]).argmin()]
        if button:
            if self.Qdot :
                cp = np.hstack([d.ap,d.bp])
            else:
                cp = d.Cp
            self.s = self.SphereArray(N=d.N,ka=d.ka,kp=d.kp,kd=d.kd,nmax=d.nmax-1,Cp=cp)
            if self.uncoupled:
                self.s0 = copy.copy(self.s)
                self.s0.solve(copt=0)
            self.id0 = id0                              #;t=time.time()
            self.ax[1].cla();self.plot_cp()             #;print('cp:',time.time()-t);t=time.time()
            self.ax[2].cla();self.plot_f()              #;print('ff:',time.time()-t);t=time.time()
            self.ax[3].cla();self.plot_fields(field=0)  #;print('nf:',time.time()-t);t=time.time()
            self.refresh()                              #;print('refresh:',time.time()-t);t=time.time()
            #update ax0
            self.h0.set_data(None,None)
            self.h1.set_data(d[self.cols0[0]],d[self.cols0[1]])#; print(self.h1)
            self.h1.set_color(col)                      #;print('misc:',time.time()-t);t=time.time()
            self.update_ax0()#self.h1)

        elif hover and not button:
            self.h0.set_data(d[self.cols0[0]],d[self.cols0[1]])
            self.h0.set_color(col)
            self.update_ax0()#self.h0)
        else:
            self.h0.set_data(None,None)
            self.update_ax0()#self.h0)
