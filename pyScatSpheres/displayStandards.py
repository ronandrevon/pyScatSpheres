'''Display utilities'''
import numpy as np
import matplotlib,os
import matplotlib.pyplot as plt
# import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
from matplotlib.collections import PatchCollection as PatchColl
from matplotlib.patches import Rectangle
from scipy import ndimage
from subprocess import check_output
from PIL import Image
# from .glob_colors import*
from . import glob_colors as colors
# Get screen info, remove toolbar
# matplotlib.rcParams['toolbar'] = 'None'
matplotlib.rcParams['backend'] = 'GTK3Agg'
matplotlib.rcParams['pcolor.shading'] = 'auto'
matplotlib.rc('text', usetex=True)
# matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
dpi = check_output("xdpyinfo | awk '/resolution/{print $2}'",shell=True).decode()
dpi = np.array(dpi.strip().split('x'),dtype=int)
screen_size = check_output("xrandr | grep \"*+\" | awk '{print $1}'",shell=True).decode().split('\n')
#choose second monitor if present
screen_size=screen_size[int(len(screen_size)>2)]
screen_size = np.array(screen_size.split('x'),dtype=int)/dpi #inches


###########################################################################
#def : plotting standard
###########################################################################
def standardDisplay(ax,labs=['','',''],name='', xylims=[], axPos=1,legOpt=None,view=None,
    fonts=[30,25,15,20], c=['k','k'], logOpt='',changeXYlims=True,is_3d=0,
    gridOn=True,gridmOn=False,ticksOn=True,title='',legLoc='best',legElt=[],
    figopt='1',xyTicks=None,xyTicksm=None,xyTickLabs=None,mg=0.05,opt='p',setPos=True,equal=False,
    pOpt=None):
    '''
    legElt : dict - {'lab1':[(1,0,0),'s-']} overrides existing handles
    opt : p(plot), s(save), c(close) '
    pOpt : p(setPos)X(changeXYlims)l(legOpt)t(ticksOn)e(equal),G(gridOn),g(gridmOn)
    view : [elev,azim] for 3D projection axes
    '''
    is_3d = "3D" in str(ax.__class__)
    if isinstance(pOpt,str):
        setPos,changeXYlims,legOpt,ticksOn,equal,gridOn,gridmOn = [s in pOpt for s in 'pXlteGg']
    axPos = get_axPos(axPos); #print(changeXYlims)
    ax.name = name;
    if isinstance(fonts,dict) : fonts = get_font(fonts,dOpt=True)
    fs,fsLeg,fsL,fsT = fonts; #print(fsLeg)
    #axis, labels, position
    if 'x' in logOpt : ax.set_xscale('log')
    if 'y' in logOpt : ax.set_yscale('log')
    ax.set_xlabel(labs[0],fontsize=fs   ,color=c[0])
    ax.set_ylabel(labs[1],fontsize=fs   ,color=c[1])
    if is_3d : ax.set_zlabel(labs[2],fontsize=fs   ,color=c[1])
    ax.axis(['off','on'][ticksOn])
    ax.set_axisbelow(True)
    if setPos : ax.set_position(axPos);
    #lims,ticks and grid
    xylims=changeAxesLim(ax,mg,xylims,is_3d,changeXYlims)#;print(changeXYlims,xylims)
    change_ticks(ax,xyTicks,xyTickLabs,xylims,is_3d,xyTicksm)
    ax.tick_params('x',labelsize=fsL    ,colors=c[0],direction='in')
    ax.tick_params('y',labelsize=fsL    ,colors=c[1],direction='in')
    if is_3d : ax.tick_params('z',labelsize=fsL,colors=c[1])
    ax.grid(False);#print(gridOn,gridmOn)
    if gridOn  : ax.grid(gridOn,which='major',color=(0.9,0.9,0.9),linestyle='-')
    if gridmOn :
        ax.minorticks_on()
        ax.grid(gridmOn,which='minor',color=(0.95,0.95,0.95),linestyle='--')
    if view : ax.view_init(elev=view[0], azim=view[1])
    #title,legend,save
    if equal :ax.axis(['equal','auto'][is_3d]);
    ax.set_title(title, {'fontsize': fsT}); #print(title)
    addLegend(ax,fsLeg,legOpt,legLoc,legElt)
    disp_quick(name,ax,opt,figopt)

def stddisp(plots=[],scat=None,texts=[],colls=[],patches=[],im=None,rot=0,
            tricont=None,contour=None,quiv=None,surfs=[],caxis=None,
            cs=None,lw=1,ms=5,marker='o',fonts={},axPos=1,imOpt='',cb_pos=[],cmap='jet',
            ax=None,fig=None,cb=None,figsize=(0.5,1),pad=None,rc='11',
            inset=None,iOpt='Gt',name='',opt='p',figopt='',pargs={},
            std=1,**kwargs):
    '''
    fonts.keys=['lab','leg','tick','text','title']
    lw,ms,marker : linewidth,markersize,marker to apply to all
    axPos   : 4-list or int - ax position (see get_axPos)
    contour : list - x,y,z,[levels]
    quiv    : functor to gradient that will be plotted on the isocontours (must return shape (2,N,N))
    imOpt   : c(colorbar) h(horizontal) v(vertical)
    cs      : C(contours),I(image),S(scatter),N(None) artist to use for colorbar
    iOpt    : inset displayStandards keys - ['axpos','xylims','ms','lw','lwp']
    kwargs  : standardDisplay
    '''
    if isinstance(fonts,dict) : fonts = get_font(fonts)
    fsT = fonts[3]; fonts=(np.array(fonts)[[0,1,2,4]]).tolist()
    if not ax : fig,ax = create_fig(figsize,pad,rc)
    if ax and not axPos : axPos=ax.get_position().bounds
    #add lines,patches,texts,scatters
    is_3d = "3D" in str(ax.__class__) # not isinstance(ax,matplotlib.axes._subplots.Subplot) and not isinstance(ax,matplotlib.axes._axes.Axes)#;print(is_3d)
    for coll in colls : ax.add_collection(coll)
    for pp in patches : ax.add_patch(pp)
    pltPlots(ax,plots,lw,ms,is_3d,pargs)
    pltTexts(ax,texts,fsT,is_3d)

    #note that at them moment only one cb per ax authorized
    # priority order here : scat,im,contour
    css = dict(zip(['I','S','C','N'],[None]*4))
    if not contour==None:css['C'] = plt_contours(ax,contour,quiv,cmap,lw)
    if not im==None:css['I']      = pltImages(ax,im,cmap,caxis,imOpt,rot)
    if not scat==None:css['S']    = pltScatter(ax,scat,'b',ms,marker,rc=='3d',cmap)
    # if 'c' in imOpt and contour==None and scat==None and im==None:cs=None
    if isinstance(cs,str) : cs=css[cs]
    add_colorbar(fig,ax,cs,cb,imOpt,cmap,is_3d,caxis=caxis,cb_pos=cb_pos)
    if is_3d :pltSurfs(ax,surfs)

    if isinstance(inset,dict):
        inset=fill_inset_dict(inset)
        ec = inset['ec']
        in_box = inset['xylims']
        xy,Wrect,Hrect = (in_box[0],in_box[2]),in_box[1]-in_box[0],in_box[3]-in_box[2]
        ax.add_patch(Rectangle(xy,Wrect,Hrect,fill=0,ec=ec,lw=inset['lwp']))
        ax2 = add_inset(fig,plots,inset['xylims'],inset['axpos'],inset['lw'],inset['ms'],iOpt)
        for axis in ['top','bottom','left','right']:
            ax2.spines[axis].set_linewidth(3)
            ax2.spines[axis].set_color(ec)
    if std:
        standardDisplay(ax,is_3d=is_3d,axPos=axPos,fonts=fonts,
            name=name,opt=opt,figopt=figopt,**kwargs)
    else:
        disp_quick(name,ax,opt,figopt)
    return fig,ax


def image_bg(im,rot=0,name='',opt='p',figopt='',**kwargs):
    ''' display an image in background to data
    fig = dsp.image_bg('image.png',xylims=[0,1,0,1])
    NOTE :
    use fix_pos(fig) after the figure has been displayed
    '''
    fig,ax = stddisp(im=im,pOpt='',opt='',rot=rot)
    ax1 = fig.add_axes([0,0,1,1],frameon=False)
    fig,ax1 = stddisp(fig=fig,ax=ax1,pOpt='GtX',opt='',**kwargs)
    fix_pos(fig)
    disp_quick(name,ax,opt,figopt)
    return fig

def fix_pos(fig):
    fig.axes[1].set_position(fig.axes[0].get_position())
    fig.canvas.draw()
    return fig

def display_solutions(disp_d,objs=[],key_pair=(),help=0,**kwargs):
    ''' display a pair data from different objects (used to compare methods)
    - disp_d : dict - {key_sol : (sol_function_str,color,display_label))
    - objs : list of tuples - [(obj,cml,leg)...] or [obj...]
    - key_pair : (str,str or list of str) where str in disp_d.keys()
    - help : for the list of available keys
    '''

    if help :print('available keys:\n',disp_d.keys());return


    kx,kys = key_pair[0],key_pair[1]
    fx_str,c,xlab = disp_d[kx]
    if isinstance(kys,str):kys=[kys] #ylab = disp_d[al[kys]][0]

    #legend for the quantities displayed
    legElt = {disp_d[ky][2]:disp_d[ky][1] for ky in kys}

    # keys legend
    plts,ylab = [],''
    for obj_t in objs:
        if isinstance(obj_t,tuple):
            obj,cml,leg = obj_t
            legElt[leg] = cml
            ml = ''.join(getCML(cml)[1:3])
        else:
            obj,ml=obj_t,'-'
        for ky in kys:
            fy_str,c = disp_d[ky][0:2]
            fx = obj.__getattribute__(fx_str)
            fy = obj.__getattribute__(fy_str)
            plts += [[fx(),fy(),[c,ml],'']]

    #display
    labs = [xlab,ylab]
    fig,ax = stddisp(plts,labs=labs,legElt=legElt,**kwargs)
    return fig,ax,labs,legElt

def add_colorbar(fig,ax,cs=None,cb=None,imOpt='c',cmap='jet',is_3d=0,l=0.03,L=0.85,caxis=None,cb_pos=[],**kwargs):
    if 'c' in imOpt :
        if not cs :
            if caxis==None:caxis=[0,1]
            cs = plt.cm.ScalarMappable(plt.Normalize(vmin=caxis[0],vmax=caxis[1]),cmap=cmap)
        if is_3d:
            cb = fig.colorbar(cs,ax=ax,vmin=caxis[0],vmax=caxis[1],cmap=cmap)
        else:
            if cb:
                cb.mappable=cs
            else:
                if not cb_pos : cb_pos = [[0.9,0.1,l,L],[0.1,0.9,L,l]]['h' in imOpt]
                orient,tickloc=[['vertical','right'],['horizontal','top']]['h' in imOpt]
                ax_cb = fig.add_axes(cb_pos)
                cb=fig.colorbar(cs,ax=ax,cax=ax_cb,orientation=orient,ticklocation=tickloc,**kwargs)
        return cb
#######################################################################
###### Old one
def stdDispPlt(plots,xlabs=['',''],name='',xlims=[],axPos=[],c=['k','k'],
               showOpt=1,ax=None,fig=None,texts=[],legOpt=1,lw=1,ms=5,fonts={},
               logOpt='',changeXYlims=True,gridOn=True,ticksOn=True,fullsize=False,title='',
               legLoc='upper left',legElt=[],xyTicks=[],colls=[],opt='',figsize=(9,9)):
    '''fonts.keys=['lab','leg','tick','text','title']'''
    if isinstance(fonts,dict) : fonts = get_font(fonts)
    fsT = fonts[3]; fonts=(np.array(fonts)[[0,1,2,4]]).tolist()
    if showOpt and not 'q' in opt : opt+='p'
    if not ax : fig,ax = plt.subplots(figsize=figsize)
    if ax and axPos==[] : axPos=ax.get_position().bounds
    if fullsize : mng = plt.get_current_fig_manager(); mng.resize(*mng.window.maxsize())
    for coll in colls : ax.add_collection(coll)
    pltPlots(ax,plots,lw,ms)
    pltTexts(ax,texts,fsT)
    standardDisplay(ax,xlabs,name,xlims,axPos ,legOpt,fonts,c,
                    logOpt=logOpt,changeXYlims=changeXYlims,gridOn=gridOn,title=title,
                    legLoc=legLoc,legElt=legElt,xyTicks=xyTicks,opt=opt)
    return fig,ax

def addyAxis(fig,ax,plots,yLab='',c='k', lw=1,ms=5,axPos=[],showOpt=1,yTicks=[],yTickLabs=[],
             **kwargs):
    axPosI = axPos if len(axPos)==4  else [0.13, 0.15, 0.75, 0.75]
    ax2 = ax.twinx()
    pltPlots(ax2,plots,lw)
    standardDisplay(ax2,['',yLab],c=['k',c],axPos=axPosI,xyTicks=[[],yTicks],xyTickLabs=[[],yTickLabs],
                    **kwargs)
    ax.set_position(axPosI)
    # if showOpt : fig.show()

def addxAxis(ax,plots,xLab='',c='k', lw=1,axPos=[],
            xTicks=[],xTickLabs=[], **kwargs):
    axPosI = axPos if len(axPos)==4  else [0.13, 0.15, 0.8, 0.675]
    ax2 = ax.twiny()
    pltPlots(ax2,plots,lw)
    standardDisplay(ax2,[xLab,''],c=[c,'k'],axPos=axPosI,xyTicks=[xTicks,[]],xyTickLabs=[xTickLabs,[]],**kwargs)
    ax2.tick_params('x',color=c)
    ax2.set_xlabel(xLab, verticalalignment='bottom')
    ax.set_position(axPosI)
    #if showOpt : plt.show()

def add_inset(fig,plots,xylims,axPosI,lw=2,ms=5,iOpt='GX'):#,**kwargs)
    ax2 = fig.add_axes(axPosI)
    pltPlots(ax2,plots,lw,ms) #; print(xylims)
    standardDisplay(ax2,labs=['',''],xylims=xylims,setPos=False,
        xyTicks=None,xyTickLabs=[[],[]],legOpt=0,pOpt=iOpt)#,**kwargs)
    return ax2

########################################################################
### def : plot calls
########################################################################
def pltPlots(ax,plots,lw0,ms0=5,is_3d=0,pargs={}):
    ''' plots a list of plots :
        - [x,y,color]([x,y,z,color] if is_3d) + [label](optional)
    '''
    if plots :
        if len(plots)>=3+is_3d and isinstance(plots[2+is_3d],str):
            plots = [plots]; #print(len(plots))
    if is_3d :
        for p in plots:
            cml = 'b' if len(p)<4 else p[3]
            lab = ''  if len(p)<5 else p[4]
            c,m,l = getCML(cml)
            ax.plot(p[0],p[1],p[2],label=lab,color=c,linestyle=l,marker=m, linewidth=lw0,markersize=ms0)
    else:
        for p in plots:
            cml = 'b' if len(p)<3 else p[2]
            lab = ''  if len(p)<4 else p[3]
            lw  = lw0 if len(p)<5 else p[4]
            c,m,l = getCML(cml)
            ax.plot(p[0],p[1],label=lab,color=c,linestyle=l,marker=m, linewidth=lw,markersize=ms0,**pargs)
def pltTexts(ax,texts,fsT,is_3d=0):
    ''''texts format : [x,y,text,color], [x,y,z,txt,color] if is_3d'''
    if texts:
         if not isinstance(texts[0],list) : texts = [texts]
    if is_3d :
        for t in texts:
            ax.text(t[0],t[1],t[2],t[3],fontsize=fsT,color=t[4])
    else:
        for t in texts:
            c_t = 'k' if len(t)<4 else t[3]
            ax.text(t[0],t[1],t[2],fontsize=fsT,color=c_t)
def pltImages(ax,im=None,cmap='viridis',caxis=None,imOpt='',rot=0):
    cs = None
    if isinstance(im,str):
        if '.npy' in im:
            im = np.load(im)
            if len(im.shape)==3:im=list(im)
        else:
            if 'g' in imOpt:
                image = Image.open(fname).convert("L")
                image = np.asarray(image)/255
                cs = ax.pcolor(im,cmap='gray',)
                im = None
            else:
                image = plt.imread(im)
                if rot:image=ndimage.rotate(image,rot)
                cs=ax.imshow(image,cmap=cmap)#,**im_args)#,origin='upper')
                im = None
    if isinstance(im,np.ndarray): im=[im]
    if isinstance(im,list):
        if len(im)>2:
            x,y,z = im[:3]
            N = int(np.sqrt(x.shape[0]))
            # dx,dy = abs(x[0,0]-x[0,1])/2, abs(y[0,0]-y[1,0])/2
            args = {}
            if caxis : args = {'vmin':caxis[0],'vmax':caxis[1]}
            if len(im)>3:args['alpha']=im[3]
            # cs=ax.pcolormesh(x-dx,y-dy,z,cmap=cmap,edgecolors='none',ec=None,**args)
            cs=ax.pcolormesh(x,y,z,shading='auto',cmap=cmap,edgecolors='none',ec=None,**args)
        else :
            # cs=ax.pcolor(im[0],cmap=cmap)
            cs=ax.imshow(im[0],cmap=cmap)
    return cs
def pltScatter(ax,scat,c='b',s=5,marker='o',proj_3D=False,cmap='jet') :
    '''
    - scat : [x,y,<z>,s,c] or [x,y,<z>,c]
    - s : int or list/np.array of int
    - c : tuple,str or list/np.array of tuple/str
    '''
    cs = None
    if len(scat) :
        if proj_3D:
            x,y,z = scat[:3]
        else :
            x,y = scat[:2]
        # color and marker size
        sc = scat[2+proj_3D:]
        if len(sc)==1 : c = sc[0]
        elif len(sc)==2: s,c = sc
        N=np.array(x).size
        if isinstance(s,int) : s = [s]*N
        if isinstance(c,tuple) or isinstance(c,str) : c = [c]*N
        if proj_3D :
            cs=ax.scatter3D(x,y,z,s=s,c=c,marker=marker,cmap=cmap)
        else :
            cs=ax.scatter(x,y,s,c,marker=marker,cmap=cmap)
    return cs

def pltSurfs(ax,surfs,c='b',a=0.2,lw=1,ec='b'):
    '''surf : [x,y,z] or [x,y,z,**kwargs] or [x,y,z,c,a,lw,ec]
    '''
    for surf in surfs:
        x,y,z = surf[:3]
        if len(surf)>3 :
            if isinstance(surf[3],dict):
                kwargs=surf[3]
                ax.plot_surface(x,y,z,**kwargs)
                return
            else:
                c,a,lw,ec = surf[3:]
        ax.plot_surface(x,y,z,color=c,alpha=a,linewidth=lw,edgecolor=c)

def plt_contours(ax,contour,quiv,cmap,lw):
    cs = ax.contour(*contour,cmap=cmap,linewidths=lw)
    if not quiv == None:
        v_cs = get_iso_contours_vertx(cs)
        x_cs,y_cs = v_cs.T
        dfx,dfy = quiv(x_cs,y_cs)
        ax.quiver(x_cs,y_cs,dfx,dfy,linewidths=lw,units='xy')
    return cs

########################################################################
# handles and properties
########################################################################
def create_fig(figsize=(0.5,1),pad=None,rc='11') :
    '''figsize :
        tuple : (width,height) normalized units
        str   : f(full),12(half),22(quarter)
        rc    : str or list - layout arrangement : '3d','11','22',...
    '''
    if isinstance(figsize,str) :
        figsize = {'f':(1,1),'12':(0.5,1),'21':(1,0.5),'22':(0.5,0.5)}[figsize]
    if isinstance(rc,str) :
        rc = {'fig':'fig','3d':'3d','11':[1,1],'21':[2,1],'12':[1,2],'22':[2,2]}[rc];#print(rc)
    figsize = tuple(np.array(figsize)*screen_size)
    if rc=='3d':
        wh = min(screen_size)
        figsize = (wh,wh);#print(figsize)
        fig = plt.figure(figsize=figsize,dpi=dpi[0])
        ax  = plt.subplot(111,projection='3d')
    elif rc=='fig':
        fig=plt.figure(figsize=figsize,dpi=dpi[0])
        return fig
    else:
        fig,ax = plt.subplots(nrows=rc[0],ncols=rc[1],figsize=figsize,dpi=dpi[0])
    ax.set_navigate(True)
    if pad : plt.tight_layout(pad)
    return fig,ax

def get_font(d_font=dict(),dOpt=False) :
    '''keys = ['lab','leg','tick','text','title']
    '''
    keys = ['lab','leg','tick','text','title']
    # vals = [30,25,15,20,30]
    vals = [25,20,15,20,20]
    font_dict = dict(zip(keys,vals))
    for k in d_font.keys() : font_dict[k] = d_font[k]
    if dOpt : keys = ['lab','leg','tick','title']
    fonts = [ font_dict[k] for k in keys]
    return  fonts

def fill_inset_dict(inset_dict=dict()):
    '''keys = ['axpos','xylims','ms','lw','lwp']'''
    keys = ['axpos','xylims','ms','lw','lwp','ec']
    vals = [[0,0,0.25,0.25],None,30,3,2,'g']
    full_dict = dict(zip(keys,vals))
    for k,v in inset_dict.items() : full_dict[k]=v
    return full_dict

def getCML(C):
    ml = ''
    if isinstance(C,list):
        c,ml = C[0], C[1]
    elif isinstance(C,tuple) or isinstance(C,np.ndarray):
        c,m,l = C, '','-'
    elif isinstance(C,str):
        c,m,l,ml = C[0],'','-',C[1:]
    else:
        c,m,l= 'b','','-'

    if not ml == '':
        l = [char for char in ml if char in '-.:']    ;l=''.join(l)
        m = [char for char in ml if char not in '-.:'];m=''.join(m)
    return c,m,l

def get_legElt(legElt):
    legE=[]
    if isinstance(legElt,dict):
        for lab,cml in legElt.items():
            c,m,l=getCML(cml)
            legE+=[Line2D([0],[0],linewidth=2,linestyle=l,color=c,marker=m,label=lab)]
    return legE

def addLegend(ax,fsLeg,legOpt=None,loc='upper left',legElt=[]):
    hs,labs = ax.get_legend_handles_labels()
    legElt = get_legElt(legElt)
    if legOpt == None : legOpt = any(labs) or any(legElt)
    if legOpt:
        out = None
        if loc[-4:]==' out' :
            loc = loc[:-4]
            if loc=='upper left' : out = (1, 1)
            elif loc=='center left' : out = (1, 0.5)

        leg = ax.legend(handles=hs+legElt,fontsize=fsLeg,loc=loc,bbox_to_anchor=out)
        leg.set_draggable(True)         # matplotlib 3.1
        #leg.draggable(True)            # matplotlib 2.2
        #if leg : leg.DraggableLegend() # matplotlib?

def getCs(name,N,cmopt=0):
    if cmopt:
        cs = matplotlib.cm.get_cmap(name)
    else:
        cmap = matplotlib.cm.get_cmap(name)
        cs = [cmap(float(i+1)/N)[0:3] for i in range(N)]
    return cs

def get_iso_contours(CS):
    '''get the list of iso-contours coordinates '''
    coords = [c.get_paths()[0].vertices for c in CS.collections]
    return coords

def get_iso_contours_vertx(cs):
    '''Get all vertices from isocontours
    returns Nx2 ndarray
    '''
    v = np.zeros((0,2))
    for coll in cs.collections:
        ps = coll.get_paths()
        viE = np.zeros((0,2))
        for p in ps:
            viE = np.vstack((viE,p.vertices))
            #print(viE.shape)
        v = np.vstack((v,viE))
    return v

def get_lims(ax,mg,xylims=None,is_3d=0):
    if not xylims:
        if is_3d :
            xmin,xmax = ax.get_xlim3d()
            ymin,ymax = ax.get_ylim3d()
            zmin,zmax = ax.get_zlim3d()
            Wx,Wy,Wz = xmax-xmin,ymax-ymin,zmax-zmin
            xmin,xmax = xmin-mg*Wx, xmax+mg*Wx
            ymin,ymax = ymin-mg*Wy, ymax+mg*Wy
            zmin,zmax = zmin-mg*Wz, zmax+mg*Wz
            xylims = [xmin,xmax,ymin,ymax,zmin,zmax]
        else:
            # data = np.array(ax.dataLim.bounds,dtype=float)
            # data[data==np.inf]=1;data[data==-np.inf]=-1;data[data==np.nan]=0;
            xmin,xmax = ax.get_xlim()
            ymin,ymax = ax.get_ylim()
            W,H = xmax-xmin,ymax-ymin
            xmin,xmax = xmin-mg*W, xmax+mg*W
            ymin,ymax = ymin-mg*H, ymax+mg*H
            xylims = [xmin,xmax,ymin,ymax]
    return xylims

def changeAxesLim(ax,mg,xylims=[],is_3d=0,changeXYlims=0):
    if isinstance(xylims,np.ndarray):xylims=xylims.tolist()
    xylims = get_lims(ax,mg,xylims,is_3d)
    if len(xylims)==4:
        xmin,xmax,ymin,ymax = xylims
    elif len(xylims)==6:
        xmin,xmax,ymin,ymax,zmin,zmax = xylims
    elif len(xylims)==3:
        xylims0 = get_lims(ax,mg,None,is_3d)
        if xylims[0]=='x':
            xmin,xmax = xylims[1:3]
            ymin,ymax = xylims0[2:]
        if xylims[0]=='y':
            xmin,xmax = xylims0[:2]
            ymin,ymax = xylims[1:3]
    if changeXYlims:
        # print([xmin,xmax,ymin,ymax])
        ax.set_xlim((xmin, xmax))
        ax.set_ylim((ymin, ymax))
    xylims = [xmin,xmax,ymin,ymax]
    if is_3d :
        if changeXYlims:ax.set_zlim3d((zmin, zmax))
        xylims +=[zmin, zmax]
    return xylims

def get_tick_array(ticks,ndim,xylims):
    if isinstance(ticks,float) or isinstance(ticks,int) : ticks = [ticks]*ndim
    for i in range(ndim):
        if isinstance(ticks[i],float) or isinstance(ticks[i],int):
            ticks[i]=np.arange(xylims[2*i+0],xylims[2*i+1],ticks[i])
    return ticks
def change_ticks(ax,ticks,tick_labs,xylims,is_3d,ticks_m):
    ndim = [2,3][is_3d]
    if ticks :
        ticks = get_tick_array(ticks,ndim,xylims)
        ax.set_xticks(ticks[0])
        ax.set_yticks(ticks[1])
        if is_3d : ax.set_zticks(ticks[2])
    if ticks_m :
        ticks_m = get_tick_array(ticks_m,ndim,xylims)
        ax.set_xticks(ticks_m[0],minor=True)
        ax.set_yticks(ticks_m[1],minor=True)
    if tick_labs:
        ax.set_xticklabels(tick_labs[0])
        ax.set_yticklabels(tick_labs[1])
        if is_3d : ax.set_yticklabels(tick_labs[2])

# def add_cb(ax):
#     sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('Greys',ns));sm.set_array([]);
#     cb=fig.colorbar(sm,boundaries=0.5+np.arange(ns+1),ticks=range(1,ns+1));
#     cb.ax.set_yticklabels(['%d' %(n) for n in Ns])

def get_axPos(axPosI):
    ''' Positions predefined
    axPos = {'T':[0.2,0.12,0.75,0.75],
        'V':[0.1, 0.1, 0.75, 0.8],
        'L':[0.02, 0.1, 0.75, 0.8],
        1:[0.15, 0.11, 0.82, 0.82],
        11:[0.1, 0.1, 0.35, 0.8],
        12:[0.6, 0.1, 0.35, 0.8],
        21:[0.07, 0.1, 0.9, 0.4],
        22:[0.07, 0.6, 0.9, 0.4],
        31:[0.15, 0.18, 0.75, 0.75]}
    '''
    axPos = {'T':[0.2,0.12,0.75,0.75],
        'V':[0.1, 0.1, 0.75, 0.8],
        'E':[0.1, 0.1, 0.75, 0.75],
        'L':[0.02, 0.1, 0.75, 0.8],
        1:[0.15, 0.11, 0.82, 0.82],
        11:[0.1, 0.1, 0.35, 0.8],
        12:[0.6, 0.1, 0.35, 0.8],
        21:[0.07, 0.1, 0.9, 0.4],
        22:[0.07, 0.6, 0.9, 0.4],
        31:[0.15, 0.18, 0.75, 0.75]}
    axPosition = axPos[1]
    if isinstance(axPosI,list) or isinstance(axPosI,tuple):
        if len(axPosI)==4 : axPosition = axPosI
    elif axPosI in axPos.keys():
        axPosition = axPos[axPosI]
    return axPosition

############
def basename(file):
    if file[-1]=='/' : file = file[:-1]
    return os.path.basename(file)
def dirname(file):
    #if file[-1]=='/' : file = file[:-1]
    return os.path.dirname(file)+'/'

def warn(msg):
    print(colors.red+msg+colors.black)

def get_figpath(file,rel_path):
    figpath=os.path.realpath(os.path.dirname(os.path.realpath(file))+rel_path)+'/'
    return figpath

def saveFig(fullpath,ax=None,png=None,fmt='',figopt='1',topt=0):
    ''' save figure with transparent option
    - figopt : t(transparent) i(quality dpi=i*96)
        - example : 't','1','2', 't1','t2',...
    - topt : bool - transparent option
    '''
    if '.' in fullpath:
        filename,fmt=fullpath,fullpath.split('.')[-1]
    else:
        if not png==None : fmt=['eps','svg','png'][png]
        filename = fullpath+'.'+fmt
    if 't' in figopt : figopt,topt=figopt[1:],1
    r_dpi= 1 if not figopt else int(figopt); #print(r_dpi)
    plt.savefig(filename, format=fmt, dpi=r_dpi*96,transparent=topt)
    print(colors.green+'Saving figure :\n'+colors.yellow+filename+colors.black)

def crop_fig(name,crop):
        cropcmd = "%dx%d+%d+%d" %tuple(crop)
        cmd = "convert %s -crop %s %s" %(name,cropcmd,name)
        check_output(cmd,shell=True)
        print(colors.yellow+name+colors.blue+' cropped to '+colors.black,cropcmd)

def disp_quick(name,ax,opt,figopt):
    if 's' in opt : saveFig(name,ax,fmt='png',figopt=figopt)
    if 'p' in opt : ax.figure.show()
    if 'c' in opt : plt.close(ax.figure)

def im2gif(figpattern,fmt='svg'):
    cmd='im2gif '+figpattern +' ' + fmt
    print(colors.magenta+cmd+'  ...'+colors.black)
    out=check_output(['/bin/bash','-i','-c',cmd]).decode()
    print(colors.green+out+colors.black)

###########################################################################
#misc
def getSymCurve(x,y,symOpt=1):
    x = np.concatenate((np.flipud(-x),x))
    y = np.concatenate((np.flipud(symOpt*y),y))
    return x,y
