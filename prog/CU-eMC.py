#!/usr/bin/python -O
# -*- coding: utf-8 -*-
#"""
#Created on Tue Jan 04 16:10:16 2011
#
#@author: haywoojr
#"""
#if a variable is declared/changed within a function it is considered local
#to that function if you need it globally then you must declare it global
#in the function that originally defines it only then it will be available
#for use but not modification in other functions
#import matplotlib
import subprocess
#matplotlib.interactive(True)
#matplotlib.use( 'Qt4Agg' )
import dicom,time
import pylab,gc
#import math as mt
import os#,xlrd
import numpy as np
import sys
#from mpl_toolkits.mplot3d import Axes3D
import Tkinter, tkFileDialog
from scipy import weave,stats,ndimage
#import weave
from scipy.interpolate import Rbf,griddata
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.widgets import CheckButtons, Button, RadioButtons,Widget
from matplotlib.mlab import dist
from matplotlib.patches import Ellipse
#butn=ndimage.imread('Button.jpg')
#from matplotlib.nxutils import points_inside_poly
#from matplotlib.path import Path
#from PyQt4.QtCore import pyqtRemoveInputHook
#from PyQt4 import QtCore#,QtGui
#import pyCudaElec
#from multiprocessing import Process, Value, Array
#Directory where the program is running from needed for linking to CUDA subroutine
progdir="/home/emcuser/eMCfiles/prog/"
cudadir='/usr/local/cuda/lib64'
#Cuda include directory
cudainc='/usr/local/cuda/include'
#change directory to progdir
os.chdir(progdir)
#os.environ['LD_LIBRARY_PATH']=cudadir+":$LD_LIBRARY_PATH"
#import pymain
#print pymain.pymain.__doc__
#directory where cuda is installed needed for linking to cuda subroutine


#used for plotting
majorLocator=MultipleLocator(20)
majorFormatter=FormatStrFormatter('%d')
minorLocator=MultipleLocator(5)
    
#where function used to find to which index a value corresponds
#the function uses weave to compile the c version of the code inline
#with the python code
def whr(x,val,n):
    whrsupcode="""
    #include <cmath>
    #include <stdlib.h>
    using namespace std;
    """
    whrcode="""
    //
        double *x1;
        int i,indx=0;
        x1=new double[n1];
        double minn;
        double t,v;
        v=val1;
        for (i=0;i<n1;i++) x1[i]=x[i];
        t=x1[0];
        minn=abs(t-v);
        for (i=1;i<n1;i++){
            t=x1[i];
            if(abs(t-v)<minn){
                minn=abs(t-v);
                indx=i;
            }
        }
        return_val=indx;
    """
    val1=np.float64(val)
    n1=int(n);
    indx=weave.inline(whrcode,['x','val1','n1'],support_code=whrsupcode)
    
    return indx
#Another where function written in python using a little different syntax
#this version might be a little faster if written in C than the above function
#because it returns immediately when the index is found I have not tested it
#too much
def wher(x,val,N):
    indx=0
    dx2=abs(x[1]-x[0])/2.0
    for i in xrange(len(x)):
        if (val >= (x[i]-dx2) and val <= (x[i]+dx2)):
            return i
    return indx
##############
#Function to define radio buttons in the plotting window
#this is a copy of the function RadioButtons in the matplotlib library
#I changed it to use Ellipse instead of Circle to get bigger buttons in the
#plot
class MyRadioButtons(Widget):
    def __init__(self, ax, labels, active=0, activecolor='blue'):
        self.activecolor = activecolor

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_navigate(False)
        txs = np.asarray([0.2,0.6])
        ys = np.asarray([0.3,0.7])
        pxs=txs-0.065
        valx,valy=np.meshgrid(txs,ys)
        valpx,valy=np.meshgrid(pxs,ys)
        cnt = 0
        axcolor = ax.get_axis_bgcolor()

        self.labels = []
        self.circles = []
        for x, y, label in zip(valx.reshape(len(labels)), valy.reshape(len(labels)), labels):
            t = ax.text(x, y, label, transform=ax.transAxes,
                        horizontalalignment='left',
                        verticalalignment='center')
            self.labels.append(t)
        
        for x,y in zip(valpx.reshape(len(labels)), valy.reshape(len(labels))):
            if cnt==active:
                facecolor = activecolor
            else:
                facecolor = axcolor
            #this is the part I changes
            p = Ellipse(xy=(x, y), width=0.06,height=0.18, facecolor=facecolor,
                       )#transform=ax.transAxes)

            self.circles.append(p)
            ax.add_patch(p)
            cnt += 1

        ax.figure.canvas.mpl_connect('button_press_event', self._clicked)
        self.ax = ax


        self.cnt = 0
        self.observers = {}

    def _clicked(self, event):
        if event.button !=1 : return
        if event.inaxes != self.ax: return
        xy = self.ax.transAxes.inverted().transform_point((event.x, event.y))
        pclicked = np.array([xy[0], xy[1]])
        def inside(p):
            pcirc = np.array([p.center[0], p.center[1]])
            return dist(pclicked, pcirc) < max(p.height,p.width)

        for p,t in zip(self.circles, self.labels):
            if t.get_window_extent().contains(event.x, event.y) or inside(p):
                inp = p
                thist = t
                break
        else: return

        for p in self.circles:
            if p==inp: color = self.activecolor
            else: color = self.ax.get_axis_bgcolor()
            p.set_facecolor(color)



        if self.drawon: self.ax.figure.canvas.draw()

        if not self.eventson: return
        for cid, func in self.observers.items():
            func(thist.get_text())


    def on_clicked(self, func):
        """
        When the button is clicked, call this func with button label

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        'remove the observer with connection id cid'
        try: del self.observers[cid]
        except KeyError: pass
###############
#function to write out an html file that reports the calculated MU, the comparison with
#the TPS and an image corresponding to the image currently being displayed
#the webpage is printed from the browser as a pdf for import to the patients chart
def saveimg(event):
    global filedir
    extent=ax.get_window_extent().transformed(fig1.dpi_scale_trans.inverted())
    print extent.extents
    fig1.savefig(filedir+'/DoseCompare.png',bbox_inches=extent.expanded(1.15,1.09))
    #save to html webpage that can be viewed printed etc
    outfile=open(filedir+'/outfile.html','w')
    outfile.write('<!DOCTYPE html>\n')
    outfile.write('<html>\n')
    outfile.write('<head>\n')
    outfile.write('<title>eMC Second Check for '+ptLastName+', '+ptFirstName+'</title>\n')
    outfile.write('<style>\n')
    outfile.write('table,tr,th {background-color:#EBF4EA;text-align:left}\n')
    outfile.write('table,tr,td {background-color:#EBF4F3;text-align:left}\n')
    outfile.write('div {\n')
    outfile.write('box-shadow: 10px 10px 5px #ECECEC;')
    outfile.write('width: 95%;')
    outfile.write('height: 95%;')
    outfile.write('position:relative;padding:0 auto;margin:0 auto;')
    outfile.write('}\n')
    outfile.write('div, img {float:left;max-width:80%; max-height:80%;left:10%}')
    outfile.write('div, table {float:left}')
    outfile.write('</style>\n')
    outfile.write('</head>\n')
    outfile.write('<body>\n')
    outfile.write('<div>')
    outfile.write('<h1>Electron Calculation Second Check</h1>\n')
    outfile.write('<table>')
    outfile.write('<tr><th>Patient Name:</th><td>'+ptLastName+', '+ptFirstName+'</td></tr>\n')
    outfile.write('<tr><th>Patient ID:</th><td>'+ptID+'</td></tr>\n')
    outfile.write('<tr><th>Plan Name:</th><td>'+planName+'</td></tr>\n')
    outfile.write('</table>')
    wd=600
    ht=wd*(extent.extents[3]-extent.extents[1])/(extent.extents[2]-extent.extents[0])
    outfile.write('<img src="DoseCompare.png">')
    outfile.write('<table>')
    outfile.write('<tr><th>Beam Name:</th><td>'+beamName+'</td><th>&nbsp</th><td>&nbsp</td><th>&nbsp</th><td>&nbsp</td></tr>\n')
    outfile.write('<tr><th>Radiation Type:</th><td>'+beamRadType+'</td><th>Cone Name:</th><td>'+beamApp+'</td><th>Beam Energy</th><td>'+"%3.2f"%beamEnergy+'</td></tr>\n')
    outfile.write('<tr><th>Gantry:</th><td>'+"%3.2f"%gAng+'</td><th>Collimator:</th><td>'+"%3.2f"%col+'</td><th>PSA:</th><td>'+"%3.2f"%psa+'</td></tr>\n')
    outfile.write('<tr><th>Second check MU:</th><td>'+"%3.2f"%MU+'</td>')
    outfile.write('<th>TPS MU:</th><td>'+"%3.2f"%tpsMU+'</td>\n')
    outfile.write('<th>Percent diff:</th><td>'+"%3.2f"%(100*(MU-tpsMU)/tpsMU)+'</td></tr>\n')
    outfile.write('</div>')
    outfile.write('</body>\n')
    outfile.close()
    filedirbk=filedir
    #filedir=filedir.replace(',','\,')
    #filedir=filedir.replace('(','\(')
    #filedir=filedir.replace(')','\)')
    filedir=filedir.replace(' ','\ ')
    print 'x-www-browser '+filedir+'/outfile.html'
    fstring='x-www-browser '+filedir+'/outfile.html'
    proc=subprocess.Popen([fstring],shell=True)
    filedir=filedirbk
    return
##############
#function to visualize the ct, dose, mu, etc.
def plot_function():
    #in python only define globals in functions if they will be modified by the function
    #otherwise python created a local version and globals are not updated
    #global curi,contColors,contPoints,rotang,ax
    #global xc,yc,zc,npix3d,xd,yd,zd,rd,blckz,blckx,blcky,conez,conex,coney
    #global isoLines,totalDose,calx,caly,calz,cDose3d
    #global isopos,ftin,refpos,contOnOff,contNames,MU,tpsMU
    global ftin,ax,fig1,axpos,npx0,npx1,npy0,npy1
    #if we have been in the function already get the axis limits
    try:
        if ftin==1:
            limits=ax.axis()
        ax.format_coord=lambda x,y : 'x=%g y=%g' % (x,y)
    except:
        ax.format_coord=lambda x,y : 'Please select the CT, RP, RD, RS data datasets first x=%g y=%g' % (x,y)
        return 0
    #if we have not been in the function get the original limits
    if ftin==0:
        (npx0,npy0),(npx1,npy1)=ax.get_position(original=True).get_points()
    fig1.delaxes(ax)
    #add the axis to the plot
    ax=fig1.add_subplot(111)
    #position of the plot
    fig1.subplots_adjust(bottom=0.16,top=0.95,left=0.3,right=0.97)     
    axpos=ax.get_position()
    #define which isodose lines to display, 10%, 30%, 50%, 75%, 90%, 95%, 100%, 105%, 110%
    isoLines=np.asarray([0.1,0.3,0.5,0.75,0.9,0.95,1.0,1.05,1.1])
    #colors for the line 10% to the left, 110% to the right
    clrs=('purple','cyan','blue','green','yellow','orange','red','magenta','white')
    #local arrays and variables for holding the data to plot
    onZ=[];dz=zd[1]-zd[0];bonZx=[];bonZy=[];conZx=[];conZy=[]
    calZx=[];calZy=[];dx=calx[1]-calx[0];
    #scale the isodose lines to cGy from % for the TPS dose
    #only used for display
    V=isoLines*totalDose;
    #print "IsoLevels= ",V
    #scale the isodose line to cCy from % for the calculated dose
    #only used for display
    mV=isoLines*totalDose;
    #title for the plot    
    tit=""
    #zdi is the Z index of the plot for the tps dose, cdi is the current index for the calculated dose
    #bdz is the distance between 
    #planes for the block, zdz is the distance between each Z index
    zdi=0;cdi=0;bdz=blckz[1]-blckz[0];zdz=abs(zc[1]-zc[0])
    try:
        #find where the Z planes in the tps dose grid are closest to the current calculated plane
        zdi=np.where(abs(zd-zc[curi])<=abs(dz/2.0))[0][0];
        #dose here is the TPS dose
        dose=np.asarray(rd.pixel_array[zdi]*rd.DoseGridScaling)#,order='fortran')
    except:
        try:
            #same as above but look over a larger region
            zdi=np.where(abs(zd-zc[curi])<=abs(dz/1.25))[0][0];
            #average between the two closest planes
            dose=(np.asarray(rd.pixel_array[zdi])+ np.asarray(rd.pixel_array[zdi+1]))*rd.DoseGridScaling/2.0
        except:
            tit=tit+' no dose on slice'		
    try:
        #find where the dose plane calculated in this simultion corresponds to the current
        cdi=np.where(abs(cdz-zc[curi])<=abs(zdz/2.0))[0][0];
        #print cdi, calx[cdi], zc[curi]
        #scale the dose on this plane to cGy
        mdose=cDose3d[:,:,cdi]*totalDose;
    except:
        try:
            #same as above but with a larger search area
            cdi=np.where(abs(cdz-zc[curi])<=abs(zdz/1.25));
            #print cdi, calx[cdi], z[curi]
            #average between planes
            mdose=(cDose3d[:,:,cdi]+cDose3d[:,:,cdi-1])*totalDose/2.0;
        except:
            tit=tit+" no Calc'd dose on slice"""
        #tit=tit+" no Calc'd dose on slice"
    #print len(blckz), zc[curi], curi
    #for each of blckz and conez find which points are on the current plane
    #append the points to the array bonZx, block x points on this Z, bonZy block y points on this Z
    #append the points to the arrat conZx, cone z points on this Z, cone y points on this Z
    for i in xrange(len(blckz)):
        if abs(blckz[i]-zc[curi])<abs(zdz/2.0):
            bonZx.append(blckx[i])
            bonZy.append(blcky[i])
    for i in xrange(len(conez)):
        if abs(conez[i]-zc[curi])<abs(zdz/2.0):
            conZx.append(conex[i])
            conZy.append(coney[i])
            
    #
    
    """calZz=[]
    for i in xrange(len(calz)):
        if abs(calz[i]-zc[curi])<abs(dz/2.0):
            calZx.append(calx[i])
                #calZx.append(calx.min())
            calZy.append(caly[i])
            #calZz.append(cDose3d[i])
                #calZy.append(calz.min())
        tx=np.arange(np.asarray(calZx).min(),np.asarray(calZx).max()+dy,dy)
        ty=np.arange(np.asarray(calZy).min(),np.asarray(calZy).max()+dy,dy)
        XI,YI=np.meshgrid(tx,ty)
        rbfi=Rbf(calZx,calZy,calZz)
        mdose=rbfi(XI,YI)
    except:
        tit=tit+'no rbf on image'"""
    #rdose=rotMatrix.rotmat(dose,np.asarray(xd,order='fortran'), \
    #np.asarray(yd,order='fortran'),rotang)
    #print dose.shape
    #find all the contours that are on this Z plane
    for i in xrange(len(contPoints)):
        if contPoints[i][2][2]==zc[curi]:
            for k in xrange(len(contNames)):
                if contNames[k]==contPoints[i][3] and contOnOff[k]:
                    onZ.append([contPoints[i][0],contPoints[i][1],contPoints[i][2]])        
    #create contours from onZ
    onZx=[];onZy=[];colors=[]
    for i in xrange(len(onZ)):
        onZx.append(np.array(onZ[i][2]).reshape(onZ[i][1],3)[:,0])
        onZy.append(np.array(onZ[i][2]).reshape(onZ[i][1],3)[:,1])
        #ronZx.append(mt.cos(rotang)*onZx[i]-mt.sin(rotang)*onZy[i])
        #ronZy.append(mt.sin(rotang)*onZx[i]+mt.cos(rotang)*onZy[i])
        #define colors to be the same as in the TPS
        for j in xrange(len(contColors)):
            if onZ[i][0]==contColors[j][0]:
                colors.append(np.array(contColors[j][1])/255.0)
    #print 'event ',event.button,' curi ',curi
    #set the limits of the plot
    xl=xc.min()
    xr=xc.max()
    yb=yc.max()
    yt=yc.min()
    ax.cla()
    tit=''
    #show the plot with the density stored in npix3d and color map of bone
    ax.imshow(npix3d[:,:,curi],cmap=pylab.cm.bone,extent=[xl,xr,yb,yt])
    #print "xl= ",xl, " xr= ",xr, " yb= ",yb," yt= ",yt
    #pylab.contour(xc,yc,npix3d[:,:,curi],cmap=pylab.cm.bone)
    #plot the contours is requested
    for i in xrange(len(onZx)):
        ax.plot(onZx[i],onZy[i],color=colors[i])
    #plot the TPS dose as contours
    try:
        #rxd=mt.cos(rotang)*xd-mt.sin(rotang)*yd
        #ryd=mt.sin(rotang)*xd+mt.cos(rotang)*yd
        ax.contour(xd,yd,dose,levels=V,linestyles='dotted',colors=clrs)
    except:
        tit=' no dose contours on slice'
    #plot the iso center point
    if abs(zc[curi]-isopos[2])<=zdz/2.0:
        ax.plot(isopos[0],isopos[1],'-o',ms=10,lw=2,alpha=0.7,mfc='orange')
    #plot the reference point
    if abs(zc[curi]-refpos[2])<=zdz/2.0:
        ax.plot(refpos[0],refpos[1],'cd',ms=10)
    #add the dose contours calcualted in this simultaion on this slice
    try:
        ax.contour(cdx,cdy,mdose,levels=mV,origin='image',colors=clrs)
    except:
        tit=" no calc'd contours on slice"
    #plot the block points    
    ax.plot(bonZx,bonZy,color='brown')
    #plot the cone points
    if len(conZx)!=0:
        ax.scatter(conZx,conZy,color='yellow',marker="x")
    #if len(calZx)!=0:
    #    ax.scatter(calZx,calZy,color='blue',marker="+")
    #adjust the plot limits "better"
    xl=min(xc.min(),calx.min())-5
    xr=max(xc.max(),calx.max())+5
    yb=max(yc.max(),caly.max())+5
    yt=min(yc.min(),caly.min())-5
    #print "xl= ",xl, " xr= ",xr, " yb= ",yb," yt= ",yt
    #set the axis title
    ax.set_title("Image Plane: "+str(zc[curi])+tit)
    #store the axis limits for replotting
    if ftin==0:
        ax.axis([xl,xr,yb,yt])
        limits=ax.axis()
        limitso=ax.axis()
        ftin=1
    else:
        ax.axis(limits)
    #add the MU comparison
    ax.text(xl*(0.95),yt*0.9,'Cal MU= %8.2f'%MU,fontsize=20,color="white",bbox=dict(facecolor='blue',alpha=0.9)) 
    ax.text(xl*(0.95),yt*0.8,'TPS MU= %8.2f'%tpsMU,fontsize=20,color="white",bbox=dict(facecolor='blue',alpha=0.9))  
    #draw the plot on the screen
    pylab.draw()
#function to plot PDD if right click in the window or profile if middle click in the window
#not working fully for profiles, need to store them and tell the proper path
def myonclick(event):
    global ax,cind,rind,indz,cDose3d,fig1,dose3
    try:
        if len(cax)==0:
            return
    except:
        return
    #print event.inaxes==ax
    if event.inaxes==ax:
        xdd6=np.asarray([0.0,0.8,1.0,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4]);
        pdd6=np.asarray([0.757,0.927,0.969,0.998,0.965,0.892,0.779,0.632,0.469,0.308,0.175,0.081,0.031,0.012]);
        xdd9=np.asarray([-1.0,-0.5,-0.09,0.0,0.8,1.0,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2]);
        pdd9=np.asarray([0.0,0.0,0.0,80.2,88.7,90.8,95.1,97.2,98.9,99.9,99.7,98,94.5,88.8,81.3,71.9,60.3,48.8,36.4,24.9,15.8,8.6,4.4,2.2,1.4,1.2])/100.;
        xdd12=np.asarray([0.0, 0.8,1.0,1.4,1.6,1.8, 2.0,2.2,2.4,2.6,2.8, 3.0,3.2,3.4, 3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0, 5.2,5.4,5.8,6.2,6.6, 7.0]);
        pdd12=np.asarray([0.86,0.917,.927,.946,.954,.963,.972,.982,.99,.997,1,.999,.993,.979,.96,.928,.885,.832,.765,.691,.608,.517,.424,.332,.178,.076,.032,.021]);
        #print event.button
        #print axpos
        dz=zd[1]-zd[0];dx=xd[1]-xd[0];
        dose3=dose3.reshape([cay.size,cax.size,caz.size])
        zdi=np.where(abs(zd-refpos[2])<=abs(dz/2.0))[0][0];
        edose=np.asarray(rd.pixel_array[zdi]*rd.DoseGridScaling)#,order='fortran')
        xdi=np.where(abs(xd-refpos[0])<=abs(dx/2.0))[0][0];
        ydi=np.where(abs(yd-refpos[1])<=abs(dx/2.0))[0][0];
        drefdose=edose[ydi,xdi]
        edose=edose/drefdose;
        print dose3.shape
        if event.button==3:
            fig1.delaxes(ax)
            ax=fig1.add_subplot(111)
            fig1.subplots_adjust(bottom=0.16,top=0.95,left=0.3,right=0.97)  
            ax.set_title('PDD')
            ax.set_autoscale_on(True)
            cind=whr(cax+isopos[0],refpos[0],cax.size)
            indz=whr(caz+isopos[2],refpos[2],caz.size)
            sxdi=whr(cdx,refpos[0],cdx.size)
            szdi=whr(cdz,refpos[2],cdz.size)
            
            print "cind= ",cind,"indz= ",indz
            xl=cay.min();xr=cay.max();yb=dose3[:,cind,indz].min();
            yt=dose3[:,cind,indz].max();
            #ax.xaxis.set_ticks(cay)
            ax.grid(True)
            majorLocator=MultipleLocator(25)
            majorFormatter=FormatStrFormatter('%d')
            minorLocator=MultipleLocator(5)
            ax.xaxis.set_major_locator(majorLocator)
            ax.xaxis.set_major_formatter(majorFormatter)
            ax.xaxis.set_minor_locator(minorLocator)
            ax.xaxis.grid(True,which='minor')
            majorLocator=MultipleLocator(0.25)
            minorLocator=MultipleLocator(0.05)
            ax.yaxis.set_major_locator(majorLocator)
            majorFormatter=FormatStrFormatter('%1.2f')
            ax.yaxis.set_major_formatter(majorFormatter)
            ax.yaxis.set_minor_locator(minorLocator)
            ax.yaxis.grid(True,which='minor')
            ax.plot(cay,dose3[:,cind,indz])
            ax.plot(yd-isopos[1],edose[:,xdi])
            ax.plot(cdy-isopos[1],cDose3d[:,sxdi,szdi])
            if abs(nomE-9)<1e-6:
                ax.plot(xdd9*10,pdd9)
            elif abs(nomE-12)<1e-6:
                ax.plot(xdd12*10,pdd12)
            elif abs(nomE-6)<1e-6:
                ax.plot(xdd6*10,pdd6)
            # xl,xr,yb,yt
           #ax.set_ylim(yb,yt)
            #ax.set_position([npx0, npy0, npx1-npx0, npy1-npy0])
            #print axpos
            ax.relim()
            ax.autoscale_view(True,True,True)
        elif event.button==2:
            profdat=np.genfromtxt('/home/Public/data/09E-10x10-20-OCR.csv',delimiter=',')
            profdat3=np.genfromtxt('/home/Public/data/09E-10x10-30-OCR.csv',delimiter=',')
            fig1.delaxes(ax)
            ax=fig1.add_subplot(111)
            fig1.subplots_adjust(bottom=0.16,top=0.95,left=0.3,right=0.97)
            ax.set_title('dMax Profile')
            ax.set_autoscale_on(True)
            rind=whr(cay,20.0,cay.size)
            cind=whr(cax,0.0,cax.size)
            indz=whr(caz,0.0,caz.size)
            sydi=whr(cay,30.0,cay.size)
            szdi=whr(caz,0.0,caz.size)
            print "rind= ",rind,"indz= ",indz
            xl=cax.min();xr=cax.max();yb=dose3[sydi,:,szdi].min();
            yt=dose3[sydi,:,szdi].max();
            ax.plot(cax,dose3[rind,:,indz]/dose3[rind,cind-5:cind+5,indz].mean())
            ax.plot(profdat[:,0],profdat[:,1]/100.0)
            ax.plot(cax,dose3[sydi,:,szdi]/dose3[sydi,cind-5:cind+5,szdi].mean())
            ax.plot(profdat3[:,0],profdat3[:,1]/100.0)
            #print xl,xr,yb,yt
           #ax.set_ylim(yb,yt)
            #ax.set_position([npx0, npy0, npx1-npx0, npy1-npy0])
            #print axpos
            ax.relim()
            ax.autoscale_view(True,True,True)
        #need machine profile here
        else:
            return 0
    else:
        return 0
    event.canvas.draw()
#function to scroll to the next plane of the dose grid    
def scroll_function(event):
    global curi,ax
    try:
        if ftin==1:
            limits=ax.axis()
        ax.format_coord=lambda x,y : 'x=%g y=%g' % (x,y)
    except:
        ax.format_coord=lambda x,y : 'Please select the CT, RP, RD, RS data datasets first x=%g y=%g' % (x,y)
        return 0

    if event.button=='up':
        curi-=1
        if curi <= 0 :
            curi=0
    else:
        curi+=1
        if curi > npix3d.shape[2]-1:
            curi=npix3d.shape[2]-1
    plot_function()
    event.canvas.draw()
    
#function to reset the plot limits to the original
def resetLimits(event):
    global ftin
    ftin=0
    try:
        plot_function()
    except:
        return 0
#get the directory and read all the files needed to to the calculations        
def getFile(event):
    global filedir,cdir,CTfiles,RSfiles,RPfiles,RDfiles,fig1
    global rd,xd,yd,zd,ftin,ds,npix3d,rp,tpsMU,tpsDose,numBeams,totalDose
    global isopos,refpos,rotang,blckx,blcky,blckz,conez,conex,coney,colAng,PSAang
    global mincone,maxcone,contColors,contPoints,bolpnts,rs
    global axContours,check,axBolus,bolus,contNames,bolNames,contOnOff,bolOnOff
    global xc,yc,zc,dy,calx,caly,calz,curi,xx,zz,MU,bolPrevConv,cDose3d,cdy,cdx,cdz
    global ptLastName,ptFirstName,ptID,planName,planEnergy,gAng,col,psa,beamName,beamApp
    global beamRadType,beamEnergy,ddy,cone,npix3do
    #do some clean up. not quite right still buggy
    try:
        cDose3d=None
    except:
        print "no cDose3d to delete"
    try:
        del bolOnOff
        del bolNames
        del rs
        del rd
        del xd,yd,zd,ds,npix3d,rp
    except:
        print "no Overrides to delete"
    gc.collect()
    z=[];x=[];y=[]
    MU=0
    ftin=0
    filedir=""
    cdir=os.getcwd()
    root=Tkinter.Tk()
    root.withdraw()
    #get patient files directory all files should be in the same directory
    filedir=tkFileDialog.askdirectory(parent=root,initialdir=cdir,\
                                          title="Please select CT image directory")
    root.destroy()
        #put code here to catch if there are no CTfiles in the folder
    try:
        os.chdir(filedir)
    except:
        return
    #get all files in directory
    files=os.listdir(".")
    CTfiles=[];RSfiles=[];RPfiles=[];RDfiles=[]
    #divide into each type of file
    #needs error checking like are these CT files the same as referenced in RT Plan etc
    for item in files:
        if item.split(".")[0]=="CT":
            CTfiles.append(item)
        elif item.split(".")[0]=="RS":
            RSfiles.append(item)
        elif item.split(".")[0]=="RP":
            RPfiles.append(item)
        elif item.split(".")[0]=="RD":
            RDfiles.append(item)
    #exit if did not find all the needed files
    if not len(CTfiles) or not len(RSfiles) or not len(RPfiles) or not len(RDfiles):
        print "All the necessary files were not found.  Please try again!"
        print "CT, RD, RP, and RS files are all needed."
        return
    #sort the files should read in each files Z plane and sort them with that but this works for now
    try:
        CTfiles.sort(key=lambda line: long(line.split(".")[12]))
    except:
        CTfiles.sort(key=lambda line: long(line.split(".")[-3]))
    #read in RT Plan
    rp=[];
    rp=dicom.read_file(RPfiles[0])
    PtSetup=rp.PatientSetups[0].PatientPosition
    print PtSetup
    #read in the dose
    rd=[];
    #rd is the TPS dose
    rd=dicom.read_file(RDfiles[0])
    #dose=np.asarray(rd.pixel_array,'float32',order='Fortran')*rd.DoseGridScaling
    #Xd, yd, zd are the x, y, z, coordinates of the TPS dose volume
    xd=np.arange(rd.Columns)*rd.PixelSpacing[1]+rd.ImagePositionPatient[0]
    yd=np.arange(rd.Rows)*rd.PixelSpacing[0]+rd.ImagePositionPatient[1]
    zd=np.asarray(rd.GridFrameOffsetVector)+rd.ImagePositionPatient[2]
    if PtSetup=='HFP':
        #for head first prone
        xd=np.arange(rd.Columns)*rd.PixelSpacing[1]-rd.ImagePositionPatient[0]
        yd=np.arange(rd.Rows)*rd.PixelSpacing[0]-rd.ImagePositionPatient[1]
        zd=np.asarray(rd.GridFrameOffsetVector)+rd.ImagePositionPatient[2]
    #read the RT Plan and define some patient relate values
    ds=dicom.read_file(RPfiles[0])
    try:
        ptLastName,ptFirstName=ds.PatientsName.split("^")
    except:
        strng=[]
        strng=ds.PatientsName.split("^")
        ptLastName=strng[0]
        ptFirstName=strng[1]
    ptID=ds.PatientID
    planName=ds.RTPlanLabel
    beamName=ds.Beams[0].BeamName
    beamRadType=ds.Beams[0].RadiationType
    beamApp=ds.Beams[0].Applicators[0].ApplicatorID
    beamEnergy=ds.Beams[0].ControlPoints[0].NominalBeamEnergy
    print ptLastName,ptFirstName,ptID
    print planName,beamName
    print beamRadType,beamApp,beamEnergy
    #read in the plan
    #take the first plan
    #Bolus is not coming over for some reason Varian verified it doesn't
    rp=[];
    rp=dicom.read_file(RPfiles[0])
    tpsMU=rp.FractionGroups[0].ReferencedBeams[0].BeamMeterset
    tpsDose=rp.FractionGroups[0].ReferencedBeams[0].BeamDose*100 #100 converts to cGy

    #needs error handling
    numBeams=rp.FractionGroups[0].NumberofBeams
    if numBeams==1:
        totalDose=rp.FractionGroups[0].NumberofFractionsPlanned* \
            rp.FractionGroups[0].ReferencedBeams[0].BeamDose
    # the iso is always with respect to the image 0 so eclipse will
    #display a different iso if you have moved the user origin
    #iso center position and reference point position
    isopos=rp.Beams[0].ControlPoints[0].IsocenterPosition
    refpos=rp.FractionGroups[0].ReferencedBeams[0].BeamDoseSpecificationPoint
    #get gantry angle, collimator andle, and Patient Support Assembly angle (or couch)
    gAng=rp.Beams[0].ControlPoints[0].GantryAngle
    colAng=rp.Beams[0].ControlPoints[0].BeamLimitingDeviceAngle
    col=colAng
    PSAang=-(rp.Beams[0].ControlPoints[0].PatientSupportAngle)+360.0
    psa=PSAang
    print gAng,colAng,PSAang
    #the coordinates are off since +y is down so have to
    #correct for that
    #if gAng <180.0 and gAng >0.0:
    rotang=(gAng)*np.pi/180.0
    #elif gAng>180.0 and gAng < 360.0:
    #    rotang=(360.0-gAng)*np.pi/180.0
    #else:
    #    rotang=0.0
    print rotang
    print isopos
    print refpos
    if PtSetup=='HFP':
         #for head first prone 
        isopos[0]=-isopos[0];isopos[1]=-isopos[1];
        refpos[0]=-refpos[0];refpos[1]=-refpos[1];
    #read in the first CT file to define arrays
    ds=dicom.read_file(CTfiles[0])
    npix3d=np.zeros((ds.Columns,ds.Rows,len(CTfiles)),dtype=np.intc)
    for i in xrange(len(CTfiles)):
        #read in all files to the array
        ds=dicom.read_file(CTfiles[i])
        npix3d[:,:,i]=ds.pixel_array*int(ds.RescaleSlope)+int(ds.RescaleIntercept)
        #Need switch for HFDR HFDL etc for now manual
        #npix3d[:,:,i]=np.rot90(npix3d[:,:,i],3)
    #read in the CT scan
    ds=[]
    #get the x, y, z corrdinates for the CT images
    for i in xrange(len(CTfiles)):
        ds.append(dicom.read_file(CTfiles[i]))
        #see if CT or VARIAN phantom file
        x0=ds[i].ImagePositionPatient[0]
        y0=ds[i].ImagePositionPatient[1]
        try:
            z.append(ds[i].SliceLocation)
        except:
            z.append(ds[i].ImagePositionPatient[2])
        #x.append((np.arange(ds[i].Columns)+0.5-(ds[i].Columns)/2.0)*ds[i].PixelSpacing[0])
        #pixel spacing is the row spacing in [0] and column spacing in [1]
        #y=rows, x=columns
        x.append(np.arange(x0,x0+ds[i].Columns*ds[i].PixelSpacing[1],ds[i].PixelSpacing[1]))
        #y.append((np.arange(ds[i].Rows)+0.5-(ds[i].Rows)/2.0)*ds[i].PixelSpacing[1])
        y.append(np.arange(y0,y0+ds[i].Rows*ds[i].PixelSpacing[0],ds[i].PixelSpacing[0]))
        #y.append((np.arange(ds[i].Rows)-(ds[i].Rows)/2.0)*ds[i].PixelSpacing[0])
        #rotx.append(mt.cos(rotang)*x[i]-mt.sin(rotang)*y[i])
        #roty.append(mt.sin(rotang)*x[i]+mt.cos(rotang)*y[i])
    npix3do=np.zeros_like(npix3d)
    npix3do[:,:,:]=npix3d[:,:,:]
    print len(z), len(set(z))
    print 'x[0]= ',x[0]
    print 'y[0]= ',y[0]
    #this does a repair on the Z coordinates needed when the sort of the CT files above
    #does not give the correct order
    if len(z) != len(set(z)):
        dzbeg=z[1]-z[0]
        dzend=z[-1]-z[-2]
        if dzend != dzbeg:
            for i in xrange(len(z)):
                if i>0 :
                    z[-1-i]=z[-1-i+1]-dzend
    print z
    #check if z has repeating 
    #create numpy arrays for the calculations these have beeter functionality than python lists
    xc=np.asarray(x[0])
    yc=np.asarray(y[0])
    ddy=abs(yc[0]-yc[1])
    print "ddy= ",ddy, " ",ds[0].PixelSpacing[0]
    zc=np.asarray(z)
    del x;del y;
   
    print isopos
    print refpos
    #colAng=50.0
    #rotang=20*np.pi/180.0#0000EE
    #PSAang=0.0

    #get block data
    if rp.Beams[0].NumberofBlocks:
        blckp=rp.Beams[0].Blocks[0].BlockNumberofPoints
        bckd=np.asarray(rp.Beams[0].Blocks[0].BlockData)
        xx=bckd.reshape(blckp,2,)[:,0]
        zz=bckd.reshape(blckp,2,)[:,1]
    else:
        #no block assume shape is the same as cone
        if beamApp=="A10":
            zz=np.asarray([50,50,-50,-50,50])
            xx=np.asarray([-50,50,50,-50,-50])
        elif beamApp=="A15":
            zz=np.asarray([75,75,-75,-75,75])
            xx=np.asarray([-75,75,75,-75,-75])        
        elif beamApp=="A20":
            zz=np.asarray([100,100,-100,-100,100])
            xx=np.asarray([-100,100,100,-100,-100])
        elif beamApp=="A06":
            zz=np.asarray([30,30,-30,-30,30])
            xx=np.asarray([-30,30,30,-30,-30])
        elif beamApp=="A25":
            zz=np.asarray([125,125,-125,-125,125])
            xx=np.asarray([-125,125,125,-125,-125])


    #zz=zz*950.0/1000.0
    #xx=xx*950.0/1000.0
    yy=np.zeros(len(xx))-50.0
    #rotate around collimator
    xxtmp=np.cos(colAng*np.pi/180)*(xx)-np.sin(colAng*np.pi/180)*(zz)
    zztmp=np.sin(colAng*np.pi/180)*(xx)+np.cos(colAng*np.pi/180)*(zz)
    #rotate around the patient using gantry
    xx=np.cos(rotang)*xxtmp-np.sin(rotang)*yy
    yytmp=np.sin(rotang)*xxtmp+np.cos(rotang)*yy
    #rotate couch
    xxtmp=np.cos(PSAang*np.pi/180.)*xx-np.sin(PSAang*np.pi/180.0)*(zztmp)
    zz=np.sin(PSAang*np.pi/180.0)*xx+np.cos(PSAang*np.pi/180.0)*(zztmp)
    xx=xxtmp.copy()
    yy=yytmp.copy()
    # if no couch rotation just rotate x's
    #according to Dicom Stnd Block data is at iso center
    #so project back to cone at 95cm
    #rotate as if it were at 0
    #convert to physical coordinates
    #zz=(zz*(950/1000))
    #iindx,iindy=np.mgrid[0:len(xx),0:len(yy)]
    #xx=xx*(950/1000.0)
    
    del xxtmp,zztmp,yytmp
    blckz=zz+isopos[2]
    
    #yy=np.linspace(-55,-50,5)
    #######rotate about isopos by gang######## translate to iso
    #iindx,iindy=np.mgrid[0:len(xx),0:len(yy)]
    blckx=(xx)+isopos[0]
    blcky=(yy)+isopos[1]
    #redefine for cudaget block data
    if rp.Beams[0].NumberofBlocks:
        blckp=rp.Beams[0].Blocks[0].BlockNumberofPoints
        bckd=np.asarray(rp.Beams[0].Blocks[0].BlockData)
        xx=bckd.reshape(blckp,2,)[:,0]
        zz=bckd.reshape(blckp,2,)[:,1]
    """sszz=30
    zz=np.asarray([sszz,sszz,-sszz,-sszz,sszz])
    xx=np.asarray([-sszz,sszz,sszz,-sszz,-sszz])"""
    #move block to end of cone for calc
    zz=zz*950.0/1000.0
    xx=xx*950.0/1000.0

    #del iindx;del iindy;
    #shift back to isopos
    #blckx=blckx+isopos[0];blcky=blcky+isopos[1];
    #blckx=blckx+50.0*mt.cos(np.pi/2.0-abs(rotang))+isopos[0]
    #have to make all more negative since eclipse calls up negative
    #blcky=blcky-50.0*mt.sin(np.pi/2.0-abs(rotang))+isopos[1]
    #blckx=blckx+isopos[0]
    #blcky=blcky+isopos[1]
    #draw define cone
    cone=rp.Beams[0].Applicators[0].ApplicatorID
    cx=[];cz=[];cy=[]
    if cone=='A15':
        cz=np.linspace(-75,75,100)*(950.0/1000.0)
        cx=np.linspace(-75,75,100)*(950.0/1000.0)
    elif cone=='A10':
        cz=np.linspace(-50,50,100)*(950.0/1000.0)
        cx=np.linspace(-50,50,100)*(950.0/1000.0)
    elif cone=='A06':
        cz=np.linspace(-30,30,100)*(950.0/1000.0)
        cx=np.linspace(-30,30,100)*(950.0/1000.0)
    elif cone=='A20':
        cz=np.linspace(-100,100,100)*(950.0/1000.0)
        cx=np.linspace(-100,100,100)*(950.0/1000.0)
    elif cone=='A25':
        cz=np.linspace(-125,125,100)*(950.0/1000.0)
        cx=np.linspace(-125,125,100)*(950.0/1000.0)
    else:
        print "no cone matches"
        return
    iindx,iindz=np.mgrid[0:len(cx),0:len(cz)]
    ccx=cx[iindx];ccz=cz[iindz]
    #cz=np.asarray([-75,75,75,-75,-75])
    #cx=np.asarray([-75,-75,75,75,-75])
    cy=np.zeros(ccz.shape)-50.0
    mincone=cz.min()
    maxcone=cz.max()
    #conez=conez+isopos[2]
    #conex=mt.cos(rotang)*cx
    #rotate collimator
    cxtmp=np.cos(colAng*np.pi/180.0)*(ccx)-np.sin(colAng*np.pi/180.0)*(ccz)
    cztmp=np.sin(colAng*np.pi/180.0)*(ccx)+np.cos(colAng*np.pi/180.0)*(ccz)
    #rotate around pat
    ccx=np.cos(rotang)*(cxtmp)-np.sin(rotang)*(cy)
    cytmp=np.sin(rotang)*(cxtmp)+np.cos(rotang)*(cy)
    #rotate patient
    cxtmp=np.cos(PSAang*np.pi/180.0)*(ccx)-np.sin(PSAang*np.pi/180.0)*(cztmp)
    ccz=np.sin(PSAang*np.pi/180.0)*(ccx)+np.cos(PSAang*np.pi/180)*(cztmp)
    #move to isopos
    conez=ccz+isopos[2]
    conex=cxtmp+isopos[0]
    coney=cytmp+isopos[1]
    
    conez=conez.reshape(conez.size)
    conex=conex.reshape(conex.size)
    coney=coney.reshape(coney.size)
    """fig2=pylab.figure()
    ax3=Axes3D(fig2)
    ax3.plot(conex,coney,conez)
    ax3.plot(blckx,blcky,blckz)
    fig2.show()"""
    #read in the structure set    
    #only take first structure set
    #need to put try except statements here set and error
    #use the error to exit safely :-);
    contColors=[];contPoints=[];bolpnts=[];
    rs=[];
    rs=dicom.read_file(RSfiles[0])
    contNames=[];bolNames=[];contOnOff=[];bolOnOff=[];bolPrevConv=[];test=[]
    
    for i in xrange(len(rs.ROIContours)):
        #if (rs.StructureSetROIs[i].ROIName=="6EBolus" or rs.StructureSetROIs[i].ROIName.lower()=="Body".lower()):
        contNames.append(rs.StructureSetROIs[i].ROIName)
        contOnOff.append(False)
        contColors.append([rs.ROIContours[i].ReferencedROINumber, \
                               rs.ROIContours[i].ROIDisplayColor])
        try:
            for j in xrange(len(rs.ROIContours[i].Contours)):
                contPoints.append([rs.ROIContours[i].ReferencedROINumber, \
                                       rs.ROIContours[i].Contours[j].NumberofContourPoints, \
                                       rs.ROIContours[i].Contours[j].ContourData,rs.StructureSetROIs[i].ROIName])
            #for ind in range(len(contPoints[i][2]))[1::3]:
            #    contPoints[i][2][ind]=-contPoints[i][2][ind]
        
        #onZx.append(np.array(onZ[i][2]).reshape(onZ[i][1],3)[:,0])        
        except:
            print "Exception taken"
    if PtSetup=='HFP':
        for i in xrange(len(contPoints)):
            for ind in range(len(contPoints[i][2]))[::3]:
                contPoints[i][2][ind]=-contPoints[i][2][ind]
                contPoints[i][2][ind+1]=-contPoints[i][2][ind+1]
    print len(contPoints),len(contNames)
    #need a way to check density overrides and to select bolus
    #Only allow 1 bolus
    #automatically add names if physical density is overridden
    for i in xrange(len(rs.ROIContours)):
        try:
                
            val=rs.RTROIObservations[i].ROIPhysicalProperties[0].ROIPhysicalPropertyValue
            val=1000*(val-1)/1
            bolNames.append(rs.StructureSetROIs[i].ROIName)
            bolPrevConv.append(False)
            bolOnOff.append(False)
            print rs.StructureSetROIs[i].ROIName+" will be given an HU of ",  val, " if desired"
            #if(rs.StructureSetROIs[i].ROIName.lower()=="Body".lower()):
            try:
                for j in xrange(len(rs.ROIContours[i].Contours)):
                    bolpnts.append([rs.ROIContours[i].ReferencedROINumber, \
                                        rs.ROIContours[i].Contours[j].NumberofContourPoints, \
                                        rs.ROIContours[i].Contours[j].ContourData,val,rs.StructureSetROIs[i].ROIName])
            except:
                print "Exception taken"
        except:
            val=-9999
    """for i in xrange(len(bolpnts)):
        for ind in range(len(bolpnts[i][2]))[::3]:
            bolpnts[i][2][ind]=-bolpnts[i][2][ind]
            bolpnts[i][2][ind+1]=-bolpnts[i][2][ind+1]"""
    fig1.delaxes(axContours)
    fig1.delaxes(axBolus)
    del check,bolus
    #define the subplots of the window that list the contour names as well as contour
    #names that have override options, all overridable contours are put in the bolus subplot
    axContours=pylab.axes([0.03,0.47,0.2,0.39])
    axContours.set_title('Contours to plot',fontsize=11)
    axContours.get_xaxis().set_visible(False)
    axContours.get_yaxis().set_visible(False)
    #check=CheckButtons(axContours,('NoBody','NoGTV'),(False,False))
    #check.on_clicked(callback.func)
    axBolus=pylab.axes([0.03,0.14,0.2,0.29])
    axBolus.set_title('Contours with override in TPS.',fontsize=11)
    axBolus.get_xaxis().set_visible(False)
    axBolus.get_yaxis().set_visible(False)


    axContours.set_title('Contours to plot',fontsize=11)
    try:
            #del check
        check=CheckButtons(axContours,contNames,contOnOff)
    except:
        check=CheckButtons(axContours,contNames,contOnOff)
    check.on_clicked(func)
    
    axBolus.set_title('Contours to override',fontsize=11)
    try:
            #del bolus
        bolus=CheckButtons(axBolus,bolNames,bolOnOff)
    except:
        bolus=CheckButtons(axBolus,bolNames,bolOnOff)
    bolus.on_clicked(bfunc)
    #for i in bolNames:
    #    print i
        #default dy=2.5
    #add code to auto calculate this i.e. do some trig here and 5cm past Rp
    #us RP~Enom/2
    #the following code sets up the grid based on the energy of the beam
    #since a 6E beam doesn't penetrate as far as a 12E beam make a smaller
    #calculation volume
    #distance between isocenter and reference point
    dist=np.sqrt((isopos[0]-refpos[0])**2+(isopos[1]-refpos[1])**2+(isopos[2]-refpos[2])**2)
    rangep=beamEnergy/.20
    dmax=beamEnergy/.40
    cay=np.arange(-50,round(dist+(rangep-dmax)+50),dy)
    #cay=np.linspace(-50,150,NROW,'float32')
    #dy=abs(cay[1]-cay[0])
    #lateral extent 3 cm beyond cone size
    caz=np.arange(mincone-30,maxcone+30,dy)
    cax=np.arange(mincone-30,maxcone+30,dy)
    NZ=caz.size
    NCOL=cax.size
    NROW=cay.size
    iindy,iindx,iindz=np.mgrid[0:NROW,0:NCOL,0:NZ]
    ccay=cay[iindy];ccax=cax[iindx];ccaz=caz[iindz]
    #print NZ, NCOL,NROW
    #########
    #print NX,NXB, NY,NYB, NZ
    #NX is the number of columns in dens3
    #NY is the number of rows in dens3
    #the width of the cone should depend on the energy will adjust later
    #cax=np.linspace(mincone-25,maxcone+25,NCOL,'float32')
    #cay=np.linspace(-50,150,NROW,'float32')
    #caz=np.linspace(mincone-25,maxcone+25,NZ,'float32')
    vol=(abs(cax[1]-cax[0])*abs(cay[1]-cay[0])*abs(caz[1]-caz[0]))
    #print vol
    #rotate to match block rotation
    #rotate to col
    caxtmp=np.cos(colAng*np.pi/180.0)*ccax-np.sin(colAng*np.pi/180.0)*ccaz
    caztmp=np.sin(colAng*np.pi/180.0)*ccax+np.cos(colAng*np.pi/180.0)*ccaz
    #rotate gantry
    ccax=np.cos(rotang)*caxtmp-np.sin(rotang)*ccay
    caytmp=np.sin(rotang)*caxtmp+np.cos(rotang)*ccay
    #rotate to psa
    caxtmp=np.cos(PSAang*np.pi/180.0)*ccax-np.sin(PSAang*np.pi/180.0)*caztmp
    ccaz=np.sin(PSAang*np.pi/180.0)*ccax+np.cos(PSAang*np.pi/180.0)*caztmp
    #x=x, y=y, z=z
    calz=ccaz.reshape(ccaz.size)
    calx=caxtmp.reshape(caxtmp.size)
    caly=caytmp.reshape(caytmp.size)
    del iindx,iindz;del iindy;
    print len(calz),len(caly),len(calx)
    #del iindx;del iindy;
    calx=calx+isopos[0]
    caly=caly+isopos[1]
    calz=calz+isopos[2]
    #convert to cm for use in cuda
    cax=cax/10.0
    cay=cay/10.0
    caz=caz/10.0    
    cdy=np.arange(caly.min(),caly.max(),dy)
    cdx=np.arange(calx.min(),calx.max(),dy)
    cdz=np.arange(calz.min(),calz.max(),dy)
    cDose3d=np.zeros([cdy.size,cdx.size,cdz.size])
    #define starting position for plot
    curi=np.where(abs(zc-isopos[2])==np.min(abs(zc-isopos[2])))[0][0]
    plot_function()
#function to turn plotting of contours on and off        
def func(label):
    global contNames,contOnOff
    for i in xrange(len(contNames)):
        #print label,contNames[i]
        if label==contNames[i]:
            contOnOff[i]= not contOnOff[i]
            #print contOnOff[i]
    plot_function()
#function to turn plotting os boluses (and contours with overrides) on and off
def bfunc(label):
    global bolOnOff,bolPrevConv
    for i in xrange(len(bolNames)):
        #print label,bolNames[i],len(bolNames)
        if label==bolNames[i]:
            bolOnOff[i]= not bolOnOff[i]
            bolPrevConv[i]=False
            #print bolOnOff[i]
#the following inc functions are used to fine tune A, E0, alp, and bnc for each energy
#these are not used after intial set up and should probably be put in a "physics" module
#that stored the values in a file read in a runtime
def incA(event):
    global A
    A+=0.5
    print 'A= ',A
def decA(event):
    global A
    A-=0.5
    print 'A= ',A
def incE(event):
    global E0
    E0+=0.1
    print 'E0= ',E0
def decE(event):
    global E0
    E0-=0.1
    print 'E0= ',E0
def incal(event):
    global alp
    alp+=0.1
    print 'alp= ',alp
def decal(event):
    global alp
    alp-=0.1
    print 'alp= ',alp
def incbnc(event):
    global bnc
    bnc+=0.01
    print 'bnc= ',bnc
def decbnc(event):
    global bnc
    bnc-=0.01
    print 'bnc= ',bnc
#Change the calculation grid size    
def setCalcGrid(label):
    global dy
    s0=1.0;s1=1.5;s2=2.5;
    dydict={'1.0 mm':s0,'1.5 mm':s1,'2.5 mm':s2}
    dy=dydict[label]
    #dy sets the resolution in all dimensions, x, y, and z for simplicity
    print "set dy to ",dy
#function to define the normalization point for the simulation, either the reference point
#or the GTV but this isn't working quite right
def setNormPoint(label):
    global npnt,cDose3d,MU,mGtvDose
    global refDose,doseWantAtdmax
    s0='GTV';s1='RefPnt';
    refDict={'RefPoint':s1,'GTV':s0}
    npnt=refDict[label]
    try:
        if npnt=='RefPoint':
            cDose3d=cDose3d*refDose
        elif npnt=='GTV' and (not np.isnan(mGtvDose)):
            cDose3d=cDose3d*mGtvDose
        else:
            return 0
        print cDose3d.max(), refDose
        if npnt=='GTV' and not(np.isnan(mGtvDose)):
                #print "Max Dose and Ref Pt dose differ my more than 40%"
            print "scaling by mean dose to gtv"
            MU=doseWantAtdmax*refDose30*vol/(mGtvDose*vol30)
            cDose3d=cDose3d/mGtvDose
        else:
            if np.isnan(mGtvDose):
                print "No GTV defined scaling to RefPoint"
            cDose3d=cDose3d/refDose
            print cDose3d.max(), refDose
            MU=doseWantAtdmax*refDose30*vol/(refDose*vol30)
        plot_function()
    except:
        print "No reference points defined yet, please Open Files and run Calc"
        return 0
        #sfile=open('refDose30.dat','w')
        #sfile.write(repr(refDose))
        #sfile.close()
    #cDose3d=cDose3d/cDose3d[rind,cind,indz]
    print "2nd Check MU= ",MU
    print "%diff with TPS= ",100*(MU-tpsMU)/tpsMU
    plot_function()
#se the number of particles per cell for the simultaion
#these givedecent results no matter which dy is chosen        
def setNpcell(label):
    global npcell,slideNp
    s0=10000;s1=20000;s2=30000;s3=50000;
    npdict={'10000':s0,'20000':s1,'30000':s2,'50000':s3}
    npcell=npdict[label]
    print npcell
#function to exit        
def quitter( event):
    sys.exit()
#############
#function to run the simultaion
def calmain( event):
    global npix3d,cDose3d,A,E0,alp,npcell
    global calx,caly,calz,cax,cay,caz,cind,rind,indz,cdy,cdx,cdz
    global MU,vol,mGtvDose,refDose,doseWantAtdmax,nomE,dose3,refDose30
#######
#get gtv needed to normalize to gtv if chosen
    gtvpnts=[];gtvs=0;
    try:
        for i in xrange(len(rs.ROIContours)):
            try:
                val=rs.RTROIObservations[i].RTROIInterpretedType
                if (val=='GTV'):
                    print rs.StructureSetROIs[i].ROIName+" is of Type "+val
                    gtvs+=1
        #if(rs.StructureSetROIs[i].ROIName.lower()=="Body".lower()):
                try:
                    for j in xrange(len(rs.ROIContours[i].Contours)):
                        gtvpnts.append([rs.ROIContours[i].ReferencedROINumber, \
                                            rs.ROIContours[i].Contours[j].NumberofContourPoints, \
                                            rs.ROIContours[i].Contours[j].ContourData,rs.StructureSetROIs[i].ROIName])
                except:
                    print "Exception taken"
            except:
                val=''
    except:
        return 0
        
    #if more than one bolus make it selectable
    # need a file that outputs the block x,y
    # blckx, blcky and blckz with correct velocity
    #how do electron "bounce" in the medium

    #modify npix3d to include manually set contour densities
    xyco=np.dstack(np.meshgrid(xc,yc)).reshape(-1,2)
    xco=xyco[:,0].astype('float32')
    yco=xyco[:,1].astype('float32')
    #speed this u
    #for i in range(npix3d.shape[2]):
    #######
    #code for calling pnpoly
    #######
    #os.chdir(cdir)
    tp0=time.time()
    supcode="""
    extern "C" {
    void pninpoly(int nvert, float * vertx, float * verty, int ntest, \
    float *testx, float *testy, int *mask,float dy);
    }//
    """
    #pninpoly is the code used to find is a point is inside a polygon Called on the GPU
    #for speed
    pcode="""//
    pninpoly(nvert,bpsx,bpsy,ntest,xco,yco,mmask,ddy);
    //23
    """
    pcneed=['nvert','bpsx','bpsy','ntest','xco','yco','mmask','ddy']
    libs=['cuda','cudart']
    dirs=['/usr/lib/',cudadir]
    rlibs=['/usr/lib/',cudadir]
    ecomparg=['-O3 -mtune=native']
    eobj=[progdir+'/pninpoly.o']
    mmask=np.zeros(xco.size).astype('intc')
    ntest=int(xco.size)
    
    #forcing recompile will take a long time in the loop should probably
    #make a function that can be recompiled once then called here
    #as written adds only three seconds versus calling the builtin
    #points_in_poly which takes a long time 200+ seconds
    name=""
    convDens=False;prevConv=False;
    npix3d[:,:,:]=npix3do[:,:,:]
    #convert the points in npix3d to over ride values if needed
    if(len(bolpnts)!=0):
        #need a way to undo this once it has been done
        bindx=[]
        for j in xrange(len(bolpnts)):
            #print "j is",j
			#dy was 1e-6
            i=np.where(abs(zc-bolpnts[j][2][2])<1.0e-6)[0][0];
            bindx.append(i)
        bindx=np.asarray(bindx)
        print bindx
        bndx=np.bincount(bindx)
        ent=0;
        for j in xrange(len(bolpnts)):
            i=np.where(abs(zc-bolpnts[j][2][2])<1.0e-6)[0][0];
            bpsverts=np.asarray(bolpnts[j][2]).reshape([bolpnts[j][1],3])[:,:-1];
            #the code below can be un commented if the contours are not overridden properly
            #for example in a scalp the bolus has two contours to define it, you only want to
            #over ride the points between the two bolus contours not all the points inside the
            #bolus contours, need to have an undo and a special button for this
            """try:
                if(bndx[i]==2 and ent==0):
                    bpsvvs=np.asarray(bpsverts)
                    ln1=len(bpsvvs)
                    ent=1
                    continue
                elif(bndx[i]==2 and ent==1):
                    ln2=len(bpsverts)
                    bpsverts=np.append([[0,0]],bpsverts)
                    bpsvvs=np.append([[0,0]],bpsvvs)
                    bpsverts=np.append(bpsverts,bpsvvs)
                    bpsverts.shape=([ln1+ln2+2,2])
                    ent=0
                    print bpsverts.shape
            except:
                continue"""
            bpsx=bpsverts[:,0].astype('float32');bpsy=bpsverts[:,1].astype('float32');
            nvert=int(bpsx.size);
            prnt=True
            #print "I=",i," J=",j
            if bolpnts[j][4]!=name:
                name=bolpnts[j][4]
                prnt=True
                for k in xrange(len(bolNames)):
                    if bolNames[k]==name:
                        convDens=bolOnOff[k]
                        prevConv=bolPrevConv[k]
                        #bolPrevConv[k]= True
            else:
                prnt= False
            if (j==0 or prnt) and convDens and (not prevConv):
                print "converting contour "+bolpnts[j][4]+" to e density of ",bolpnts[j][3]
        #print bpsx.min(),bpsx.max(),bpsy.min(),bpsy.max()
        #type_converters=weave.converters.blitz,
            #print convDens,j
            #weave compiles the gpu pninpoly code to be used in python
            if convDens and (not prevConv):
                weave.inline(pcode,pcneed,support_code=supcode,\
                                         libraries=libs,library_dirs=dirs,extra_compile_args=ecomparg,extra_objects=eobj,\
                                         runtime_library_dirs=rlibs,compiler='gcc')
                mask=mmask==1
                #p=Path(bpsverts)
                #mask=p.contains_points(xyco)
                #mask=points_inside_poly(xyco,bpsverts)
                #override the points in mask with the override value in bolpnts[j][3]
                npix3d[:,:,i].reshape(xco.size)[mask]=bolpnts[j][3];
       
    print "time to convert ",time.time()-tp0
    ################
    try:
        del bpsverts;del bpsx;del bpsy;del tp0;
    except:
        print "no bolus points"
        #NROW=75;NCOL=75;NZ=75;
    ########### 
    #auto calculate the grid sixe, may have changed if dy changed
    dist=np.sqrt((isopos[0]-refpos[0])**2+(isopos[1]-refpos[1])**2+(isopos[2]-refpos[2])**2)
    rangep=beamEnergy/.20
    dmax=beamEnergy/.40
    cay=np.arange(-50,round(dist+(rangep-dmax)+50),dy)
    #cay=np.linspace(-50,150,NROW,'float32')
    #dy=abs(cay[1]-cay[0])
    #lateral extent
    caz=np.arange(mincone-30,maxcone+30,dy)
    cax=np.arange(mincone-30,maxcone+30,dy)
    NZ=caz.size
    NCOL=cax.size
    NROW=cay.size
    iindy,iindx,iindz=np.mgrid[0:NROW,0:NCOL,0:NZ]
    ccay=cay[iindy];ccax=cax[iindx];ccaz=caz[iindz]
    #print NZ, NCOL,NROW
    #########
    #print NX,NXB, NY,NYB, NZ
    #NX is the number of columns in dens3
    #NY is the number of rows in dens3
    #the width of the cone should depend on the energy will adjust later
    #cax=np.linspace(mincone-25,maxcone+25,NCOL,'float32')
    #cay=np.linspace(-50,150,NROW,'float32')
    #caz=np.linspace(mincone-25,maxcone+25,NZ,'float32')
    vol=(abs(cax[1]-cax[0])*abs(cay[1]-cay[0])*abs(caz[1]-caz[0]))
    print "vol= ",vol
    #rotate to match block rotation
    #rotate to col
    caxtmp=np.cos(colAng*np.pi/180.0)*ccax-np.sin(colAng*np.pi/180.0)*ccaz
    caztmp=np.sin(colAng*np.pi/180.0)*ccax+np.cos(colAng*np.pi/180.0)*ccaz
    #rotate gantry
    ccax=np.cos(rotang)*caxtmp-np.sin(rotang)*ccay
    caytmp=np.sin(rotang)*caxtmp+np.cos(rotang)*ccay
    #rotate to psa
    caxtmp=np.cos(PSAang*np.pi/180.0)*ccax-np.sin(PSAang*np.pi/180.0)*caztmp
    ccaz=np.sin(PSAang*np.pi/180.0)*ccax+np.cos(PSAang*np.pi/180.0)*caztmp
    #x=x, y=y, z=z
    rrcaz=(ccaz[0,0,:]+isopos[2]).astype('float32')
    rrcax=(caxtmp[0,:,0]+isopos[0]).astype('float32')
    rrcay=(caytmp[:,0,0]+isopos[1]).astype('float32')
    calz=ccaz.reshape(ccaz.size)
    calx=caxtmp.reshape(caxtmp.size)
    caly=caytmp.reshape(caytmp.size)
    del iindx,iindz;del iindy;
    print len(calz),len(caly),len(calx)
    #del iindx;del iindy;
    calx=calx+isopos[0]
    caly=caly+isopos[1]
    calz=calz+isopos[2]
    #convert to cm for use in cuda
    cax=cax/10.0
    cay=cay/10.0
    caz=caz/10.0
    #get dens array from npix3d dens3 will be used in the simultion and is rotated so that
    #y is down with respect to the beams eye view
    try:
        del dens3;
        dens3=np.zeros([NROW,NCOL,NZ],dtype=np.intc);
    except: 
        dens3=np.zeros([NROW,NCOL,NZ],dtype=np.intc);
    #dens3[:,:,:]=npix3d[ndxmnt:ndxmxt,ndymnt:ndymxt,ndzmn:ndzmx]
    #assuming that all X,Y are the same in each image
    #Weave this section so it will actually finish
    support_code="""
    #include <cmath>
    #include <iostream>
    #include <omp.h>
    using namespace std;
    int wher(double x[], double val, int n){
        int i,indx=0;
        double dx2=abs(x[1]-x[0])/2.0;
        for (i=0;i<n;i++){
            if(val>(x[i]-dx2) && val <=(x[i]+dx2)){
                return i;
            }
        }
        return indx;
    }
    //utility to find an index in an array 1D on the GPU
    //binary search
    int wher0(const double x[],const double val, const int n){
    int nabove=n+1,nbelow=0,mid;  
    //pass in dlo/2
    while(nabove -nbelow >1){
        mid=(nabove+nbelow)/2;
        if(val==x[mid-1])return mid-1;
        if(val< x[mid-1]){ nabove=mid;}
        else{ nbelow=mid;}
        }
        return nbelow-1;
}

    """
    #the code below needs a check that will ensure the density gets to be -1000
    #if calx etc are outside of the xc,yc,zc.  Right now just chooses the npix3d
    #value that is closest to the calx...since there is a border aroung npix3d
    #of air then it works out.
    #uses openmp to do this in a parallel manner can also be done on the gpu
    gcode="""//4
    int i,j,k;
    int rind,cind,indz;
    int rindm,cindm,indzm;
    int rindp,cindp,indzp;
    double dy2=dy;
    //double *x, *y, *z;
    //x=new double[NXC];
    //y=new double[NYC];
    //z=new double[NZC];
    //for(i=0;i<NXC;i++) x[i]=xc[i];
    //for(i=0;i<NYC;i++) y[i]=yc(i);
    //for(i=0;i<NZC;i++) z[i]=zc(i);
    #pragma omp parallel for private(i,j,k,indz,rind,cind,indzm,indzp,rindm,rindp,cindm,cindp)
    for(i=0;i<NROW;i++){
        for(j=0;j<NCOL;j++){
            for(k=0;k<NZ;k++){
                indzm=wher(zc,calz[k+j*NZ+i*NCOL*NZ]-dy2,NZC);
                indz=wher(zc,calz[k+j*NZ+i*NCOL*NZ],NZC);
                indzp=wher(zc,calz[k+j*NZ+i*NCOL*NZ]+dy2,NZC);
                rindm=wher(yc,caly[k+j*NZ+i*NCOL*NZ]-dy2,NYC);
                rind=wher(yc,caly[k+j*NZ+i*NCOL*NZ],NYC);
                rindp=wher(yc,caly[k+j*NZ+i*NCOL*NZ]+dy2,NYC);
                cindm=wher(xc,calx[k+j*NZ+i*NCOL*NZ]-dy2,NXC);
                cind=wher(xc,calx[k+j*NZ+i*NCOL*NZ],NXC);
                cindp=wher(xc,calx[k+j*NZ+i*NCOL*NZ]+dy2,NXC);
                if(indz==0 || indz==NZC || rind==0 || rind==NYC || cind==0 || cind==NXC){
                //npix3d[indz+cind*NZC+rind*NXC*NZC];
                    dens3[k+j*NZ+i*NCOL*NZ]=npix3d[indz+cind*NZC+rind*NXC*NZC];
                //    npix3d[indzm+cind*NZC+rind*NXC*NZC]+npix3d[indz+cind*NZC+rind*NXC*NZC]+npix3d[indzp+cind*NZC+rind*NXC*NZC]+\
                //    npix3d[indz+cindm*NZC+rind*NXC*NZC]+npix3d[indz+cind*NZC+rind*NXC*NZC]+npix3d[indz+cindp*NZC+rind*NXC*NZC]+\
                //    npix3d[indz+cind*NZC+rindm*NXC*NZC]+npix3d[indz+cind*NZC+rind*NXC*NZC]+npix3d[indz+cind*NZC+rindp*NXC*NZC])/9.0;
                } else {
                    dens3[k+j*NZ+i*NCOL*NZ]=(\
                    npix3d[indzm+cind*NZC+rind*NXC*NZC]+npix3d[indz+cind*NZC+rind*NXC*NZC]+npix3d[indzp+cind*NZC+rind*NXC*NZC]+\
                    npix3d[indz+cindm*NZC+rind*NXC*NZC]+npix3d[indz+cind*NZC+rind*NXC*NZC]+npix3d[indz+cindp*NZC+rind*NXC*NZC]+\
                    npix3d[indz+cind*NZC+rindm*NXC*NZC]+npix3d[indz+cind*NZC+rind*NXC*NZC]+npix3d[indz+cind*NZC+rindp*NXC*NZC])/9.0;
                }
            }
        }
    }
    ////delete x;delete y; delete z;
    """
    cneed=['NROW','NCOL','NZ','dens3','NXC','NYC','NZC','xc','yc','zc','npix3d','calx','caly','calz','dy']
    NXC=int(xc.size);NYC=int(yc.size);NZC=int(zc.size);
    ecomparg=['-O3 -mtune=native -fopenmp']
    #type_converters=weave.converters.blitz,
    weave.inline(gcode,cneed,support_code=support_code, libraries=['gomp'],\
                     extra_compile_args=ecomparg,compiler='gcc')
#rewrite this to execute on the GPU compare to weave, and f2py fortran
#add f2py call for calcmain see if that is better and will work on laptop
#changed lead to bone for testing 06-26-13 11.4
#convert HU to physical density using openmp could probably do all of this in one function call
#but this works
    code="""
    int i,j,k;
    int mat;//,*cdens;
    #pragma omp parallel for private(i,j,k,mat)
    for (i=0;i<NROW;i++){
        for (j=0;j<NCOL;j++){
            for (k=0;k<NZ;k++){
                mat=dens3[k+j*NZ+i*NCOL*NZ];
                //air
                if (mat<=-875) {dens3[k+j*NZ+i*NCOL*NZ]=0;
                //ilung
                } else if(mat<=-650){dens3[k+j*NZ+i*NCOL*NZ]=1;
                //elung
                } else if(mat<=-350){dens3[k+j*NZ+i*NCOL*NZ]=2;
                //fat
                } else if(mat<=-20){dens3[k+j*NZ+i*NCOL*NZ]=3;
                //water
                } else if(mat<=20){dens3[k+j*NZ+i*NCOL*NZ]=4;
                //muscle
                } else if(mat<=185){dens3[k+j*NZ+i*NCOL*NZ]=5;
                //cbone
                //changed for expanders or wires
                } else if(mat<=350){dens3[k+j*NZ+i*NCOL*NZ]=6;
                //sbone
                } else if(mat<=1500){dens3[k+j*NZ+i*NCOL*NZ]=7;
                //lead but called sbone for now
                } else if(mat<=5000){dens3[k+j*NZ+i*NCOL*NZ]=7;
                } else {
                //make it air if not found
                dens3[k+j*NZ+i*NCOL*NZ]=0;
                }
            }
        }
    }
    """
    cneed=['NROW','NCOL','NZ','dens3']
    ecomparg=['-O3 -mtune=native -fopenmp']
    #type_converters=weave.converters.blitz,
    weave.inline(code,cneed,support_code=support_code, libraries=['gomp'],\
                     extra_compile_args=ecomparg,compiler='gcc')
    
    dose3=np.zeros(dens3.shape).astype(np.float32)
    fl=open('dens3.dat','wb')
    #store the density file for later
    dbin=dens3.tostring()
    fl.write('%d\n%s\n'%(len(dbin),str(dens3.shape)))
    fl.write(dbin)
    fl.close()
    ####
    calcscode="""
    #include <cmath>
    #include <iostream>
    #include <cuda.h>
    
    using namespace std;

    extern "C" {
    void pymain (float E0,int NROW,int NCOL,int NZ, int npcell, int nMP, int nthreads, float alp, \
            float A, float bnc, int dens[],	\
            float x[], float y[], float z[], float dose[], float strt, \
            float stop, float xmin, float xmax,float blkx[], float blky[], \
            int nblk,float ymin,float ymax,float zmin,float zmax,float phper,float felec[],\
            float sx, float ix, float sy, float iy, float sz, float iz);
    }
    float maxxval(float x[],long n){
        float maxx;
        long i;
        maxx=x[0];
        for (i=0;i<n;i++){
            if(x[i]>=maxx){
            maxx=x[i];
            }
        }
        return(maxx);
    }
    /*int wher(double x[], double val,int n){
        int i,indx=0;
        double dx2=abs(x[1]-x[0])/2.0;
        for (i=0;i<n;i++){
            if(val > (x[i]-dx2) && val < (x[i]+dx2)){
                return indx=i;
            }
        }
        return indx;
    }*/
    """
    calccode="""
    //345
    //cout<<"testing some c++ code "<<!((int)(1.0f-0.2f)*(1.0f+1.0f))<<endl;
    cout<<"before pymain"<<endl;
    pymain(E0,NROW,NCOL,NZ,npcell,nMP,nthreads,alp,A,bnc,dens3,cax,cay,caz,dose3,\
    strt,sstop,cxmin,cxmax,bxx,byy,nblk,cymin,cymax,czmin,czmax,phper,felec,sx,ix,sy,iy,sz,iz);
    cout<<maxxval(dose3,NROW*NCOL*NZ)<<endl;
    ////cout<<"after pymain"<<endl;01

    """
    t0=time.time()
    #npcell=100000
    nomE=rp.Beams[0].ControlPoints[0].NominalBeamEnergy
    #A=20;E0=10.7;alp=1.6;
    #A=35.0;E0=10.7;alp=1.6;phper=0.02;
    #A=A*12.0/nomE;E0=E0*nomE/12.0;alp=alp*nomE/12.0;phper=phper*nomE/12.0;
    #cone factors to scale each cone to 1 cGy at it's dmax for each energy 0.984
    #need permachine factors
    coneFacs={'A066':0.96,'A106':1.00,'A156':0.997,'A206':1.015,'A256':1.015,\
              'A069':.984,'A109':1.00,'A159':0.993,'A209':0.985,'A259':0.964,\
              'A0612':0.973,'A1012':1.00,'A1512':0.989,'A2012':0.979,'A2512':0.952,\
              'A0615':0.973,'A1015':1.00,'A1515':0.989,'A2015':0.979,'A2515':0.952}
    """coneFacs={'A066':1.0,'A106':1.00,'A156':1.0,'A206':1.0,'A256':1.0,\
              'A069':1.0,'A109':1.00,'A159':1.0,'A209':1.0,'A259':1.0,\
              'A0612':1.0,'A1012':1.00,'A1512':1.0,'A2012':1.0,'A2512':1.0}"""
    #used to scale simultaion dose based on energy and grid size selected
    #this should probably be a single value for the highest resolution but if one had a
    #slow computer they might calibrate using a lower resolution need to put these in a file
    #with the A, E0, Alp, phper, bnc, etc
    refFac={'6E2.5':3.6778505588E-011,'6E1.5':3.59170158892e-11,'6E1.0':3.59936e-11,\
            '9E2.5':3.3765505438E-011,'9E1.5':3.25399637791e-11,'9E1.0':3.27292E-11,\
            '12E2.5':3.3467284499E-011,'12E1.5':3.17686224549e-11,'12E1.0':3.13e-11,\
            '15E2.5':3.3467284499E-011,'15E1.5':3.17686224549e-11,'15E1.0':2.9e-11}
    """refFac={'6E5.0':3.5E-011,'6E2.5':3.5E-011,'6E1.0':3.5E-011,\
            '9E5.0':3.2E-011,'9E2.5':3.2E-011,'9E1.0':3.2E-011,\
            '12E5.0':3.2E-011,'12E2.5':3.2E-011,'12E1.0':3.2E-011}"""
    #A is the scattering amplitude, E0 is the nominal energy, and alp is the spread of a 
    #gaussian rand
    #        om distribution each one is defined for each energy during the calibration
    #phase and not changed after
    #A=A*9.0/nomE; E0=E0*nomE/9.0;alp=alp*nomE/9.0        
    print A, E0, alp
    
    if (abs(nomE-12.0)<1e-6):
        #06-26-13;
        A=33.0;E0=10.5;alp=1.4;
        #A=1;E0=11.5;alp=1;
        
        refDose30=refFac[str(int(nomE))+'E'+str(dy)];vol30=15.625;phper=0.025;#refDose30=3.11711e-13;
    elif (abs(nomE-9.0)<1e-6):
        #06-26-13
        A=42.166;E0=7.825;alp=1.07
        refDose30=refFac[str(int(nomE))+'E'+str(dy)];vol30=15.625;phper=0.01
    elif (abs(nomE-6.0)<1e-6):
		#05-13-14
        A=60.0;E0=5.7;alp=0.7;
        refDose30=refFac[str(int(nomE))+'E'+str(dy)];phper=0.005;
    elif (abs(nomE-15.0)<1e-6):
        A=26.4;E0=12.75;alp=1.85
        refDose30=refFac[str(int(nomE))+'E'+str(dy)];phper=0.035
    elif (abs(nomE-18.0)<1e-6):
        refDose30=3.0e-13;E0=15.0;phper=0.05;
    
    print A, E0, alp
    print "TPS Energy ",nomE
    #number of blocks and threads for CUDA need to make these interactive so user can fine 
    #tune it
    nMP=np.intc(1024);nthreads=np.intc(1024);
    #nMP=int(480);nthreads=int(576);
    doseWantAtdmax=(tpsDose)
    print doseWantAtdmax, tpsDose
    felec=np.zeros(1).astype('float32');
    #convert to 1D array for use in CUDA
    dose3=dose3.reshape(NROW*NCOL*NZ).astype('float32');
    dens3=dens3.reshape(NROW*NCOL*NZ).astype('intc');
        #cax=cax.astype('float32')
        #cay=cay.astype('float32')
        #caz=caz.astype('float32')
    #needs to be defined consistently between here and cuda
    mx=min(abs(cax[1]-cax[0]),abs(cay[1]-cay[0]),abs(caz[1]-caz[0]));
    E0=np.float32(E0);alp=np.float32(alp);A=np.float32(A);bnc=np.float32(0.5/dy);
    cax=cax.astype('float32');cay=cay.astype('float32');caz=caz.astype('float32')
    sx,ix,p,q,r=stats.linregress(range(cax.size),cax);sy,iy,p,q,r=stats.linregress(range(cay.size),cay)
    sz,iz,p,q,r=stats.linregress(range(caz.size),caz);
    sx=np.float32(sx);ix=np.float32(ix);sy=np.float32(sy);iy=np.float32(iy);sz=np.float32(sz);iz=np.float32(iz)
    """#Where is the refpoint in the simulation index rotate from psa
    caxtmp=np.cos(-PSAang*np.pi/180.0)*(refpos[0]-isopos[0])-np.sin(-PSAang*np.pi/180.0)*(refpos[2]-isopos[2])
    caztmp=np.sin(-PSAang*np.pi/180.0)*(refpos[0]-isopos[0])+np.cos(-PSAang*np.pi/180.0)*(refpos[2]-isopos[2])
    #rotate from gantry
    ccax=np.cos(-rotang)*caxtmp-np.sin(-rotang)*(refpos[1]-isopos[1])
    caytmp=np.sin(-rotang)*caxtmp+np.cos(-rotang)*(refpos[1]-isopos[1])
    #rotate from col
    caxtmp=np.cos(-colAng*np.pi/180.0)*ccax-np.sin(-colAng*np.pi/180.0)*caztmp
    caztmp=np.sin(-colAng*np.pi/180.0)*ccax+np.cos(-colAng*np.pi/180.0)*caztmp

    refposc=np.zeros(3).astype('float32')
    refposc[0]=caxtmp/10.0;refposc[1]=caytmp/10.0;refposc[2]=caztmp/10.0
    print "refposc= ", refposc"""
    strt=np.float32(mincone/10.0);sstop=np.float(maxcone/10.0);
    #cxmin=grid min in x and y same with cxmax
    cxmin=np.float32(cax.min());cxmax=np.float32(cax.max());
    cymin=np.float32(cay.min());cymax=np.float32(cay.max());
    czmin=np.float32(caz.min());czmax=np.float32(caz.max());
    #in Dicom block is at iso which is centered at CAX,CAY,CAZ properly by definition
    #move xx back for input into cuda, zz wasn't changed
    #block x and z positions
    bxx=((xx).astype('float32'))/10.0;byy=(zz).astype('float32')/10.0;nblk=np.intc(len(bxx));
    #print "Blocks: ",bxx,byy
    MU=(0);       
    clibs=['cuda','cudart','gsl','gslcblas']
    cdirs=['/usr/lib/',cudadir]
    crlibs=['/usr/lib/',cudadir,progdir]
    incdirs=[cudainc]
    cecomparg=['-O3 -mtune=native']
    ceobj=[progdir+'/PYelec8-7.5.o']
    #npcell=100000
    #python defined variables needed in the GPU code
    calcneed=['E0','NROW','NCOL','NZ','npcell','nMP','nthreads','alp','A','bnc','felec',\
                  'dens3','cax','cay','caz','dose3','strt','sstop','cxmin','cxmax','cymin',\
                  'cymax','czmin','czmax','bxx','byy','nblk','phper','sx','ix','sy','iy','sz','iz']
    #call CUDA calculation
    weave.inline(calccode,calcneed,support_code=calcscode, \
                     libraries=clibs,library_dirs=cdirs,
                     extra_compile_args=cecomparg,
                     extra_objects=ceobj,\
                     runtime_library_dirs=crlibs,include_dirs=incdirs,compiler='gcc',verbose=2)
    #print "NROW= ",NROW, " NCOL= ",NCOL, " NZ= ",NZ, " nblk= ",nblk
    #pymain.pymain(E0,np.intc(npcell),nMP,nthreads,alp,A,bnc,dens3,cax,cay,caz,dose3,strt,sstop,cxmin,cxmax,\
    #                          bxx,byy,cymin,cymax,czmin,czmax,phper,felec,sx,ix,sy,iy,sz,iz)#,nrow=np.intc(NROW),ncol=np.intc(NCOL),nz=np.intc(NZ),nblk=np.intc(nblk))
    t1=time.time()
        #pyCudaElec.pymain(E0,npcell,nMP,nthreads,alp,A,bnc,dens3,cax,cay,caz,dose3,strt,\
        #sstop,cxmin,cxmax,bxx,byy,cymin,cymax,czmin,czmax)
      #dose3/=totelec
    print npcell
    #dose3 is the simulation dose rotated so that y down is beams eye view
    dose3=dose3*1.602e-10*coneFacs[cone+str(int(nomE))]/float(dy*npcell)#*coneFacs[cone+str(int(nomE))]
    #dose3=dose3*coneFacs[cone+str(int(nomE))]*(abs(sstop-strt)**2)/float(dy**3)#/float(felec[0])
    """if dy==5:
        cDose3d=dose3.reshape(NROW,NCOL,NZ)/1.29;#*np.float64((abs(sstop-strt)**2))/vol;
    elif dy==2.5:
        cDose3d=dose3.reshape(NROW,NCOL,NZ)/1.11;#*np.float64((abs(sstop-strt)**2))/vol;
    elif dy==1:"""
    cDose3d=dose3.reshape(NROW,NCOL,NZ)#*np.float64((abs(sstop-strt)**2))/vol;
        
    print "Max dose= ",dose3.max()
    print "Max cdose3d= ",cDose3d.max()
    #rind is the y index that corresponds to 0, cind is the x index that corresponds to 0
    #indz is the z index that corresponds to 0
    #these are consistent with Eclipse z,x are planes in the beam's eye view
    #y is increasing down away from the beam
    rind=whr(cay,0,cay.size)
    cind=whr(cax,0,cax.size)
    indz=whr(caz,0,caz.size)
    dose3=dose3.reshape(NROW*NCOL*NZ)
    print "dens at isopos=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at isopos=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,.25,cay.size)
    print "dens at 2.5=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 2.5=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,.5,cay.size)
    print "dens at 5=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 5=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,.6,cay.size)
    print "dens at 6=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 6=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,.7,cay.size)
    print "dens at 7=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 7=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,.8,cay.size)
    print "dens at 8=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 8=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,.9,cay.size)
    print "dens at 9=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 9=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.0,cay.size)
    print "dens at 10.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 10.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.1,cay.size)
    print "dens at 11=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 11=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.2,cay.size)
    print "dens at 12=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 12=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.3,cay.size)
    print "dens at 13=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 13=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.4,cay.size)
    print "dens at 14=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 14=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.5,cay.size)
    print "dens at 15=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 15=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.6,cay.size)
    print "dens at 16=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 16=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.7,cay.size)
    print "dens at 17=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 17=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.8,cay.size)
    print "dens at 18=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 18=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,1.9,cay.size)
    print "dens at 19=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 19=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.0,cay.size)
    print "dens at 20.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 20.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.1,cay.size)
    print "dens at 21.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 21.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.2,cay.size)
    print "dens at 22.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 22.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.3,cay.size)
    print "dens at 23.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 23.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.4,cay.size)
    print "dens at 24.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 24.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.5,cay.size)
    print "dens at 25.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 25.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.6,cay.size)
    print "dens at 26.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 26.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.7,cay.size)
    print "dens at 27.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 27.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.8,cay.size)
    print "dens at 28.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 28.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,2.9,cay.size)
    print "dens at 29.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 29.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,3.0,cay.size)
    print "dens at 30.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 30.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,3.1,cay.size)
    print "dens at 31.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 31.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,3.2,cay.size)
    print "dens at 32.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 32.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,3.3,cay.size)
    print "dens at 33.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 33.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,3.4,cay.size)
    print "dens at 34.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 34.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,3.5,cay.size)
    print "dens at 35.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 35.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,5.0,cay.size)
    print "dens at 50.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 50.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,5.5,cay.size)
    print "dens at 55.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 55.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,6.0,cay.size)
    print "dens at 60.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 60.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,6.6,cay.size)
    print "dens at 66.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 66.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    rind=whr(cay,7.1,cay.size)
    print "dens at 71.0=",dens3[indz+cind*NZ+rind*NZ*NCOL]
    print "dose at 71.0=",dose3[indz+cind*NZ+rind*NZ*NCOL]
    #convert back to mm
    #del tdose;   
    
    cax*=10.0;cay*=10.0;caz*=10.0

    #x=x, y=y, z=z
    # find points in gtv
        #cind=np.where(abs((cax+isopos[0])-refpos[0])<=min(abs((cax+isopos[0])-refpos[0])))[0][0];
    cind=whr(calx,refpos[0],calx.size)
        #rind=np.where(abs((cay+isopos[1])-refpos[1])<=min(abs((cay+isopos[1])-refpos[1])))[0][0];
    rind=whr(caly,refpos[1],caly.size)
        #indz=np.where(abs((caz+isopos[2])-refpos[2])<=min(abs((caz+isopos[2])-refpos[2])))[0][0];
    indz=whr(calz,refpos[2],calz.size)
    print calx.size
    print "xs ",cind,calx[cind],refpos[0]
    print "ys ",rind,caly[rind],refpos[1]
    print "zs ",indz,calz[indz],refpos[2]
    
    #rotate from psa
    caxtmp=np.cos(-PSAang*np.pi/180.0)*(calx[cind]-isopos[0])-np.sin(-PSAang*np.pi/180.0)*(calz[indz]-isopos[2])
    caztmp=np.sin(-PSAang*np.pi/180.0)*(calx[cind]-isopos[0])+np.cos(-PSAang*np.pi/180.0)*(calz[indz]-isopos[2])
    #rotate from gantry
    ccax=np.cos(-rotang)*caxtmp-np.sin(-rotang)*(caly[rind]-isopos[1])
    caytmp=np.sin(-rotang)*caxtmp+np.cos(-rotang)*(caly[rind]-isopos[1])
    #rotate from col
    caxtmp=np.cos(-colAng*np.pi/180.0)*ccax-np.sin(-colAng*np.pi/180.0)*caztmp
    caztmp=np.sin(-colAng*np.pi/180.0)*ccax+np.cos(-colAng*np.pi/180.0)*caztmp
    #x=x, y=y, z=z
    cind=wher(cax,caxtmp,cax.size)
    rind=wher(cay,caytmp,cay.size)
    indz=wher(caz,caztmp,caz.size)
    print "xs ",cind,cax[cind],caxtmp
    print "ys ",rind,cay[rind],caytmp
    print "zs ",indz,caz[indz],caztmp

    print "Checks ",cind, rind, indz, cDose3d.max()
    """Put code here to adjust normalization method based on cone size?"""
    """if dy==2.5:
        refDose=cDose3d[rind,cind,indz]
    elif dy==1.5:
        refDose=cDose3d[rind-1:rind+1,cind-1:cind+1,indz-1:indz+1].mean()
    elif dy==1.0:"""
    refDose=cDose3d[rind-1:rind+1,cind-1:cind+1,indz-1:indz+1].mean()
    print "**********"
    print "%diff refdose= ",(refDose-cDose3d[rind,cind,indz])*200.0/(refDose+cDose3d[rind,cind,indz])
    print "**********"
    print "refdose= ",refDose,"cDose3d= ",cDose3d[rind,cind,indz]
    ##############
    dose3=dose3.reshape(NROW,NCOL,NZ)/refDose;
    #get new y index at 0
    print "init indices= ",rind,cind,indz    
    #get statistics on dose3
    minDose3=dose3[rind-2:rind+2,cind-1:cind+1,indz-1:indz+1].min()
    aveDose3=dose3[rind-2:rind+2,cind-1:cind+1,indz-1:indz+1].mean()
    maxDose3=dose3[rind-2:rind+2,cind-1:cind+1,indz-1:indz+1].max()
    stdDev=dose3[rind-2:rind+2,cind-1:cind+1,indz-1:indz+1].std()
    varDose3=dose3[rind-2:rind+2,cind-1:cind+1,indz-1:indz+1].var()
    print "Dose3 stats"
    print "Min= ",minDose3," Mean= ",aveDose3," Max= ",maxDose3
    print "StdDev= ",stdDev, " Variance= ",varDose3
    #plots dose plane at reference point
    #define which isodose lines to display, 10%, 30%, 50%, 75%, 90%, 95%, 100%, 105%, 110%
    isoLines=np.asarray([0.1,0.3,0.5,0.75,0.9,0.95,1.0,1.05,1.1])
    #colors for the line 10% to the left, 110% to the right
    clrs=('purple','cyan','blue','green','yellow','orange','red','magenta','white')
    V=isoLines*totalDose
    #new figure for contours at depth of refpoint
    """figP=pylab.figure()
    axP=figP.add_subplot(111)
    axP.contour(cax,caz,dose3[rind,:,:],levels=V,colors=clrs)
    figP.show()
    #get indices for 0 0 plot
    rind=wher(cay,0,cay.size)
    cind=wher(cax,0,cax.size)
    indz=wher(caz,0,caz.size)
    print "adj indices= ",rind,cind,indz
    #pylab.plot(cay,dose3[:,cind,indz])
    #will plot a line of the Depth dose
    fig2=pylab.figure()
    ax2=fig2.add_subplot(111)
    ax2.plot(cay[rind:],dose3[rind:,cind,indz])
    fig2.show()"""
    ######

    maxDose=cDose3d.max()
    cDose3d=cDose3d/refDose
    #mu is dose*refdose10x10*areacone/(refdosecone*area10x10)
    print "caz= ",abs(caz[1]-caz[2]),"cax= ",abs(cax[1]-cax[2]),"sstop= ",sstop,"strt= ",strt
    print "felec= ",felec
    print "Dose want Dmax= ",doseWantAtdmax
    MU=doseWantAtdmax*refDose30/(refDose)
    npnt='RefPnt'
    sfile=open('refDose30.dat','w')
    sfile.write(repr(refDose)+'\n')
    sfile.write(repr(vol)+'\n')
    sfile.write(repr(A)+'\n')
    sfile.write(repr(E0)+'\n')
    sfile.write(repr(dy)+'\n')
    sfile.write(repr(npcell)+'\n')
    sfile.write(repr(alp)+'\n')
    sfile.close()
    #cDose3d=cDose3d/cDose3d[rind,cind,indz]
    print "2nd Check MU= ",MU
    print "%diff with TPS= ",100*(MU-tpsMU)/tpsMU
    print "refdose= ",refDose,"vol= ",vol#, mGtvDose
    mask=np.isnan(cDose3d)
    cDose3d[mask]=0.0
    mask=np.isinf(cDose3d)
    cDose3d[mask]=1.0e6;
    
    #interpolate cDose3d here to higher res
    ddy=1.0
    caax=np.arange(cax.min(),cax.max()+ddy,ddy)
    caay=np.arange(cay.min(),cay.max()+ddy,ddy)
    caaz=np.arange(caz.min(),caz.max()+ddy,ddy)
    cDose3dr=np.zeros([caay.size,caax.size,caaz.size])
    print "interpolating Dose to Visualization"
    isupcode="""
#include <iostream>
#include <cuda.h>
#include <cmath>
#include <stdlib.h>
using namespace std;

extern "C" {
  void interps(int npnts, float *caax, float *caay, float *caaz,float *cDose3dr, \
	       int NX, float *cax, int NY, float *cay,			\
	       int NZ, float *caz, float * cDose3d,float * mnmxs);
} 

"""
    icode="""
//test007
cout<<"before interps"<<endl;
interps(NPNTS,ccax,ccay,ccaz,cDose3dr,NCOL,cax,NROW,cay,NZ,caz,cDose3d,mnmxs);
"""
    ######
    #will be used below also
    NNZ=int(caaz.size)
    NNCOL=int(caax.size)
    NNROW=int(caay.size)
    iindy,iindx,iindz=np.mgrid[0:NNROW,0:NNCOL,0:NNZ]
    try:
        del ccay,ccax,ccaz
    except MemoryError:
        print 'Memory Error'
    else:
        ccay=caay[iindy];ccax=caax[iindx];ccaz=caaz[iindz]
    ######
    NPNTS=int(ccay.size)
    ccay=ccay.reshape(iindy.size).astype('float32')
    ccax=ccax.reshape(iindx.size).astype('float32')
    ccaz=ccaz.reshape(iindz.size).astype('float32')
    cDose3dr=cDose3dr.reshape(NNZ*NNCOL*NNROW).astype('float32')
    mnmxs=np.asarray([cay.min(),cay.max(),cax.min(),cax.max(),caz.min(),caz.max()]).astype('float32')
    print NPNTS, mnmxs
    print ccay.min(),ccay.max(),ccax.min(),ccax.max(),ccaz.min(),ccaz.max()
    icneed=['cax','ccax','cay','ccay','caz','ccaz','cDose3d','cDose3dr','mnmxs','NPNTS','NCOL','NROW','NZ']
    #weave.inline(icode,icneed,support_code=isupcode)
    clibs=['cuda','cudart']
    cdirs=['/usr/lib/',cudadir]
    crlibs=['/usr/lib/',cudadir,progdir]
    incdirs=[cudainc]
    cecomparg=['-O3 -mtune=native']
    ieobj=[progdir+'/interp.o']
    weave.inline(icode,icneed,support_code=isupcode, \
                     libraries=clibs,library_dirs=cdirs,extra_compile_args=cecomparg,extra_objects=ieobj,\
                     runtime_library_dirs=crlibs,include_dirs=incdirs,compiler='gcc')
    print "Finished interp"
 #rotate cDose3dr coordinates to visualization
    ccay=ccay.reshape(iindy.shape)
    ccax=ccax.reshape(iindx.shape)
    ccaz=ccaz.reshape(iindz.shape)
    #rotate to match block rotation
    #rotate to col
    caxtmp=np.cos(colAng*np.pi/180.0)*ccax-np.sin(colAng*np.pi/180.0)*ccaz
    caztmp=np.sin(colAng*np.pi/180.0)*ccax+np.cos(colAng*np.pi/180.0)*ccaz
    #rotate gantry
    ccax=np.cos(rotang)*caxtmp-np.sin(rotang)*ccay
    caytmp=np.sin(rotang)*caxtmp+np.cos(rotang)*ccay
    #rotate to psa
    caxtmp=np.cos(PSAang*np.pi/180.0)*ccax-np.sin(PSAang*np.pi/180.0)*caztmp
    ccaz=np.sin(PSAang*np.pi/180.0)*ccax+np.cos(PSAang*np.pi/180.0)*caztmp
    #print caytmp[:,0,0]
    callz=ccaz.reshape(ccaz.size)+isopos[2]
    callx=caxtmp.reshape(caxtmp.size)+isopos[0]
    cally=caytmp.reshape(caytmp.size)+isopos[1]
    #cDose3d=cDose3d.reshape(cDose3d.size)
    print cally.min(),cally.max(),callx.min(),callx.max(),callz.min(),callz.max()
    dyc=abs(yc[0]-yc[1])
    #NCDY=abs(whr(yc,cally.min(),yc.size)-whr(yc,cally.max(),yc.size))
    #NCDX=abs(whr(xc,callx.min(),xc.size)-whr(xc,callx.max(),xc.size))
    #NCDZ=abs(whr(zc,callz.min(),zc.size)-whr(zc,callz.max(),zc.size))
    #was dyc*2
    cdy=np.arange(cally.min(),cally.max(),dyc*2.5).astype('float32')
    cdx=np.arange(callx.min(),callx.max(),dyc*2.5).astype('float32')
    cdz=np.arange(callz.min(),callz.max(),dyc*2.5).astype('float32')
    NCDY=cdy.size;NCDX=cdx.size;NCDZ=cdz.size
    cdose2d=np.zeros([NCDY,NCDX,NCDZ]).astype('float32')-np.float32(9999999)
    #iindy,iindx,iindz=np.mgrid[0:NCDY,0:NCDX,0:NCDZ]
    #ccdy=cdy[iindy].reshape(iindy.size).astype('float32');
    #ccdx=cdx[iindx].reshape(iindx.size).astype('float32');
    #ccdz=cdz[iindz].reshape(iindz.size).astype('float32');
    #print "counts ",counts.shape, cdose2d.shape
    print "put cDose3dr to cdose2d" 
    copyCDsupcode="""
#include <cuda.h>
#include <iostream>
using namespace std;
extern "C" {
  void copyCDs(int nx, float *cdx,int ny, float *cdy,int nz, float *cdz, float * cdose2d, \
	       int npcd3, float *rax, float *ray, float *raz, float * cDose3d,
	       float dx, float dy, float dz);
} 
int whr(float * x, float val,int n){
        int i,indx;
        float minn;
        minn=abs(x[0]-val);
        indx=0;
        for (i=1;i<n;i++){
            if(abs(x[i]-val)<minn){
                minn=abs(x[i]-val);
                indx=i;
            }
        }
        return indx;
}
"""
    copyCDcode="""
//I need to adjust these params to make them make sense I am just lazy
//test025
cout<<"before copyCDs2"<<endl;
copyCDs(NCDX, cdx,NCDY, cdy,NCDZ, cdz,cdose2d,NPNTScD3, callx,cally,callz,cDose3dr,dcx,dcy,dcz);
//
"""
    #NCCDX=int(ccdx.size)
    NPNTScD3=int(cally.size)
    print rrcay.min(),rrcay.max(),rrcax.min(),rrcax.max(),rrcaz.min(),rrcaz.max()
    print cally.min(),cally.max(),callx.min(),callx.max(),callz.min(),callz.max()
    dcx=float(abs(cdx[0]-cdx[1]));dcy=float(abs(cdy[0]-cdy[1]));dcz=float(abs(cdz[0]-cdz[1]));
    copyCDneeds=['callx','cally','callz','NPNTScD3','cdx','cdy','cdz','NCDX','NCDY','NCDZ','cdose2d','cDose3dr','dcx','dcy','dcz']
    #cdose2d=cdose2d.reshape(NCDY*NCDX*NCDZ)

    cplibs=['cuda','cudart']
    cpdirs=['/usr/lib/',cudadir]
    cprlibs=['/usr/lib64/',cudadir,progdir]
    cpincdirs=[cudainc]
    cpecomparg=['-O3 -mtune=native']
    cpeobj=[progdir+'/copyCD.o']
    print "weave copyCD"
    weave.inline(copyCDcode,copyCDneeds,support_code=copyCDsupcode, \
                     libraries=cplibs,library_dirs=cpdirs,extra_compile_args=cecomparg,extra_objects=cpeobj,\
                     runtime_library_dirs=cprlibs,include_dirs=cpincdirs,compiler='gcc')
    #print weave.inline(copyCDcode,copyCDneeds,support_code=copyCDsupcode,extra_compile_args=cpecomparg)
    cdose2d=cdose2d.reshape([NCDY,NCDX,NCDZ]).astype('float64')
    print cdose2d.min(),cdose2d.max()
    #min dose outside of body
    #print NNROW,NNCOL,NNZ
    #print rrcay,rrcax,rrcaz
    k=[]
    minCd=cdose2d.min()
    xyco=np.dstack(np.meshgrid(cdx,cdy)).reshape(-1,2)
    xco=xyco[:,0].astype('float32')
    yco=xyco[:,1].astype('float32')
    mmask=np.zeros(xco.size).astype('int32')
    ntest=int(xco.size)
    print "Sweep dose outside body"
    dz=abs(zc[0]-zc[1])
    print 'dz= ',dz
    
    if(len(contPoints)!=0):
        name="BODY"
        for i in xrange(len(contPoints)):
            if contPoints[i][3].upper()[0:4]==name:
                k.append(contPoints[i])
        #print k
        for j in k:
            try:
                i=np.where(abs(cdz-j[2][2])<=dz/2.0);
                #print i
                bpsverts=np.asarray(j[2]).reshape([j[1],3])[:,:-1];
                bpsx=bpsverts[:,0].astype('float32');bpsy=bpsverts[:,1].astype('float32');
                nvert=int(bpsx.size);
                weave.inline(pcode,pcneed,support_code=supcode,\
                                         libraries=libs,library_dirs=dirs,extra_compile_args=ecomparg,extra_objects=eobj,\
                                         runtime_library_dirs=rlibs,compiler='gcc')
                mask=mmask==1
                cdose2d[:,:,i[0][0]].reshape(xco.size)[~mask]=minCd;
            except:
                #print i
                continue
    
    print cdose2d.min(),cdose2d.max()        
    mask=np.isnan(cdose2d)
    cdose2d[mask]=0.0;
    mask=np.isinf(cdose2d)
    cdose2d[mask]=1.0e6;
    mask=cdose2d==-9999999
    cdose2d[mask]=0.0;
    dose3=cDose3d.copy()
    cDose3d=cdose2d
    #ave dose over chamber
    """k=[]
    avechamdose=[]
    if(len(contPoints)!=0):
        name="7002"
        for i in xrange(len(contPoints)):
            if contPoints[i][3].upper()==name.upper():
                k.append(contPoints[i])
        #print k
        for j in k:
            try:
                i=np.where(abs(cdz-j[2][2])<=dz/2.0);
                #print i
                bpsverts=np.asarray(j[2]).reshape([j[1],3])[:,:-1];
                bpsx=bpsverts[:,0].astype('float32');bpsy=bpsverts[:,1].astype('float32');
                nvert=int(bpsx.size);
                weave.inline(pcode,pcneed,support_code=supcode,\
                                         libraries=libs,library_dirs=dirs,extra_compile_args=ecomparg,extra_objects=eobj,\
                                         runtime_library_dirs=rlibs,compiler='gcc')
                mask=mmask==1
                avechamdose.append(cDose3d[:,:,i[0][0]].reshape(xco.size)[mask].mean());
            except:
                #print i
                continue
    vals=np.asarray(avechamdose)
    print "ave dose to chamber 7002", vals.min(), " ",vals.mean()," ",vals.max()
    """#ave dose over chamber
    """k=[]
    avechamdose=[]
    if(len(contPoints)!=0):
        name="7001"
        for i in xrange(len(contPoints)):
            if contPoints[i][3].upper()==name.upper():
                k.append(contPoints[i])
        #print k
        for j in k:
            try:
                i=np.where(abs(cdz-j[2][2])<=dz/2.0);
                #print i
                bpsverts=np.asarray(j[2]).reshape([j[1],3])[:,:-1];
                bpsx=bpsverts[:,0].astype('float32');bpsy=bpsverts[:,1].astype('float32');
                nvert=int(bpsx.size);
                weave.inline(pcode,pcneed,support_code=supcode,\
                                         libraries=libs,library_dirs=dirs,extra_compile_args=ecomparg,extra_objects=eobj,\
                                         runtime_library_dirs=rlibs,compiler='gcc')
                mask=mmask==1
                avechamdose.append(cDose3d[:,:,i[0][0]].reshape(xco.size)[mask].mean());
            except:
                #print i
                continue
    vals=np.asarray(avechamdose)
    print "ave dose to chamber 7001", vals.min(), " ",vals.mean()," ",vals.max()
    #ave dose over chamber
    k=[]
    avechamdose=[]
    if(len(contPoints)!=0):
        name="Old Chamber"
        for i in xrange(len(contPoints)):
            if contPoints[i][3].upper()==name.upper():
                k.append(contPoints[i])
        #print k
        for j in k:
            try:
                i=np.where(abs(cdz-j[2][2])<=dz/2.0);
                #print i
                bpsverts=np.asarray(j[2]).reshape([j[1],3])[:,:-1];
                bpsx=bpsverts[:,0].astype('float32');bpsy=bpsverts[:,1].astype('float32');
                nvert=int(bpsx.size);
                weave.inline(pcode,pcneed,support_code=supcode,\
                                         libraries=libs,library_dirs=dirs,extra_compile_args=ecomparg,extra_objects=eobj,\
                                         runtime_library_dirs=rlibs,compiler='gcc')
                mask=mmask==1
                avechamdose.append(cDose3d[:,:,i[0][0]].reshape(xco.size)[mask].mean());
            except:
                #print i
                continue
    vals=np.asarray(avechamdose)
    print "ave dose to chamber Old", vals.min(), " ",vals.mean()," ",vals.max()
    ####"""
    
    del dens3;del cdose2d; #del dose3;
    #del dens3;
    print 'cdose max',cDose3d.max()
    print 'time to produce image %g' %(t1-t0)
    print 'refdos= ',refDose
    print "felec= ",felec
    print "dy= ",dy
    plot_function()
    
def myGrid(pnts2,dat,cdxg,cdyg):
    return griddata(pnts2,dat,(cdxg,cdyg),method='cubic',fill_value=0.0);

def main():
    global ax,axContours,check,axBolus,bolus,npcell,dy,npnt,refDose30,vol30,slideNp,tpsMU
    global cax,cay,caz,cDose3d,refpos,fig1,bar,A,alp,E0,bnc
    #set some defaults
    #A=16;E0=5.4;alp=0.8;bnc=0.2;
    #for gaussian energy simulation
    #need to add gaussian profile to positions of x,z particles to match in air profile?!?
    #A=33.0;E0=10.8;alp=1.4;
    #print butn.shape
    #butn[:,0:10,0]=153
    #butn[:,0:10,1]=204
    #butn[:,0:10,2]=255
    A=42.166;E0=7.825;alp=1.07
    bnc=float(0.2);
    #for uniform energy simultaion
    #A=35.0;E0=10;alp=2.0;bnc=1.0;
    npcell=int(10000);dy=1.0;npnt='RefPnt';refDose30=3.0e-13;vol30=15.625
    tpsMU=0

    fig1 = pylab.figure(figsize=(11,8.5))
    #fig1.patch.set_facecolor('#99ccff')
    #thismanager=pylab.get_current_fig_manager()
    #thismanager.window.setGeometry(0,0,936,720)

    fig1.canvas.mpl_connect('scroll_event', scroll_function)
    ax=pylab.axes([0.3,0.16,0.97-0.3,0.95-0.16])
    #fig1.subplots_adjust(bottom=0.16,top=0.95,left=0.3,right=0.97)
    ax.set_title('Calculation')
    fig1.canvas.mpl_connect('button_press_event', myonclick)
    #place buttons
    axGetFile=pylab.axes([0.03,0.9,0.2,0.075])
    bGetFile=Button(axGetFile, 'Open File Directory')#,color='#99ccff',hovercolor='#3399ff')
    bGetFile.on_clicked(getFile)
    axQuit=pylab.axes([0.87,0.015,0.1,0.05])
    bQuit=Button(axQuit,'Quit')#,color='#9999ff',hovercolor='#ff0000')
    bQuit.on_clicked(quitter)
    axIA=pylab.axes([0.47,0.06,0.07,0.035])
    bIA=Button(axIA,'Inc A')#,color='#99ccff',hovercolor='#3399ff')
    bIA.on_clicked(incA)
    axDA=pylab.axes([0.47,0.015,0.07,0.035])
    bDA=Button(axDA,'Dec A')#,color='#99ccff',hovercolor='#3399ff')
    bDA.on_clicked(decA)
    axIE=pylab.axes([0.55,0.06,0.07,0.035])
    bIE=Button(axIE,'Inc E')#,color='#99ccff',hovercolor='#3399ff')
    bIE.on_clicked(incE)
    axDE=pylab.axes([0.55,0.015,0.07,0.035])
    bDE=Button(axDE,'Dec E')#,color='#99ccff',hovercolor='#3399ff')
    bDE.on_clicked(decE)
    axIal=pylab.axes([0.63,0.06,0.07,0.035])
    bIal=Button(axIal,'Inc al')#,color='#99ccff',hovercolor='#3399ff')
    bIal.on_clicked(incal)
    axDal=pylab.axes([0.63,0.015,0.07,0.035])
    bDal=Button(axDal,'Dec al')#,color='#99ccff',hovercolor='#3399ff')
    bDal.on_clicked(decal)
    axIbnc=pylab.axes([0.71,0.06,0.07,0.035])
    bIbnc=Button(axIbnc,'Inc bnc')#,color='#99ccff',hovercolor='#3399ff')
    bIbnc.on_clicked(incbnc)
    axDbnc=pylab.axes([0.71,0.015,0.07,0.035])
    bDbnc=Button(axDbnc,'Dec bnc')#,color='#99ccff',hovercolor='#3399ff')
    bDbnc.on_clicked(decbnc)
    axCalc=pylab.axes([0.795,0.015,0.07,0.035])
    bCalc=Button(axCalc,'Calc')#,color='#66ffcc',hovercolor='#00ff99')
    bCalc.on_clicked(calmain)
    axSave=pylab.axes([0.795,0.06,0.07,0.035])
    bSave=Button(axSave,'Save')#,color='#99ccff',hovercolor='#3399ff')
    bSave.on_clicked(saveimg)
    axReset=pylab.axes([0.87,0.075,0.1,0.05])
    bReset=Button(axReset,'Reset Zoom')#,color='#99ccff',hovercolor='#3399ff')
    bReset.on_clicked(resetLimits)
    #place radio button
    rax=pylab.axes([0.35,0.015,0.1,0.08])
    #rax.text(0.0,1.01,'Calculation Grid Size')
    rax.set_title('Calculation Grid Size',fontsize=11)
    radio=RadioButtons(rax,('1.0 mm','1.5 mm','2.5 mm'),active=0)
    radio.on_clicked(setCalcGrid)
    nax=pylab.axes([0.03,0.015,0.1,0.08])
    nax.set_title('Normalization',fontsize=11)
    nradio=RadioButtons(nax,('RefPoint','GTV'))
    nradio.on_clicked(setNormPoint)
    npax=pylab.axes([0.15,0.015,0.18,0.08])
    npax.set_title('Number of e per cell', fontsize=11)
    nprad=MyRadioButtons(npax,('30000','50000','10000','20000'),active=2)
    nprad.on_clicked(setNpcell)
    #place check boxes
    axContours=pylab.axes([0.03,0.47,0.2,0.39])
    axContours.set_title('Contours to plot',fontsize=11)
    axContours.get_xaxis().set_visible(False)
    axContours.get_yaxis().set_visible(False)
    check=CheckButtons(axContours,(''),(False))
    check.on_clicked(func)
    axBolus=pylab.axes([0.03,0.14,0.2,0.29])
    axBolus.set_title('Contours with override in TPS.',fontsize=11)
    axBolus.get_xaxis().set_visible(False)
    axBolus.get_yaxis().set_visible(False)
    bolus=CheckButtons(axBolus,(''),(False))
    bolus.on_clicked(bfunc)

    #axNpcell=pylab.axes([0.2,0.025,0.3,0.05])
    #axNpcell.set_title('Number of e per cell', fontsize=11)
    #slideNp=Slider(axNpcell,'',100,10000,valinit=500,dragging=True,valfmt='%d')
    #slideNp.on_changed(callback.setNpcell)
    pylab.show()
    
if __name__ == '__main__':

    main()
