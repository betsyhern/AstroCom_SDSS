'''Initialized 7/13/15'''

import os as os
import numpy as np
import astropy.io as aio
from astropy.table import Table,join,Column
import astropy.constants as aconst
import astropy.cosmology as acosmo
import matplotlib.pyplot as plt
#import spherematch as sp
 
fnS1='Simard_Table3.fit' #1D fits 
fnS2='Simard_Table1.fit' #2D fits with n_b = 4
fnM ='Mendel_dusty.fit'  #dusty masses combines Tables 3 and 4

def create_S11_catalog():
    hdulist1=aio.fits.open(fnS1)
    hdulist2=aio.fits.open(fnS2)
    Simard_table1=hdulist1[1].data
    Simard_table3=hdulist2[1].data
    shared=['objID','z','Sp','Scale']
    t=join(Simard_table3,Simard_table1,keys=shared)
    #remove galaxies with no absolute magnitudes
    t=t[np.isfinite(t['ggMag_1'])]
    w=np.where(t['objID']==587741532231565557)
    t.remove_row(int(w[0]))
    #remove galaxies which are pure bulge in g and pure disk in r
    w=np.where(t['objID']==587738410589552933)
    t.remove_row(int(w[0]))
    w=np.where(t['objID']==587739811572809954)
    t.remove_row(int(w[0]))
    col_cosi=Column(name='cosi',data=np.cos(np.radians(t['i'])))
    t.add_column(col_cosi)
    return t

def create_bestfit():
    #this creates a best fit catalog where each galaxy is given a 2D or
    #1D fit based on the F-test and other statistics. If 1D is a better fit
    #then the values for disk or bulge are replaced with the 1D fit and the
    #other values are set to zero
    S11=create_S11_catalog()
    
def create_mass_catalog():
    S_table=create_S11_catalog()
    hdulist=aio.fits.open(fnM)
    M_table=hdulist[1].data
    shared=['objID']
    t=join(M_table,S_table,keys=shared)
    t.remove_column('Sp') #they are all 2
    t.rename_column('z_1','z')
    t.remove_column('z_2') #they should be the same
    t.rename_column('PpS_2','PpS') # these should be the same but are slightly
    t.remove_column('PpS_1')       # different
    t.rename_column('__B_T_g','B/T_g')
    t.rename_column('__B_T_r','B/T_r')
    #remove error ranges on masses
    t.remove_column('b_logM')
    t.remove_column('B_logM1')
    t.remove_column('b_logMt')
    t.remove_column('B_logMt1')
    t.remove_column('b_logMb')
    t.remove_column('B_logMb1')
    t.remove_column('b_logMd')
    t.remove_column('B_logMd1')
    #replace values where 2D fit failed with 1D values
    type1=(t['Type']==1)
    type2=(t['Type']==2)
    pure_bulge=np.logical_or(t['B/T_g']==1,t['B/T_r']==1)
    pure_bulge=np.logical_or(type1,pure_bulge)
    pure_bulge=np.logical_and(~type2,pure_bulge)
    t['logMd'][pure_bulge]=0.0
    t['logMb'][pure_bulge]=t['logM'][pure_bulge]
    t['logMt'][pure_bulge]=t['logM'][pure_bulge]
    t['B/T_g'][pure_bulge]=1.0
    t['B/T_r'][pure_bulge]=1.0
    t['i'][pure_bulge]=-90.
    t['cosi'][pure_bulge]=-1.0
    t['gdMag'][pure_bulge]=-99.9 #luminosity =0
    t['rdMag'][pure_bulge]=-99.9 #luminosity =0
    pure_disk=np.logical_or(t['B/T_g']==0,t['B/T_r']==0)
    pure_disk=np.logical_or(type2,pure_disk)
    pure_disk=np.logical_and(~type1,pure_disk)
    t['logMb'][pure_disk]=0.0
    t['logMd'][pure_disk]=t['logM'][pure_disk]
    t['logMt'][pure_disk]=t['logM'][pure_disk]
    t['B/T_g'][pure_disk]=0.0
    t['B/T_r'][pure_disk]=0.0
    t['gbMag'][pure_disk]=-99.9
    t['rbMag'][pure_disk]=-99.9
    return t

def create_disk_catalog():
    t=create_mass_catalog()
    disks=np.logical_and(t['B/T_r'] < 0.35,t['PpS']<0.32) #PpS is low then you need disk and bulge, otherwise if high you only need the disk or the bulge
    return t[disks]

def create_full_catalog():
    S11_table=create_S11_catalog()
    hdulist=aio.fits.open(fnM)
    M_table=hdulist[1].data
    shared=['objID','z']
    t=join(S11_table,M_table,keys=shared,join_type='left')
    return t

def weird_stuff(t):
    #This should only be run on the full catalog as read from fits file
    print "There are some cases with m and z but no M"
    w=np.where(np.logical_or(np.isfinite(t['ggMag_1'])==False,
                             np.isfinite(t['ggMag_2'])==False))
    print 'objID','z','gg2d_1D','ggMag_1D','rg2d_1D','rgMag_1D','logM','logMt'
    for i in w[0]:
        print t['objID'][i],t['z'][i],t['gg2d_1'][i],t['ggMag_1'][i] \
            ,t['rg2d_1'][i],t['rgMag_1'][i],t['logM'][i],t['logMt'][i] \
            ,t['gg2d_2'][i],t['ggMag_2'][i],t['rgMag_2'][i]
    print "There are cases that are 100% in one band and 100% disk in the other"
    w1=np.where(np.logical_and(t['__B_T_r']==0.0,t['__B_T_g']==1.0))
    w2=np.where(np.logical_and(t['__B_T_r']==1.0,t['__B_T_g']==0.0))

def check_column(t):
#this routine goes through the columns and checks the range of each column
#and that their values are finite
    names=t.colnames
    for name in names:
        max=np.amax(t[name])
        min=np.amin(t[name])
        fin=len(t)-np.sum(np.isfinite(t[name]))
        print name+" range of {0} to {1} and {2} NaN".format(max,min,fin)

def tplot(t,name1,name2):
    plt.plot(t[name1],t[name2],'.')
    plt.xlabel('$'+name1+'$')
    plt.ylabel('$'+name2+'$')
    plt.show()

def iplot(t):
    plt.figure()
    cosi1=t['cosi'][t['Type']==2]
    g3=np.logical_and(t['Type']==3,t['PpS'] < 0.32)
    cosi2=t['cosi'][g3]
    cosi3=t['cosi'][np.logical_and(g3,t['B/T_r'] < 0.35)]
    cosi4=t['cosi'][np.logical_and(t['B/T_r'] < 0.35,t['PpS'] < 0.32)]
    names=['Type2','Type3','Good 3','Disks']
    n,bins,patches=plt.hist([cosi1,cosi2,cosi3,cosi4],18,histtype='step'
                            ,label=names,range=[0.1,1.0])
    plt.xlabel('b/a')
    plt.legend(loc=2)
    plt.show()

def colorplot(t):
    color1=t['ggMag_1']-t['rgMag_1']
    color2=t['ggMag_2']-t['rgMag_2']
    edgeon=t['cosi'] < 0.3
    faceon=t['cosi'] > 0.85
    plt.figure()
    n,b,p=plt.hist([color2[faceon],color2[edgeon]],30,histtype='step'
         ,range=[-0.5,1.5],color=['blue','red'],label=['face-on','edge-on'])
    plt.xlabel('g-r')
    plt.legend()
    plt.show()
    
def check_ds(ds):
    names=ds.field_list()
    for name in names:
        max=np.amax(t[name])
        min=np.amin(t[name])
        fin=len(t)-np.sum(np.isfinite(t[name]))
        print name+" range of {0} to {1} and {2} NaN".format(max,min,fin)
