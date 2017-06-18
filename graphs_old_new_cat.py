'''Initialized 6/11/15'''

import os as os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.constants as const
import astropy.cosmology as cosmo
import matplotlib.pyplot as plt
from decimal import *
import long_catalog_sdss as lcs
import catalog_sdss_condensed as csc


#reduce catalog sample to only view galaxies that have prominent bulges.
B_To=np.logical_and(lcs.cat['B_T_g']<=0.3,lcs.cat['B_T_r']<=0.3)
B_Tn=np.logical_and(csc.cat['B/T_g']<=0.3,csc.cat['B/T_r']<=0.3)
B_Toldcut=np.logical_and(np.logical_and(lcs.cat['B_T_g']>0.,lcs.cat['B_T_g']<=0.3),np.logical_and(lcs.cat['B_T_r']!=0.,lcs.cat['B_T_r']<=0.3))
B_Tnewcut=np.logical_and(np.logical_and(csc.cat['B/T_g']>0.,csc.cat['B/T_g']<=0.3),np.logical_and(csc.cat['B/T_r']!=0.,csc.cat['B/T_r']<=0.3))

#command for cut of edge-on galaxies with B/T 0-0.3
edge_old=lcs.cat['i_cos'][B_To]<0.30
edge_new=csc.cat['i'][B_Tn]<0.30
edge_oldcut=lcs.cat['i_cos'][B_Toldcut]<0.30
edge_newcut=csc.cat['i'][B_Tnewcut]<0.30

print ' edge old, edge new, edge oldcut, edge newcut'
print len(edge_old)
print len(edge_new)
print len(edge_oldcut)
print len(edge_newcut)
print ''

#command of cut for face-on galaxies with B/T 0-0.3
face_old=lcs.cat['i_cos'][B_To]>0.83
face_new=csc.cat['i'][B_Tn]>0.83
face_oldcut=lcs.cat['i_cos'][B_Toldcut]>0.84
face_newcut=csc.cat['i'][B_Tnewcut]>0.84

print ' face old, face new, face oldcut, face newcut'
print len(face_old)
print len(face_new)
print len(face_oldcut)
print len(face_newcut)
print ''
print 'old dust'
print len(lcs.cat['logM'][edge_old])
print len(lcs.cat['logM'][face_old])
print ''
print 'new dust'
print len(csc.cat['logM_dt'][edge_new])
print len(csc.cat['logM_dt'][face_new])
print ''
print 'old dust cut'
print len(lcs.cat['logM'][edge_oldcut])
print len(lcs.cat['logM'][face_oldcut])
print ''
print 'new dust cut'
print len(csc.cat['logM_dt'][edge_newcut])
print len(csc.cat['logM_dt'][face_newcut])
print ''
print 'new dust free'
print len(csc.cat['logM_df'][edge_new])
print len(csc.cat['logM_df'][face_new])
print ''
print 'new dust free cut'
print len(csc.cat['logM_df'][edge_newcut])
print len(csc.cat['logM_df'][face_newcut])
print ''

#command of cut for galaxies removing edge-on and face-on samples with B/T 0-0.3
middle_old=np.logical_and(lcs.cat['i_cos'][B_To]>=0.30,lcs.cat['i_cos'][B_To]<=0.83)
middle_new=np.logical_and(csc.cat['i'][B_Tn]>=0.30,csc.cat['i'][B_Tn]<=0.83)
middle_oldcut=np.logical_and(lcs.cat['i_cos'][B_Toldcut]>=0.30,lcs.cat['i_cos'][B_Toldcut]<=0.84)
middle_newcut=np.logical_and(csc.cat['i'][B_Tnewcut]>=0.30,csc.cat['i'][B_Tnewcut]<=0.84)

print len(lcs.cat['ggMag'])
print len(lcs.cat['rgMag'])
print len(B_To)

#command to create the g-r magnitude array of B/T 0-0.3
g_r_old=np.subtract(lcs.cat['ggMag'],lcs.cat['rgMag'])[B_To]
g_r_new=np.subtract(csc.cat['ggMag'],csc.cat['rgMag'])[B_Tn]
g_r_oldcut=np.subtract(lcs.cat['ggMag'],lcs.cat['rgMag'])[B_Toldcut]
g_r_newcut=np.subtract(csc.cat['ggMag'],csc.cat['rgMag'])[B_Tnewcut]

#Histogram of inclination data set after using cosine for 2011 catalog
def old_hist(data,title=None):
	plt.hist(data, bins=10, range=(data.min(),data.max()), color='g')
	plt.axis([0,1.1,0,180000],fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.show()

#old_hist(lcs.icos_use, title='Inclination cosine angle frequecy 2011 catalog\nOLD CATALOG')	
#old_hist(csc.cut1['i'], title='Inclination cosine angle frequecy 2011 catalog\nNEW CATALOG')

#Histogram of our catalog cosine inclinations 
def cos_angle_hist(data,color=None, title=None):
	plt.hist(data, bins=10, range=(data.min(),data.max()), color=color, alpha=0.6)	
	plt.axis([0,1.1,0,100000],fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.show()

#cos_angle_hist(lcs.cat['i_cos'], color= 'y', title='Inclination Angle Frequency of Our\nCatalog Cosine Angle Min to Max\nOLD CATALOG')
#cos_angle_hist(csc.cat['i'], color= 'm',title='Inclination Angle Frequency of Our\nCatalog Cosine Angle Min to Max\nNEW CATALOG')

#Histogram of our catalog cosine inclinations of old and new catalogs together 
def old_new_angle_hist(old, new):
	plt.hist(old, bins=10, range=(old.min(),old.max()), color='y',alpha=0.6,label='old catalog')
	plt.hist(new, bins=10, range=(new.min(),new.max()), color='m',alpha=0.6,label='new catalog')	
	plt.axis([0,1.1,0,100000],fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	plt.title('Inclination Angle Frequency of Our\nCatalog Cosine Angle Min to Max\nOLD and NEW CATALOG', y=1.05,fontsize=17, fontweight='bold')
	plt.legend(loc=0)
	plt.show()

#old_new_angle_hist(lcs.cat['i_cos'],csc.cat['i'])

#Normalized histogram of 2011 and Our Catalog
def cos_hist_old_new(orig,cat,color1=None, color2=None, title=None):
	plt.hist(orig, bins=10, range=(orig.min(), orig.max()),color=color1, normed=1, label='2011 catalog')
	plt.hist(cat,bins=10, range=(cat.min(), cat.max()), color=color2, alpha=0.6, normed=1, label='our catalog')
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Angle', fontsize=16, fontweight='bold')
	plt.axis([0,1.1,0,1.8],fontsize=30, fontweight='bold')
	plt.legend(loc=0)
	plt.show()

#cos_hist_old_new(lcs.icos_use,lcs.cat['i_cos'], color1='g', color2='y', title='Normalized Inclination Frequency of\n2011 Catalog and Our Catalog\nOLD CATALOG')
#cos_hist_old_new(csc.cut1['i'],csc.cat['i'], color1='g', color2='m', title='Normalized Inclination Frequency of\n2011 Catalog and Our Catalog\nNEW CATALOG')

# frequency of inclination cosine angles
def cat_hist(all,cut,color1=None, color2=None,title=None):
	plt.hist(all, bins=10, range=(all.min(),all.max()), color=color1, alpha=0.6, label='Our catalog galaxies')
	plt.hist(all[cut], bins=10, range=(all.min(),all.max()), color=color2, alpha=0.6, label='Our catalog galaxies with a\nsmall bulge and large disk')
	plt.axis([0,1.1,0,100000],fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.legend(loc=0)
	plt.show()

#cat_hist(lcs.cat['i_cos'],B_To, color1='y', color2='c', title='Inclination of Cosine Angle Frequency of All\nGalaxies in Our Catalog and Galaxies with B/T<0.3\n in Our OLD CATALOG')
#cat_hist(csc.cat['i'],B_Tn, color1='m', color2='b', title='Inclination of Cosine Angle Frequency of All\nGalaxies in Our Catalog and Galaxies with B/T<0.3\nin Our NEW CATALOG')

#histogram of g-r magnitudes
def mag_hist(data, edge, face, A=25, title=None):
	plt.hist(data[face], bins=A, range=(0,1), color='b', histtype='stepfilled', alpha=0.7, label='face-on',normed=1)
	plt.hist(data[edge], bins=A, range=(0,1), color='r', histtype='stepfilled', alpha=0.8, label='edge-on',normed=1)
	plt.hist(data, bins=A, range=(0,1), color='c', histtype='stepfilled', alpha=0.7, label='all our galaxies\nwith small bulge\nand large disk',normed=1)
	plt.axvline(data[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.axvline(data.mean(), color='k', linestyle='-', label='mean all our galaxies\nwith small bulge\nand large disk: %.2f'% np.mean(data), linewidth=2)
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel('g-r (color)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.axis([-0.2,1.1,0,4.5],fontsize=30, fontweight='bold')
	plt.legend(loc=2)
	plt.show()

#mag_hist(g_r_old,edge_old,face_old,title='Normalized Frequency of green-band Magnitudes\nMinus red-band Magnitudes with B/T<0.3\nOLD CATALOG')
#mag_hist(g_r_new,edge_new,face_new,title='Normalized Frequency of green-band Magnitudes\nMinus red-band Magnitudes with B/T<0.3\nNEW CATALOG')
#mag_hist(g_r_oldcut,edge_oldcut,face_oldcut,title='Normalized Frequency of green-band Magnitudes\nMinus red-band Magnitudes with B/T<0.3 and not =0\nOLD CATALOG')
#mag_hist(g_r_newcut,edge_newcut,face_newcut,title='Normalized Frequency of green-band Magnitudes\nMinus red-band Magnitudes with B/T<0.3 and not =0\nNEW CATALOG')

#histogram of stellar masses of disc 0-0.3 cut, edge-on(cosine<0.3), and face-on(cosine>0.83)
def mass_hist(data, edge, face, middle, A=50, title=None, xlabel=None):
	n=len(data[middle])
	B=data.min()
	C=data.max()
	plt.hist(data[face],bins=A, range=(B,C),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	plt.hist(data[edge], bins=A, range=(B,C), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	plt.hist(data, bins=A, range=(B,C), color='c', histtype='stepfilled', alpha=0.7, label='all our galaxies with a\nsmall bulge and large disk', normed=1)
	plt.axvline(data[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.axvline(data.mean(), color='k', linestyle='-', label='mean all our galaxies with a\nsmall bulge and large disk: %.2f'% np.mean(data), linewidth=2)
	plt.legend(loc=0)
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel('logM (solar mass)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.show()

#mass_hist(lcs.cat['logM'][B_To],edge_old,face_old,middle_old,title='Normalized Frequency of logM (solar mass) with B/T<0.3\nOLD CATALOG WITH DUST')
#mass_hist(csc.cat['logM_dt'][B_Tn],edge_new,face_new,middle_new,title='Normalized Frequency of logM (solar mass) with B/T<0.3\nNEW CATALOG WITH DUST')
#mass_hist(csc.cat['logM_df'][B_Tn],edge_new,face_new,middle_new,title='Normalized Frequency of logM (solar mass) with B/T<0.3\nNEW CATALOG DUST FREE')
#mass_hist(lcs.cat['logM'][B_Toldcut],edge_oldcut,face_oldcut,middle_oldcut,title='Normalized Frequency of logM (solar mass) with\nB/T<0.3 and not=0 OLD CATALOG WITH DUST')
#mass_hist(csc.cat['logM_dt'][B_Tnewcut],edge_newcut,face_newcut,middle_newcut,title='Normalized Frequency of logM (solar mass) with\nB/T<0.3 and not=0 NEW CATALOG WITH DUST')
#mass_hist(csc.cat['logM_df'][B_Tnewcut],edge_newcut,face_newcut,middle_newcut,title='Normalized Frequency of logM (solar mass) with\nB/T<0.3 and not=0 NEW CATALOG DUST FREE')
