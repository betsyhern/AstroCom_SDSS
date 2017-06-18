'''initialized 7/18/14'''

import os as os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.constants as const
import astropy.cosmology as cosmo
import matplotlib.pyplot as plt
import matplotlib as mpl
import long_catalog_sdss as cs

#import catalog from catalog_sdss file
cat=cs.cat

#reduce catalog sample to only view galaxies that have prominent bulges.
B_T=np.logical_and(cat['B_T_g']<=0.3,cat['B_T_r']<=0.3)

#B_T index reduces arrays and leaves those with disk 0.3 or less.
objID_B_T=cat['objID'][B_T]
z_pap1_B_T=cat['z_pap1'][B_T]
deg_B_T=cat['i_deg'][B_T]
cos_B_T=cat['i_cos'][B_T]
Ft_pap1_B_T=cat['Ft_pap1'][B_T]
ggMag_B_T=cat['ggMag'][B_T]
rgMag_B_T=cat['rgMag'][B_T]
gg2d_B_T=cat['gg2d'][B_T]
rg2d_B_T=cat['rg2d'][B_T]
Rhlg_B_T=cat['Rhlg'][B_T]
Rhlr_B_T=cat['Rhlr'][B_T]
B_T_g_B_T=cat['B_T_g'][B_T]
B_T_r_B_T=cat['B_T_r'][B_T]
z_pap2_B_T=cat['z_pap2'][B_T]
logM_B_T=cat['logM'][B_T]
b_logM_B_T=cat['b_logM'][B_T]
B_logM1_B_T=cat['B_logM1'][B_T]
Type_B_T=cat['Type'][B_T]
logMt_B_T=cat['logMt'][B_T]
b_logMt_B_T=cat['b_logMt'][B_T]
B_logMt1_B_T=cat['B_logMt1'][B_T]
logMb_B_T=cat['logMb'][B_T]
b_logMb_B_T=cat['b_logMb'][B_T]
B_logMb1_B_T=cat['B_logMb1'][B_T]
logMd_B_T=cat['logMd'][B_T]
b_logMd_B_T=cat['b_logMd'][B_T]
B_logMd1_B_T=cat['B_logMd1'][B_T]
Ft_pap2_B_T=cat['Ft_pap2'][B_T]
dBD_B_T=cat['dBD'][B_T]

#print catalog length and type; print B_T index length, sample, and array type; and objID_B_T length
debug=10
if debug==1:
	print 'catalog length'
	print len(cat)
	print 'catalog array type'
	print type(cat)
	print ' '
	print 'length of B_T index of B/T 0.3 or less'
	print len(B_T)
	print 'B_T index sample'
	print B_T[0:10]
	print 'B_T array type'
	print type(B_T)
	print ' '
	print 'objID_B_T length'
	print len(objID_B_T)
	print ' '
	print 'length of galaxies with cosine<0.30'
	a=cos_B_T<0.30
	b=cos_B_T[a]
	print len(b)
	print ' '
	print 'length of galaxies with cosine>0.83'
	c=cos_B_T>0.83
	d=cos_B_T[c]
	print len(d)
	print ' '
	print 'length of galaxies with 0.30<=cosine<=0.83'
	e=np.logical_and(cos_B_T>=0.30,cos_B_T<=0.83)
	n=len(cos_B_T[e])
	f=cos_B_T[e][0:n:4]
	print len(f)
	print ' '

#command for cut of edge-on galaxies with B/T 0-0.3
edge=cos_B_T<0.30

#command of cut for face-on galaxies with B/T 0-0.3
face=cos_B_T>0.83

#command of cut for galaxies removing edge-on and face-on samples with B/T 0-0.3
middle=np.logical_and(cos_B_T>=0.30,cos_B_T<=0.83)

#command to create the g-r magnitude array
g_rMag=cat['gg2d']-cat['rgMag']
g_rMag_B_T=ggMag_B_T-rgMag_B_T

#Histogram of inclination data set after using cosine for 2011 catalog
def old_hist(data):
	plt.hist(data, bins=10, range=(data.min(),data.max()), color='g')
	plt.axis([0,1.1,0,180000],fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	plt.title('Inclination cosine angle frequecy 2011 catalog', y=1.05,fontsize=17, fontweight='bold')
	plt.show()
	
old_hist(cs.icos_use)

#Histogram of inclination data set after using cosine with B/T 0.3 or less.
#A is number of bins, B is min range, C is max range
#angle=0: min thru max, angle=1: min  thru 0.3, degree=2: 0.83 thru max
#cosine min-0.30 has 42551 galaxies and cosine 0.83-max has 42193 galaxies 0.30<=cosine<=0.83 has 40484 galaxies with step 4
def angle_hist(data, B_T, A1=10, A2=5, A3=84, angle=0, title=None, xlabel=None, axis=None):
	B=data.min()
	C=data.max()
	if angle==0: #cosine
		plt.hist(data, bins=A1, range=(B,C), color='g')	
	if angle==1: #cosine
		plt.hist(data, bins=A2, range=(B,0.30), color='g')
	if angle==2: #cosine
		plt.hist(data, bins=A2, range=(0.83,C), color='g')
	if angle==3: #cosine
		plt.hist(data, bins=A2, range=(0.30,0.83), color='g')
	if angle==4: #degrees
		plt.hist(data, bins=A3, range=(B,C), color='g')
	if angle==5: #degrees
		plt.hist(data, bins=A3, range=(B,84.0), color='g')
	if angle==6: #degrees
		plt.hist(data, bins=A1, range=(84.0,C), color='g')
	if angle==7:
		plt.hist(data, bins=A1, range=(B,C), color='y', alpha=0.5, label='all our catalog galaxies')
		plt.hist(data[B_T], bins=A1, range=(B,C), color='c', label='our catalog galaxies with a\nsmall bulge and large disk')
	if angle==8:
		plt.hist(data, bins=A1, range=(B,C), color='y', alpha=0.5, label='all catalog galaxies', normed=1)
		plt.hist(data[B_T], bins=A1, range=(B,C), color='g', alpha=0.5, label='galaxies B/T<0.30', normed=1)
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel(xlabel,fontsize=16, fontweight='bold')
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.legend(loc=0)
	plt.show()

#angle_hist(cat['i_deg'],B_T,angle=0,title='Inclination Angle Frequency in Degrees Min to Max', xlabel='Inclination Angle (degrees)',axis=[-5,90,0,110000])
#angle_hist(cat['i_cos'],B_T,angle=0,title='Inclination Angle Frequency of Cosine Angle Min to Max', xlabel='Inclination Cosine Angle', axis=[0,1.1,0,100000])
#angle_hist(deg_B_T,B_T,angle=4,title='Inclination Angle Frequency in Degrees Min to Max with B/T 0-0.3', xlabel='Inclination Angle (degrees)')
#angle_hist(cos_B_T,B_T,angle=0,title='Inclination Angle Frequency of Cosine\nAngle Min to Max with B/T 0-0.3', xlabel='Inclination Cosine Angle',axis=[0,1.1,0,32000])
#angle_hist(cat['i_cos'],B_T,angle=7,title='Inclination of Cosine Angle Frequency of All\nGalaxies in Catalog and Galaxies with B/T<0.3', xlabel='Inclination Cosine Angle',axis=[0,1.1,0,100000])
#angle_hist(cat['i_cos'],B_T,angle=8,title='Normalized Inclination of Cosine Angle Frequency of\nAll Galaxies in Catalog and Galaxies with B/T<0.3', xlabel='Inclination Cosine Angle',axis=[0,1.1,0,1.8])

#histogram of stellar masses of disc 0-0.3 cut, edge-on(cosine<0.3), and face-on(cosine>0.83)
def mass_hist(data, edge, face, middle, A=50, title=None, xlabel=None):
	n=len(data[middle])
	B=data.min()
	C=data.max()
	plt.hist(data[face],bins=A, range=(B,C),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	plt.hist(data[edge], bins=A, range=(B,C), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	plt.hist(data, bins=A, range=(B,C), color='c', histtype='stepfilled', alpha=0.7, label='all our galaxies with a\nsmall bulge and large disk', normed=1)
	#plt.hist(data[middle][0:n:4], bins=A, range=(B,C), color='c',  histtype='stepfilled', alpha=0.7, label='reduced middle')
	plt.axvline(data[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.axvline(data.mean(), color='k', linestyle='-', label='mean all our galaxies with a\nsmall bulge and large disk: %.2f'% np.mean(data), linewidth=2)
	#plt.axvline(data[middle][0:n:4].mean(), color='k', linestyle='-', linewidth=2, label='mean reduced middle %.2f'% np.mean(data[middle][0:n:4]))
	plt.legend(loc=0)
	#plt.axis([6,13,0,0.85],fontsize=30, fontweight='bold')
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel(xlabel,fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.show()

#mass_hist(logM_B_T,edge,face,middle,title='Normalized Frequency of logM (solar mass) with B/T<0.3',xlabel='logM (solar mass)')
#mass_hist(logMt_B_T,edge,face,middle,title='Normalized Frequency Sum (logM bulge + logM disk)\nin solar masses with B/T<0.3',xlabel='logM bulge + logM disk (solar mass)')
#mass_hist(logMb_B_T,edge,face,middle,title='Normalized Frequency of logM bulge (solar mass)\nwith B/T<0.3',xlabel='logM bulge (solar mass)')
#mass_hist(logMd_B_T,edge,face,middle,title='Normalized Frequency of logM disk (solar mass)\nwith B/T<0.3',xlabel='logM disk (solar mass)')

#histogram of g-r magnitudes
def mag_hist(data, edge, face, A=25, title=None, xlabel=None, axis=None):
	plt.hist(data[face], bins=A, range=(0,1), color='b', histtype='stepfilled', alpha=0.7, label='face-on',normed=1)
	plt.hist(data[edge], bins=A, range=(0,1), color='r', histtype='stepfilled', alpha=0.8, label='edge-on',normed=1)
	plt.hist(data, bins=A, range=(0,1), color='c', histtype='stepfilled', alpha=0.7, label='all our galaxies\nwith small bulge\nand large disk',normed=1)
	plt.axvline(data[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.axvline(data.mean(), color='k', linestyle='-', label='mean all our galaxies\nwith small bulge\nand large disk: %.2f'% np.mean(data), linewidth=2)
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel(xlabel,fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.legend(loc=2)
	plt.show()

#mag_hist(g_rMag_B_T,edge, face,title='Normalized Frequency of green-band Magnitudes\nMinus red-band Magnitudes with B/T<0.3',xlabel='g-r (color)',axis=[-0.2,1.1,0,4.5])

#UGRC
def cos_hist_old_new(icos_use,use_sort,A=10):
	B=icos_use.min()
	C=icos_use.max()
	plt.hist(icos_use, bins=A, range=(B, C),color='g', normed=1, label='2011 catalog')
	plt.hist(icos_use[use_sort],bins=A, range=(B, C),color='y', alpha=0.5, normed=1, label='our catalog')
	plt.title('Normalized Inclination of\n2011 Catalog and Our Catalog',fontsize=17, fontweight='bold')
	plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle', fontsize=16, fontweight='bold')
	plt.axis([0,1.1,0,1.8],fontsize=30, fontweight='bold')
	plt.legend(loc=0)
	plt.show()

#cos_hist_old_new(cs.icos_use,cs.use_sort)
#angle_hist(cat['i_cos'],B_T,angle=7,title='Inclination Frequency of Galaxies', xlabel='Inclination Cosine Angle',axis=[0,1.1,0,100000])
#mag_hist(g_rMag_B_T,edge, face,title='Normalized Frequency of Green Luminosity\nminus Red Luminosity for Galaxies',xlabel='G-R (color)',axis=[-0.2,1.1,0,4.5])
#mass_hist(logM_B_T,edge,face,middle,title='Normalized Frequency of logM (solar mass) Galaxies',xlabel='logM (solar mass)')
