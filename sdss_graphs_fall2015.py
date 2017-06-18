'''Initialized 12/30/15'''
import os as os
import numpy as np
import astropy.io as aio
from astropy.table import Table,join,Column,unique
import astropy.constants as aconst
import astropy.cosmology as acosmo
import matplotlib.pyplot as plt
from sdss import *

'''
cat=open_flux_mass_catalog()
#cat=open_zoo_mass_line_catalog()
disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
'''
#Histogram of our catalog cosine inclinations of old and new catalogs together 
def angle_hist(data,cut):
	plt.hist(data, bins=6, range=(0,1), color='y',alpha=0.6,label='Excluding Schawinski Full Catalog')
	plt.hist(data[cut], bins=6, range=(0,1), color='m',alpha=0.6,label='Excluding Schawinski Disk Sample')	
	plt.axis([-0.1,1.1,0,110000],fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	plt.title('Inclination Cosine Angles', y=1.05,fontsize=17, fontweight='bold')
	plt.legend(loc=0)
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.12, .5)
	plt.show()

#angle_hist(cat['cosi'],disk)

'''
cat=open_flux_mass_catalog()
disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
disk_cat=cat[disk]
HaHb=disk_cat['H_alpha']/disk_cat['H_beta']
cut=HaHb<12
HaHb12=HaHb[cut]
cutHa_gr=np.logical_and(disk_cat['HA_SN']>3,cut)
cutHb_gr=np.logical_and(disk_cat['HB_SN']>3,cut)
cutHa_le=np.logical_and(disk_cat['HA_SN']<3,cut)
cutHb_le=np.logical_and(disk_cat['HB_SN']<3,cut)
HaHb_Ha_gr=HaHb[cutHa_gr]
HaHb_Hb_gr=HaHb[cutHb_gr]
HaHb_Ha_le=HaHb[cutHa_le]
HaHb_Hb_le=HaHb[cutHb_le]
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.86
edge_Ha_gr=np.logical_and(cutHa_gr,edge)
face_Ha_gr=np.logical_and(cutHa_gr,face)
edge_Hb_gr=np.logical_and(cutHb_gr,edge)
face_Hb_gr=np.logical_and(cutHb_gr,face)
edge_Ha_le=np.logical_and(cutHa_le,edge)
face_Ha_le=np.logical_and(cutHa_le,face)
edge_Hb_le=np.logical_and(cutHb_le,edge)
face_Hb_le=np.logical_and(cutHb_le,face)

print 'cat length: '+str(len(cat))
print 'disk cat length: '+str(len(disk_cat))
print 'HaHb length: '+str(len(HaHb))
print 'HaHb < 12 length: '+str(len(HaHb12))
print 'edge-on length: '+str(len(cat[edge]))
print 'face-on length: '+str(len(cat[face]))
print ''
print 'Halpha SN>3'
print 'HaHb length: '+str(len(HaHb_Ha_gr))
print 'edge-on length: '+str(len(cat[edge_Ha_gr]))
print 'face-on length: '+str(len(cat[face_Ha_gr]))
print ''
print 'Halpha SN<3'
print 'HaHb length: '+str(len(HaHb_Ha_le))
print 'edge-on length: '+str(len(cat[edge_Ha_le]))
print 'face-on length: '+str(len(cat[face_Ha_le]))
print ''
print 'Hbeta SN>3'
print 'HaHb length: '+str(len(HaHb_Hb_gr))
print 'edge-on length: '+str(len(cat[edge_Hb_gr]))
print 'face-on length: '+str(len(cat[face_Hb_gr]))
print ''
print 'Hbeta SN<3'
print 'HaHb length: '+str(len(HaHb_Hb_le))
print 'edge-on length: '+str(len(cat[edge_Hb_le]))
print 'face-on length: '+str(len(cat[face_Hb_le]))
'''
# creates flux histogram from Brinchmann, Simard and Mendel
def flux_hist(data, edge, face, A=20,axis=None, title=None):
	plt.hist(np.log10(data[edge]), bins=A, range=(0,2), align='right', color='r',histtype='stepfilled', alpha=0.8,label='edge-on')
	plt.hist(np.log10(data[face]), bins=A, range=(0,2), align='right', color='b',histtype='stepfilled', alpha=0.7,label='face-on')
	plt.axvline(np.log10(data[face]).mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.3f'% np.mean(np.log10(data[face])))
	plt.axvline(np.log10(data[edge]).mean(), color='k', linestyle='-', linewidth=2, label='mean edge-on: %.3f'% np.mean(np.log10(data[edge])))
	plt.legend(loc=2)
	plt.axis(axis,fontsize=40, fontweight='bold')
	plt.title(title, y=1.05,fontsize=18, fontweight='bold')
	plt.xlabel(r'log H$\alpha$/H$\beta$ (erg/s/cm^2)', fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.1)
	plt.gca().yaxis.set_label_coords(-.1, .5)
	plt.show()

#flux_hist(HaHb, edge, face, axis=[0,1.2,0,5000], title=r'H$\alpha$/H$\beta$ Flux Ratio < 12')
#flux_hist(HaHb, edge_Ha_gr, face_Ha_gr,axis=[0,1.2,0,5000], title=r'H$\alpha$/H$\beta$ Flux Ratio < 12 H$\alpha$ SN > 3')
#flux_hist(HaHb, edge_Ha_le, face_Ha_le,axis=[0,1.2,0,80], title=r'H$\alpha$/H$\beta$ Flux Ratio < 12 H$\alpha$ SN < 3')

# creates flux ratio vs inclination scatter plot from Brinchmann, Simard and Mendel
def flux_scatter(X, Y, axis=None, title=None):
	plt.scatter(X, Y, color='g',marker='.', label=r'H$\alpha$/H$\beta$ ratio')
	coeffs=np.polyfit(X,Y,1)
	trendline=np.poly1d(coeffs)
	correlation=np.corrcoef(X,Y)[0,1]
	correlationSq=correlation**2
	results={'polynomial':coeffs.tolist(),'correlation':correlation,'determination':correlationSq}
	plt.plot(X,trendline(X),'k-',label='trendline equation'+str(trendline)+'\n$R^2$= '+str(round(results['determination'],3)))
	plt.axis(axis,fontsize=40, fontweight='bold')
	plt.title(title, y=1.05,fontsize=18, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle', fontsize=16, fontweight='bold')
	plt.ylabel(r'H$\alpha$/H$\beta$ (erg/s/cm^2)',fontsize=16, fontweight='bold')
	plt.legend(loc=0)
	plt.gca().xaxis.set_label_coords(.5, -.06)
	plt.gca().yaxis.set_label_coords(-.04, .5)
	plt.show()

#flux_scatter(disk_cat['cosi'][cut],HaHb12, axis=[0,1.05,0,18], title=r'H$\alpha$/H$\beta$ Flux Ratio < 12 versus Inclination')
#flux_scatter(disk_cat['cosi'][cutHa_gr],HaHb_Ha_gr, axis=[0,1.05,0,18], title=r'H$\alpha$/H$\beta$ Flux Ratio < 12\n'+r'versus Inclination H$\alpha$ SN > 3')   
#flux_scatter(disk_cat['cosi'][cutHa_le],HaHb_Ha_le, axis=[0,1.05,0,18], title=r'H$\alpha$/H$\beta$ Flux Ratio < 12\n'+r'versus Inclination H$\alpha$ SN < 3')    
 
'''
cat=open_zoo_mass_line_catalog()
cut=cat['L_O3']>0 #zoo lum
#cut=np.logical_and(cat['L_O3']>0,cat['OIII_tot']>0) #Brinchmann zoo
cat=cat[cut]
#check_column(cat)
disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
disk_cat=cat[disk]
print 'zoo disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.27
face=disk_cat['cosi']>0.80
print 'zoo edge-on length: '+str(len(cat[edge]))
print 'zoo face-on length: '+str(len(cat[face]))
'''
'''
cat=open_line_mass_catalog()
cut=cat['OIII_tot']>0
cat=cat[cut]
disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
disk_cat=cat[disk]
print 'line mass disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.87
print 'line mass edge-on length: '+str(len(cat[edge]))
print 'line mass face-on length: '+str(len(cat[face]))
'''
'''
cat=open_line_mass_catalog()
cut=cat['H_alpha']>0
cat=cat[cut]
#check_column(cat)
disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
disk_cat=cat[disk]
print 'line mass disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.87
print 'line mass edge-on length: '+str(len(cat[edge]))
print 'line mass face-on length: '+str(len(cat[face]))
'''
# creates luminosity histogram from zoo
def lumlog_hist(data, edge, face, A=20,title=None, axis=None, xlabel=None):
	plt.rc('text', usetex=False)
	plt.hist(np.log10(data[edge]), bins=A, align='right', color='r',histtype='stepfilled', alpha=0.8,label='edge-on')
	plt.hist(np.log10(data[face]), bins=A, align='right', color='b',histtype='stepfilled', alpha=0.7,label='face-on')
	plt.axvline(np.log10(data[face]).mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(np.log10(data[face])))
	plt.axvline(np.log10(data[edge]).mean(), color='k', linestyle='-', linewidth=2, label='mean edge-on: %.2f'% np.mean(np.log10(data[edge])))
	plt.title(title, y=1.05,fontsize=18, fontweight='bold')
	plt.xlabel(xlabel, fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
	plt.axis(axis,fontsize=40, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.1)
	plt.gca().yaxis.set_label_coords(-.1, .5)
	plt.legend(loc=1)
	plt.show()

#lumlog_hist(disk_cat['L_O3'], edge, face, A=55,title='OIII Luminosity',xlabel='log$_{10}$ OIII luminosity (ergs/s)',axis=[37.5,42.5,0,900])
#lumlog_hist(disk_cat['OIII_tot'], edge, face, A=30,title='Brinchmann OIII Flux',xlabel='log$_{10}$ OIII Flux (ergs/s/cm$^2$)', axis=[0,2.5,0,5000])
#lumlog_hist(disk_cat['OIII_tot'], edge, face, A=18,title='Brinchmann and Zoo OIII Flux',xlabel='log$_{10}$ OIII Flux (ergs/s/cm$^2$)', axis=[0,3.5,0,1200])
#lumlog_hist(disk_cat['H_alpha'], edge, face, A=25,title=r'Brinchmann H$\alpha$ Flux',xlabel=r'log$_{10}$ H$\alpha$ (ergs/s/cm$^2$)', axis=[0.5,4.5,0,4000])


# creates mass histogram
def mass_hist(data, edge, face, A=50, axis=None):
	B=data.min()
	C=data.max()
	plt.hist(data[face],bins=A, range=(B,C),color='b',  histtype='stepfilled', alpha=0.6, label='face-on')
	plt.hist(data[edge], bins=A, range=(B,C), color='r',  histtype='stepfilled', alpha=0.6, label='edge-on')
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=2, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=3, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.legend(loc=0)
	plt.title('Stellar Mass Estimates', y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel('logM$\odot$',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.1)
	plt.gca().yaxis.set_label_coords(-.1, .5)
	plt.show()

# mass histograms for zoo and mendel
#mass_hist(disk_cat['LOG_MSTELLAR'], edge, face, A=10, axis=[9,12.5,0,820]) #zoo
#mass_hist(disk_cat['logM'], edge, face, A=15, axis=[9,12.5,0,600])	#mendel and zoo
#mass_hist(disk_cat['logM'], edge, face, A=15, axis=[8,12,0,6000])	#mendel
'''
cat=open_flux_mass_catalog()
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['sfr_flag']==0)
cut3=np.logical_and(cat['OIII_tot_SN']>3,cat['HA_SN']>3)
cut10=np.logical_and(cat['OIII_tot_SN']>10,cat['HA_SN']>10)
cut20=np.logical_and(cat['OIII_tot_SN']>20,cat['HA_SN']>20)
disk=cat[cut]
disk_SN3=cat[cut3]
disk_SN10=cat[cut10]
disk_SN20=cat[cut20]
edge=cat['cosi'][cut]<0.30
face=cat['cosi'][cut]>0.86
edge_SN3=cat['cosi'][cut3]<0.30
face_SN3=cat['cosi'][cut3]>0.81
edge_SN10=cat['cosi'][cut10]<0.30
face_SN10=cat['cosi'][cut10]>0.82
edge_SN20=cat['cosi'][cut20]<0.30
face_SN20=cat['cosi'][cut20]>0.82
print 'cat length: '+str(len(cat))
print ''
print 'disk cat length: '+str(len(disk))
print 'edge-on length: '+str(len(disk[edge]))
print 'face-on length: '+str(len(disk[face]))
print 'SFR min: '+str(np.amin(disk['sfr']))
print 'SFR max: '+str(np.amax(disk['sfr']))
print 'SFR mean: '+str(np.mean(disk['sfr']))
print ''
print 'Halpha and OIII total SN>3'
print 'disk cat length: '+str(len(disk_SN3))
print 'edge-on length: '+str(len(disk_SN3[edge_SN3]))
print 'face-on length: '+str(len(disk_SN3[face_SN3]))
print 'SFR min: '+str(np.amin(disk_SN3['sfr']))
print 'SFR max: '+str(np.amax(disk_SN3['sfr']))
print 'SFR mean: '+str(np.mean(disk_SN3['sfr']))
print ''
print 'disk cat length SN>10: '+str(len(disk_SN10))
print 'edge-on length SN>10: '+str(len(disk_SN10[edge_SN10]))
print 'face-on length SN>10: '+str(len(disk_SN10[face_SN10]))
print 'SFR min SN>10: '+str(np.amin(disk_SN10['sfr']))
print 'SFR max SN>10: '+str(np.amax(disk_SN10['sfr']))
print 'SFR mean SN>10: '+str(np.mean(disk_SN10['sfr']))
print ''
print 'disk cat length SN>20: '+str(len(disk_SN20))
print 'edge-on length SN>20: '+str(len(disk_SN20[edge_SN20]))
print 'face-on length SN>20: '+str(len(disk_SN20[face_SN20]))
print 'SFR min SN>20: '+str(np.amin(disk_SN20['sfr']))
print 'SFR max SN>20: '+str(np.amax(disk_SN20['sfr']))
print 'SFR mean SN>20: '+str(np.mean(disk_SN20['sfr']))

def SFR_hist(data,edge,face, A=20, title=None):
	plt.hist(data[edge],bins=A, color='r',  histtype='stepfilled', alpha=0.7, label='edge-on galaxies')
	plt.hist(data[face], bins=A, color='b',  histtype='stepfilled', alpha=0.7, label='face-on galaxies')
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean edge-on: %.3f'% np.mean(data[edge]))
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=3, label='mean face-on: %.3f'% np.mean(data[face]))
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	#plt.axis([-2,2,0,2200],fontsize=30, fontweight='bold')
	plt.legend(loc=2)
	plt.xlabel('log SFR (M$\odot$/yr)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.09, .5)
	plt.show()

SFR_hist(disk['sfr'],edge,face,title='Frequency of Galaxy Star Formation Rates All SN')
SFR_hist(disk_SN3['sfr'],edge_SN3,face_SN3, A=50,title='Frequency of Galaxy Star Formation Rates\n'+r'H$\alpha$ and OIII total SN > 3')
SFR_hist(disk_SN10['sfr'],edge_SN10,face_SN10,title='Frequency of Galaxy Star Formation Rates SN > 10')
SFR_hist(disk_SN20['sfr'],edge_SN20,face_SN20,A=10,title='Frequency of Galaxy Star Formation Rates SN > 20')
'''
'''
cat=open_flux_mass_catalog()
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['sfr_flag']==0)
HaO3_g=np.logical_and(cat['HA_SN']>10,cat['OIII_tot_SN']>3)
HaHbO3_g=np.logical_and(HaO3_g,cat['HB_SN']>10)
cutHa_g=np.logical_and(cat['HA_SN']>10,cut)
cutHa_l=np.logical_and(cat['HA_SN']<10,cut)
cutHb_g=np.logical_and(cat['HB_SN']>10,cut)
cutHb_l=np.logical_and(cat['HB_SN']<10,cut)
cutO3_g=np.logical_and(cat['OIII_tot_SN']>10,cut)
cutO3_l=np.logical_and(cat['OIII_tot_SN']<10,cut)
cutHaO3_g=np.logical_and(HaO3_g,cut)
cutHaHbO3_g=np.logical_and(HaHbO3_g,cut)
disk=cat[cut]
disk_Ha_g=cat[cutHa_g]
disk_Ha_l=cat[cutHa_l]
disk_Hb_g=cat[cutHb_g]
disk_Hb_l=cat[cutHb_l]
disk_O3_g=cat[cutO3_g]
disk_O3_l=cat[cutO3_l]
disk_HaO3_g=cat[cutHaO3_g]
disk_HaHbO3_g=cat[cutHaHbO3_g]
edge=cat['cosi'][cut]<0.30
face=cat['cosi'][cut]>0.86
edge_Ha_g=cat['cosi'][cutHa_g]<0.30
face_Ha_g=cat['cosi'][cutHa_g]>0.86
edge_Ha_l=cat['cosi'][cutHa_l]<0.30
face_Ha_l=cat['cosi'][cutHa_l]>0.86
edge_Hb_g=cat['cosi'][cutHb_g]<0.30
face_Hb_g=cat['cosi'][cutHb_g]>0.86
edge_Hb_l=cat['cosi'][cutHb_l]<0.30
face_Hb_l=cat['cosi'][cutHb_l]>0.86
edge_O3_g=cat['cosi'][cutO3_g]<0.30
face_O3_g=cat['cosi'][cutO3_g]>0.85
edge_O3_l=cat['cosi'][cutO3_l]<0.30
face_O3_l=cat['cosi'][cutO3_l]>0.90
edge_HaO3_g=cat['cosi'][cutHaO3_g]<0.30
face_HaO3_g=cat['cosi'][cutHaO3_g]>0.85
edge_HaHbO3_g=cat['cosi'][cutHaHbO3_g]<0.30
face_HaHbO3_g=cat['cosi'][cutHaHbO3_g]>0.85

print len(disk[edge])
print disk['sfr'][edge][0:20]
print ''
print len(disk_Ha_g[edge_Ha_g])
print disk_Ha_g['sfr'][edge_Ha_g][0:20]
print len(disk_Ha_l[edge_Ha_l])
print disk_Ha_l['sfr'][edge_Ha_l][0:20]
print ''
print len(disk_Hb_g[edge_Hb_g])
print disk_Hb_g['sfr'][edge_Hb_g][0:20]
print len(disk_Hb_l[edge_Hb_l])
print disk_Hb_l['sfr'][edge_Hb_l][0:20]
print ''
print len(disk_O3_g[edge_O3_g])
print disk_O3_g['sfr'][edge_O3_g][0:20]
print len(disk_O3_l[edge_O3_l])
print disk_O3_l['sfr'][edge_O3_l][0:20]
print ''
print len(disk_HaO3_g[edge_HaO3_g])
print disk_HaO3_g['sfr'][edge_HaO3_g][0:20]
print ''
print len(disk_HaHbO3_g[edge_HaHbO3_g])
print disk_HaHbO3_g['sfr'][edge_HaHbO3_g][0:20]
print ''

print 'cat length: '+str(len(cat))
print 'disk length: '+str(len(disk))
print 'disk edge-on length: '+str(len(disk[edge]))
print 'disk face-on length: '+str(len(disk[face]))
print 'SFR min: '+str(np.amin(disk['sfr']))
print 'SFR max: '+str(np.amax(disk['sfr']))
print 'SFR mean: '+str(np.mean(disk['sfr']))
print ''
print 'Halpha SN>10'
print 'disk cat length: '+str(len(disk_Ha_g))
print 'disk edge-on length: '+str(len(disk_Ha_g[edge_Ha_g]))
print 'disk face-on length: '+str(len(disk_Ha_g[face_Ha_g]))
print 'SFR min: '+str(np.amin(disk_Ha_g['sfr']))
print 'SFR max: '+str(np.amax(disk_Ha_g['sfr']))
print 'SFR mean: '+str(np.mean(disk_Ha_g['sfr']))
print ''
print 'Halpha SN<3'
print 'disk cat length: '+str(len(disk_Ha_l))
print 'disk edge-on length: '+str(len(disk_Ha_l[edge_Ha_l]))
print 'disk face-on length: '+str(len(disk_Ha_l[face_Ha_l]))
print 'SFR min: '+str(np.amin(disk_Ha_l['sfr']))
print 'SFR max: '+str(np.amax(disk_Ha_l['sfr']))
print 'SFR mean: '+str(np.mean(disk_Ha_l['sfr']))
print ''
print 'OIII total SN>10'
print 'disk cat length: '+str(len(disk_O3_g))
print 'disk edge-on length: '+str(len(disk_O3_g[edge_O3_g]))
print 'disk face-on length: '+str(len(disk_O3_g[face_O3_g]))
print 'SFR min: '+str(np.amin(disk_O3_g['sfr']))
print 'SFR max: '+str(np.amax(disk_O3_g['sfr']))
print 'SFR mean: '+str(np.mean(disk_O3_g['sfr']))
print ''
print 'OIII total SN<10'
print 'disk cat length: '+str(len(disk_O3_tot_l))
print 'disk edge-on length: '+str(len(disk_O3_tot_l[edge_O3_tot_l]))
print 'disk face-on length: '+str(len(disk_O3_tot_l[face_O3_tot_l]))
print 'SFR min: '+str(np.amin(disk_O3_tot_l['sfr']))
print 'SFR max: '+str(np.amax(disk_O3_tot_l['sfr']))
print 'SFR mean: '+str(np.mean(disk_O3_tot_l['sfr']))
print ''
print ' Halpha and OIII total SN>10'
print 'disk cat length: '+str(len(disk_HaO3_g))
print 'disk edge-on length: '+str(len(disk_HaO3_g[edge_HaO3_g]))
print 'disk face-on length: '+str(len(disk_HaO3_g[face_HaO3_g]))
print 'SFR min: '+str(np.amin(disk_HaO3_g['sfr']))
print 'SFR max: '+str(np.amax(disk_HaO3_g['sfr']))
print 'SFR mean: '+str(np.mean(disk_HaO3_g['sfr']))
print ''
print ' Halpha, Hbeta and OIII total SN>10'
print 'disk cat length: '+str(len(disk_HaHbO3_g))
print 'disk edge-on length: '+str(len(disk_HaHbO3_g[edge_HaHbO3_g]))
print 'disk face-on length: '+str(len(disk_HaHbO3_g[face_HaHbO3_g]))
print 'SFR min: '+str(np.amin(disk_HaHbO3_g['sfr']))
print 'SFR max: '+str(np.amax(disk_HaHbO3_g['sfr']))
print 'SFR mean: '+str(np.mean(disk_HaHbO3_g['sfr']))

def SFR_hist(data,edge,face, A=10, axis=None, title=None):
	plt.hist(data[edge],bins=A, color='r',  histtype='stepfilled', alpha=0.7, label='edge-on galaxies')
	plt.hist(data[face], bins=A, color='b',  histtype='stepfilled', alpha=0.7, label='face-on galaxies')
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean edge-on: %.3f'% np.mean(data[edge]))
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=3, label='mean face-on: %.3f'% np.mean(data[face]))
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.legend(loc=2)
	plt.xlabel('log SFR (M$\odot$/yr)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.09, .5)
	plt.show()

#SFR_hist(disk['sfr'],edge,face,axis=[-2.5,2,0,4000],title=r'Frequency of Galaxy Star Formation Rates')
#SFR_hist(disk_Ha_g['sfr'],edge_Ha_g,face_Ha_g,axis=[-2.5,2,0,3000],title=r'Frequency of Galaxy Star Formation Rates H$\alpha$ SN > 10')
#SFR_hist(disk_Ha_l['sfr'],edge_Ha_l,face_Ha_l,A=7, axis=[-1.5,2,0,50],title=r'Frequency of Galaxy Star Formation Rates H$\alpha$ SN < 3')
#SFR_hist(disk_O3_g['sfr'],edge_O3_g,face_O3_g,axis=[-2.5,2,0,3000],title='Frequency of Galaxy Star Formation Rates OIII total SN > 10')
#SFR_hist(disk_O3_tot_l['sfr'],edge_O3_tot_l,face_O3_tot_l,axis=[-2.5,2,0,600],title='Frequency of Galaxy Star Formation Rates OIII total SN < 3')
#SFR_hist(disk_HaO3_g['sfr'],edge_HaO3_g,face_HaO3_g,axis=[-2.5,2,0,3000],title='Frequency of Galaxy Star Formation Rates\n'+r'H$\alpha$ and OIII total SN > 3')
#SFR_hist(disk_HaHbO3_g['sfr'],edge_HaHbO3_g,face_HaHbO3_g,axis=[-2.5,2,0,3000],title='Frequency of Galaxy Star Formation Rates\n'+r'H$\alpha$, H$\beta$ and OIII total SN > 3')
'''
#HERE
'''
cat=open_line_mass_catalog() 
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['sfr_flag']==0)
cutHa=np.logical_and(cut,cat['H_alpha']>0)
cutHa_e=np.logical_and(cat['H_alpha']>0,cat['H_alpha_e']>0)
cutHa_e=np.logical_and(cut,cutHa_e)
disk=cat[cut]
disk_Ha=cat[cutHa]
disk_Ha_e=cat[cutHa_e]
Ha_SN=disk_Ha_e['H_alpha']/disk_Ha_e['H_alpha_e']
cutSN3=Ha_SN>3
disk_Ha_SN3=disk_Ha_e[cutSN3]
edge=cat['cosi'][cut]<0.30
face=cat['cosi'][cut]>0.87
edge_Ha=cat['cosi'][cutHa]<0.30
face_Ha=cat['cosi'][cutHa]>0.87
edge_Ha_SN3=disk_Ha_e['cosi'][cutSN3]<0.30
face_Ha_SN3=disk_Ha_e['cosi'][cutSN3]>0.86
'''

cat=open_zoo_mass_line_catalog()
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['sfr_flag']==0)
cutHa=np.logical_and(cut,cat['H_alpha']>0)
cutHa_e=np.logical_and(cat['H_alpha']>0,cat['H_alpha_e']>0)
cutHa_e=np.logical_and(cut,cutHa_e)
disk=cat[cut]
disk_Ha=cat[cutHa]
disk_Ha_e=cat[cutHa_e]
Ha_SN=disk_Ha_e['H_alpha']/disk_Ha_e['H_alpha_e']
cutSN3=Ha_SN>3
disk_Ha_SN3=disk_Ha_e[cutSN3]
edge=cat['cosi'][cut]<0.27
face=cat['cosi'][cut]>0.80
edge_Ha=cat['cosi'][cutHa]<0.27
face_Ha=cat['cosi'][cutHa]>0.80
edge_Ha_SN3=disk_Ha_e['cosi'][cutSN3]<0.27
face_Ha_SN3=disk_Ha_e['cosi'][cutSN3]>0.80
'''
print 'cat length: '+str(len(cat))
print 'disk length: '+str(len(disk))
print 'disk edge-on length: '+str(len(disk[edge]))
print 'disk face-on length: '+str(len(disk[face]))
print 'SFR min: '+str(np.amin(disk['sfr']))
print 'SFR max: '+str(np.amax(disk['sfr']))
print 'SFR mean: '+str(np.mean(disk['sfr']))
print ''
print 'Halpha>0'
print 'disk cat length: '+str(len(disk_Ha))
print 'disk edge-on length: '+str(len(disk_Ha[edge_Ha]))
print 'disk face-on length: '+str(len(disk_Ha[face_Ha]))
print 'SFR min: '+str(np.amin(disk_Ha['sfr']))
print 'SFR max: '+str(np.amax(disk_Ha['sfr']))
print 'SFR mean: '+str(np.mean(disk_Ha['sfr']))
print ''
print 'Halpha SN>3'
print 'disk cat length: '+str(len(disk_Ha_SN3))
print 'disk edge-on length: '+str(len(disk_Ha_SN3[edge_Ha_SN3]))
print 'disk face-on length: '+str(len(disk_Ha_SN3[face_Ha_SN3]))
print 'SFR min: '+str(np.amin(disk_Ha_SN3['sfr']))
print 'SFR max: '+str(np.amax(disk_Ha_SN3['sfr']))
print 'SFR mean: '+str(np.mean(disk_Ha_SN3['sfr']))
print ''
'''
def SFR_hist(data,edge,face, A=10, axis=None, title=None):
	plt.hist(data[edge],bins=A, color='r',  histtype='stepfilled', alpha=0.7, label='edge-on galaxies')
	plt.hist(data[face], bins=A, color='b',  histtype='stepfilled', alpha=0.7, label='face-on galaxies')
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean edge-on: %.3f'% np.mean(data[edge]))
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=3, label='mean face-on: %.3f'% np.mean(data[face]))
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.legend(loc=2)
	plt.xlabel('log SFR (M$\odot$/yr)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.09, .5)
	plt.show()

#SFR_hist(disk['sfr'],edge,face,axis=[-2,2.5,0,8000],title=r'Star Formation Rates')
#SFR_hist(disk_Ha['sfr'],edge_Ha,face_Ha,axis=[-2,2,0,8000],title=r'Brinchmann Galaxy Star Formation Rates H$\alpha$ > 0')
#SFR_hist(disk_Ha_SN3['sfr'],edge_Ha_SN3,face_Ha_SN3, axis=[-2,2.5,0,8000],title=r'Star Formation Rates H$\alpha$ Signal to Noise > 3')

#SFR_hist(disk['sfr'],edge,face,axis=[-2,1,0,1000],title=r'Star Formation Rates')
#SFR_hist(disk_Ha['sfr'],edge_Ha,face_Ha,axis=[-2,1,0,800],title=r'Zoo Galaxy Star Formation Rates H$\alpha$ > 0')
#SFR_hist(disk_Ha_SN3['sfr'],edge_Ha_SN3,face_Ha_SN3, axis=[-2,1,0,800],title=r'Zoo Galaxy Star Formation Rates H$\alpha$ > 0 SN > 3')

'''
cat=open_flux_mass_catalog()
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['sfr_flag']==0)
cutHa=np.logical_and(cat['HA_SN']>3,cut)
cutHa10=np.logical_and(cat['HA_SN']>10,cut)
disk=cat[cut]
disk_Ha=cat[cutHa]
disk_Ha10=cat[cutHa10]
edge=cat['cosi'][cut]<0.30
face=cat['cosi'][cut]>0.86
edge_Ha=cat['cosi'][cutHa]<0.30
face_Ha=cat['cosi'][cutHa]>0.86
edge_Ha10=cat['cosi'][cutHa10]<0.30
face_Ha10=cat['cosi'][cutHa10]>0.86

print 'cat length: '+str(len(cat))
print ''
print 'disk cat length: '+str(len(disk))
print 'edge-on length: '+str(len(disk[edge]))
print 'face-on length: '+str(len(disk[face]))
print 'Halpha min: '+str(np.amin(np.log10(disk['H_alpha'])))
print 'Halpha max: '+str(np.amax(np.log10(disk['H_alpha'])))
print 'Halpha mean: '+str(np.mean(np.log10(disk['H_alpha'])))
print ''
print 'Halpha SN>3'
print 'disk cat length: '+str(len(disk_Ha))
print 'edge-on length: '+str(len(disk[edge_Ha]))
print 'face-on length: '+str(len(disk[face_Ha]))
print 'Halpha min: '+str(np.amin(np.log10(disk_Ha['H_alpha'])))
print 'Halpha max: '+str(np.amax(np.log10(disk_Ha['H_alpha'])))
print 'Halpha mean: '+str(np.log10(np.mean(disk_Ha['H_alpha'])))
print ''
print 'Halpha SN>10'
print 'disk cat length: '+str(len(disk_Ha10))
print 'edge-on length: '+str(len(disk[edge_Ha10]))
print 'face-on length: '+str(len(disk[face_Ha10]))
print 'Halpha min: '+str(np.amin(np.log10(disk_Ha10['H_alpha'])))
print 'Halpha max: '+str(np.amax(np.log10(disk_Ha10['H_alpha'])))
print 'Halpha mean: '+str(np.log10(np.mean(disk_Ha10['H_alpha'])))

def flux_hist(data,edge,face, A=7, title=None):
	plt.hist(data[edge],bins=A, color='r',  histtype='stepfilled', alpha=0.7, label='edge-on galaxies')
	plt.hist(data[face], bins=A, color='b',  histtype='stepfilled', alpha=0.7, label='face-on galaxies')
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean edge-on: %.3f'% np.mean(data[edge]))
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=3, label='mean face-on: %.3f'% np.mean(data[face]))
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	#plt.axis([-2,2,0,2200],fontsize=30, fontweight='bold')
	plt.legend(loc=1)
	plt.xlabel('log SFR (M$\odot$/yr)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.09, .5)
	plt.show()

#flux_hist(np.log10(disk['H_alpha']),edge,face,title=r'Frequency of H$\alpha$ Flux')
#flux_hist(np.log10(disk_Ha['H_alpha']),edge_Ha,face_Ha,title=r'Frequency of H$\alpha$ Flux SN > 3')
#flux_hist(np.log10(disk_Ha10['H_alpha']),edge_Ha10,face_Ha10,title=r'Frequency of H$\alpha$ Flux SN > 10')
'''
'''
def test_zoo_cat():
	hdulist1=aio.fits.open(fnzphot)
	hdulist2=aio.fits.open(fnzspec)
	hdulist3=aio.fits.open(fnzTable2)
	hdulist4=aio.fits.open(fnzTable3)
	hdulist5=aio.fits.open(fnzGZ)    
	data1=hdulist1[1].data
	data2=hdulist2[1].data
	data3=hdulist3[1].data
	data4=hdulist4[1].data
	data5=hdulist5[1].data

	print 'fnzphot'
	print 'number of columns list: '+str(len(data1[0]))
	print 'length of column in list: '+str(len(data1['dr8objid']))
	print 'fnzphot names'
	print ''
	print 'fnzspec'
	print 'number of columns list: '+str(len(data2[0]))
	print 'length of column in list: '+str(len(data2['dr8objid']))
	print 'fnzspec names'
	print ''
	print 'fnzTable2'
	print 'number of columns list: '+str(len(data3[0]))
	print 'length of column in list: '+str(len(data3['OBJID']))
	print 'fnzTable2 names'
	print data3.names
	print ''
	print 'fnzTable3'
	print 'number of columns list: '+str(len(data4[0]))
	print 'length of column in list: '+str(len(data4['OBJID']))
	print 'fnzTable3 names'
	print data4.names
	print ''
	print 'fnzGZ'
	print 'number of lists in hdulist: '+str(len(data5[0]))
	print 'length of list in hdulist: '+str(len(data5['OBJID']))
	print 'length of column in list in hdulist: '+str(len(data5['OBJID'][0]))
	print 'fnzGZ names'
	print data5.names
	print ''
'''	