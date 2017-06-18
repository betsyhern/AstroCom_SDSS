'''Initialized 7/24/15'''

import os as os
import numpy as np
import astropy.io as aio
from astropy.table import Table,join,Column,unique
import astropy.constants as aconst
import astropy.cosmology as acosmo
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u 
import matplotlib.pyplot as plt
from astropy.table import Table

def	read_fits_table(filename):
	hdulist=aio.fits.open(filename)
	tbdata=hdulist[1].data
	return tbdata
	
'''	 
agnBr='Brinchmann2014.fit' 
agnBe='Best2005.fit' 
agnSa ='Saintonge2012.fit' 

agn_table1=read_fits_table(agnBr)
agn1_cut=agn_table1['Plate']!=0
agn_table1=agn_table1[agn1_cut]

agn_table2=read_fits_table(agnBe)
agn_table3=read_fits_table(agnSa)
'''

fnS1='Simard_Table3.fit' #1D fits 
fnS2='Simard_Table1.fit' #2D fits with n_b = 4
fnM ='Mendel_dusty.fit'  #dusty masses combines Tables 3 and 4
fnCl='gal_iclass_table_dr7_v5_2.fits'	#objID, class 
fnLi='gal_line_dr7_v5_2.fit'			#OIII 5007, NII6584, H_alpha, H_beta
fnSFR='gal_totsfr_dr7_v5_2.fits'			##MPA DR7 galaxy masses (AVG, flag)
fnsSFR='gal_totspecsfr_dr7_v5_2.fits'		#AVG, flag


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

def create_MPA_catalog():
	galCl_table=read_fits_table(fnCl)
	galLi_table=read_fits_table(fnLi)
	SFR_table=read_fits_table(fnSFR)
	sSFR_table=read_fits_table(fnsSFR)
	S_M_table=create_mass_catalog()
	SDSS_table=Table([galCl_table['objID'],		
		galCl_table['I_CLASS'],
		galLi_table['OIII_5007_FLUX'],
		galLi_table['NII_6584_FLUX'],
		galLi_table['SII_6717_FLUX'],
		galLi_table['OI_6300_FLUX'],
		galLi_table['H_ALPHA_FLUX'],
		galLi_table['H_BETA_FLUX'],
		SFR_table['AVG'],
		sSFR_table['AVG'],
		SFR_table['FLAG'],
		sSFR_table['FLAG']],
		names=('objID','class','OIII','NII','SII','OI','H_a','H_b','SFR avg','sSFR avg','SFR flag','sSFR flag'), 
		dtype=(np.int64,np.int32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.float32,np.int16,np.int16))
	SDSS_table=unique(SDSS_table,keys='objID')
	shared=['objID']
	t=join(S_M_table,SDSS_table,keys=shared)
	t.rename_column('_RA','RA')
	t.rename_column('_DE','DE')
	return t

cat=create_MPA_catalog()
'''
disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
edge=cat['cosi'][disk]<0.30
face=cat['cosi'][disk]>0.87

# frequency of inclination cosine angles
def cat_hist(all,cut,color1=None, color2=None,title=None):
	plt.hist(all, bins=12, range=(all.min(),all.max()), color=color1, alpha=0.6, label='All Catalog Galaxies')
	plt.hist(all[cut], bins=12, range=(all.min(),all.max()), color=color2, alpha=0.6, label='Catalog Galaxies with B/T < 0.35')
	plt.axis([-0.1,1.1,0,160000],fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.legend(loc=2)
	plt.show()

#cat_hist(cat['cosi'],disk,color1='y',color2='m',title='Inclination of Cosine Angle Frequency')

#histogram of stellar masses of disc 0-0.3 cut, edge-on(cosine<0.3), and face-on(cosine>0.83)
def mass_hist(data, edge, face, A=50, title=None, xlabel=None):
	B=data.min()
	C=data.max()
	plt.hist(data[face],bins=A, range=(B,C),color='b',  histtype='stepfilled', alpha=0.6, label='face-on')
	plt.hist(data[edge], bins=A, range=(B,C), color='r',  histtype='stepfilled', alpha=0.6, label='edge-on')
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=2, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=3, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.legend(loc=0)
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel('logM$\odot$',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.axis([8,12,0,2000],fontsize=30, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.1)
	plt.gca().yaxis.set_label_coords(-.1, .5)
	plt.show()

mass_hist(cat['logM'][disk],edge,face,title='Frequency of logM$\odot$')

disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
disk_cat=cat[disk]
cut=disk_cat['H_a']>0
disk_cat=disk_cat[cut]
print 'disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.87
print 'edge-on length index: '+str(len(edge))
print 'face-on length index: '+str(len(face))
print 'edge lenght: '+str(len(disk_cat['H_a'][edge]))
print 'face lenght: '+str(len(disk_cat['H_a'][face]))
print''
print 'H_a'
print 'disk cat lenght: '+str(len(disk_cat['H_a']))
print 'edge lenght: '+str(len(disk_cat['H_a'][edge]))
print 'face lenght: '+str(len(disk_cat['H_a'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['H_a']))
print 'disk cat max: '+str(np.amax(disk_cat['H_a']))
print 'disk cat mean: '+str(np.mean(disk_cat['H_a']))
print 'edge mean: '+str(np.mean(disk_cat['H_a'][edge]))
print 'face mean: '+str(np.mean(disk_cat['H_a'][face]))
print''

plt.rc('text', usetex=False)
plt.hist(np.log10(disk_cat['H_a'][edge]), bins=20, align='right', color='r',histtype='stepfilled', alpha=0.8,label=r'H$\alpha$ edge-on')
plt.hist(np.log10(disk_cat['H_a'][face]), bins=20, align='right', color='b',histtype='stepfilled', alpha=0.7,label=r'H$\alpha$ face-on')
plt.axvline(np.log10(disk_cat['H_a'][face]).mean(), color='k', linestyle='--', linewidth=2, label=r'mean H$\alpha$ face-on: %.2f'% np.mean(np.log10(disk_cat['H_a'][face])))
plt.axvline(np.log10(disk_cat['H_a'][edge]).mean(), color='k', linestyle='-', linewidth=2, label=r'mean H$\alpha$ edge-on: %.2f'% np.mean(np.log10(disk_cat['H_a'][edge])))
plt.title(r'H$\alpha$ Flux Frequency', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel(r'log H$\alpha$ Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.axis([0,4,0,5000],fontsize=40, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=2)
plt.show()

cut=disk_cat['H_b']>0
disk_cat=disk_cat[cut]
print 'disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.87
print 'edge-on length index: '+str(len(edge))
print 'face-on length index: '+str(len(face))
print 'edge lenght: '+str(len(disk_cat['H_b'][edge]))
print 'face lenght: '+str(len(disk_cat['H_b'][face]))
print''
print 'H_b'
print 'disk cat lenght: '+str(len(disk_cat['H_b']))
print 'edge lenght: '+str(len(disk_cat['H_b'][edge]))
print 'face lenght: '+str(len(disk_cat['H_b'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['H_b']))
print 'disk cat max: '+str(np.amax(disk_cat['H_b']))
print 'disk cat mean: '+str(np.mean(disk_cat['H_b']))
print 'edge mean: '+str(np.mean(disk_cat['H_b'][edge]))
print 'face mean: '+str(np.mean(disk_cat['H_b'][face]))
print''

plt.rc('text', usetex=False)
plt.hist(np.log10(disk_cat['H_b'][edge]), bins=20, align='right', color='r',histtype='stepfilled', alpha=0.8,label=r'H$\beta$ edge-on')
plt.hist(np.log10(disk_cat['H_b'][face]), bins=20, align='right', color='b',histtype='stepfilled', alpha=0.7,label=r'H$\beta$ face-on')
plt.axvline(np.log10(disk_cat['H_b'][face]).mean(), color='k', linestyle='--', linewidth=2, label=r'mean H$\beta$ face-on: %.2f'% np.mean(np.log10(disk_cat['H_b'][face])))
plt.axvline(np.log10(disk_cat['H_b'][edge]).mean(), color='k', linestyle='-', linewidth=2, label=r'mean H$\beta$ edge-on: %.2f'% np.mean(np.log10(disk_cat['H_b'][edge])))
plt.title(r'H$\beta$ Flux Frequency', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel(r'log H$\beta$ Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.axis([-0.5,3,0,5000],fontsize=40, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=2)
plt.show()

cut=disk_cat['OIII']>0
disk_cat=disk_cat[cut]
print 'disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.87
print 'edge-on length index: '+str(len(edge))
print 'face-on length index: '+str(len(face))
print 'edge lenght: '+str(len(disk_cat['OIII'][edge]))
print 'face lenght: '+str(len(disk_cat['OIII'][face]))
print''
print 'OIII'
print 'disk cat lenght: '+str(len(disk_cat['OIII']))
print 'edge lenght: '+str(len(disk_cat['OIII'][edge]))
print 'face lenght: '+str(len(disk_cat['OIII'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['OIII']))
print 'disk cat max: '+str(np.amax(disk_cat['OIII']))
print 'disk cat mean: '+str(np.mean(disk_cat['OIII']))
print 'edge mean: '+str(np.mean(disk_cat['OIII'][edge]))
print 'face mean: '+str(np.mean(disk_cat['OIII'][face]))
print''

plt.rc('text', usetex=False)
plt.hist(np.log10(disk_cat['OIII'][edge]), bins=20, align='right', color='r',histtype='stepfilled', alpha=0.8,label='OIII edge-on')
plt.hist(np.log10(disk_cat['OIII'][face]), bins=20, align='right', color='b',histtype='stepfilled', alpha=0.7,label='OIII face-on')
plt.axvline(np.log10(disk_cat['OIII'][face]).mean(), color='k', linestyle='--', linewidth=2, label='mean OIII face-on: %.2f'% np.mean(np.log10(disk_cat['OIII'][face])))
plt.axvline(np.log10(disk_cat['OIII'][edge]).mean(), color='k', linestyle='-', linewidth=2, label='mean OIII edge-on: %.2f'% np.mean(np.log10(disk_cat['OIII'][edge])))
plt.title('OIII Flux Frequency', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel('OIII Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.axis([-0.5,3,0,6200],fontsize=40, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=2)
plt.show()

cut=disk_cat['NII']>0
disk_cat=disk_cat[cut]
print 'disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.87
print 'edge-on length index: '+str(len(edge))
print 'face-on length index: '+str(len(face))
print 'edge lenght: '+str(len(disk_cat['NII'][edge]))
print 'face lenght: '+str(len(disk_cat['NII'][face]))
print''
print 'NII'
print 'disk cat lenght: '+str(len(disk_cat['NII']))
print 'edge lenght: '+str(len(disk_cat['NII'][edge]))
print 'face lenght: '+str(len(disk_cat['NII'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['NII']))
print 'disk cat max: '+str(np.amax(disk_cat['NII']))
print 'disk cat mean: '+str(np.mean(disk_cat['NII']))
print 'edge mean: '+str(np.mean(disk_cat['NII'][edge]))
print 'face mean: '+str(np.mean(disk_cat['NII'][face]))
print''

plt.rc('text', usetex=False)
plt.hist(np.log10(disk_cat['NII'][edge]), bins=20, align='right', color='r',histtype='stepfilled', alpha=0.8,label='NII edge-on')
plt.hist(np.log10(disk_cat['NII'][face]), bins=20, align='right', color='b',histtype='stepfilled', alpha=0.7,label='NII face-on')
plt.axvline(np.log10(disk_cat['NII'][face]).mean(), color='k', linestyle='--', linewidth=2, label='mean NII face-on: %.2f'% np.mean(np.log10(disk_cat['NII'][face])))
plt.axvline(np.log10(disk_cat['NII'][edge]).mean(), color='k', linestyle='-', linewidth=2, label='mean NII edge-on: %.2f'% np.mean(np.log10(disk_cat['NII'][edge])))
plt.title('NII Flux Frequency', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel('NII Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.axis([0,3.5,0,8000],fontsize=40, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=2)
plt.show()

cut=disk_cat['SII']>0
disk_cat=disk_cat[cut]
print 'disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.86
print 'edge-on length index: '+str(len(edge))
print 'face-on length index: '+str(len(face))
print 'edge lenght: '+str(len(disk_cat['SII'][edge]))
print 'face lenght: '+str(len(disk_cat['SII'][face]))
print''
print 'NII'
print 'disk cat lenght: '+str(len(disk_cat['SII']))
print 'edge lenght: '+str(len(disk_cat['SII'][edge]))
print 'face lenght: '+str(len(disk_cat['SII'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['SII']))
print 'disk cat max: '+str(np.amax(disk_cat['SII']))
print 'disk cat mean: '+str(np.mean(disk_cat['SII']))
print 'edge mean: '+str(np.mean(disk_cat['SII'][edge]))
print 'face mean: '+str(np.mean(disk_cat['SII'][face]))
print''

plt.rc('text', usetex=False)
plt.hist(np.log10(disk_cat['SII'][edge]), bins=20, align='right', color='r',histtype='stepfilled', alpha=0.8,label='SII edge-on')
plt.hist(np.log10(disk_cat['SII'][face]), bins=20, align='right', color='b',histtype='stepfilled', alpha=0.7,label='SII face-on')
plt.axvline(np.log10(disk_cat['SII'][face]).mean(), color='k', linestyle='--', linewidth=2, label='mean SII face-on: %.2f'% np.mean(np.log10(disk_cat['SII'][face])))
plt.axvline(np.log10(disk_cat['SII'][edge]).mean(), color='k', linestyle='-', linewidth=2, label='mean SII edge-on: %.2f'% np.mean(np.log10(disk_cat['SII'][edge])))
plt.title('SII Flux Frequency', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel('SII Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.axis([-0.5,3,0,7000],fontsize=40, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=2)
plt.show()

cut=disk_cat['OI']>0
disk_cat=disk_cat[cut]
print 'disk cat length: '+str(len(disk_cat))
edge=disk_cat['cosi']<0.30
face=disk_cat['cosi']>0.86
print 'edge-on length index: '+str(len(edge))
print 'face-on length index: '+str(len(face))
print 'edge lenght: '+str(len(disk_cat['OI'][edge]))
print 'face lenght: '+str(len(disk_cat['OI'][face]))
print''
print 'OI'
print 'disk cat lenght: '+str(len(disk_cat['OI']))
print 'edge lenght: '+str(len(disk_cat['OI'][edge]))
print 'face lenght: '+str(len(disk_cat['OI'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['OI']))
print 'disk cat max: '+str(np.amax(disk_cat['OI']))
print 'disk cat mean: '+str(np.mean(disk_cat['OI']))
print 'edge mean: '+str(np.mean(disk_cat['OI'][edge]))
print 'face mean: '+str(np.mean(disk_cat['OI'][face]))
print''

plt.rc('text', usetex=False)
plt.hist(np.log10(disk_cat['OI'][edge]), bins=20, align='right', color='r',histtype='stepfilled', alpha=0.8,label='OI edge-on')
plt.hist(np.log10(disk_cat['OI'][face]), bins=20, align='right', color='b',histtype='stepfilled', alpha=0.7,label='OI face-on')
plt.axvline(np.log10(disk_cat['OI'][face]).mean(), color='k', linestyle='--', linewidth=2, label='mean OI face-on: %.2f'% np.mean(np.log10(disk_cat['OI'][face])))
plt.axvline(np.log10(disk_cat['OI'][edge]).mean(), color='k', linestyle='-', linewidth=2, label='mean OI edge-on: %.2f'% np.mean(np.log10(disk_cat['OI'][edge])))
plt.title('OI Flux Frequency', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel('OI Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.axis([-0.5,2.5,0,9000],fontsize=40, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=2)
plt.show()

disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
disk_SFR=np.logical_and(disk,cat['SFR flag']==0)
disk_SFR_cat=cat[disk_SFR]
print 'disk cat length: '+str(len(disk_SFR_cat))
edge_SFR=cat['cosi'][disk_SFR]<0.30
face_SFR=cat['cosi'][disk_SFR]>0.87
print 'edge-on length: '+str(len(disk_SFR_cat[edge_SFR]))
print 'face-on length: '+str(len(disk_SFR_cat[face_SFR]))
print 'SFR min: '+str(np.amin(disk_SFR_cat['SFR avg']))
print 'SFR max: '+str(np.amax(disk_SFR_cat['SFR avg']))
print 'SFR mean: '+str(np.mean(disk_SFR_cat['SFR avg']))

def SFR_hist(data,edge,face, A=20, title=None):
	plt.hist(data[edge],bins=A, color='r',  histtype='stepfilled', alpha=0.7, label='edge-on galaxies')
	plt.hist(data[face], bins=A, color='b',  histtype='stepfilled', alpha=0.7, label='face-on galaxies')
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=3, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.axis([-2,2,0,3000],fontsize=30, fontweight='bold')
	plt.legend(loc=2)
	plt.xlabel('log SFR (M$\odot$/yr)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.09, .5)
	plt.show()

SFR_hist(disk_SFR_cat['SFR avg'],edge_SFR,face_SFR,title='Frequency of Galaxy Star Formation Rates')

'''
cut=cat['OIII']>0
cat=cat[cut]
print 'cat length OIII cut: '+str(len(cat))
cut=cat['NII']>0
cat=cat[cut]
print 'cat length NII cut: '+str(len(cat))
cut=cat['SII']>0
cat=cat[cut]
print 'cat length SII cut: '+str(len(cat))
cut=cat['OI']>0
cat=cat[cut]
print 'cat length OI cut: '+str(len(cat))
cut=cat['H_a']>0
cat=cat[cut]
print 'cat length H_a cut: '+str(len(cat))
cut=cat['H_b']>0
cat=cat[cut]
print 'cat length H_b cut: '+str(len(cat))
print ''

disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
disk_cat=cat[disk]
print 'disk cat length: '+str(len(disk_cat))
edge=cat['cosi'][disk]<0.30
face=cat['cosi'][disk]>0.85

print 'edge-on length index: '+str(len(edge))
print 'face-on length index: '+str(len(face))
print ''

print 'OIII'
print 'disk cat lenght: '+str(len(disk_cat['OIII']))
print 'edge lenght: '+str(len(disk_cat['OIII'][edge]))
print 'face lenght: '+str(len(disk_cat['OIII'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['OIII']))
print 'disk cat max: '+str(np.amax(disk_cat['OIII']))
print 'disk cat mean: '+str(np.mean(disk_cat['OIII']))
print 'edge mean: '+str(np.mean(disk_cat['OIII'][edge]))
print 'face mean: '+str(np.mean(disk_cat['OIII'][face]))
print ''
print 'NII'
print 'disk cat lenght: '+str(len(disk_cat['NII']))
print 'edge lenght: '+str(len(disk_cat['NII'][edge]))
print 'face lenght: '+str(len(disk_cat['NII'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['NII']))
print 'disk cat max: '+str(np.amax(disk_cat['NII']))
print 'disk cat mean: '+str(np.mean(disk_cat['NII']))
print 'edge mean: '+str(np.mean(disk_cat['NII'][edge]))
print 'face mean: '+str(np.mean(disk_cat['NII'][face]))
print''
print 'H_a'
print 'disk cat lenght: '+str(len(disk_cat['H_a']))
print 'edge lenght: '+str(len(disk_cat['H_a'][edge]))
print 'face lenght: '+str(len(disk_cat['H_a'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['H_a']))
print 'disk cat max: '+str(np.amax(disk_cat['H_a']))
print 'disk cat mean: '+str(np.mean(disk_cat['H_a']))
print 'edge mean: '+str(np.mean(disk_cat['H_a'][edge]))
print 'face mean: '+str(np.mean(disk_cat['H_a'][face]))
print''
print 'H_b'
print 'disk cat lenght: '+str(len(disk_cat['H_b']))
print 'edge lenght: '+str(len(disk_cat['H_b'][edge]))
print 'face lenght: '+str(len(disk_cat['H_b'][face]))
print 'disk cat min: '+str(np.amin(disk_cat['H_b']))
print 'disk cat max: '+str(np.amax(disk_cat['H_b']))
print 'disk cat mean: '+str(np.mean(disk_cat['H_b']))
print 'edge mean: '+str(np.mean(disk_cat['H_b'][edge]))
print 'face mean: '+str(np.mean(disk_cat['H_b'][face]))


log_NII_Ha=np.log10(disk_cat['NII']/disk_cat['H_a'])
print 'NII_Ha min: '+str(np.amin(log_NII_Ha))
print 'NII_Ha max: '+str(np.amax(log_NII_Ha))
print ''

log_SII_Ha=np.log10(disk_cat['SII']/disk_cat['H_a'])
print 'SII_Ha min: '+str(np.amin(log_SII_Ha))
print 'SII_Ha max: '+str(np.amax(log_SII_Ha))
print ''

log_OI_Ha=np.log10(disk_cat['OI']/disk_cat['H_a'])
print 'OI_Ha min: '+str(np.amin(log_OI_Ha))
print 'OI_Ha max: '+str(np.amax(log_OI_Ha))
print ''

log_OIII_Hb=np.log10(disk_cat['OIII']/disk_cat['H_b'])
print 'OIII_Hb min: '+str(np.amin(log_OIII_Hb))
print 'OIII_Hb max: '+str(np.amax(log_OIII_Hb))
print ''

eqn_NII_Ha=((0.61/(log_NII_Ha-0.47))+1.19)
eqn_NII_Ha_K=((0.61/(log_NII_Ha-0.05))+1.3)
eqn_SII_Ha=((0.72/(log_SII_Ha-0.32))+1.30)
eqn_OI_Ha=((0.73/(log_OI_Ha+0.59))+1.33)

V = [np.log10(0.6) for x in range(0,len(log_OIII_Hb))]
H_NII_Ha=log_NII_Ha>=np.log10(0.6)
H_NII_Ha=log_NII_Ha[H_NII_Ha]
H = [np.log10(3) for x in range(0,len(H_NII_Ha))]

# Two subplots, unpack the axes array immediately
#f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
#ax1.plot(x, y)
#ax1.set_title('Sharing Y axis')
#ax2.scatter(x, y)
#f,(ax1)=plt.subplots(1,3,sharey=True)
plt.rc('text', usetex=False)
#plt.scatter(log_NII_Ha,log_OIII_Hb,c='k',marker='.', label='all B/T < 0.35 galaxies')
plt.scatter(log_NII_Ha[face][0:len(log_NII_Ha[face]):15],log_OIII_Hb[face][0:len(log_OIII_Hb[face]):15],c='b',marker='o',label='face-on')
plt.scatter(log_NII_Ha[edge][0:len(log_NII_Ha[edge]):15],log_OIII_Hb[edge][0:len(log_OIII_Hb[edge]):15],c='r',marker='o',label='edge-on')
plt.scatter(log_NII_Ha[edge],eqn_NII_Ha[edge],color='r',marker='.')
plt.scatter(log_NII_Ha[face],eqn_NII_Ha[face],color='b',marker='.')
plt.scatter(log_NII_Ha,eqn_NII_Ha,color='y',alpha=0.7,marker='.')
plt.scatter(log_NII_Ha[edge],eqn_NII_Ha_K[edge],alpha=0.4,color='r')
plt.scatter(log_NII_Ha[face],eqn_NII_Ha_K[face],alpha=0.4,color='b',marker='.')
plt.scatter(log_NII_Ha,eqn_NII_Ha_K,color='c',alpha=0.4,marker='.')
plt.plot(H_NII_Ha,H,color='g',linestyle='--',linewidth=2.5)
plt.plot(V,log_OIII_Hb,color='g',linestyle='--',linewidth=2.5)
plt.title(r'OIII/H$\beta$ versus NII/H$\alpha$ BPT Diagram', y=1.05,fontsize=24, fontweight='bold')
plt.axis([-1.5,0.68,-2,2.3],fontsize=40, fontweight='bold')
plt.xlabel(r'log NII/H$\alpha$', fontsize=20, fontweight='bold')
plt.ylabel(r'log OIII/H$\beta$', fontsize=20, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.09)
plt.gca().yaxis.set_label_coords(-.07, .5)
plt.legend(loc=2)
plt.show()
'''
ax2.scatter(log_SII_Ha,log_OIII_Hb,c='k',marker='.',label='disk galaxy sample')
ax2.scatter(log_SII_Ha,eqn_SII_Ha,color='g',marker='o',label='disk galaxy sample\nseparation')
ax2.scatter(log_SII_Ha[edge],eqn_SII_Ha[edge],color='m',marker='^',label='edge-on separation')
ax2.scatter(log_SII_Ha[face],eqn_SII_Ha[face],color='c',marker='.',label='face-on separation')
ax2.title('OIII/H_beta versus SII/H_alpha\nOptical Diagnostic Diagram (Flux > 0)', y=1.05,fontsize=18, fontweight='bold')
ax2.axis([-2,2,-4,4],fontsize=40, fontweight='bold')
ax2.xlabel('NII/H_alpha', fontsize=16, fontweight='bold')
ax2.ylabel('OIII/H_beta', fontsize=16, fontweight='bold')
ax2.gca().xaxis.set_label_coords(.5, -.07)
ax2.gca().yaxis.set_label_coords(-.07, .5)
ax2.legend(loc=0)

plt.scatter(log_OI_Ha,log_OIII_Hb,c='k',marker='.',label='disk galaxy sample')
plt.scatter(log_OI_Ha,eqn_OI_Ha,color='g',marker='o',label='disk galaxy sample\nseparation')
plt.scatter(log_OI_Ha[edge],eqn_OI_Ha[edge],color='m',marker='^',label='edge-on separation')
plt.scatter(log_OI_Ha[face],eqn_OI_Ha[face],color='c',marker='.',label='face-on separation')
plt.title('OIII/H_beta versus SII/H_alpha\nOptical Diagnostic Diagram (Flux > 0)', y=1.05,fontsize=18, fontweight='bold')
plt.axis([-2,2,-4,4],fontsize=40, fontweight='bold')
plt.xlabel('NII/H_alpha', fontsize=16, fontweight='bold')
plt.ylabel('OIII/H_beta', fontsize=16, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.07)
plt.gca().yaxis.set_label_coords(-.07, .5)
#plt.legend(loc=0)
plt.show()



print len(cat['z'])
cut=cat['OIII']>0
cat=cat[cut]
print len(cat['z'])
cut=cat['NII']>0
cat=cat[cut]
print len(cat['z'])
'''
def element(Type=None):
	if Type=='OIII':
		plt.hist(cat['OIII'], 10, align='right', color='g')
		plt.title('OIII frequency', y=1.05,fontsize=17, fontweight='bold')
		plt.xlabel('OIII 5007', fontsize=16, fontweight='bold')
		plt.axvline(cat['OIII'].mean(), color='k', linestyle='-', label="mean cat['OIII']: %.2f"% np.mean(cat['OIII']), linewidth=1)
		#plt.axis([0,2e9,0,4],fontsize=30, fontweight='bold')
	if Type=='NII':
		plt.hist(cat['NII'], 10, align='right', color='b')
		plt.title('NII frequency', y=1.05,fontsize=17, fontweight='bold')
		plt.xlabel('NII 6584', fontsize=16, fontweight='bold')
		#plt.axis([0,1.6e11,0,4],fontsize=30, fontweight='bold')
	if Type=='H_a':
		plt.hist(cat['H_a'], 60)#, range=(0,6), align='right', color='g')
		plt.title('H_a $\alpha$ frequency', y=1.05,fontsize=17, fontweight='bold')
		plt.xlabel('H_a $\alpha$', fontsize=16, fontweight='bold')
	if Type=='H_b':
		plt.hist(cat['H_b'], 60)#, range=(0,6), align='right', color='g')
		plt.title('H_b $\beta$ frequency', y=1.05,fontsize=17, fontweight='bold')
		plt.xlabel('H_b $\beta$', fontsize=16, fontweight='bold')
	plt.ylabel('Frequency', fontsize=16, fontweight='bold')
	plt.legend(loc=0)
	plt.show()
	
'''
print np.amin(cat['OIII'])
print np.amax(cat['OIII'])
print np.amin(cat['NII'])
print np.amax(cat['NII'])
print np.amin(cat['H_a'])
print np.amax(cat['H_a'])
print np.amin(cat['H_b'])
print np.amax(cat['H_b'])

element(Type='OIII')
#element(Type='NII')


print len(cat)
cut=np.logical_and(cat['H_a']>0,cat['H_b']>0)
cut=np.logical_and(cut,cat['OIII']>0)
cut=np.logical_and(cut,cat['NII']>0)
cat=cat[cut]
print len(cat)
print np.amin(cat['OIII'])
print np.amax(cat['OIII'])
print np.amin(cat['NII'])
print np.amax(cat['NII'])

zero_a=np.where(cat['H_a']==0)
zero_b=np.where(cat['H_b']==0)
#print zero_a
#print zero_b

print cat['H_a'][10192]
print ''
print cat['H_b'][7243]

print 'OIII_5007_FLUX'
print cat['OIII'][0:16]
print cat['OIII'][16:31]
print ''
print 'NII_6584_FLUX'
print cat['NII'][0:16]
print cat['NII'][16:31]
print ''
print 'H_GAMMA_FLUX'
print cat['H_a'][0:16]
print cat['H_a'][16:31]
print ''
print 'H_BETA_FLUX'
print cat['H_b'][0:16]
print cat['H_b'][16:31]
print ''
print 'SFR flag'
print cat['SFR flag'][0:16]
print cat['SFR flag'][16:31]
print ''
print 'sSFR flag'
print cat['sSFR flag'][0:16]
print cat['sSFR flag'][16:31]
print ''

X=cat['NII']/cat['H_a']
Y=cat['OIII']/cat['H_b']
print np.mean(X)
print np.mean(Y)

plt.scatter(X,Y,color='k', marker='.', label='all galaxies')
plt.title('OIII5007/Hbeta versus NII6584/Halpha', fontsize=17, fontweight='bold')
plt.xlabel('NII6584/Halpha', fontsize=16, fontweight='bold')
plt.ylabel('OIII5007/Hbeta', fontsize=16, fontweight='bold')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.001,10000)
plt.ylim(0.001,10000)
#plt.axis([-1,100,-1,1000],fontsize=30, fontweight='bold')
#plt.xticks(range(100)[0:100:10])
#plt.yticks(range(24)[0:24:2])
plt.legend(loc=0)
plt.show()

c1 = SkyCoord(ra=agn_table1['RAJ2000']*u.degree, dec=agn_table1['DEJ2000']*u.degree)  
catalog1 = SkyCoord(ra=cat['_RA']*u.degree, dec=cat['_DE']*u.degree)  
idx1, d2d1, d3d1 = c1.match_to_catalog_sky(catalog1) 
#cat['objID'].format='18d'
t1 = Table([cat['objID'][idx1],agn_table1['Plate'],agn_table1['MJD'],agn_table1['Fib'],agn_table1['RAJ2000'],cat['_RA'][idx1],
	agn_table1['DEJ2000'],cat['_DE'][idx1]], names=('objID','Plate','MJD','Fib','agnRA', 'catRA', 'agnDE','catDE'), 
	dtype=(np.int64,np.int16,np.int32,np.int16,np.float32,np.float32,np.float32,np.float32))


print 'table 1'	
print t1[0:16]
print t1[16:31]
print t1[31:46]
print t1[46:60]
print ''

print ''
print len(agn_table1['RAJ2000'])
print len(idx1)
print idx1[0:2]
print 'agn1 and cat'
print 'RA'
print agn_table1['RAJ2000']
print cat['_RA'][idx1]
print ''
print 'DEC'
print agn_table1['DEJ2000']
print cat['_DE'][idx1]
print ''

c2 = SkyCoord(ra=agn_table2['RAJ2000']*u.degree, dec=agn_table2['DEJ2000']*u.degree)  
catalog2 = SkyCoord(ra=cat['_RA']*u.degree, dec=cat['_DE']*u.degree)  
idx2, d2d2, d3d2 = c2.match_to_catalog_sky(catalog2) 

t2 = Table([agn_table2['RAJ2000'], cat['_RA'][idx2], agn_table2['DEJ2000'], cat['_DE'][idx2]], 
	names=('agnRA', 'catRA', 'agnDE','catDE'), dtype=('f7', 'f7', 'f7','f7'))
#2712
print 'table 2'
print t2[0:16]
print t2[16:31]
print t2[31:46]
print t2[46:60]
print ''


print len(agn_table2['RAJ2000'])
print len(idx2)
print idx2[0:2]
print 'agn2 and cat'
print 'RA'
print agn_table2['RAJ2000'][0]
print cat['_RA'][18558]
print agn_table2['RAJ2000'][1]
print cat['_RA'][18323]
print ''
print 'DEC'
print agn_table2['DEJ2000'][0]
print cat['_DE'][18558]
print agn_table2['DEJ2000'][1]
print cat['_DE'][18323]
print ''


c3 = SkyCoord(ra=agn_table3['_RA']*u.degree, dec=agn_table3['_DE']*u.degree)  
catalog3 = SkyCoord(ra=cat['_RA']*u.degree, dec=cat['_DE']*u.degree)  
idx3, d2d3, d3d3 = c3.match_to_catalog_sky(catalog3) 

t3 = Table([agn_table3['_RA'], cat['_RA'][idx3], agn_table3['_DE'], cat['_DE'][idx3]], 
	names=('agnRA', 'catRA', 'agnDE','catDE'), dtype=('f7', 'f7', 'f7','f7'))
#308
print 'table 3'	
print t3[0:16]
print t3[16:31]
print t3[31:46]
print t3[46:60]
print ''

print len(agn_table3['_RA'])
print len(idx3)
print idx3[0:2]
print 'agn3 and cat'
print 'RA'
print agn_table3['_RA'][0]
print cat['_RA'][130432]
print agn_table3['_RA'][1]
print cat['_RA'][637053]
print ''
print 'DEC'
print agn_table3['_DE'][0]
print cat['_DE'][130432]
print agn_table3['_DE'][1]
print cat['_DE'][637053]
print ''

def create_bestfit():
    #this creates a best fit catalog where each galaxy is given a 2D or
    #1D fit based on the F-test and other statistics. If 1D is a better fit
    #then the values for disk or bulge are replaced with the 1D fit and the
    #other values are set to zero
    S11=create_S11_catalog()

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
'''
'''
print 'I_CLASS'
print type(galCl_table['I_CLASS'][0])
print ''
print 'objID'
print type(galCl_table['objID'][0])
print ''
print 'RA'
print type(galIn_table['RA'][0])
print ''
print 'DEC'
print type(galIn_table['DEC'][0])
print ''
print 'Z'
print type(galIn_table['Z'][0])
print ''
print 'OIII_FLUX'
print type(galLi_table['OIII_FLUX'][0])
print ''
print 'OIII_4363_FLUX'
print type(galLi_table['OIII_4363_FLUX'][0])
print ''
print 'OIII_4959_FLUX'
print type(galLi_table['OIII_4959_FLUX'][0])
print ''
print 'OIII_5007_FLUX'
print type(galLi_table['OIII_5007_FLUX'][0])
print ''
print 'NII_6548_FLUX'
print type(galLi_table['NII_6548_FLUX'][0])
print ''
print 'NII_6584_FLUX'
print type(galLi_table['NII_6584_FLUX'][0])
print ''
print 'H_GAMMA_FLUX'
print type(galLi_table['H_GAMMA_FLUX'][0])
print ''
print 'H_BETA_FLUX'
print type(galLi_table['H_BETA_FLUX'][0])
print ''
print 'gal AVG'
print type(galTo_table['AVG'][0])
print ''
print 'gal FLAG'
print type(galTo_table['FLAG'][0])
print ''
print 'SFR AVG'
print type(SFR_table['AVG'][0])
print ''
print 'SFR AVG'
print type(SFR_table['FLAG'][0])
print ''

print 'OIII_FLUX'
print galLi_table['OIII_FLUX'][0:20]
print ''
print 'OIII_4363_FLUX'
print galLi_table['OIII_4363_FLUX'][0:20]
print ''
print 'OIII_4959_FLUX'
print galLi_table['OIII_4959_FLUX'][0:20]
print ''
print 'OIII_5007_FLUX'
print galLi_table['OIII_5007_FLUX'][0:20]
print ''
print 'NII_6548_FLUX'
print galLi_table['NII_6548_FLUX'][0:20]
print ''
print 'NII_6584_FLUX'
print galLi_table['NII_6584_FLUX'][0:20]
print ''
print 'H_GAMMA_FLUX'
print galLi_table['H_GAMMA_FLUX'][0:20]
print ''
print 'H_BETA_FLUX'
print galLi_table['H_BETA_FLUX'][0:20]
print ''
print 'gal AVG'
print galTo_table['AVG'][0:20]
print ''
print 'gal FLAG'
print galTo_table['FLAG'][0:20]
print ''
print 'SFR AVG'
print SFR_table['AVG'][0:20]
print ''
print 'SFR AVG'
print SFR_table['FLAG'][0:20]
'''
#graphs
'''plt.hist(disk_cat['cosi'], bins=10, color='c', alpha=0.7,label='disk cat')
plt.title('Cosine Angle Frequnecy', y=1.05,fontsize=17, fontweight='bold')
plt.xlabel('Cosine Angle', fontsize=16, fontweight='bold')
plt.ylabel('Frequency', fontsize=16, fontweight='bold')
plt.axis([0,1.1,0,14000],fontsize=30, fontweight='bold')
plt.show()

plt.hist(disk_cat['OIII'], bins=60, color='k',linestyle='dotted',linewidth=3,histtype='step',normed=1,label='OIII values 0-100 range')
plt.hist(disk_cat['OIII'][edge], bins=60, align='right', color='r',histtype='stepfilled', alpha=0.8,normed=1,label='OIII edge-on')
plt.hist(disk_cat['OIII'][face], bins=60, align='right', color='b',histtype='stepfilled', alpha=0.7,normed=1,label='OIII face-on')
plt.axvline(disk_cat['OIII'][face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(disk_cat['OIII'][face]))
plt.axvline(disk_cat['OIII'][edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(disk_cat['OIII'][edge]))
plt.axvline(disk_cat['OIII'].mean(), color='k', linestyle='-', label='mean 0 < OIII < 100 values: %.2f'% np.mean(disk_cat['OIII']), linewidth=2)
plt.axis([-5,105,0,0.08],fontsize=40, fontweight='bold')
plt.title('Normalized OIII $\lambda$5007 Flux Frequency \n OIII value 0-100 range, B/T_r<0.35 and PpS<0.32', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel('OIII $\lambda$5007 Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.1, .5)
plt.legend(loc=0)
plt.show()

plt.hist(disk_cat['NII'], bins=60, color='k',linestyle='dotted',linewidth=3,histtype='step',normed=1,label='NII values 0-100 range')
plt.hist(disk_cat['NII'][edge], bins=60, align='right', color='r',histtype='stepfilled', alpha=0.8,normed=1,label='NII edge-on')
plt.hist(disk_cat['NII'][face], bins=60, align='right', color='b',histtype='stepfilled', alpha=0.7,normed=1,label='NII face-on')
plt.axvline(disk_cat['NII'][face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(disk_cat['NII'][face]))
plt.axvline(disk_cat['NII'][edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(disk_cat['NII'][edge]))
plt.axvline(disk_cat['NII'].mean(), color='k', linestyle='-', label='mean 0 < NII < 100 values: %.2f'% np.mean(disk_cat['NII']), linewidth=2)
plt.axis([-5,105,0,0.03],fontsize=40, fontweight='bold')
plt.title('NII $\lambda$6584 Flux Frequency \n NII value 0-100 range, B/T_r<0.35 and PpS<0.32', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel('NII $\lambda$6584 Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Normalized Frequency (counts)', fontsize=16, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=0)
plt.show()

plt.hist(disk_cat['H_a'], bins=20, color='k',linestyle='dotted',linewidth=3,histtype='step',normed=1,label='Halpha values 0-30 range')
plt.hist(disk_cat['H_a'][edge], bins=20, align='right', color='r',histtype='stepfilled', alpha=0.8,normed=1,label='Halpha edge-on')
plt.hist(disk_cat['H_a'][face], bins=20, align='right', color='b',histtype='stepfilled', alpha=0.7,normed=1,label='Halpha face-on')
plt.axvline(disk_cat['H_a'][face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(disk_cat['H_a'][face]))
plt.axvline(disk_cat['H_a'][edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(disk_cat['H_a'][edge]))
plt.axvline(disk_cat['H_a'].mean(), color='k', linestyle='-', label='mean 0 < Halpha < 30 values: %.2f'% np.mean(disk_cat['H_a']), linewidth=2)
#plt.axis([-1,31,0,0.1],fontsize=40, fontweight='bold')
plt.title('Normalized Halpha Flux Frequency \n Halpha value 0-100 range', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel('H alpha Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=0)
plt.show()

plt.hist(disk_cat['H_b'], bins=20, color='k',linestyle='dotted',linewidth=3,histtype='step',normed=1,label='Hbeta values 0-60 range')
plt.hist(disk_cat['H_b'][edge], bins=20, align='right', color='r',histtype='stepfilled', alpha=0.8,normed=1,label='Hbeta edge-on')
plt.hist(disk_cat['H_b'][face], bins=20, align='right', color='b',histtype='stepfilled', alpha=0.7,normed=1,label='Hbeta face-on')
plt.axvline(disk_cat['H_b'][face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(disk_cat['H_b'][face]))
plt.axvline(disk_cat['H_b'][edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(disk_cat['H_b'][edge]))
plt.axvline(disk_cat['H_b'].mean(), color='k', linestyle='-', label='mean 0 < Hbeta < 60 values: %.2f'% np.mean(disk_cat['H_b']), linewidth=2)
plt.axis([-1,62,0,0.05],fontsize=40, fontweight='bold')
plt.title('Normalized Hbeta Flux Frequency \n Hbeta value 0-60 range, B/T_r<0.35 and PpS<0.32', y=1.05,fontsize=18, fontweight='bold')
plt.xlabel('H beta Flux (erg/s/cm^2)', fontsize=16, fontweight='bold')
plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
plt.gca().xaxis.set_label_coords(.5, -.1)
plt.gca().yaxis.set_label_coords(-.12, .5)
plt.legend(loc=0)
plt.show()
'''
