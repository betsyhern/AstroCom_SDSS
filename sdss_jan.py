'''Inintialized 9/17/15'''

import os as os
import numpy as np
import astropy.io as aio
from astropy.table import Table,join,Column,unique
import astropy.constants as aconst
import astropy.cosmology as acosmo
import matplotlib.pyplot as plt

#The tables from Simard 2011 and Mendel 2014
fnS1='Simard_Table3.fit' #1D fits 
fnS2='Simard_Table1.fit' #2D fits with n_b = 4
fnM ='Mendel_dusty.fit'  #dusty masses combines Tables 3 and 4

#The following come from www.mpa-garching.mpg.de/SDSS/DR7/raw_data.html
fngal='gal_info_dr7_v5_2.fit'
fnline='gal_line_dr7_v5_2.fit'

#The following files come from home.strw.leidenuniv.nl/~jarle/SDSS/
fninfo='gal_info_dr7_v5_2.fit' #photometric properties
fnline='gal_line_dr7_v5_2.fit' #spectral line properties
fnkcor='gal_kcorrect_dr7_v5_2.fits'
fnsfr='gal_totsfr_dr7_v5_2.fits' #MPA DR7 galaxy sfr
fnclass='gal_iclass_table_dr7_v5_2.fits' #galaxy class and ids

#Tables from galaxyzoo.org
fnzphot='zoo2MainPhotoz.fits'
fnzspec='zoo2MainSpecz.fits'
fnzTable2='GalaxyZoo1_DR_table2.fits'
fnzTable3='GalaxyZoo1_DR_table3.fits'
fnzGZ='schawinski_GZ_2010_catalogue.fits' #has sfr and oIII

#Tables from code below merging, Brinchmann, Simard, Mendel and Schawinski (zoo)
fnfluxmass='flux_mass_catalog.fits' 
fnlinemass='line_mass_catalog.fits'
fnzoomassline='zoo_mass_line_catalog.fits'
 
#creates catalog from SDSS DR7 tables (Brinchmann method) 
# chooses the first objID is there are multiple ones
def create_line_catalog():
    hdulist1=aio.fits.open(fninfo)
    hdulist2=aio.fits.open(fnline)
    hdulist3=aio.fits.open(fnkcor)
    hdulist4=aio.fits.open(fnsfr)
    hdulist5=aio.fits.open(fnclass)    
    data1=hdulist1[1].data
    data2=hdulist2[1].data
    data3=hdulist3[1].data
    data4=hdulist4[1].data
    data5=hdulist5[1].data
    SDSS_table=Table([data5['OBJID'],
    data5['I_CLASS'],
    data4['AVG'],
    data4['FLAG'],
    data3['MODEL'],	#nans present in values
    data3['MODEL_ABSMAG'],	#nans present in values
    data3['MODEL_MASS'],
    data2['SII_6717_FLUX'],
    data2['SII_6717_FLUX_ERR'],
    data2['SII_6717_CONT'],
    data2['SII_6717_CONT_ERR'],
    data2['NII_6584_FLUX'],
    data2['NII_6584_FLUX_ERR'],
    data2['NII_6584_CONT'],
    data2['NII_6584_CONT_ERR'],
    data2['OI_6300_FLUX'],
    data2['OI_6300_FLUX_ERR'],
    data2['OI_6300_CONT'],
    data2['OI_6300_CONT_ERR'],
    data2['OII_FLUX'],
    data2['OII_FLUX_ERR'],
    data2['OIII_FLUX'],
    data2['OIII_FLUX_ERR'],
    data2['OIII_5007_FLUX'],
    data2['OIII_5007_FLUX_ERR'],
    data2['OIII_5007_CONT'],
    data2['OIII_5007_CONT_ERR'],
    data2['H_ALPHA_FLUX'],
    data2['H_ALPHA_FLUX_ERR'],
    data2['H_ALPHA_CONT'],
    data2['H_ALPHA_CONT_ERR'],
    data2['H_BETA_FLUX'],
    data2['H_BETA_FLUX_ERR'],
    data2['H_BETA_CONT'],
    data2['H_BETA_CONT_ERR'],
    data2['H_GAMMA_FLUX'],
    data2['H_GAMMA_FLUX_ERR'],
    data2['H_GAMMA_CONT'],
    data2['H_GAMMA_CONT_ERR'],
    data2['H_DELTA_FLUX'],
    data2['H_DELTA_FLUX_ERR'],
    data2['H_DELTA_CONT'],
    data2['H_DELTA_CONT_ERR'],
    data1['RA'],
    data1['DEC'],
    data1['z']],
    names=('objID','iclass','sfr','sfr_flag','model','model_absmag','mass','SII','SII_e',
    	'SII_cont','SII_econt','NII','NII_e','NII_cont','NII_econt','OI','OI_e','OI_cont',
    	'OI_econt','OII','OII_e','OIII_tot','OIII_tot_e','OIII','OIII_e','OIII_cont',
    	'OIII_econt','H_alpha','H_alpha_e','H_alpha_cont','H_alpha_econt','H_beta','H_beta_e',
    	'H_beta_cont','H_beta_econt','H_gamma','H_gamma_e','H_gamma_cont','H_gamma_econt',
    	'H_delta','H_delta_e','H_delta_cont','H_delta_econt','RA','DEC','z'))
    t=unique(SDSS_table,keys='objID')
    #print data3.names
    #print t.columns
    return t

# adds to the SDSS DR7 (Brinchmann method) catalog. 
# selects positive flux and adds signal to noise variable
def create_line_SN_catalog():
	t=create_line_catalog()
	#select only positive flux values
	cut=np.logical_and(t['H_alpha']>0,t['H_beta']>0)
	t=t[cut]
	cut=np.logical_and(t['H_gamma']>0,t['H_delta']>0)
	t=t[cut]
	cut=np.logical_and(t['SII']>0,t['NII']>0)
	t=t[cut]
	cut=np.logical_and(t['OI']>0,t['OIII']>0)
	t=t[cut]
	cut=np.logical_and(t['SII_econt']>0,t['NII_econt']>0)
	t=t[cut]
	cut=np.logical_and(t['OI_econt']>0,t['H_alpha_econt']>0)
	t=t[cut]
	cut=t['OII']>0
	t=t[cut]
	#create signal to noise table and merge with positive flux table
	SII_SN=Column(name='SII_SN', data=t['SII_cont']/t['SII_econt'])
	NII_SN=Column(name='NII_SN', data=t['NII_cont']/t['NII_econt'])
	OI_SN=Column(name='OI_SN', data=t['OI_cont']/t['OI_econt'])
	OIII_SN=Column(name='OIII_SN', data=t['OIII_cont']/t['OIII_econt'])
	OII_SN=Column(name='OII_SN', data=t['OII']/t['OII_e'])
	OIII_tot_SN=Column(name='OIII_tot_SN', data=t['OIII_tot']/t['OIII_tot_e'])
	HA_SN=Column(name='HA_SN', data=t['H_alpha_cont']/t['H_alpha_econt'])
	HB_SN=Column(name='HB_SN', data=t['H_beta_cont']/t['H_beta_econt'])
	HC_SN=Column(name='HC_SN', data=t['H_gamma_cont']/t['H_gamma_econt'])
	HD_SN=Column(name='HD_SN', data=t['H_delta_cont']/t['H_delta_econt'])
	t.add_column(SII_SN)
	t.add_column(NII_SN)
	t.add_column(OI_SN)
	t.add_column(OIII_SN)
	t.add_column(OII_SN)
	t.add_column(OIII_tot_SN)
	t.add_column(HA_SN)
	t.add_column(HB_SN)
	t.add_column(HC_SN)
	t.add_column(HD_SN)
	return t

# creates Simard catalog from table 1 and table 3 data
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
    #include cosine angle ('cosi'). rename inclination 'i' to 'degi' and 'e_i' to 'e_degi'
    col_cosi=Column(name='cosi',data=np.cos(np.radians(t['i'])))
    t.add_column(col_cosi)
    t.rename_column('i','degi')
    t.rename_column('e_i','e_degi')
    #rename columns with same name
    # Simard_table1 columns will have S1 at the end of column name (replaces _2)
    # Simard_table3 columns will have S3 at the end of column name (replaces _1)
    t.rename_column('gg2d_2','gg2d_S1')
    t.rename_column('gg2d_1','gg2d_S3')
    t.rename_column('rg2d_2','rg2d_S1')
    t.rename_column('rg2d_1','rg2d_S3')
    t.rename_column('Rhlg_2','Rhlg_S1')
    t.rename_column('Rhlg_1','Rhlg_S3')
    t.rename_column('e_2','e_S1')
    t.rename_column('e_1','e_S3')
    t.rename_column('S2r_2','S2r_S1')
    t.rename_column('S2r_1','S2r_S3')
    t.rename_column('ggMag_2','ggMag_S1')
    t.rename_column('ggMag_1','ggMag_S3')
    t.rename_column('rgMag_2','rgMag_S1')
    t.rename_column('rgMag_1','rgMag_S3')
    t.rename_column('S2g_2','S2g_S1')
    t.rename_column('S2g_1','S2g_S3')
    t.rename_column('__B_T_g','B/T_g')
    t.rename_column('__B_T_r','B/T_r')
    #convert to solar luminosities ugriz_{0.1}=6.80,5.45,4.76,4.58,4.51
    #Blanton 03, but from http://mips.as.arizona.edu/~cnaw/sun.html
    Lsun_g=5.45
    Lsun_r=4.76
    col_Ldg=Column(name='L_dg',data=10**(-0.4*(t['gdMag']-Lsun_g)))
    col_Lbg=Column(name='L_bg',data=10**(-0.4*(t['gbMag']-Lsun_g)))
    col_Lg_S3=Column(name='L_g_S3',data=10**(-0.4*(t['ggMag_S3']-Lsun_g)))
    col_Lg_S1=Column(name='L_g_S1',data=10**(-0.4*(t['ggMag_S1']-Lsun_g)))
    t.add_column(col_Ldg)
    t.add_column(col_Lbg)
    t.add_column(col_Lg_S3)
    t.add_column(col_Lg_S1)
    col_Ldr=Column(name='L_dr',data=10**(-0.4*(t['rdMag']-Lsun_r)))
    col_Lbr=Column(name='L_br',data=10**(-0.4*(t['rbMag']-Lsun_r)))
    col_Lr_S3=Column(name='L_r_S3',data=10**(-0.4*(t['rgMag_S3']-Lsun_r)))
    col_Lr_S1=Column(name='L_r_S1',data=10**(-0.4*(t['rgMag_S1']-Lsun_r)))
    t.add_column(col_Ldr)
    t.add_column(col_Lbr)
    t.add_column(col_Lr_S3)
    t.add_column(col_Lr_S1)
    return t

def create_bestfit():
    #this creates a best fit catalog where each galaxy is given a 2D or
    #1D fit based on the F-test and other statistics. If 1D is a better fit
    #then the values for disk or bulge are replaced with the 1D fit and the
    #other values are set to zero
    S11=create_S11_catalog()

# prints 'GalaxyZoo1_DR_table2.fits' column names 
def zoo_tables():
    hdulist1=aio.fits.open(fnzTable2)
    hdulist2=aio.fits.open(fnzTable3)
    hdulist3=aio.fits.open(fnzphot)
    hdulist4=aio.fits.open(fnzspec)
    hdulist5=aio.fits.open(fnzGZ)
    table1=hdulist1[1].data
    table2=hdulist2[1].data
    table3=hdulist3[1].data
    table4=hdulist4[1].data
    table5=hdulist5[1].data
    print table1.names
    print table1[0]
    return "ok"

# merges the Simard and Mendel catalogs
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
    #remove error ranges on masses
    t.remove_column('b_logM')
    t.remove_column('B_logM1')
    t.remove_column('b_logMt')
    t.remove_column('B_logMt1')
    t.remove_column('b_logMb')
    t.remove_column('B_logMb1')
    t.remove_column('b_logMd')
    t.remove_column('B_logMd1')
    #mass to light ratios
    MLg_S3=Column(name='M/L_g_S3',data=(10**t['logM'])/t['L_g_S3'])
    MLg_S1=Column(name='M/L_g_S1',data=(10**t['logMt'])/t['L_g_S1'])
    MLgd=Column(name='M/L_dg',data=(10**t['logMd'])/t['L_dg'])
    MLgb=Column(name='M/L_bg',data=(10**t['logMb'])/t['L_bg'])
    t.add_column(MLg_S3)
    t.add_column(MLg_S1)
    t.add_column(MLgd)
    t.add_column(MLgb)
    MLr_S3=Column(name='M/L_r_S3',data=(10**t['logM'])/t['L_r_S3'])
    MLr_S1=Column(name='M/L_r_S1',data=(10**t['logMt'])/t['L_r_S1'])
    MLrd=Column(name='M/L_dr',data=(10**t['logMd'])/t['L_dr'])
    MLrb=Column(name='M/L_br',data=(10**t['logMb'])/t['L_br'])
    t.add_column(MLr_S3)
    t.add_column(MLr_S1)
    t.add_column(MLrd)
    t.add_column(MLrb)
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
    t['degi'][pure_bulge]=-90.
    t['cosi'][pure_bulge]=-1.0
    t['gdMag'][pure_bulge]=-99.9
    t['rdMag'][pure_bulge]=-99.9
    t['L_dg'][pure_bulge]=0.0
    t['L_dr'][pure_bulge]=0.0
    t['M/L_dg'][pure_bulge]=0.0
    t['M/L_dr'][pure_bulge]=0.0
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
    t['L_bg'][pure_disk]=0.0
    t['L_br'][pure_disk]=0.0
    t['M/L_bg'][pure_disk]=0.0
    t['M/L_br'][pure_disk]=0.0
    #print t['z'][0:10]
    return t

# merges the Simard/Mendel catalog with SDSS DR7 (Brinchmann method) catalog
'''creates a fits file for this catalog'''
def create_line_mass_catalog():
    S_M_table=create_mass_catalog()
    B_table=create_line_catalog()
    shared=['objID']
    t=join(S_M_table,B_table,keys=shared)
    t.remove_column('z_1') #nearly the same as Brinchmann z
    t.rename_column('z_2','z')
    t.remove_column('_RA') #same as RA from Brinchmann
    t.remove_column('_DE') #same as RA from Brinchmann
    tab=Table(t)
    tab.write('line_mass_catalog.fits',format='fits')
    return t

# opens the Simard/Mendel/Brinchmann catalog fits file ...
# (no positive flux nor signal to noise variable)
def open_line_mass_catalog():
	hdulist=aio.fits.open(fnlinemass)
	t=hdulist[1].data
	return t

# merges the Simard/Mendel catalog with ...
# SDSS DR7 (Brinchmann method) catalog possessing positive flux and signal to noise
'''creates a fits file for this catalog'''
def create_flux_mass_catalog():
    S_M_table=create_mass_catalog()
    B_table=create_line_SN_catalog()
    shared=['objID']
    t=join(S_M_table,B_table,keys=shared)
    t.remove_column('z_1') #nearly the same as Brinchmann z
    t.rename_column('z_2','z')
    t.remove_column('_RA') #same as RA from Brinchmann
    t.remove_column('_DE') #same as RA from Brinchmann
    tab=Table(t)
    tab.write('flux_mass_catalog.fits',format='fits')
    return t

# opens the Simard/Mendel/Brinchmann catalog fits file ...
# that includes positive flux and signal to noise
def open_flux_mass_catalog():
	hdulist=aio.fits.open(fnfluxmass)
	t=hdulist[1].data
	return t

# merges the Simard/Mendel catalog zoo table 2
def create_zoo_catalog():
    mass_table=create_mass_catalog()
    mass_table.rename_column('objID','OBJID')
    hdulist=aio.fits.open(fnzTable2)
    Z_table= hdulist[1].data
    shared=['OBJID']
    t=join(mass_table,Z_table,keys=shared)
    t.rename_column('OBJID','objID')
    t.remove_column('RA')
    t.remove_column('DEC')
    return t

# merges the Simard/Mendel catalog with ...
# Schawinski (zoo) which possess AGN mass and OIII luminosities
def create_zoospec_catalog():
    mass_table=create_mass_catalog()
    mass_table.rename_column('objID','OBJID')
    hdulist=aio.fits.open(fnzGZ)
    Z_table= (hdulist[1].data)
    #OBJID is string must convert to ints
    ids=np.array(Z_table['OBJID'][0],dtype='int')
    print len(ids)
    bptclass=np.array(Z_table['BPT_CLASS'][0])
    logmstellar=np.array(Z_table['LOG_MSTELLAR'][0])
    LO3=np.array(Z_table['L_O3'][0])
    tmpT=Table([ids,bptclass,logmstellar,LO3]
               ,names=['OBJID','BPT_CLASS','LOG_MSTELLAR','L_O3'])
    shared=['OBJID']
    t=join(mass_table,tmpT,keys=shared)
    t.rename_column('OBJID','objID')
    return t
create_zoospec_catalog()
# merges the Simard/Mendel/Schawinski (zoo) catalog with ...
# SDSS DR7 tables (Brinchmann method) that has no positive flux nor signal to noise
'''creates a fits file for this catalog'''
def create_zoo_mass_line_catalog():
    zoo_table=create_zoospec_catalog()
    line_table=create_line_catalog()
    shared=['objID']
    t=join(zoo_table,line_table,keys=shared)
    t.remove_column('z_1') #nearly the same as Brinchmann z
    t.rename_column('z_2','z')
    t.remove_column('_RA') #same as RA from Brinchmann
    t.remove_column('_DE') #same as RA from Brinchmann
    tab=Table(t)
    tab.write('zoo_mass_line_catalog.fits',format='fits')
    return t

# opens the Simard/Mendel/Brinchmann/Schawinski (zoo) catalog ...
# fits file (no positive flux nor signal to noise variable)
def open_zoo_mass_line_catalog():
	hdulist=aio.fits.open(fnzoomassline)
	t=hdulist[1].data
	return t

def create_disk_catalog():
    t=create_mass_catalog()
    disks=np.logical_and(t['B/T_r'] < 0.35,t['PpS']<0.32)
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
    #This just prints out the weird stuff for identification
    print "There are some cases with m and z but no M"
    w=np.where(np.logical_or(np.isfinite(t['ggMag_S3'])==False,
                             np.isfinite(t['ggMag_S1'])==False))
    print 'objID','z','gg2d_1D','ggMag_1D','rg2d_1D','rgMag_1D','logM','logMt'
    for i in w[0]:
        print t['objID'][i],t['z'][i],t['gg2d_S3'][i],t['ggMag_S3'][i] \
            ,t['rg2d_S3'][i],t['rgMag_S3'][i],t['logM'][i],t['logMt'][i] \
            ,t['gg2d_S1'][i],t['ggMag_S1'][i],t['rgMag_S1'][i]
    print "There are cases of 100% bulge in one band and 100% disk in the other"
    w1=np.where(np.logical_and(t['__B_T_r']==0.0,t['__B_T_g']==1.0))
    w2=np.where(np.logical_and(t['__B_T_r']==1.0,t['__B_T_g']==0.0))

# prints min, max and number of nans for tables by name 
def check_names(t):
    names=t.names
    for name in names:
        max=np.amax(t[name])
        min=np.amin(t[name])
        nans=np.isnan(t[name])
        sum=np.nansum(nans)
        print name+" Min: {0}   Max: {1}    Nans found: {2}".format(min,max,sum)

# prints min, max and number of nans for tables by column header
def check_column(t):
    names=t.colnames
    for name in names:
        max=np.amax(t[name])
        min=np.amin(t[name])
        nans=np.isnan(t[name])
        sum=np.nansum(nans)
        print name+" Min: {0}   Max: {1}    Nans found: {2}".format(min,max,sum)

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
    color1=t['ggMag_S3']-t['rgMag_S3']
    color2=t['ggMag_S1']-t['rgMag_S1']
    edgeon=t['cosi'] < 0.3
    faceon=t['cosi'] > 0.85
    plt.figure()
    n,b,p=plt.hist([color2[faceon],color2[edgeon]],30,histtype='step'
         ,range=[-0.5,1.5],color=['blue','red'],label=['face-on','edge-on'])
    plt.xlabel('g-r')
    plt.legend()
    plt.show()
