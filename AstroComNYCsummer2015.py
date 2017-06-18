'''Initialized 6/14/15'''

import os as os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.constants as const
import astropy.cosmology as cosmo
import matplotlib.pyplot as plt
from decimal import *
from catalog_sdss_condensed import *

#returns new dictionary with applied True/False index
def slice_dict(dict,index):
	sort_dict={}
	for key in range(0,len(dict)):
		sort_key=dict.keys()[key]
		sort_value=dict.values()[key][index]	
		sort_dict[sort_key]=sort_value
	return sort_dict

# creates arrays that exclues 2 exact values
def exclude_val(item,val1,val2):
	index=[]
	for i in range(0,len(item)):
		if item[i]==val1 or item[i]==val2:
			index.append(False)
		else:
			index.append(True)
	arr=np.array(index)
	return arr

# creates dictionary using exclude_val function to exclues 2 exact values
def no_extremes(dict,first,second,val1,val2):
	index1=exclude_val(first,val1,val2)
	index2=exclude_val(second,val1,val2)
	cut=np.logical_and(index1,index2)
	cut_dict={}
	for key in range(0,len(dict)):
		cut_key=dict.keys()[key]
		cut_value=dict.values()[key][cut]	
		cut_dict[cut_key]=cut_value
	return cut_dict
	
# shows the min and max values in the keys of dictionaries. if nans are present it would show up as either min and or nan
def len_min_max(dict):
	for i in range(0,len(dict)):
		print dict.keys()[i]
		print 'length: '+str(len(dict.values()[i]))
		print 'minimum: '+str(np.amin(dict.values()[i]))
		print 'maximum: '+str(np.amax(dict.values()[i]))
		print ' '

#displays sum of inf in keys of dictionary
def inf_total(dict):
	for i in range(0,len(dict)):
		print dict.keys()[i]
		print 'total inf'
		print np.nansum(np.isinf(dict.values()[i]))
		print ' '

#displays sum of nans in keys of dictionary
def nan_total(dict):
	for i in range(0,len(dict)):
		print dict.keys()[i]
		print 'total nans'
		print np.nansum(np.isnan(dict.values()[i]))
		print ' '

# creates a scatter plot with a trendline
def trendplot(X,Y,color=None,marker=None,labelplot=None,labeltrend=None,title=None,xtitle=None,ytitle=None,loc=0,axis=None):
	coeffs=np.polyfit(X,Y,1)
	trendline=np.poly1d(coeffs)
	correlation=np.corrcoef(X,Y)[0,1]
	correlationSq=correlation**2
	results={'polynomial':coeffs.tolist(),'correlation':correlation,'determination':correlationSq}
	
	plt.scatter(X,Y,color=color, marker=marker, label=labelplot)
	plt.plot(X,trendline(X),'k-',label=labeltrend+'\ny='+str(trendline)+'\n$R^2$= '+str(round(results['determination'],3)))
	plt.title(title, fontsize=22, fontweight='bold', y=1.05)
	plt.ylabel(ytitle, fontsize=19, fontweight='bold')
	plt.xlabel(xtitle, fontsize=19, fontweight='bold')
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.legend(loc=loc)
	plt.show()
		
#create index to focus on Type 3 galaxies
Type=cat['Type']
Type3_index=Type==3 

# dictionaries of inclination, magnitude, B/T ratio & log masses from the catalog in the catalog_sdss_condensed module
i={'i':cat['i']}
Mag={'ggMag':cat['ggMag'],'rgMag':cat['rgMag']}
B_T={'B/T_g':cat['B/T_g'],'B/T_r':cat['B/T_r']}
logMass={'logM_dt':cat['logM_dt'],'logMb+d_dt':cat['logMb+d_dt'],'logMb_dt':cat['logMb_dt'],'logMd_dt':cat['logMd_dt'], 
		 'logM_df':cat['logM_df'],'logMb+d_df':cat['logMb+d_df'],'logMb_df':cat['logMb_df'],'logMd_df':cat['logMd_df']} 

# new dictionaries with only Type 3 galaxies
i_type3=slice_dict(i,Type3_index)
Mag_type3=slice_dict(Mag,Type3_index)
B_T_type3=slice_dict(B_T,Type3_index)
logMass_type3=slice_dict(logMass,Type3_index)

# dictionaries of magnitude, B/T ratio & log masses from the catalog using index cut removing B/T=0 or 1.
cut_i=no_extremes(i_type3,B_T_type3['B/T_g'],B_T_type3['B/T_r'],0,1)
cut_Mag=no_extremes(Mag_type3,B_T_type3['B/T_g'],B_T_type3['B/T_r'],0,1)
cut_B_T=no_extremes(B_T_type3,B_T_type3['B/T_g'],B_T_type3['B/T_r'],0,1)
cut_logMass=no_extremes(logMass_type3,B_T_type3['B/T_g'],B_T_type3['B/T_r'],0,1)

# Dictionary solar luminosities of green: ggSM (total), ggSM_b (bulge), ggSM_d (disk)
# and red: rgSM (total), rgSM_b (bulge), rgSM_d (disk)
SM_lum={'ggSM':10**(-0.4*(cut_Mag['ggMag']-5.45)),'ggSM_b':(10**(-0.4*(cut_Mag['ggMag']-5.45)))*cut_B_T['B/T_g'],
	    'ggSM_d':(10**(-0.4*(cut_Mag['ggMag']-5.45)))*(1-cut_B_T['B/T_g']),'rgSM':10**(-0.4*(cut_Mag['rgMag']-4.76)),
	    'rgSM_b':(10**(-0.4*(cut_Mag['rgMag']-4.76)))*cut_B_T['B/T_r'],
	    'rgSM_d':(10**(-0.4*(cut_Mag['rgMag']-4.76)))*(1-cut_B_T['B/T_r'])}

# List of solar luminosity values
SM_lumg_list=[SM_lum['ggSM'],SM_lum['ggSM'],SM_lum['ggSM_b'],SM_lum['ggSM_d']]
SM_lumr_list=[SM_lum['rgSM'],SM_lum['rgSM'],SM_lum['rgSM_b'],SM_lum['rgSM_d']]

# List of modified solar luminosity values
SM_lumg_list_mod=[SM_lum['ggSM']*1.5,SM_lum['ggSM']*1.77,SM_lum['ggSM_b']*1.5,SM_lum['ggSM_d']*1.5]
SM_lumr_list_mod=[SM_lum['rgSM']*1.32,SM_lum['rgSM']*1.55,SM_lum['rgSM_b']*1.5,SM_lum['rgSM_d']*1.5]

# List of values in the log mass dictionary
logMass_dt_list=[cut_logMass['logM_dt'],cut_logMass['logMb+d_dt'],cut_logMass['logMb_dt'],cut_logMass['logMd_dt']]
logMass_df_list=[cut_logMass['logM_df'],cut_logMass['logMb+d_df'],cut_logMass['logMb_df'],cut_logMass['logMd_df']]
 
# List of mass values converted from log
masscon_dt_list=[]
for i in range(0,len(logMass_dt_list)):
	value=10**logMass_dt_list[i]
	masscon_dt_list.append(value)

masscon_df_list=[]
for i in range(0,len(logMass_df_list)):
	value=10**logMass_df_list[i]
	masscon_df_list.append(value)

# List of mass-light ratios values
m_l_val=[]
for i in range(0,4):
	m_l_val.append(masscon_dt_list[i]/SM_lumg_list[i])	
for i in range(0,4):
	m_l_val.append(masscon_df_list[i]/SM_lumg_list[i])
for i in range(0,4):
	m_l_val.append(masscon_dt_list[i]/SM_lumr_list[i])
for i in range(0,4):
	m_l_val.append(masscon_df_list[i]/SM_lumr_list[i])

# List of modified mass-light ratios values
m_l_val_mod=[]
for i in range(0,4):
	m_l_val_mod.append(masscon_dt_list[i]/SM_lumg_list_mod[i])	
for i in range(0,4):
	m_l_val_mod.append(masscon_df_list[i]/SM_lumg_list_mod[i])
for i in range(0,4):
	m_l_val_mod.append(masscon_dt_list[i]/SM_lumr_list_mod[i])
for i in range(0,4):
	m_l_val_mod.append(masscon_df_list[i]/SM_lumr_list_mod[i])
	

# List of mass-light ratio names
m_l_name=['M_dt/ggSM','Mb+d_dt/ggSM','Mb_dt/ggSM_b','Md_dt/ggSM_d','M_df/ggSM','Mb+d_df/ggSM','Mb_df/ggSM_b','Md_df/ggSM_d',
	'M_dt/rgSM','Mb+d_dt/rgSM','Mb_dt/rgSM_b','Md_dt/rgSM_d','M_df/rgSM','Mb+d_df/rgSM','Mb_df/rgSM_b','Md_df/rgSM_d']

# Dictionary composed of mass-light ratios values with mass-light ratio names
m_l={}
if len(m_l_val)==len(m_l_name):
	for i in range(0,len(m_l_name)):
		key=m_l_name[i]
		value=m_l_val[i]
		m_l[key]=value

# Dictionary composed of mass-light ratios values with mass-light ratio names
m_l_mod={}
if len(m_l_val_mod)==len(m_l_name):
	for i in range(0,len(m_l_name)):
		key=m_l_name[i]
		value=m_l_val_mod[i]
		m_l_mod[key]=value

# cuts to focus on disk and bulge galaxies
disk_r=cut_B_T['B/T_r'] < 0.12
bulge_r=cut_B_T['B/T_r'] > 0.71
disk_g=cut_B_T['B/T_g'] < 0.05
bulge_g=cut_B_T['B/T_g'] > 0.67

#cut to choose disk galaxies in green and red band
BT_disk_type3=np.logical_and(cut_B_T['B/T_g']<=0.3,cut_B_T['B/T_r']<=0.3)

# command for cut of edge-on galaxies with B/T 0.01-0.3 and face-on 0.83-0.99
edge=cut_i['i'][BT_disk_type3]<0.3
middle=np.logical_and(cut_i['i'][BT_disk_type3]>=0.3,cut_i['i'][BT_disk_type3]<=0.81)
face=cut_i['i'][BT_disk_type3]>0.81

#command to create the g-r magnitude array of B/T 0-0.3
g_r=np.subtract(cut_Mag['ggMag'],cut_Mag['rgMag'])[BT_disk_type3]

# Four subplots sharing both x/y axes
def multi_hist(data1, data1_edge, data2, data2_edge, data3, data3_edge, data4, data4_edge, edge, face, A=20, title1=None, title2=None, title3=None, title4=None):
	f, ax = plt.subplots(2, 2,sharex=True,sharey=True)
	ax[0,0].hist(data1[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	ax[0,0].hist(data1_edge[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	ax[0,0].hist(data1, bins=A, range=(0,4), color='k', histtype='step', linestyle='dotted', linewidth=3, label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	ax[0,0].axvline(data1[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data1[face]))
	ax[0,0].axvline(data1_edge[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data1_edge[edge]))
	ax[0,0].axvline(data1.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data1), linewidth=2)
	ax[0,0].legend(loc=0)
	
	ax[0,1].hist(data2[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	ax[0,1].hist(data2_edge[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	ax[0,1].hist(data2, bins=A, range=(0,4), color='k', histtype='step', linestyle='dotted', linewidth=3, label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	ax[0,1].axvline(data2[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data2[face]))
	ax[0,1].axvline(data2_edge[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data2_edge[edge]))
	ax[0,1].axvline(data2.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data2), linewidth=2)
	ax[0,1].legend(loc=0)
	
	ax[1,0].hist(data3[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	ax[1,0].hist(data3_edge[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	ax[1,0].hist(data3, bins=A, range=(0,4), color='k', histtype='step', linestyle='dotted', linewidth=3, label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	ax[1,0].axvline(data3[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data3[face]))
	ax[1,0].axvline(data3_edge[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data3_edge[edge]))
	ax[1,0].axvline(data3.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data3), linewidth=2)
	ax[1,0].legend(loc=0)
	
	ax[1,1].hist(data4[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	ax[1,1].hist(data4_edge[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	ax[1,1].hist(data4, bins=A, range=(0,4), color='k', histtype='step', linestyle='dotted', linewidth=3, label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	ax[1,1].axvline(data4[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data4[face]))
	ax[1,1].axvline(data4_edge[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data4_edge[edge]))
	ax[1,1].axvline(data4.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data4), linewidth=2)
	ax[1,1].legend(loc=0)
	
	ax[0,0].set_title(title1, fontsize=14)
	ax[0,1].set_title(title2, fontsize=14)
	ax[1,0].set_title(title3, fontsize=14)
	ax[1,1].set_title(title4, fontsize=14)
	# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
	plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
	ax[1,0].set_xlabel('Mass/Light Ratio (M$\odot$/L$\odot$)',fontsize=16, x=1)
	ax[0,0].set_ylabel('Frequency (counts)',fontsize=16)
	ax[1,0].set_ylabel('Frequency (counts)',fontsize=16)
	plt.show()

'''multi_hist(m_l['M_dt/rgSM'][BT_disk_type3],m_l['Mb+d_dt/rgSM'][BT_disk_type3],m_l['M_dt/ggSM'][BT_disk_type3],
	m_l['Mb+d_dt/ggSM'][BT_disk_type3],edge,face, title1='M_dt/rgSM Normalized Frequency of M/L ratio (M$\odot$/L$\odot$) with B/T<0.3',
	title2='Mb+d_dt/rgSM Normalized Frequency of M/L ratio (M$\odot$/L$\odot$) with B/T<0.3', title3='M_dt/ggSM Normalized Frequency of M/L ratio (M$\odot$/L$\odot$) with B/T<0.3',
	title4='Mb+d_dt/ggSM Normalized Frequency of M/L ratio (M$\odot$/L$\odot$) with B/T<0.3')

multi_hist(m_l['M_dt/rgSM'][BT_disk_type3],m_l_mod['M_dt/rgSM'][BT_disk_type3],m_l['Mb+d_dt/rgSM'][BT_disk_type3],m_l_mod['Mb+d_dt/rgSM'][BT_disk_type3],m_l['M_dt/ggSM'][BT_disk_type3],
	m_l_mod['M_dt/ggSM'][BT_disk_type3],m_l['Mb+d_dt/ggSM'][BT_disk_type3],m_l_mod['Mb+d_dt/ggSM'][BT_disk_type3],edge,face, title1='Modified M_dt/rgSM Normalized Frequency\nof M/L ratio (M$\odot$/L$\odot$) with B/T<0.3',
	title2='Modified Mb+d_dt/rgSM Normalized Frequency\nof M/L ratio (M$\odot$/L$\odot$) with B/T<0.3', title3='Modified M_dt/ggSM Normalized Frequency\nof M/L ratio (M$\odot$/L$\odot$) with B/T<0.3',
	title4='Modified Mb+d_dt/ggSM Normalized Frequency\nof M/L ratio (M$\odot$/L$\odot$) with B/T<0.3')'''

#histogram of mass/light ratio  and modified mass/light ratio of disc 0.01-0.3 cut, edge-on(cosine<0.3) 
def m_l_hist_edge(data1, data2, edge, A=20, title=None):
	plt.hist(data1[edge],bins=A, range=(0,4),color='y',  histtype='stepfilled', alpha=0.7, label='edge-on M/L ratio', normed=1)
	plt.hist(data2[edge], bins=A, range=(0,4), color='m',  histtype='stepfilled', alpha=0.8, label='edge-on modified M/L ratio', normed=1)
	plt.axvline(data1[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data1[edge]))
	plt.axvline(data2[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data2[edge]))
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.legend(loc=0)
	plt.xlabel('Mass/Light Ratio (M$\odot$/L$\odot$)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.show()

m_l_hist_edge(m_l['M_dt/rgSM'][BT_disk_type3],m_l_mod['M_dt/rgSM'][BT_disk_type3],edge,title='M_dt/rgSM Normalized Frequency of M/L ratio and\nmodified M/L ratio (M$\odot$/L$\odot$) with B/T<0.3')