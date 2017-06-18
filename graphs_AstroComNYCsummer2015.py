'''Initialized 6/24/15'''
import os as os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.constants as const
import astropy.cosmology as cosmo
import matplotlib.pyplot as plt
from decimal import *
from AstroComNYCsummer2015 import *

#python 2.7.10 anaconda 2.1.0 matplotlib 1.4.3

#histograms of mass-light ratios for dust and dust free with green and red bands
def m_l_hist(Type=None):
	if Type=='M_dt_g':
		plt.hist(m_l['M_dt/ggSM'], 60, range=(0,6), align='right', color='g')
		plt.title('M_dt/ggSM frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_dt_g':
		plt.hist(m_l['Mb+d_dt/ggSM'], 60, range=(0,6), align='right', color='g')
		plt.title('Mb+d_dt/ggSM frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_dt_g':
		plt.hist(m_l['Mb_dt/ggSM_b'], 60, range=(0,6), align='right', color='g')
		plt.title('Mb_dt/ggSM_b frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_g':
		plt.hist(m_l['Md_dt/ggSM_d'], 60, range=(0,6), align='right', color='g')
		plt.title('Md_dt/ggSM_d frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='M_dt_r':
		plt.hist(m_l['M_dt/rgSM'], 60, range=(0,6), align='right', color='r')
		plt.title('M_dt/rgSM frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_dt_r':
		plt.hist(m_l['Mb+d_dt/rgSM'], 60, range=(0,6), align='right', color='r')
		plt.title('Mb+d_dt/rgSM frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_dt_r':
		plt.hist(m_l['Mb_dt/rgSM_b'], 60, range=(0,6), align='right', color='r')
		plt.title('Mb_dt/rgSM_b frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_r':
		plt.hist(m_l['Md_dt/rgSM_d'], 60, range=(0,6), align='right', color='r')
		plt.title('Md_dt/rgSM_d frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='M_df_g':
		plt.hist(m_l['M_df/ggSM'], 60, range=(0,6), align='right', color='g')
		plt.title('M_df/ggSM frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_df_g':
		plt.hist(m_l['Mb+d_df/ggSM'], 60, range=(0,6), align='right', color='g')
		plt.title('Mb+d_df/ggSM frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_df_g':	
		plt.hist(m_l['Mb_df/ggSM_b'], 60, range=(0,6), align='right', color='g')
		plt.title('Mb_df/ggSM_b frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_df_g':	
		plt.hist(m_l['Md_df/ggSM_d'], 60, range=(0,6), align='right', color='g')
		plt.title('Md_df/ggSM_d frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='M_df_r':	
		plt.hist(m_l['M_df/rgSM'], 60, range=(0,6), align='right', color='r')
		plt.title('M_df/rgSM frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_df_r':
		plt.hist(m_l['Mb+d_df/rgSM'], 60, range=(0,6), align='right', color='r')
		plt.title('Mb+d_df/rgSM frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_df_r':	
		plt.hist(m_l['Mb_df/rgSM_b'], 60, range=(0,6), align='right', color='r')
		plt.title('Mb_df/rgSM_b frequency', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_df_r':	
		plt.hist(m_l['Md_df/rgSM_d'], 60, range=(0,6), align='right', color='r')
		plt.title('Md_df/rgSM_d frequency', y=1.05,fontsize=17, fontweight='bold')
	plt.ylabel('Frequency', fontsize=16, fontweight='bold')
	plt.xlabel('mass light ratio', fontsize=16, fontweight='bold')
	plt.show()


# histograms mass-light ratio of dust and dust free using red and green bands focusing on galaxies with prominent disks and prominent bulges
def m_l_bd(Type=None):
	if Type=='M_dt_r':
		plt.hist(m_l['M_dt/rgSM'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['M_dt/rgSM'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['M_dt/rgSM'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.2: %.2f'% np.mean(m_l['M_dt/rgSM'][disk_r]))
		plt.axvline(m_l['M_dt/rgSM'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['M_dt/rgSM'][bulge_r]))
		plt.title('M_dt/rgSM frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_dt_r':	
		plt.hist(m_l['Mb+d_dt/rgSM'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['Mb+d_dt/rgSM'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Mb+d_dt/rgSM'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.2: %.2f'% np.mean(m_l['Mb+d_dt/rgSM'][disk_r]))
		plt.axvline(m_l['Mb+d_dt/rgSM'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Mb+d_dt/rgSM'][bulge_r]))
		plt.title('Mb+d_dt/rgSM frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_df_r_nomean':	
		plt.hist(m_l['Mb_dt/rgSM_b'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['Mb_dt/rgSM_b'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.title('Mb_dt/rgSM_b frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_df_r':		
		plt.hist(m_l['Mb_dt/rgSM_b'][disk_r], 100, range=(0,10), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['Mb_dt/rgSM_b'][bulge_r], 100, range=(0,10), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Mb_dt/rgSM_b'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.2: %.2f'% np.mean(m_l['Mb_dt/rgSM_b'][disk_r]))
		plt.axvline(m_l['Mb_dt/rgSM_b'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Mb_dt/rgSM_b'][bulge_r]))
		plt.title('Mb_dt/rgSM_b frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_r_nomean':	
		plt.hist(m_l['Md_dt/rgSM_d'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['Md_dt/rgSM_d'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.title('Md_dt/rgSM_d frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_r':	
		plt.hist(m_l['Md_dt/rgSM_d'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['Md_dt/rgSM_d'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Md_dt/rgSM_d'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.2: %.2f'% np.mean(m_l['Md_dt/rgSM_d'][disk_r]))
		plt.axvline(m_l['Md_dt/rgSM_d'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Md_dt/rgSM_d'][bulge_r]))
		plt.title('Md_dt/rgSM_d frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='M_dt_g':
		plt.hist(m_l['M_dt/ggSM'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['M_dt/ggSM'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.axvline(m_l['M_dt/ggSM'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['M_dt/ggSM'][disk_g]))
		plt.axvline(m_l['M_dt/ggSM'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['M_dt/ggSM'][bulge_g]))
		plt.title('M_dt/ggSM frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_dt_g':	
		plt.hist(m_l['Mb+d_dt/ggSM'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Mb+d_dt/ggSM'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.axvline(m_l['Mb+d_dt/ggSM'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['Mb+d_dt/ggSM'][disk_g]))
		plt.axvline(m_l['Mb+d_dt/ggSM'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['Mb+d_dt/ggSM'][bulge_g]))
		plt.title('Mb+d_dt/ggSM frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_dt_g_nomean':	
		plt.hist(m_l['Mb_dt/ggSM_b'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Mb_dt/ggSM_b'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.title('Mb_dt/ggSM_b frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_dt_g':		
		plt.hist(m_l['Mb_dt/ggSM_b'][disk_g], 300, range=(0,30), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Mb_dt/ggSM_b'][bulge_g], 300, range=(0,30), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.axvline(m_l['Mb_dt/ggSM_b'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['Mb_dt/ggSM_b'][disk_g]))
		plt.axvline(m_l['Mb_dt/ggSM_b'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['Mb_dt/ggSM_b'][bulge_g]))
		plt.title('Mb_dt/ggSM_b frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_g_nomean':	
		plt.hist(m_l['Md_dt/ggSM_d'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Md_dt/ggSM_d'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.title('Md_dt/ggSM_d frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_g':	
		plt.hist(m_l['Md_dt/ggSM_d'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Md_dt/ggSM_d'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.axvline(m_l['Md_dt/ggSM_d'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['Md_dt/ggSM_d'][disk_g]))
		plt.axvline(m_l['Md_dt/ggSM_d'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['Md_dt/ggSM_d'][bulge_g]))
		plt.title('Md_dt/ggSM_d frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='M_df_r':
		plt.hist(m_l['M_df/rgSM'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['M_df/rgSM'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['M_df/rgSM'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.2: %.2f'% np.mean(m_l['M_df/rgSM'][disk_r]))
		plt.axvline(m_l['M_df/rgSM'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['M_df/rgSM'][bulge_r]))
		plt.title('M_df/rgSM frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_df_r':	
		plt.hist(m_l['Mb+d_df/rgSM'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['Mb+d_df/rgSM'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Mb+d_df/rgSM'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.2: %.2f'% np.mean(m_l['Mb+d_df/rgSM'][disk_r]))
		plt.axvline(m_l['Mb+d_df/rgSM'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Mb+d_df/rgSM'][bulge_r]))
		plt.title('Mb+d_df/rgSM frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_df_r':
		plt.hist(m_l['Mb_df/rgSM_b'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['Mb_df/rgSM_b'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Mb_df/rgSM_b'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.2: %.2f'% np.mean(m_l['Mb_df/rgSM_b'][disk_r]))
		plt.axvline(m_l['Mb_df/rgSM_b'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Mb_df/rgSM_b'][bulge_r]))
		plt.title('Mb_df/rgSM_b frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_df_r':	
		plt.hist(m_l['Md_df/rgSM_d'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.2')
		plt.hist(m_l['Md_df/rgSM_d'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Md_df/rgSM_d'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.2: %.2f'% np.mean(m_l['Md_df/rgSM_d'][disk_r]))
		plt.axvline(m_l['Md_df/rgSM_d'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Md_df/rgSM_d'][bulge_r]))
		plt.title('Md_df/rgSM_d frequency with B/T < 0.2 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='M_df_g':	
		plt.hist(m_l['M_df/ggSM'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['M_df/ggSM'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.axvline(m_l['M_df/ggSM'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['M_df/ggSM'][disk_g]))
		plt.axvline(m_l['M_df/ggSM'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['M_df/ggSM'][bulge_g]))
		plt.title('M_df/ggSM frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_df_g':		
		plt.hist(m_l['Mb+d_df/ggSM'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Mb+d_df/ggSM'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.axvline(m_l['Mb+d_df/ggSM'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['Mb+d_df/ggSM'][disk_g]))
		plt.axvline(m_l['Mb+d_df/ggSM'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['Mb+d_df/ggSM'][bulge_g]))
		plt.title('Mb+d_df/ggSM frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_df_g_nomean':
		plt.hist(m_l['Mb_df/ggSM_b'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Mb_df/ggSM_b'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.title('Mb_df/ggSM_b frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_df_g':
		plt.hist(m_l['Mb_df/ggSM_b'][disk_g], 120, range=(0,12), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Mb_df/ggSM_b'][bulge_g], 120, range=(0,12), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.axvline(m_l['Mb_df/ggSM_b'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['Mb_df/ggSM_b'][disk_g]))
		plt.axvline(m_l['Mb_df/ggSM_b'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['Mb_df/ggSM_b'][bulge_g]))
		plt.title('Mb_df/ggSM_b frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_df_g':	
		plt.hist(m_l['Md_df/ggSM_d'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.1')
		plt.hist(m_l['Md_df/ggSM_d'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.7')
		plt.axvline(m_l['Md_df/ggSM_d'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['Md_df/ggSM_d'][disk_g]))
		plt.axvline(m_l['Md_df/ggSM_d'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['Md_df/ggSM_d'][bulge_g]))
		plt.title('Md_df/ggSM_d frequency with B/T < 0.1 and B/T > 0.7', y=1.05,fontsize=17, fontweight='bold')
	plt.ylabel('Frequency', fontsize=16, fontweight='bold')
	plt.xlabel('mass light ratio', fontsize=16, fontweight='bold')
	plt.legend(loc=0)
	plt.show()
		
def m_l_bd_type3(Type=None):
	if Type=='M_dt_r':
		plt.hist(m_l['M_dt/rgSM'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.12')
		plt.hist(m_l['M_dt/rgSM'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['M_dt/rgSM'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.12: %.2f'% np.mean(m_l['M_dt/ggSM'][disk_g]))
		plt.axvline(m_l['M_dt/rgSM'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['M_dt/ggSM'][bulge_g]))
		plt.title('M_dt/rgSM frequency with B/T < 0.12 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_dt_r':	
		plt.hist(m_l['Mb+d_dt/rgSM'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.12')
		plt.hist(m_l['Mb+d_dt/rgSM'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Mb+d_dt/rgSM'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.12: %.2f'% np.mean(m_l['Mb+d_dt/ggSM'][disk_g]))
		plt.axvline(m_l['Mb+d_dt/rgSM'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Mb+d_dt/ggSM'][bulge_g]))
		plt.title('Mb+d_dt/rgSM frequency with B/T < 0.12 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')		
	if Type=='B_dt_r_nomean':
		plt.hist(m_l['Mb_dt/rgSM_b'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.12')
		plt.hist(m_l['Mb_dt/rgSM_b'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.title('Mb_dt/rgSM_b frequency with B/T < 0.12 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_dt_r':
		plt.hist(m_l['Mb_dt/rgSM_b'][disk_r], 100, range=(0,10), align='right', color='r',alpha=0.5, label='B/T < 0.12')
		plt.hist(m_l['Mb_dt/rgSM_b'][bulge_r], 100, range=(0,10), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Mb_dt/rgSM_b'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.12: %.2f'% np.mean(m_l['Mb_dt/rgSM_b'][disk_r]))
		plt.axvline(m_l['Mb_dt/rgSM_b'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Mb_dt/rgSM_b'][bulge_r]))
		plt.title('Mb_dt/rgSM_b frequency with B/T < 0.12 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_r_nomean':	
		plt.hist(m_l['Md_dt/rgSM_d'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.12')
		plt.hist(m_l['Md_dt/rgSM_d'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.title('Md_dt/rgSM_d frequency with B/T < 0.12 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_r':
		plt.hist(m_l['Md_dt/rgSM_d'][disk_r], 60, range=(0,6), align='right', color='r',alpha=0.5, label='B/T < 0.12')
		plt.hist(m_l['Md_dt/rgSM_d'][bulge_r], 60, range=(0,6), align='right', color='r',alpha=0.8,label='B/T > 0.71')
		plt.axvline(m_l['Md_dt/rgSM_d'][disk_r].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.12: %.2f'% np.mean(m_l['Md_dt/rgSM_d'][disk_r]))
		plt.axvline(m_l['Md_dt/rgSM_d'][bulge_r].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.71: %.2f'% np.mean(m_l['Md_dt/rgSM_d'][bulge_r]))
		plt.title('Md_dt/rgSM_d frequency with B/T < 0.12 and B/T > 0.71', y=1.05,fontsize=17, fontweight='bold')
	if Type=='M_dt_g':
		plt.hist(m_l['M_dt/ggSM'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.05')
		plt.hist(m_l['M_dt/ggSM'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.67')
		plt.axvline(m_l['M_dt/ggSM'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.05: %.2f'% np.mean(m_l['M_dt/ggSM'][disk_g]))
		plt.axvline(m_l['M_dt/ggSM'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.67: %.2f'% np.mean(m_l['M_dt/ggSM'][bulge_g]))
		plt.title('M_dt/ggSM frequency with B/T < 0.05 and B/T > 0.67', y=1.05,fontsize=17, fontweight='bold')
	if Type=='BD_dt_g':	
		plt.hist(m_l['Mb+d_dt/ggSM'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.05')
		plt.hist(m_l['Mb+d_dt/ggSM'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.67')
		plt.axvline(m_l['Mb+d_dt/ggSM'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.05: %.2f'% np.mean(m_l['Mb+d_dt/ggSM'][disk_g]))
		plt.axvline(m_l['Mb+d_dt/ggSM'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.67: %.2f'% np.mean(m_l['Mb+d_dt/ggSM'][bulge_g]))
		plt.title('Mb+d_dt/ggSM frequency with B/T < 0.05 and B/T > 0.67', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_dt_g_nomean':
		plt.hist(m_l['Mb_dt/ggSM_b'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.05')
		plt.hist(m_l['Mb_dt/ggSM_b'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.67')
		plt.title('Mb_dt/ggSM_b frequency with B/T < 0.05 and B/T > 0.67', y=1.05,fontsize=17, fontweight='bold')
	if Type=='B_dt_g':	
		plt.hist(m_l['Mb_dt/ggSM_b'][disk_g], 500, align='right', color='g',alpha=0.5, label='B/T < 0.05')
		plt.hist(m_l['Mb_dt/ggSM_b'][bulge_g], 500, align='right', color='g',alpha=0.8,label='B/T > 0.67')
		plt.axvline(m_l['Mb_dt/ggSM_b'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.05: %.2f'% np.mean(m_l['Mb_dt/ggSM_b'][disk_g]))
		plt.axvline(m_l['Mb_dt/ggSM_b'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.67: %.2f'% np.mean(m_l['Mb_dt/ggSM_b'][bulge_g]))
		plt.title('Mb_dt/ggSM_b frequency with B/T < 0.05 and B/T > 0.67', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_g_nomean':	
		plt.hist(m_l['Md_dt/ggSM_d'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.05')
		plt.hist(m_l['Md_dt/ggSM_d'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.67')
		plt.title('Md_dt/ggSM_d frequency with B/T < 0.05 and B/T > 0.67', y=1.05,fontsize=17, fontweight='bold')
	if Type=='D_dt_g':	
		plt.hist(m_l['Md_dt/ggSM_d'][disk_g], 60, range=(0,6), align='right', color='g',alpha=0.5, label='B/T < 0.05')
		plt.hist(m_l['Md_dt/ggSM_d'][bulge_g], 60, range=(0,6), align='right', color='g',alpha=0.8,label='B/T > 0.67')
		plt.axvline(m_l['Md_dt/ggSM_d'][disk_g].mean(), color='k', linestyle='--', linewidth=2, label='mean B/T < 0.1: %.2f'% np.mean(m_l['Md_dt/ggSM_d'][disk_g]))
		plt.axvline(m_l['Md_dt/ggSM_d'][bulge_g].mean(), color='k', linestyle='-', linewidth=2, label='mean B/T > 0.7: %.2f'% np.mean(m_l['Md_dt/ggSM_d'][bulge_g]))
		plt.title('Md_dt/ggSM_d frequency with B/T < 0.05 and B/T > 0.67', y=1.05,fontsize=17, fontweight='bold')
	plt.ylabel('Frequency', fontsize=16, fontweight='bold')
	plt.xlabel('mass light ratio', fontsize=16, fontweight='bold')
	plt.legend(loc=0)
	plt.show()

'''
trendplot(cut_B_T['B/T_r'][0:len(cut_B_T['B/T_r'])],m_l['M_dt/rgSM'][0:len(m_l['M_dt/rgSM'])], color='r', marker='.',
	labelplot='$\mathtt{B/T_r}$', labeltrend='$\mathtt{TRENDLINE}$', title='M$_{dt}$/rgSM verus B/T$_r$', 
	xtitle='B/T$_r$', ytitle='M$_{dt}$/rgSM', axis=[-0.05,1.1,-0.5,10])

trendplot(cut_B_T['B/T_r'][0:len(cut_B_T['B/T_r'])],m_l['Mb+d_dt/rgSM'][0:len(m_l['Mb+d_dt/rgSM'])], color='r', marker='.',
	labelplot='$\mathtt{B/T_r}$', labeltrend='$\mathtt{TRENDLINE}$', title='Mb+d$_{dt}$/rgSM verus B/T$_r$', 
	xtitle='B/T$_r$', ytitle='Mb+d$_{dt}$/rgSM',axis=[-0.05,1.1,-0.5,12])

trendplot(cut_B_T['B/T_r'][0:len(cut_B_T['B/T_r'])],m_l['Mb_dt/rgSM_b'][0:len(m_l['Mb_dt/rgSM_b'])], color='r', marker='.',
	labelplot='$\mathtt{B/T_r}$', labeltrend='$\mathtt{TRENDLINE}$', title='Mb$_{dt}$/rgSM$_b$ verus B/T$_r$', 
	xtitle='B/T$_r$', ytitle='Mb$_{dt}$/rgSM$_b$',axis=[-0.05,1.1,-10,250])
	
trendplot(cut_B_T['B/T_r'][0:len(cut_B_T['B/T_r'])],m_l['Mb_dt/rgSM_b'][0:len(m_l['Mb_dt/rgSM_b'])], color='r', marker='.',
	labelplot='$\mathtt{B/T_r}$', labeltrend='$\mathtt{TRENDLINE}$', title='Mb$_{dt}$/rgSM$_b$ verus B/T$_r$ ZOOM IN', 
	xtitle='B/T$_r$', ytitle='Mb$_{dt}$/rgSM$_b$',axis=[-0.05,1.1,-5,50])

trendplot(cut_B_T['B/T_r'][0:len(cut_B_T['B/T_r'])],m_l['Md_dt/rgSM_d'][0:len(m_l['Md_dt/rgSM_d'])], color='r', marker='.',
	labelplot='$\mathtt{B/T_r}$', labeltrend='$\mathtt{TRENDLINE}$', title='Md$_{dt}$/rgSM$_d$ verus B/T$_r$', 
	xtitle='B/T$_r$', ytitle='Md$_{dt}$/rgSM$_d$', loc=2, axis=[-0.05,1.1,-5,140])

trendplot(cut_B_T['B/T_r'][0:len(cut_B_T['B/T_r'])],m_l['Md_dt/rgSM_d'][0:len(m_l['Md_dt/rgSM_d'])], color='r', marker='.',
	labelplot='$\mathtt{B/T_r}$', labeltrend='$\mathtt{TRENDLINE}$', title='Md$_{dt}$/rgSM$_d$ verus B/T$_r$ ZOOM IN', 
	xtitle='B/T$_r$', ytitle='Md$_{dt}$/rgSM$_d$', loc=2, axis=[-0.05,1.1,-1,12])

trendplot(cut_B_T['B/T_g'][0:len(cut_B_T['B/T_g'])],m_l['M_dt/ggSM'][0:len(m_l['M_dt/ggSM'])], color='g', marker='.',
	labelplot='$\mathtt{B/T_g}$', labeltrend='$\mathtt{TRENDLINE}$', title='M$_{dt}$/ggSM verus B/T$_g$', 
	xtitle='B/T$_g$', ytitle='M$_{dt}$/ggSM',axis=[-0.05,1.1,-1,30])

trendplot(cut_B_T['B/T_g'][0:len(cut_B_T['B/T_g'])],m_l['Mb+d_dt/ggSM'][0:len(m_l['Mb+d_dt/ggSM'])], color='g', marker='.',
	labelplot='$\mathtt{B/T_g}$', labeltrend='$\mathtt{TRENDLINE}$', title='Mb+d$_{dt}$/ggSM verus B/T$_g$', 
	xtitle='B/T$_g$', ytitle='Mb+d$_{dt}$/ggSM',axis=[-0.05,1.1,-1,35])

trendplot(cut_B_T['B/T_g'][0:len(cut_B_T['B/T_g'])],m_l['Mb_dt/ggSM_b'][0:len(m_l['Mb_dt/ggSM_b'])], color='g', marker='.',
	labelplot='$\mathtt{B/T_g}$', labeltrend='$\mathtt{TRENDLINE}$', title='Mb$_{dt}$/ggSM$_b$ verus B/T$_g$', 
	xtitle='B/T$_g$', ytitle='Mb$_{dt}$/ggSM$_b$',axis=[-0.05,1.1,-100,2500])

trendplot(cut_B_T['B/T_g'][0:len(cut_B_T['B/T_g'])],m_l['Mb_dt/ggSM_b'][0:len(m_l['Mb_dt/ggSM_b'])], color='g', marker='.',
	labelplot='$\mathtt{B/T_g}$', labeltrend='$\mathtt{TRENDLINE}$', title='Mb$_{dt}$/ggSM$_b$ verus B/T$_g$ ZOOM IN', 
	xtitle='B/T$_g$', ytitle='Mb$_{dt}$/ggSM$_b$',axis=[-0.05,1.1,-10,250])

trendplot(cut_B_T['B/T_g'][0:len(cut_B_T['B/T_g'])],m_l['Md_dt/ggSM_d'][0:len(m_l['Md_dt/ggSM_d'])], color='g', marker='.',
	labelplot='$\mathtt{B/T_g}$', labeltrend='$\mathtt{TRENDLINE}$', title='Md$_{dt}$/ggSM$_d$ verus B/T$_g$', 
	xtitle='B/T$_g$', ytitle='Md$_{dt}$/ggSM$_d$', loc=2,axis=[-0.05,1.1,-10,800])

trendplot(cut_B_T['B/T_g'][0:len(cut_B_T['B/T_g'])],m_l['Md_dt/ggSM_d'][0:len(m_l['Md_dt/ggSM_d'])], color='g', marker='.',
	labelplot='$\mathtt{B/T_g}$', labeltrend='$\mathtt{TRENDLINE}$', title='Md$_{dt}$/ggSM$_d$ verus B/T$_g$ ZOOM IN', 
	xtitle='B/T$_g$', ytitle='Md$_{dt}$/ggSM$_d$', loc=2,axis=[-0.05,1.1,-1,40])

'''
#histogram of stellar masses of disc 0-0.3 cut, edge-on(cosine<0.3), and face-on(cosine>0.91)
def mass_hist(data, edge, face, A=20, title=None):
	B=data.min()
	C=data.max()
	plt.hist(data[face],bins=A, range=(B,C),color='b',  histtype='stepfilled', alpha=0.6, label='face-on type 3 galaxies', normed=1)
	plt.hist(data[edge], bins=A, range=(B,C), color='r',  histtype='stepfilled', alpha=0.6, label='edge-on type 3 galaxies', normed=1)
	plt.hist(data, bins=A, range=(B,C), histtype='step', label='Type 3 galaxies with a\nsmall bulge and large disk', normed=1)
	plt.axvline(data[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on type 3 galaxies: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on type 3 galaxies: %.2f'% np.mean(data[edge]))
	plt.axvline(data.mean(), color='k', linestyle='-', label='mean type 3 galaxies with a\nsmall bulge and large disk: %.2f'% np.mean(data), linewidth=2)
	plt.legend(loc=0)
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel('Mass/Light Ratio',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.show()

#histogram of mass/light ratio of disc 0.01-0.3 cut, edge-on(cosine<0.3), and face-on(cosine>0.81)
def m_l_hist(data, edge, face, A=20, title=None):
	plt.hist(data[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	plt.hist(data[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	plt.hist(data, bins=A, range=(0,4), color='k', histtype='step', linestyle='dotted', linewidth=3,label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	plt.axvline(data[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.axvline(data.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data), linewidth=2)
	plt.legend(loc=0)
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel('Mass/Light Ratio (M$\odot$/L$\odot$)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.show()

#m_l_hist(m_l['M_dt/rgSM'][BT_disk_type3],edge,face, title='Normalized Frequency of M/L ratio for\nM$_{dt}$/rgSM (M$\odot$/L$\odot$) with B/T<0.3')
#m_l_hist(m_l['Mb+d_dt/rgSM'][BT_disk_type3],edge,face, title='Normalized Frequency of M/L ratio for\nMb+d$_{dt}$/rgSM (M$\odot$/L$\odot$) with B/T<0.3')
#m_l_hist(m_l['M_dt/ggSM'][BT_disk_type3],edge,face,title='Normalized Frequency of M/L ratio for\nM$_{dt}$/ggSM (M$\odot$/L$\odot$) with B/T<0.3')
#m_l_hist(m_l['Mb+d_dt/ggSM'][BT_disk_type3],edge,face, title='Normalized Frequency of M/L ratio for\nMb+d$_{dt}$/ggSM (M$\odot$/L$\odot$) with B/T<0.3')

#command to create the g-r magnitude array of B/T 0-0.3
g_r=np.subtract(cut_Mag['ggMag'],cut_Mag['rgMag'])[BT_disk_type3]

#histogram of g-r magnitudes
def mag_hist(data, edge, face, A=25, title=None):
	plt.hist(data[face], bins=A, range=(0,1), color='b', histtype='stepfilled', alpha=0.7, label='face-on',normed=1)
	plt.hist(data[edge], bins=A, range=(0,1), color='r', histtype='stepfilled', alpha=0.8, label='edge-on',normed=1)
	plt.hist(data, bins=A, range=(0,1), color='c', histtype='stepfilled', alpha=0.7, label='all type  galaxies\nwith B/T<0.3',normed=1)
	plt.axvline(data[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data[face]))
	plt.axvline(data[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data[edge]))
	plt.axvline(data.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data), linewidth=2)
	plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.xlabel('g-r (color)',fontsize=16, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.axis([-0.2,1.1,0,4.5],fontsize=30, fontweight='bold')
	plt.legend(loc=2)
	plt.show()

#mag_hist(g_r,edge,face,title='Normalized Frequency of green-band Magnitudes Minus\nred-band Magnitudes for type 3 galaxies with B/T<0.3')

# Four subplots sharing both x/y axes
def multi_hist(data1, data2, data3, data4, edge, face, A=20, title=None):
	#plt.suptitle(title, fontsize=17, fontweight='bold')
	f, ax = plt.subplots(2, 2,sharex=True,sharey=True)
	ax[0,0].hist(data1[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	ax[0,0].hist(data1[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	ax[0,0].hist(data1, bins=A, range=(0,4), histtype='step', label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	ax[0,0].axvline(data1[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data1[face]))
	ax[0,0].axvline(data1[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data1[edge]))
	ax[0,0].axvline(data1.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data1), linewidth=2)
	ax[0,0].legend(loc=0)
	
	ax[0,1].hist(data2[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	ax[0,1].hist(data2[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	ax[0,1].hist(data2, bins=A, range=(0,4), histtype='step', label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	ax[0,1].axvline(data2[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data2[face]))
	ax[0,1].axvline(data2[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data2[edge]))
	ax[0,1].axvline(data2.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data2), linewidth=2)
	ax[0,1].legend(loc=0)
	
	ax[1,0].hist(data3[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	ax[1,0].hist(data3[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	ax[1,0].hist(data3, bins=A, range=(0,4), histtype='step', label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	ax[1,0].axvline(data3[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data3[face]))
	ax[1,0].axvline(data3[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data3[edge]))
	ax[1,0].axvline(data3.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data3), linewidth=2)
	ax[1,0].legend(loc=0)
	
	ax[1,1].hist(data4[face],bins=A, range=(0,4),color='b',  histtype='stepfilled', alpha=0.7, label='face-on', normed=1)
	ax[1,1].hist(data4[edge], bins=A, range=(0,4), color='r',  histtype='stepfilled', alpha=0.8, label='edge-on', normed=1)
	ax[1,1].hist(data4, bins=A, range=(0,4), histtype='step', label='all type 3 galaxies\nwith B/T<0.3', normed=1)
	ax[1,1].axvline(data4[face].mean(), color='k', linestyle='--', linewidth=2, label='mean face-on: %.2f'% np.mean(data4[face]))
	ax[1,1].axvline(data4[edge].mean(), color='k', linestyle='-.', linewidth=3, label='mean edge-on: %.2f'% np.mean(data4[edge]))
	ax[1,1].axvline(data4.mean(), color='k', linestyle='-', label='mean type 3 galaxies\nwith B/T<0.3: %.2f'% np.mean(data4), linewidth=2)
	ax[1,1].legend(loc=0)
	
	ax[0,0].set_title('M_dt/rgSM', fontsize=14)
	ax[0,1].set_title('Mb+d_dt/rgSM', fontsize=14)
	ax[1,0].set_title('M_dt/ggSM', fontsize=14)
	ax[1,1].set_title('Mb+d_dt/ggSM', fontsize=14)
	# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
	plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
	plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
	ax[1,0].set_xlabel('Mass/Light Ratio (M$\odot$/L$\odot$)',fontsize=16, x=1)
	ax[0,0].set_ylabel('Frequency (counts)',fontsize=16)
	ax[1,0].set_ylabel('Frequency (counts)',fontsize=16)
	plt.show()

multi_hist(m_l['M_dt/rgSM'][BT_disk_type3],m_l['Mb+d_dt/rgSM'][BT_disk_type3],m_l['M_dt/ggSM'][BT_disk_type3],
	m_l['Mb+d_dt/ggSM'][BT_disk_type3],edge,face,title='Normalized Frequency of M/L ratio (M$\odot$/L$\odot$) with B/T<0.3')