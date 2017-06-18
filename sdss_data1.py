'''Initialized 6/6/14'''
import os as os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.constants as const
import astropy.cosmology as cosmo
import matplotlib.pyplot as plt

def	read_fits_table(filename):
	hdulist=fits.open(filename)
	tbdata=hdulist[1].data
	return tbdata

def read_sdss(filename):
    tbdata=read_fits_table(filename)
    print len(tbdata)
    return tbdata[luc_sample(tbdata)]
    
#function for name, length, sample and type of array.
def len_sample(data,name):
	print name
	print 'length'
	print len(data)
	print 'sample'
	print data[0:10]
	print 'array type'
	print type(data)
	print ' '

#start here
data1=read_fits_table('SDSS_1D.fit')
data2=read_fits_table('SDSS_2D.fit') 	
data3=read_fits_table('SDSS_3.fit')		

#function for name, length, sample and type of array.
#debug command allows use of following function when it's greater than 10 otherwise won't show.				
debug=9
if debug>=10:
	len_sample(data2['objID'],"data2['objID']")
	len_sample(data2['z'],"data2['z']")
	len_sample(data2['i'],"data2['i']")
	len_sample(data2['Pps'],"data2['Pps'] (Ft_pap1)")		
	len_sample(data2['ggMag'],"data2['ggMag']")
	len_sample(data2['rgMag'],"data2['rgMag']")
	len_sample(data2['gg2d'],"data2['gg2d']")
	len_sample(data2['rg2d'],"data2['rg2d']")
	len_sample(data2['Rhlg'],"data2['Rhlg']")
	len_sample(data2['Rhlr'],"data2['Rhlr']")
	len_sample(data2['__B_T_g'],"data2['__B_T_g']")
	len_sample(data2['__B_T_r'],"data2['__B_T_r']")
	len_sample(data3['objID'],"data3['objID']")
	len_sample(data3['z'],"data3['z']")
	len_sample(data3['logM'],"data3['logM']")
	len_sample(data3['b_logM'], "data3['b_logM']")
	len_sample(data3['B_logM1'],"data3['B_logM1']")
	len_sample(data3['Type'],"data3['Type']")
	len_sample(data3['logMt'],"data3['logMt']")
	len_sample(data3['b_logMt'],"data3['b_logMt']")
	len_sample(data3['B_logMt1'],"data3['B_logMt1']")
	len_sample(data3['logMb'],"data3['logMb']")
	len_sample(data3['b_logMb'],"data3['b_logMb']")
	len_sample(data3['B_logMb1'],"data3['B_logMb1']")
	len_sample(data3['logMd'],"data3['logMd']")
	len_sample(data3['b_logMd'],"data3['b_logMd']")
	len_sample(data3['B_logMd1'],"data3['B_logMd1']")
	len_sample(data3['PpS'],"data3['PpS'] (Ft_pap2)")
	len_sample(data3['dBD'],"data3['dBD']")

#change data2['objID'] from character array to integer array
data2_obj_int=np.array(data2['objID'],dtype=np.int)

debug=9
if debug>=10:
	print "new data2['objID'] array type"
	print type(data2_obj_int)
	print ' '

#print data2 and data3 array lengths
print 'length of data2 arrays'
print len(data2_obj_int)
print len(data2['z'])
print len(data2['i'])
print len(data2['Pps'])		
print len(data2['ggMag'])
print len(data2['rgMag'])
print len(data2['gg2d'])
print len(data2['rg2d'])
print len(data2['Rhlg'])
print len(data2['Rhlr'])
print len(data2['__B_T_g'])
print len(data2['__B_T_r'])
print ' '
print 'length of data3 arrays'
print len(data3['objID'])
print len(data3['z'])
print len(data3['logM'])
print len(data3['b_logM'])
print len(data3['B_logM1'])
print len(data3['Type'])
print len(data3['logMt'])
print len(data3['b_logMt'])
print len(data3['B_logMt1'])
print len(data3['logMb'])
print len(data3['b_logMb'])
print len(data3['B_logMb1'])
print len(data3['logMd'])
print len(data3['b_logMd'])
print len(data3['B_logMd1'])
print len(data3['PpS'])
print len(data3['dBD'])
print ' '

#function to display number of nans in array
def nan_total(data,name):
	print name
	print 'total nans'
	nans=np.isnan(data)
	sum=np.nansum(nans)
	print sum
	print ' '

#function to display number of nans in array
#debug command allows use of following function when it's greater than 10 otherwise won't show.	
debug=9																					
if debug>=10:
	nan_total(data2['z'],"data2['z']")
	nan_total(data2['i'],"data2['i']")
	nan_total(data2['Pps'],"data2['Pps']")		
	nan_total(data2['ggMag'],"data2['ggMag']")
	nan_total(data2['rgMag'],"data2['rgMag']")
	nan_total(data2['gg2d'],"data2['gg2d']")
	nan_total(data2['rg2d'],"data2['rg2d']")
	nan_total(data2['Rhlg'],"data2['Rhlg']")
	nan_total(data2['Rhlr'],"data2['Rhlr']")
	nan_total(data2['__B_T_g'],"data2['__B_T_g']")
	nan_total(data2['__B_T_r'],"data2['__B_T_r']")
	nan_total(data3['z'],"data3['z']")
	nan_total(data3['logM'],"data3['logM']")
	nan_total(data3['b_logM'],"data3['b_logM']")
	nan_total(data3['B_logM1'],"data3['B_logM1']")
	nan_total(data3['Type'],"data3['Type']")
	nan_total(data3['logMt'],"data3['logMt']")
	nan_total(data3['b_logMt'],"data3['b_logMt']")
	nan_total(data3['B_logMt1'],"data3['B_logMt1']")
	nan_total(data3['logMb'],"data3['logMb']")
	nan_total(data3['b_logMb'],"data3['b_logMb']")
	nan_total(data3['B_logMb1'],"data3['B_logMb1']")
	nan_total(data3['logMd'],"data3['logMd']")
	nan_total(data3['b_logMd'],"data3['b_logMd']")
	nan_total(data3['B_logMd1'],"data3['B_logMd1']")
	nan_total(data3['PpS'],"data3['PpS']")
	nan_total(data3['dBD'],"data3['dBD']")


# command to exclude nans usisng the data with the largest number of nans														
finite1=np.isfinite(data2['ggMag'])
finite2=np.isfinite(data2['Pps'])
finite3=np.isfinite(data3['logMb'])
finite4=np.isfinite(data3['logMd'])
use=np.logical_and(finite1,finite2) 		
keep=np.logical_and(finite3,finite4)

debug=0																			
if debug==1:
	print "data2 length"
	print len(data2['i'])
	print "'finite' command length for data2"
	print len(use)
	print "'finite' command sample for data2"
	print use[0:5]
	print ' ' 
	print "data3 length"
	print len(data3['z'])
	print "'finite' command length for data3"
	print len(keep)
	print "'finite' command sample for data3"
	print use[0:5]
	print ' '

# exclusion of rows with nans 													
objID_use = data2_obj_int[use]
z_pap1_use = data2['z'][use]
ideg_use = data2['i'][use]
Ft_pap1_use = data2['Pps'][use]				#changed data2 Pps to Ft_pap1
ggMag_use = data2['ggMag'][use]
rgMag_use = data2['rgMag'][use]
gg2d_use = data2['gg2d'][use]
rg2d_use = data2['rg2d'][use]
Rhlg_use = data2['Rhlg'][use]
Rhlr_use = data2['Rhlr'][use]
B_T_r_use = data2['__B_T_g'][use]
B_T_g_use = data2['__B_T_r'][use]
objID_keep = data3['objID'][keep]
z_pap2_keep = data3['z'][keep]
logM_keep = data3['logM'][keep]
b_logM_keep = data3['b_logM'][keep]
B_logM1_keep = data3['B_logM1'][keep]
Type_keep = data3['Type'][keep]
logMt_keep = data3['logMt'][keep]
b_logMt_keep = data3['b_logMt'][keep]
B_logMt1_keep = data3['B_logMt1'][keep]
logMb_keep = data3['logMb'][keep]
b_logMb_keep = data3['b_logMb'][keep]
B_logMb1_keep = data3['B_logMb1'][keep]
logMd_keep = data3['logMd'][keep]
b_logMd_keep = data3['b_logMd'][keep]
B_logMd1_keep = data3['B_logMd1'][keep]
Ft_pap2_keep = data3['PpS'][keep]			#changed data3 PpS to Ft_pap2
dBD_keep = data3['dBD'][keep]

#verify lengths of arrays are the same for data2 and data3 sample				
print "length of data2 arrays after command to exclude nans"
print len(objID_use)
print len(z_pap1_use)
print len(ideg_use)
print len(Ft_pap1_use)
print len(ggMag_use)
print len(rgMag_use)
print len(gg2d_use)
print len(rg2d_use)
print len(Rhlg_use)
print len(Rhlr_use)
print len(B_T_g_use)
print len(B_T_r_use)
print ' '
print "length of data3 arrays after command to exclude nans"
print len(objID_keep)
print len(z_pap2_keep)
print len(logM_keep)
print len(b_logM_keep)
print len(B_logM1_keep)
print len(Type_keep)
print len(logMt_keep)
print len(b_logMt_keep)
print len(B_logMt1_keep)
print len(logMb_keep)
print len(b_logMb_keep)
print len(B_logMb1_keep)
print len(logMd_keep)
print len(b_logMd_keep)
print len(B_logMd1_keep)
print len(Ft_pap2_keep)
print len(dBD_keep)
print ' '


#check for nans																	
print "maximums of data2 arrays after command to exclude nans"
print np.amax(z_pap1_use)
print np.amax(ideg_use)
print np.amax(Ft_pap1_use)
print np.amax(ggMag_use)
print np.amax(rgMag_use)
print np.amax(gg2d_use)
print np.amax(rg2d_use)
print np.amax(Rhlg_use)
print np.amax(Rhlr_use)
print np.amax(B_T_g_use)
print np.amax(B_T_r_use)
print ' '
print "maximums of data3 arrays after command to exclude nans"
print np.amax(z_pap2_keep)
print np.amax(logM_keep)
print np.amax(b_logM_keep)
print np.amax(B_logM1_keep)
print np.amax(Type_keep)
print np.amax(logMt_keep)
print np.amax(b_logMt_keep)
print np.amax(B_logMt1_keep)
print np.amax(logMb_keep)
print np.amax(b_logMb_keep)
print np.amax(B_logMb1_keep)
print np.amax(logMd_keep)
print np.amax(b_logMd_keep)
print np.amax(B_logMd1_keep)
print np.amax(Ft_pap2_keep)
print np.amax(dBD_keep)
print ' '

#change from degrees to radians then use cosine to make angles equidistant ranging from 0-1
irad=np.radians(ideg_use)
icos_use=np.cos(irad)

debug=4
if debug>5:
	print 'data2 samples of inclination in degrees'
	print ideg_use[0:10]
	print ' '
	print 'data2 samples of inclination in radians'
	print irad[0:10]
	print ' '
	print 'data2 samples of inclination of cosine angle'
	print icos_use[0:10]
	print ' '

#check length of arrays and presence of nans									
def len_min_max(data,name):
	print name
	print 'length'
	print len(data)
	print 'minimum'
	print np.amin(data)
	print 'maximum'
	print np.amax(data)
	print ' '

#function checks length of arrays and presence of nans	
#debug command allows use of following function when it's equal to 1 otherwise won't show.		
debug=0
if debug==1:
	print 'lengths, minimums, and maximums of data2 arrays'
	print 'length of objID_use'
	print len(objID_use)
	print ' '
	len_min_max(z_pap1_use,'z_pap1_use')
	len_min_max(ideg_use,'ideg_use')
	len_min_max(icos_use,'icos_use')
	len_min_max(Ft_pap1_use,'Ft_pap1_use')
	len_min_max(ggMag_use,'ggMag_use')
	len_min_max(rgMag_use,'rgMag_use')
	len_min_max(gg2d_use,'gg2d_use')
	len_min_max(rg2d_use,'rg2d_use')
	len_min_max(Rhlg_use,'Rhlg_use')
	len_min_max(Rhlr_use,'Rhlr_use')
	len_min_max(B_T_g_use,'B_T_g_use')
	len_min_max(B_T_r_use,'B_T_r_use')
	print 'lengths, minimums, and maximums of data3 arrays'
	print 'length of objID_keep'
	len(objID_keep)
	print ' '
	len_min_max(z_pap2_keep,'z_pap2_keep')
	len_min_max(logM_keep,'logM_keep')
	len_min_max(b_logM_keep,'b_logM_keep')
	len_min_max(B_logM1_keep,'B_logM1_keep')
	len_min_max(Type_keep,'Type_keep')
	len_min_max(logMt_keep,'logMt_keep')
	len_min_max(b_logMt_keep,'b_logMt_keep')
	len_min_max(B_logMt1_keep,'B_logMt1_keep')
	len_min_max(logMb_keep,'logMb_keep')
	len_min_max(b_logMb_keep,'b_logMb_keep')
	len_min_max(B_logMb1_keep,'B_logMb1_keep')
	len_min_max(logMd_keep,'logMd_keep')
	len_min_max(b_logMd_keep,'b_logMd_keep')
	len_min_max(B_logMd1_keep,'B_logMd1_keep')
	len_min_max(Ft_pap2_keep,'Ft_pap2_keep')
	len_min_max(dBD_keep,'dBD_keep')
	
# m total data2 samples, n steps between samples, 1 being following value
m =len(ideg_use)
n = 1
data_objID = objID_use[0:m:n]
data_z_pap1 = z_pap1_use[0:m:n]
data_ideg = ideg_use[0:m:n]		#inclination in degrees
data_icos = icos_use[0:m:n]		#inclination after cosine
data_Ft_pap1 = Ft_pap1_use[0:m:n]
data_ggMag = ggMag_use[0:m:n]
data_rgMag = rgMag_use[0:m:n]
data_gg2d = gg2d_use[0:m:n]
data_rg2d = rg2d_use[0:m:n]
data_Rhlg = Rhlg_use[0:m:n]
data_Rhlr = Rhlr_use[0:m:n]
data_B_T_g = B_T_g_use[0:m:n]
data_B_T_r = B_T_r_use[0:m:n]

#Histogram of entire inclination data set in degrees
#degree=0: min thru max, degree=1: min  thru 84.0, degree=2: 84.0 thru max
#A is number of bins, B is min range, C is max range
def deg_hist(data_ideg, A1=84, A2=10, degrees=0, B=data_ideg.min(), C=data_ideg.max()):
	if degrees==0:
		plt.hist(data_ideg, bins=A2, range=(B, C), color='g')
		#plt.axis([-5,90,0,27000],fontsize=30, fontweight='bold')
		plt.title('Inclination Angle Frequency In Degrees 2011 Catalog', y=1.05,fontsize=17, fontweight='bold')
	if degrees==1:
		plt.hist(data_ideg, bins=A1, range=(data_ideg.min(), 84.0))
		plt.title('Inclination Angle Frequency min to 84 (degrees)')
	if degrees==2:
		plt.hist(data_ideg, bins=A2, range=(84.0, data_ideg.max()))
		plt.title('Inclination Angle Frequency 84 (degrees) to max')
	plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Angle (degrees)', fontsize=16, fontweight='bold')
	plt.show()

#deg_hist(data_ideg, degrees=0)

#Histogram of entire inclination data set after using cosine.
#A is number of bins, B is min range, C is max range
def cos_hist(data_icos,A=10, B=data_icos.min(), C=data_icos.max()):
	plt.axis([0,1.1,0,180000],fontsize=30, fontweight='bold')
	plt.hist(data_icos, bins=A, range=(B, C),color='g')
	plt.title('Inclination Cosine Angle Frequency 2011 Catalog', y=1.05,fontsize=17, fontweight='bold')
	plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle', fontsize=16, fontweight='bold')
	plt.show()

#cos_hist(data_icos)

#Histogram of entire probability data set. A is number of bins.
def Ft_hist(dataname, A=80):
	plt.hist(data_Ft_pap1, bins=A)
	plt.title('Probability Frequency')
	plt.ylabel('Frequency (counts)')
	plt.xlabel('Fit Probability (decimal)')
	plt.show()

#Ft_hist(data_Ft_pap1)

#Histogram comparing inclination in degrees versus fits probability.
#A is bin numbers, B is min inclination, C is max inclination, D is min Fits and F is max Fits.
def deg_Ft_hist(data_ideg, data_Ft_pap1,degrees=0,A=100,B=data_ideg.min(),C=data_ideg.max(),D=data_Ft_pap1.min(),F=data_Ft_pap1.max()):
	if degrees==0:
		plt.hist2d(data_ideg, data_Ft_pap1, bins=A, range=([B,C],[D,F]))
		plt.title('Fit Probability versus Inclination')
	if degrees==1:
		plt.hist2d(data_ideg, data_Ft_pap1, bins=A, range=([B,84.0],[D,F]))
		plt.title('Fit Probability versus Inclination min to 84 (degree)')
	if degrees==2:
		plt.hist2d(data_ideg, data_Ft_pap1, bins=A, range=([84.0,C],[D,F]))
		plt.title('Fit Probability versus Inclination 84 to max (degree)')
	plt.ylabel('Fit Probability (decimal)')
	plt.xlabel('Inclination (degrees)')
	plt.show() 

#deg_Ft_hist(data_ideg, data_Ft_pap1, degrees=0)

#Histogram comparing inclination of cos angle versus fits probability.
#A is bin numbers, B is min inclination, C is max inclination, D is min Fits and F is max Fits.
def cos_Ft_hist(data_icos,data_Ft_pap1,A=100,B=data_icos.min(),C=data_icos.max(),D=data_Ft_pap1.min(),F=data_Ft_pap1.max()):
	plt.hist2d(data_icos, data_Ft_pap1, bins=A, range=([B,C],[D,F]))
	plt.title('Fit Probability versus Inclination')
	plt.ylabel('Fit Probability (decimal)')
	plt.xlabel('Inclination (degrees)')
	plt.show() 

#cos_Ft_hist(data_icos, data_Ft_pap1)

#Plot comparing inclination and fit probablity
def i_Ft_plot(angle=0):
	if angle==0:
		plt.plot(data_icos, data_Ft_pap1t, 'b^')
		plt.title('Fit Probability versus Inclination Cosine Angle')
		plt.xlabel('Inclination cos angle')
		plt.axis([0, 1.1, -.1, 1.1])
	if angle==1:
		plt.plot(data_ideg, data_Ft_pap1, 'y^')
		plt.title('Fit Probability versus Inclination in Degrees')
		plt.xlabel('Inclination angle (degrees)')
		plt.axis([-5, 90, -.1, 1.1])
	plt.ylabel('Fit Probability (decimal)')
	plt.show()

#i_Ft_plot(angle=0)
#i_Ft_plot(angle=1)

#Plot comparing angle and magnitude of green and red bands
def i_Mag_plot(angle=0):
	if angle==0:
		plt.plot(data_icos, data_ggMag, 'go', label='green band')
		plt.plot(data_icos, data_rgMag, 'r^', label='red band')
		plt.title('Magnitude versus Inclination of Cosine Angle')
		plt.xlabel('inclination of cos angle')
		#plt.axis([0, 1.1, -40, 20])
	if angle==1:
		plt.plot(data_ideg, data_ggMag, 'go', label='green band')
		plt.plot(data_ideg, data_rgMag, 'r^', label='red band')
		plt.title('Magnitude versus Inclination of Angle in Degrees')
		plt.xlabel('inclination of angle (degrees)')
		#plt.axis([-5, 90, -40, 20])
	plt.legend(loc=1)
	plt.ylabel('g magnitude/r magnitude')
	plt.show()

#i_Mag_plot(angle=0)	
#i_Mag_plot(angle=1)

#command to limit data2 arrays to include population that have larger disks (0.3 or less)
use2=np.logical_and(B_T_g_use<=0.3,B_T_r_use<=0.3)

debug=0
if debug==1:
	print "length of 'use2' command and sample"
	print len(use2)
	print use2[0]
	print ' '

#limit arrays to include population that have larger disks (0.3 or less)
objID_use2=objID_use[use2]
z_pap1_use2=z_pap1_use[use2]
icos_use2=icos_use[use2]
ideg_use2=ideg_use[use2]
Ft_pap1_use2=Ft_pap1_use[use2]
ggMag_use2=ggMag_use[use2]
rgMag_use2=rgMag_use[use2]
gg2d_use2=gg2d_use[use2]
rg2d_use2=rg2d_use[use2]
Rhlg_use2=Rhlg_use[use2]
Rhlr_use2=Rhlr_use[use2]
B_T_g_use2=B_T_g_use[use2]
B_T_r_use2=B_T_r_use[use2]

#check lengths for arrays are the same
print "length of data2 arrays after command to view samples with more disk than bulge"
print len(objID_use2)
print len(z_pap1_use2)
print len(ideg_use2)
print len(icos_use2)
print len(Ft_pap1_use2)
print len(ggMag_use2)
print len(rgMag_use2)
print len(gg2d_use2)
print len(rg2d_use2)
print len(Rhlg_use2)
print len(Rhlr_use2)
print len(B_T_g_use2)
print len(B_T_r_use2)
print ' '

# d samples with disk 0.3 or less, e steps between samples, 1 being following value
d =len(ideg_use2)
e = 1
data_objID2 = objID_use2[0:d:e]
data_z_pap1_2 = z_pap1_use2[0:d:e]
data_ideg2 = ideg_use2[0:d:e]		
data_icos2 = icos_use2[0:d:e]
data_Ft_pap1_2 = Ft_pap1_use2[0:d:e]
data_ggMag2 = ggMag_use2[0:d:e]
data_rgMag2 = rgMag_use2[0:d:e]
data_gg2d2 = gg2d_use2[0:d:e]
data_rg2d2 = rg2d_use2[0:d:e]
data_Rhlg2 = Rhlg_use2[0:d:e]
data_Rhlr2 = Rhlr_use2[0:d:e]
data_B_T_g2 = B_T_g_use2[0:d:e]
data_B_T_r2 = B_T_r_use2[0:d:e]

#Histogram of inclination data set in degrees with disk 0.3 or less
#degree=0: min thru max, degree=1: min  thru 84.0, degree=2: 84.0 thru max
#A is number of bins, B is min range, C is max range
def deg_hist_d(data_ideg2, A1=84, A2=10, degrees=0, B=data_ideg2.min(), C=data_ideg2.max()):
	if degrees==0:
		plt.hist(data_ideg2, bins=A1, range=(B, C))
		plt.title('Inclination Angle Frequency with Disk Total 0-0.3')
	if degrees==1:
		plt.hist(data_ideg2, bins=A1, range=(data_ideg.min(), 84.0))
		plt.title('Inclination Angle Frequency with Disk Total 0-0.3 min to 84 (degrees)')
	if degrees==2:
		plt.hist(data_ideg2, bins=A2, range=(84.0, data_ideg.max()))
		plt.title('Inclination Angle Frequency with Disk Total 0-0.3 84 (degrees) to max')
	plt.ylabel('Frequency (counts)')
	plt.xlabel('Inclination Angle (degrees)')
	plt.show()

#deg_hist_d(data_ideg2, degrees=2)

#Histogram of inclination data set after using cosine with disk 0.3 or less.
#A is number of bins, B is min range, C is max range
def cos_hist_d(data_icos2,A=10, B=data_icos2.min(), C=data_icos2.max()):
	plt.hist(data_icos2, bins=A, range=(B, C))
	plt.title('Inclination Cosine Angle Frequency with Disk Total 0-0.3')
	plt.ylabel('Frequency (counts)')
	plt.xlabel('Inclination Cosine Angle')
	plt.show()

#cos_hist_d(data_icos2)

#Histogram of entire probability data set. A is number of bins.
def Ft_hist_d(dataname, A=80):
	plt.hist(data_Ft_pap1_2, bins=A)
	plt.title('Probability Frequency with Disk Total 0-0.3')
	plt.ylabel('Frequency (counts)')
	plt.xlabel('Fit Probability (decimal)')
	plt.show()

#Ft_hist_d(data_Ft_pap1_2)

#Histogram comparing inclination in degrees versus fits probability with Disk Total 0-0.3.
#A is bin numbers, B is min inclination, C is max inclination, D is min Fits and F is max Fits.
def deg_Ft_hist_d(data_ideg2, data_Ft_pap1_2,degrees=0,A=100,B=data_ideg2.min(),C=data_ideg2.max(),D=data_Ft_pap1_2.min(),F=data_Ft_pap1_2.max()):
	if degrees==0:
		plt.hist2d(data_ideg2, data_Ft_pap1_2, bins=A, range=([B,C],[D,F]))
		plt.title('Fit Probability versus Inclination with Disk Total 0-0.3')
	if degrees==1:
		plt.hist2d(data_ideg2, data_Ft_pap1_2, bins=A, range=([B,84.0],[D,F]))
		plt.title('Fit Probability versus Inclination min to 84 (degree) with Disk Total 0-0.3')
	if degrees==2:
		plt.hist2d(data_ideg2, data_Ft_pap1_2, bins=A, range=([84.0,C],[D,F]))
		plt.title('Fit Probability versus Inclination 84 to max (degree) with Disk Total 0-0.3')
	plt.ylabel('Fit Probability (decimal)')
	plt.xlabel('Inclination (degrees)')
	plt.show() 

#deg_Ft_hist_d(data_ideg2, data_Ft_pap1_2, degrees=0)

#Histogram comparing inclination of cos angle versus fits probability with Disk Total 0-0.3.
#A is bin numbers, B is min inclination, C is max inclination, D is min Fits and F is max Fits.
def cos_Ft_hist_d(data_icos2, data_Ft_pap1_2,A=100,B=data_icos2.min(),C=data_icos2.max(),D=data_Ft_pap1_2.min(),F=data_Ft_pap1_2.max()):
	plt.hist2d(data_icos2, data_Ft_pap1_2, bins=A, range=([B,C],[D,F]))
	plt.title('Fit Probability versus Inclination with Disk Total 0-0.3')
	plt.ylabel('Fit Probability (decimal)')
	plt.xlabel('Inclination (degrees)')
	plt.show() 

#cos_Ft_hist_d(data_icos2, data_Ft2)

#Plot comparing inclination and fit probablity with Disk Total 0-0.3
def i_Ft_plot_d(angle=0):
	if angle==0:
		plt.plot(data_icos2, data_Ft_pap1_2, 'b^')
		plt.title('Fit Probability versus Inclination Cosine Angle with Disk Total 0-0.3')
		plt.xlabel('Inclination cos angle')
		plt.axis([0, 1.1, -.1, 1.1])
	if angle==1:
		plt.plot(data_ideg2, data_Ft_pap1_2, 'y^')
		plt.title('Fit Probability versus Inclination in Degrees with Disk Total 0-0.3')
		plt.xlabel('Inclination angle (degrees)')
		plt.axis([-5, 90, -.1, 1.1])
	plt.ylabel('Fit Probability (decimal)')
	plt.show()

#i_Ft_plot_d(angle=1)

#Plot comparing angle and magnitude of green and red bands with Disk Total 0-0.3
def i_Mag_plot_d(angle=0):
	if angle==0:
		plt.plot(data_icos2, data_ggMag2, 'go', label='green band')
		plt.plot(data_icos2, data_rgMag2, 'r^', label='red band')
		plt.title('Magnitude versus Inclination of Cosine Angle with Disk Total 0-0.3')
		plt.xlabel('inclination of cos angle')
		plt.axis([0, 1.1, -40, 20])
	if angle==1:
		plt.plot(data_ideg2, data_ggMag2, 'go', label='green band')
		plt.plot(data_ideg2, data_rgMag2, 'r^', label='red band')
		plt.title('Magnitude versus Inclination of Angle in Degrees with Disk Total 0-0.3')
		plt.xlabel('inclination of angle (degrees)')
		plt.axis([-5, 90, -40, 20])
	plt.legend(loc=1)
	plt.ylabel('g magnitude/r magnitude')
	plt.show()
	
#i_Mag_plot_d(angle=0)

#check if object ID elements are ascending											
test1=np.all(objID_use[1:]>objID_use[:-1])
test2=np.all(objID_keep[1:]>objID_keep[:-1])

debug=0
if debug==1:
	print 'Check if elements in array are ascending'
	print 'objID_use'
	print test1
	print ' '
	print 'objID_keep'
	print test2
	print ' '
	
#sorting index for ascending data3[objID] array										
keep_sort=np.argsort(objID_keep)

debug=0
if debug==1:
	print 'keep_sort'
	print keep_sort[0:20]
	print len(keep_sort)
	print ' '

# use ascending sorting index and check if data3[objID] is ascending				
objID_keep2=objID_keep[keep_sort]
test3=np.all(objID_keep2[1:]>objID_keep2[:-1])

debug=0
if debug==1:
	print "after sorting data3 objID"
	print objID_keep2[0:10]
	print test3
	print ' '

#use_sort command creates an array of indices to reduce sample from first paper to only use those featured in second paper
use_sort=np.searchsorted(objID_use,objID_keep2)

#use use_sort command to decrease paper 1 arrays that correspond to paper 2 objID's						
objID_sort=objID_use[use_sort]
z_pap1_sort=z_pap1_use[use_sort]
ideg_sort=ideg_use[use_sort]	
icos_sort=icos_use[use_sort]	
Ft_pap1_sort=Ft_pap1_use[use_sort]
ggMag_sort=ggMag_use[use_sort]
rgMag_sort=rgMag_use[use_sort]
gg2d_sort=gg2d_use[use_sort]
rg2d_sort=rg2d_use[use_sort]
Rhlg_sort=Rhlg_use[use_sort]
Rhlr_sort=Rhlr_use[use_sort]
B_T_g_sort=B_T_g_use[use_sort]
B_T_r_sort=B_T_r_use[use_sort]

#use keep_sort command so that paper 2 arrays match with the new order of objID's from paper 2			
z_pap2_sort=z_pap2_keep[keep_sort]
logM_sort=logM_keep[keep_sort]
b_logM_sort=b_logM_keep[keep_sort]
B_logM1_sort=B_logM1_keep[keep_sort]
Type_sort=Type_keep[keep_sort]
logMt_sort=logMt_keep[keep_sort]
b_logMt_sort=b_logMt_keep[keep_sort]
B_logMt1_sort=B_logMt1_keep[keep_sort]
logMb_sort=logMb_keep[keep_sort]
b_logMb_sort=b_logMb_keep[keep_sort]
B_logMb1_sort=B_logMb1_keep[keep_sort]
logMd_sort=logMd_keep[keep_sort]
b_logMd_sort=b_logMd_keep[keep_sort]
B_logMd1_sort=B_logMd1_keep[keep_sort]
Ft_pap2_sort=Ft_pap2_keep[keep_sort]
dBD_sort=dBD_keep[keep_sort]

#check if red shifts match in both papers match
debug=0
if debug ==1:
	plt.plot(z_pap1_sort, z_pap2_sort, 'b^')
	plt.title('z_pap2_sort versus z_pap1_sort')
	plt.xlabel('z_pap1_sort')
	plt.ylabel('z_pap2_sort')
	plt.show()
	
	plt.plot(Ft_pap1_sort, Ft_pap2_sort, 'b^')
	plt.title('Ft_pap2_sort versus Ft_pap1_sort')
	plt.xlabel('Ft_pap1_sort')
	plt.ylabel('Ft_pap2_sort')
	plt.show()

#function that displays the array name a sample and the length to verify all sorted arrays are equal
def test_sort(data,name):
	print name
	print len(data)
	print data[0:10]
	print ' '

#check use_sort command works as expected 											
debug=0
if debug==1:
	test_sort(objID_use,'objID paper1 not sorted')
	test_sort(objID_keep,'objID paper2 not sorted')
	test_sort(use_sort,'sort index for paper 1')
	test_sort(keep_sort,'sort index for paper 1')
	test_sort(objID_keep2,'sorted paper 2 objID')
	test_sort(objID_sort,'objID_sort (sorted objID for paper 1 and 2)')
	test_sort(z_pap1_sort,'z_pap1_sort')
	test_sort(ideg_sort,'ideg_sort')
	test_sort(icos_sort,'icos_sort')
	test_sort(Ft_pap1_sort,'Ft_pap1_sort')
	test_sort(ggMag_sort,'ggMag_sort')
	test_sort(rgMag_sort,'rgMag_sort')
	test_sort(gg2d_sort,'gg2d_sort')
	test_sort(rg2d_sort,'rg2d_sort')
	test_sort(Rhlg_sort,'Rhlg_sort')
	test_sort(Rhlr_sort,'Rhlr_sort')
	test_sort(B_T_g_sort,'B_T_g_sort')
	test_sort(B_T_r_sort,'B_T_r_sort')
	test_sort(z_pap2_sort,'z_pap2_sort')
	test_sort(logM_sort,'logM_sort')
	test_sort(b_logM_sort,'b_logM_sort')
	test_sort(B_logM1_sort,'B_logM1_sort')
	test_sort(Type_sort,'Type_sort')
	test_sort(logMt_sort,'logMt_sort')
	test_sort(b_logMt_sort,'b_logMt_sort')
	test_sort(B_logMt1_sort,'B_logMt1_sort')
	test_sort(logMb_sort,'logMb_sort')
	test_sort(b_logMb_sort,'b_logMb_sort')
	test_sort(B_logMb1_sort,'B_logMb1_sort')
	test_sort(logMd_sort,'logMd_sort')
	test_sort(b_logMd_sort,'b_logMd_sort')
	test_sort(B_logMd1_sort,'B_logMd1_sort')
	test_sort(Ft_pap2_sort,'Ft_pap2_sort')
	test_sort(dBD_sort,'dBD_sort')
	

#create empty array to later put names for the catalog arrays											
def sdss_data(size):
    schema=[('objID',int),			#objID that are the same in paper 1 and 2 with no nans
    		('z_pap1',float),		#redshift from paper 1
    		('i_deg',float),		#inclination in degrees from paper 1
    		('i_cos',float),		#inclination of cosine angle from paper 1
    		('Ft_pap1',float),		#fits probability from paper 1
			('ggMag',float),
			('rgMag',float),
			('gg2d',float),
			('rg2d',float),
			('Rhlg',float),
			('Rhlr',float),
			('B_T_g',float),
			('B_T_r',float),
			('z_pap2',float),
 			('logM',float),
			('b_logM',float),
 			('B_logM1',float),
 			('Type',float),
 			('logMt',float),
 			('b_logMt',float),
 			('B_logMt1',float),
 			('logMb',float),
 			('b_logMb',float),
 			('B_logMb1',float),
 			('logMd',float),
 			('b_logMd',float),
 			('B_logMd1',float),
 			('Ft_pap2',float), 	#fits probability from paper 2
 			('dBD',float)]

    sp=np.zeros(size,dtype=schema)
    return sp

#creating the catalog by assigning names and values of the arrays							\
cat_len=len(objID_sort)
cat=sdss_data(cat_len)
cat['objID']=objID_sort
cat['z_pap1']=z_pap1_sort
cat['i_deg']=ideg_sort	
cat['i_cos']=icos_sort
cat['Ft_pap1']=Ft_pap1_sort
cat['ggMag']=ggMag_sort
cat['rgMag']=rgMag_sort
cat['gg2d']=gg2d_sort
cat['rg2d']=rg2d_sort
cat['Rhlg']=Rhlg_sort
cat['Rhlr']=Rhlr_sort
cat['B_T_g']=B_T_g_sort
cat['B_T_r']=B_T_r_sort
cat['z_pap2']=z_pap2_sort
cat['logM']=logM_sort
cat['b_logM']=b_logM_sort
cat['B_logM1']=B_logM1_sort
cat['Type']=Type_sort
cat['logMt']=logMt_sort
cat['b_logMt']=b_logMt_sort
cat['B_logMt1']=B_logMt1_sort
cat['logMb']=logMb_sort
cat['b_logMb']=b_logMb_sort
cat['B_logMb1']=B_logMb1_sort
cat['logMd']=logMd_sort
cat['b_logMd']=b_logMd_sort
cat['B_logMd1']=B_logMd1_sort
cat['Ft_pap2']=Ft_pap2_sort
cat['dBD']=dBD_sort

def check_cat(data,name):
	print name
	print 'length'
	print len(data)
	print 'minimum'
	print np.amin(data)
	print 'maximum'
	print np.amax(data)
	print 'sample'
	print data[0:10]
	print 'array type'
	print type(data)
	print ' '

#check catalog arrays for length, min and max
debug=9
if debug>=10:
	print 'length, minimum, and maximum of catalog arrays'
	print "cat['objID']"
	print 'length'
	print len(cat['objID'])
	print 'sample'
	print cat['objID'][0:10]
	print 'array type'
	print type(cat['objID'])
	print ' '
	print check_cat(cat['z_pap1'],"cat['z_pap1']")
	print check_cat(cat['i_deg'],"cat['i_deg']")
	print check_cat(cat['i_cos'],"cat['i_cos']")
	print check_cat(cat['Ft_pap1'],"cat['Ft_pap1']")
	print check_cat(cat['ggMag'],"cat['ggMag']")
	print check_cat(cat['rgMag'],"cat['rgMag']")
	print check_cat(cat['gg2d'],"cat['gg2d']")
	print check_cat(cat['rg2d'],"cat['rg2d']")
	print check_cat(cat['Rhlg'],"cat['Rhlg']")
	print check_cat(cat['Rhlr'],"cat['Rhlr']")
	print check_cat(cat['B_T_g'],"cat['B_T_g']")
	print check_cat(cat['B_T_r'],"cat['B_T_r']")
	print check_cat(cat['z_pap2'],"cat['z_pap2']")
	print check_cat(cat['logM'],"cat['logM']")
	print check_cat(cat['b_logM'],"cat['b_logM']")
	print check_cat(cat['B_logM1'],"cat['B_logM1']")
	print check_cat(cat['Type'],"cat['Type']")
	print check_cat(cat['logMt'],"cat['logMt']")
	print check_cat(cat['b_logMt'],"cat['b_logMt']")
	print check_cat(cat['B_logMt1'],"cat['B_logMt1']")
	print check_cat(cat['logMb'],"cat['logMb']")
	print check_cat(cat['b_logMb'],"cat['b_logMb']")
	print check_cat(cat['B_logMb1'],"cat['B_logMb1']")
	print check_cat(cat['logMd'],"cat['logMd']")
	print check_cat(cat['b_logMd'],"cat['b_logMd']")
	print check_cat(cat['B_logMd1'],"cat['B_logMd1']")
	print check_cat(cat['Ft_pap2'],"cat['Ft _pap2']")
	print check_cat(cat['dBD'],"cat['dBD']")

def cos_hist_old_new(icos_use,use_sort,A=10):
	B=icos_use.min()
	C=icos_use.max()
	plt.hist(icos_use, bins=A, range=(B, C),color='g', normed=1, label='2011 catalog')
	plt.hist(icos_use[use_sort],bins=A, range=(B, C),color='y', alpha=0.5, normed=1, label='our catalog')
	plt.title('Normalized Inclination Cosine Angle Frequency\nof 2011 Catalog and Our Catalog', y=1.05,fontsize=17, fontweight='bold')
	plt.ylabel('Frequency (counts)', fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle', fontsize=16, fontweight='bold')
	plt.axis([0,1.1,0,1.8],fontsize=30, fontweight='bold')
	plt.legend(loc=0)
	plt.show()

#cos_hist_old_new(icos_use,use_sort)

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

cos_hist_old_new(icos_use,use_sort)