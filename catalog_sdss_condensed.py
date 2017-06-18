'''Initialized 6/3/15'''

import os as os
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.constants as const
import astropy.cosmology as cosmo
import matplotlib.pyplot as plt
from decimal import *

'''#NOT USED IN CODE BELLOW: opens and reads fits file and returns a hdulist
def	read_fits_table(filename):
	hdulist=fits.open(filename)
	tbdata=hdulist[1].data
	return tbdata'''

#FUNCTIONS USED TO CREATE CATALOG
#opens and reads a fits file then returns a dictionary
def dict_fits(filename):
	dict={}
	hdulist=fits.open(filename)[1].data
	for i in range(0,len(hdulist.columns)):
		name=hdulist.columns[i].name
		dict[name]=hdulist[name]
	return dict

#returns dictionary with less or no nans. Combination of np.isfinite and np.logical_and form the True/False index 
#to decrease the array. It is best to use the arrays with greatest number of nans in the dictionary.
def cut(dict,first,second):
	fin1=np.isfinite(first)
	fin2=np.isfinite(second)
	cut=np.logical_and(fin1,fin2)
	cut_dict={}
	for key in range(0,len(dict)):
		cut_key=dict.keys()[key]
		cut_value=dict.values()[key][cut]	
		cut_dict[cut_key]=cut_value
	return cut_dict

#returns new dictionary with applied True/False index
def slice_dict(dict,index):
	sort_dict={}
	for key in range(0,len(dict)):
		sort_key=dict.keys()[key]
		sort_value=dict.values()[key][index]	
		sort_dict[sort_key]=sort_value
	return sort_dict

#turns array of floats to a list of strings with finite number of decimal places (n)
def trunc(arr,n):
	trun=[]
	for i in range(0,len(arr)):
		s='{}'.format(arr[i])
		j,l,m=s.partition('.')
		trun.append('.'.join([j,(m+'0'*n)[:n]]))
	return trun

#returns an array index of True/False based on whether the values in the list of the same index are equal
def equal_arr(list1,list2):
	equal_index=[]
	equal=len(list1)==len(list2)
	for i in range(0,len(list1)):
		if list1[i]==list2[i]:
			equal_index.append(True)
		else:
			equal_index.append(False)
	equal_arr=np.array(equal_index)
	return equal_arr


#DEBUG FUNCTIONS USED TO TEST 
#shows the min and max values in the keys of dictionaries. if nans are present it would show up as either min and or nan
def len_min_max(dict):
	for i in range(0,len(dict)):
		print dict.keys()[i]
		print 'length: '+str(len(dict.values()[i]))
		print 'minimum: '+str(np.amin(dict.values()[i]))
		print 'maximum: '+str(np.amax(dict.values()[i]))
		print ' '
		
#displays sum of nans in keys of dictionary
def nan_total(dict):
	for i in range(0,len(dict)):
		print dict.keys()[i]
		print 'total nans'
		print np.nansum(np.isnan(dict.values()[i]))
		print ' '
		
#prints in each type, length, minimum, maximum, number of nans, and sample
def the_works(dict):
	print 'item: '+str(type(dict))
	for i in range(0,len(dict)):
		print dict.keys()[i]
		print 'length: '+str(len(dict.values()[i]))
		print 'minimum: '+str(np.amin(dict.values()[i]))
		print 'maximum: '+str(np.amax(dict.values()[i]))
		print 'total nans: '+str(np.nansum(np.isnan(dict.values()[i])))
		print ' value type: '
		print type(dict.values()[i])
		print type(dict.values()[i][0])
		print 'sample:'
		print dict.values()[i][0:50]
		print ' '

#START HERE
#Binary fits files are turned to dictionaries
fits1=dict_fits('SDSS_dust_pap1.fit')
fits2=dict_fits('SDSS_dust_pap2.fit')
fits3=dict_fits('SDSS_dustfree.fit')
print "fits1['i']: " +str(len(fits1['i']))+' original 2011 catalog'
#Removed keys and corresponding values in the dictionary.
#S2g & S2g TOO MANY NANS. zmin & zmax for paper 2 dust and dustfree too different to use as reference in merging
fits1.pop('S2g')
fits1.pop('S2r')
fits2.pop('zmin')
fits2.pop('zmax')
fits3.pop('zmin')
fits3.pop('zmax')

#check if object ID elements are ascending	
#creates index such that fits2['objID'] array is ascending
#new dictionary to follow ascending objID values in paper 2 for dust and test to confirm			
test_obj1=np.all(fits1['objID'][1:]>fits1['objID'][:-1])
test_obj2=np.all(fits2['objID'][1:]>fits2['objID'][:-1])
test_obj3=np.all(fits3['objID'][1:]>fits3['objID'][:-1])

index_sort_ID2=np.argsort(fits2['objID'])	
fits2_new=slice_dict(fits2,index_sort_ID2)								#new dictionary for fits2

test_again_obj2=np.all(fits2_new['objID'][1:]>fits2_new['objID'][:-1])

#change fits1['objID'] from character array to integer array
objID_fits1=np.array(fits1['objID'],dtype=np.int)
fits1['objID']=objID_fits1

#change inclination from degrees to radians. Then use cosine to make angles equidistant ranging from 0-1
icos1=np.cos(np.radians(fits1['i']))
fits1['i']=icos1

#new dictionaries removing almost all nans, exception z in paper 1 that still has nans. 
cut1=cut(fits1,fits1['ggMag'],fits1['PpS'])
cut2=cut(fits2_new,fits2_new['logMb'],fits2_new['logMd'])
cut3=cut(fits3,fits3['logMb'],fits3['logMd'])
print "cut1['i']: "+str(len(cut1['i']))+' after nan removal'
#creates index to reduce objIDs in paper 1 so the new dictionary only has objIDs found in paper 2
index_ID_12=np.searchsorted(cut1['objID'],cut2['objID'])
cut1_new=slice_dict(cut1,index_ID_12)

#checks if arrays of objIDs are equal
test_equal_ID_12=np.array_equal(cut1_new['objID'],cut2['objID'])
test_equal_ID_23=np.array_equal(cut2['objID'],cut3['objID'])

index_ID_cut_12=np.equal(cut1_new['objID'],cut2['objID'])			#index to obtain slice where objID in paper 1 equal objID in paper 2

#new dictionaries where the objIDs for both papers are equal size and order and test
sort_ID1=slice_dict(cut1_new,index_ID_cut_12)
sort_ID2=slice_dict(cut2,index_ID_cut_12)
sort_ID3=slice_dict(cut3,index_ID_cut_12)

test_ID_12_equal=np.array_equal(sort_ID1['objID'],sort_ID2['objID'])
test_ID_23_equal=np.array_equal(sort_ID2['objID'],sort_ID3['objID'])
print "sort_ID1['i']: "+str(len(sort_ID1['i']))+' after object IDs match in paper 1 and 2'
#tests if redshifts are equal in paper 1 and 2 with dust. tests if redshifts are equal in paper 2 with dust and dustfree
test_equal_z12=np.array_equal(sort_ID1['z'],sort_ID2['z'])
test_equal_z23=np.array_equal(sort_ID2['z'],sort_ID3['z'])

#turns redshift floats to a list of strings with finite number of decimal places
#returns an array index of True and False based on whether the values in the z list of the same index are equal
list_z1=trunc(sort_ID1['z'],3)
list_z2=trunc(sort_ID2['z'],3)
index_z_cut=equal_arr(list_z1,list_z2)

#new dictionaries where the objIDs and redshift for both papers are equal size and order, as well as the test for z
zID1=slice_dict(sort_ID1,index_z_cut)
zID2=slice_dict(sort_ID2,index_z_cut)
zID3=slice_dict(sort_ID3,index_z_cut)

test_z_12_equal=set(trunc(zID1['z'],3))==set(trunc(zID2['z'],3))
test_z_23_equal=np.array_equal(zID2['z'],zID3['z'])
print "zID1['i']: "+str(len(zID1['i']))+' after redshifts are equal'
#turns RA floats to a list of strings with finite number of decimal places
#returns an array index of True and False based on whether the values in the RA list of the same index are equal
list_RA1=trunc(zID1['_RA'],2)
list_RA2=trunc(zID2['_RA'],2)
index_RA_cut=equal_arr(list_RA1,list_RA2)

#new dictionaries where the objIDs, z, RA, and DE for both papers are equal size and order, as well as the test for RA and DE
RA1=slice_dict(zID1,index_RA_cut)
RA2=slice_dict(zID2,index_RA_cut)
RA3=slice_dict(zID3,index_RA_cut)

test_RA_12_equal=set(trunc(RA1['_RA'],2))==set(trunc(RA2['_RA'],2))
test_DE_12_equal=set(trunc(RA1['_DE'],2))==set(trunc(RA2['_DE'],2))
print "RA1['i']: "+str(len(RA1['i']))+' after RA are equal'
# catalog dictionaries with attributes from all three fits files
# Excluded: Paper 1 S2g & S2g (band image smoothness parameter) TOO MANY NANS. 
#			Paper 2 zmin & zmax dust and dustfree too different. 

cat={'objID':RA1['objID'], 'z':RA1['z'], 'i':RA1['i'], 'Scale':RA1['Scale'], 'e':RA1['e'], 'nb':RA1['nb'],'Sp':RA1['Sp'],
	'Vmax':RA1['Vmax'], 'Type':RA2['Type'], 'Ft1':RA1['PpS'], 'Ft2':RA2['PpS'], 'dBD_dust':RA2['dBD'], 'dBD_df':RA3['dBD'],
	'ggMag':RA1['ggMag'], 'rgMag':RA1['rgMag'], 'e_ggMag':RA1['e_ggMag'],'e_rgMag':RA1['e_rgMag'], 'gg2d':RA1['gg2d'], 
	'rg2d':RA1['rg2d'], 'Rhlg':RA1['Rhlg'], 'Rhlr':RA1['Rhlr'],'B/T_g':RA1['__B_T_g'], 'B/T_r':RA1['__B_T_r'], 
	'RA':RA1['_RA'], 'DE':RA1['_DE'], 'Rd':RA1['Rd'], 'Re':RA1['Re'],'phib':RA1['phib'], 'phid':RA1['phid'],
	'logM_dt':RA2['logM'], 'logMb+d_dt':RA2['logMt'], 'logMb_dt':RA2['logMb'],'logMd_dt':RA2['logMd'], 
	'b_logM_dt':RA2['b_logM'], 'b_logMb+d_dt':RA2['b_logMt'], 'b_logMb_dt':RA2['b_logMb'], 'b_logMd_dt':RA2['b_logMd'],
	'B_logM_dt':RA2['B_logM1'], 'B_logMb+d_dt':RA2['B_logMt1'], 'B_logMb_dt':RA2['B_logMb1'],'B_logMd_dt':RA2['B_logMd1'],
	'logM_df':RA3['logM'], 'logMb+d_df':RA3['logMt'], 'logMb_df':RA3['logMb'],'logMd_df':RA3['logMd'],
	'b_logM_df':RA3['b_logM'], 'b_logMb+d_df':RA3['b_logMt'], 'b_logMb_df':RA3['b_logMb'],'b_logMd_df':RA3['b_logMd'],
	'B_logM_df':RA3['B_logM1'], 'B_logMb+d_df':RA3['B_logMt1'], 'B_logMb_df':RA3['B_logMb1'], 'B_logMd_df':RA3['B_logMd1']}