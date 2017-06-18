'''Initialized 4/9/16'''

from sdss_spring2016 import *


'''cat=open_line_mass_catalog()
disk=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)'''

#Histogram of our catalog cosine inclinations of Simard, Mendel, Brinchmann
def angle_hist(data,cut):
	plt.hist(data, bins=6, range=(0,1), color='y',alpha=0.6,label='Full Catalog')
	plt.hist(data[cut], bins=6, range=(0,1), color='m',alpha=0.6,label='Disk Sample')	
	plt.axis([-0.1,1.1,0,160000],fontsize=30, fontweight='bold')
	plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	plt.title('Inclination Cosine Angles', y=1.05,fontsize=17, fontweight='bold')
	plt.legend(loc=2)
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.12, .5)
	plt.show()

#angle_hist(cat['cosi'],disk)

'''
cat=open_line_mass_catalog()
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['sfr_flag']==0)
disk=cat[cut]
edge=disk['cosi']<0.30	
face=disk['cosi']>0.87
print len(disk[edge])
print len(disk[face])
mass_half1=disk['logM']<=10.6753
mass_half2=disk['logM']>10.6753
disk_h1=disk[mass_half1]
disk_h2=disk[mass_half2]
mass_q1=disk_h1['logM']<=10.2025
mass_q2=disk_h1['logM']>10.2025
mass_q3=disk_h2['logM']<=11.0263
mass_q4=disk_h2['logM']>11.0263
disk_q1=disk_h1[mass_q1]
disk_q2=disk_h1[mass_q2]
disk_q3=disk_h2[mass_q3]
disk_q4=disk_h2[mass_q4]
mass_ei1=disk_q1['logM']<=9.76605
mass_ei2=disk_q1['logM']>9.76605
mass_ei3=disk_q2['logM']<=10.4737
mass_ei4=disk_q2['logM']>10.4737
mass_ei5=disk_q3['logM']<=10.8534
mass_ei6=disk_q3['logM']>10.8534
mass_ei7=disk_q4['logM']<=11.2354
mass_ei8=disk_q4['logM']>11.2354
disk_ei1=disk_q1[mass_ei1]
disk_ei2=disk_q1[mass_ei2]
disk_ei3=disk_q2[mass_ei3]
disk_ei4=disk_q2[mass_ei4]
disk_ei5=disk_q3[mass_ei5]
disk_ei6=disk_q3[mass_ei6]
disk_ei7=disk_q4[mass_ei7]
disk_ei8=disk_q4[mass_ei8]
edge_h1=disk_h1['cosi']<0.30	
face_h1=disk_h1['cosi']>0.83
edge_h2=disk_h2['cosi']<0.29	
face_h2=disk_h2['cosi']>0.90

edge_q1=disk_q1['cosi']<0.30	
face_q1=disk_q1['cosi']>0.81
edge_q2=disk_q2['cosi']<0.30	
face_q2=disk_q2['cosi']>0.84
edge_q3=disk_q3['cosi']<0.30	
face_q3=disk_q3['cosi']>0.87
edge_q4=disk_q4['cosi']<0.30	
face_q4=disk_q4['cosi']>0.93

edge_ei1=disk_ei1['cosi']<0.23	
face_ei1=disk_ei1['cosi']>0.80
edge_ei2=disk_ei2['cosi']<0.30	
face_ei2=disk_ei2['cosi']>0.84
edge_ei3=disk_ei3['cosi']<0.25	
face_ei3=disk_ei3['cosi']>0.87
edge_ei4=disk_ei4['cosi']<0.25	
face_ei4=disk_ei4['cosi']>0.87
edge_ei5=disk_ei5['cosi']<0.23	
face_ei5=disk_ei5['cosi']>0.90
edge_ei6=disk_ei6['cosi']<0.27	
face_ei6=disk_ei6['cosi']>0.90
edge_ei7=disk_ei7['cosi']<0.30	
face_ei7=disk_ei7['cosi']>0.90
edge_ei8=disk_ei8['cosi']<0.30	
face_ei8=disk_ei8['cosi']>0.96

print len(disk_h1[edge_h1])	
print len(disk_h1[face_h1])
print len(disk_h2[edge_h2])
print len(disk_h2[face_h2])
print ''
print len(disk_q1[edge_q1])	
print len(disk_q1[face_q1])
print len(disk_q2[edge_q2])
print len(disk_q2[face_q2])
print len(disk_q3[edge_q3])
print len(disk_q3[face_q3])
print len(disk_q4[edge_q4])
print len(disk_q4[face_q4])

print len(disk_ei1[edge_ei1])	
print len(disk_ei1[face_ei1])
print len(disk_ei2[edge_ei2])
print len(disk_ei2[face_ei2])
print len(disk_ei3[edge_ei3])
print len(disk_ei3[face_ei3])
print len(disk_ei4[edge_ei4])
print len(disk_ei4[face_ei4])
print len(disk_ei5[edge_ei5])	
print len(disk_ei5[face_ei5])
print len(disk_ei6[edge_ei6])
print len(disk_ei6[face_ei6])
print len(disk_ei7[edge_ei7])
print len(disk_ei7[face_ei7])
print len(disk_ei8[edge_ei8])
print len(disk_ei8[face_ei8])

print "half"
print "first half"
print np.amin(disk_h1['logM'])
print np.amax(disk_h1['logM'])
print ''
print "second half"
print np.amin(disk_h2['logM'])
print np.amax(disk_h2['logM'])
print ''
print "quarters"
print "first quarter"
print np.amin(disk_q1['logM'])
print np.amax(disk_q1['logM'])
print ''
print "second quarter"
print np.amin(disk_q2['logM'])
print np.amax(disk_q2['logM'])
print ''
print "third quarter"
print np.amin(disk_q3['logM'])
print np.amax(disk_q3['logM'])
print ''
print "fourth quarter"
print np.amin(disk_q4['logM'])
print np.amax(disk_q4['logM'])
print ''
print "eights"
print "first eighth"
print np.amin(disk_ei1['logM'])
print np.amax(disk_ei1['logM'])
print ''
print "second eighth"
print np.amin(disk_ei2['logM'])
print np.amax(disk_ei2['logM'])
print ''
print "third eighth"
print np.amin(disk_ei3['logM'])
print np.amax(disk_ei3['logM'])
print ''
print "fourth eighth"
print np.amin(disk_ei4['logM'])
print np.amax(disk_ei4['logM'])
print ''
print "fifth eighth"
print np.amin(disk_ei5['logM'])
print np.amax(disk_ei5['logM'])
print ''
print "sixth eighth"
print np.amin(disk_ei6['logM'])
print np.amax(disk_ei6['logM'])
print ''
print "seventh eighth"
print np.amin(disk_ei7['logM'])
print np.amax(disk_ei7['logM'])
print ''
print "final eighth"
print np.amin(disk_ei8['logM'])
print np.amax(disk_ei8['logM'])
print ''

def SFR_hist(data,edge,face, A=None, axis=None,title=None):
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

SFR_hist(disk['sfr'],edge,face,A=10,axis=[-1.7,1.4,0,6000],title='Star Formation Rates')
SFR_hist(disk_h1['sfr'],edge_h1,face_h1,A=8,axis=[-2.5,1.3,0,3500],title='Star Formation Rates lower half mass')
SFR_hist(disk_h2['sfr'],edge_h2,face_h2,A=8,axis=[-1.5,1.5,0,3000],title='Star Formation Rates upper half mass')

SFR_hist(disk_q1['sfr'],edge_q1,face_q1,A=8,axis=[-2,1,0,1600],title='Star Formation Rates first quarter mass')
SFR_hist(disk_q2['sfr'],edge_q2,face_q2,A=8,axis=[-2,1.5,0,2000],title='Star Formation Rates second quater mass')
SFR_hist(disk_q3['sfr'],edge_q3,face_q3,A=8,axis=[-1.5,1.5,0,1600],title='Star Formation Rates third quarter mass')
SFR_hist(disk_q4['sfr'],edge_q4,face_q4,A=8,axis=[-1.3,1.5,0,1000],title='Star Formation Rates fourth quarter mass')

SFR_hist(disk_ei1['sfr'],edge_ei1,face_ei1,A=6,axis=[-2.3,1,0,700],title='Star Formation Rates first eighth mass')
SFR_hist(disk_ei2['sfr'],edge_ei2,face_ei2,A=6,axis=[-1.7,1,0,1000],title='Star Formation Rates second eight mass')
SFR_hist(disk_ei3['sfr'],edge_ei3,face_ei3,A=6,axis=[-2.3,1.2,0,750],title='Star Formation Rates third eighth mass')
SFR_hist(disk_ei4['sfr'],edge_ei4,face_ei4,A=6,axis=[-2.1,1.5,0,1000],title='Star Formation Rates fourth eighth mass')
SFR_hist(disk_ei5['sfr'],edge_ei5,face_ei5,A=6,axis=[-2,1.5,0,700],title='Star Formation Rates fifth eighth mass')
SFR_hist(disk_ei6['sfr'],edge_ei6,face_ei6,A=6,axis=[-1.7,1.5,0,700],title='Star Formation Rates sixth eighth mass')
SFR_hist(disk_ei7['sfr'],edge_ei7,face_ei7,A=6,axis=[-1.6,1.5,0,900],title='Star Formation Rates seventh eighth mass')
SFR_hist(disk_ei8['sfr'],edge_ei8,face_ei8,A=6,axis=[-1.5,1.5,0,450],title='Star Formation Rates final eighth mass')


cat=open_line_mass_catalog()
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['sfr_flag']==0)
disk=cat[cut]
edge=disk['cosi']<0.30	
face=disk['cosi']>0.87
#print np.mean(disk['z'])
z_half1=disk['z']<=0.0978091
z_half2=disk['z']>0.0978091
disk_h1=disk[z_half1]
disk_h2=disk[z_half2]
#print np.mean(disk_h1['z'])
#print np.mean(disk_h2['z'])
z_q1=disk_h1['z']<=0.0618372
z_q2=disk_h1['z']>0.0618372
z_q3=disk_h2['z']<=0.141237
z_q4=disk_h2['z']>0.141237
disk_q1=disk_h1[z_q1]
disk_q2=disk_h1[z_q2]
disk_q3=disk_h2[z_q3]
disk_q4=disk_h2[z_q4]
#print np.mean(disk_q1['z'])
#print np.mean(disk_q2['z'])
#print np.mean(disk_q3['z'])
#print np.mean(disk_q4['z'])
z_ei1=disk_q1['z']<=0.0414743
z_ei2=disk_q1['z']>0.0414743
z_ei3=disk_q2['z']<=0.0790339
z_ei4=disk_q2['z']>0.0790339
z_ei5=disk_q3['z']<=0.118754
z_ei6=disk_q3['z']>0.118754
z_ei7=disk_q4['z']<=0.172588
z_ei8=disk_q4['z']>0.172588
disk_ei1=disk_q1[z_ei1]
disk_ei2=disk_q1[z_ei2]
disk_ei3=disk_q2[z_ei3]
disk_ei4=disk_q2[z_ei4]
disk_ei5=disk_q3[z_ei5]
disk_ei6=disk_q3[z_ei6]
disk_ei7=disk_q4[z_ei7]
disk_ei8=disk_q4[z_ei8]
edge_h1=disk_h1['cosi']<0.30	
face_h1=disk_h1['cosi']>0.80
edge_h2=disk_h2['cosi']<0.30	
face_h2=disk_h2['cosi']>0.93

edge_q1=disk_q1['cosi']<0.27	
face_q1=disk_q1['cosi']>0.80
edge_q2=disk_q2['cosi']<0.30	
face_q2=disk_q2['cosi']>0.83
edge_q3=disk_q3['cosi']<0.30	
face_q3=disk_q3['cosi']>0.90
edge_q4=disk_q4['cosi']<0.30	
face_q4=disk_q4['cosi']>0.97

edge_ei1=disk_ei1['cosi']<0.25	
face_ei1=disk_ei1['cosi']>0.80
edge_ei2=disk_ei2['cosi']<0.30	
face_ei2=disk_ei2['cosi']>0.80
edge_ei3=disk_ei3['cosi']<0.25	
face_ei3=disk_ei3['cosi']>0.86
edge_ei4=disk_ei4['cosi']<0.26	
face_ei4=disk_ei4['cosi']>0.87
edge_ei5=disk_ei5['cosi']<0.26	
face_ei5=disk_ei5['cosi']>0.90
edge_ei6=disk_ei6['cosi']<0.30	
face_ei6=disk_ei6['cosi']>0.92
edge_ei7=disk_ei7['cosi']<0.30	
face_ei7=disk_ei7['cosi']>0.96
edge_ei8=disk_ei8['cosi']<0.30	
face_ei8=disk_ei8['cosi']>0.99

print len(disk_h1[edge_h1])	
print len(disk_h1[face_h1])
print len(disk_h2[edge_h2])
print len(disk_h2[face_h2])
print ''
print len(disk_q1[edge_q1])	
print len(disk_q1[face_q1])
print len(disk_q2[edge_q2])
print len(disk_q2[face_q2])
print len(disk_q3[edge_q3])
print len(disk_q3[face_q3])
print len(disk_q4[edge_q4])
print len(disk_q4[face_q4])
print ''

print len(disk_ei1[edge_ei1])	
print len(disk_ei1[face_ei1])
print len(disk_ei2[edge_ei2])
print len(disk_ei2[face_ei2])
print len(disk_ei3[edge_ei3])
print len(disk_ei3[face_ei3])
print len(disk_ei4[edge_ei4])
print len(disk_ei4[face_ei4])
print len(disk_ei5[edge_ei5])	
print len(disk_ei5[face_ei5])
print len(disk_ei6[edge_ei6])
print len(disk_ei6[face_ei6])
print len(disk_ei7[edge_ei7])
print len(disk_ei7[face_ei7])
print len(disk_ei8[edge_ei8])
print len(disk_ei8[face_ei8])

print "half"
print "first half"
print np.amin(disk_h1['z'])
print np.amax(disk_h1['z'])
print ''
print "second half"
print np.amin(disk_h2['z'])
print np.amax(disk_h2['z'])
print ''

print "quarters"
print "first quarter"
print np.amin(disk_q1['z'])
print np.amax(disk_q1['z'])
print ''
print "second quarter"
print np.amin(disk_q2['z'])
print np.amax(disk_q2['z'])
print ''
print "third quarter"
print np.amin(disk_q3['z'])
print np.amax(disk_q3['z'])
print ''
print "fourth quarter"
print np.amin(disk_q4['z'])
print np.amax(disk_q4['z'])
print ''

print "eights"
print "first eighth"
print np.amin(disk_ei1['z'])
print np.amax(disk_ei1['z'])
print ''
print "second eighth"
print np.amin(disk_ei2['z'])
print np.amax(disk_ei2['z'])
print ''
print "third eighth"
print np.amin(disk_ei3['z'])
print np.amax(disk_ei3['z'])
print ''
print "fourth eighth"
print np.amin(disk_ei4['z'])
print np.amax(disk_ei4['z'])
print ''
print "fifth eighth"
print np.amin(disk_ei5['z'])
print np.amax(disk_ei5['z'])
print ''
print "sixth eighth"
print np.amin(disk_ei6['z'])
print np.amax(disk_ei6['z'])
print ''
print "seventh eighth"
print np.amin(disk_ei7['z'])
print np.amax(disk_ei7['z'])
print ''
print "final eighth"
print np.amin(disk_ei8['z'])
print np.amax(disk_ei8['z'])
print ''

def SFR_hist(data,edge,face, A=None, axis=None,title=None):
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

#SFR_hist(disk['sfr'],edge,face,A=10,axis=[-1.7,1.4,0,6000],title='Star Formation Rates')
SFR_hist(disk_h1['sfr'],edge_h1,face_h1,A=8,axis=[-2.5,1.5,0,4300],title='Star Formation Rates lower half mass')
SFR_hist(disk_h2['sfr'],edge_h2,face_h2,A=8,axis=[-1.6,1.5,0,2000],title='Star Formation Rates upper half mass')

SFR_hist(disk_q1['sfr'],edge_q1,face_q1,A=8,axis=[-2,1,0,2000],title='Star Formation Rates first quarter mass')
SFR_hist(disk_q2['sfr'],edge_q2,face_q2,A=8,axis=[-1.8,1.2,0,2200],title='Star Formation Rates second quater mass')
SFR_hist(disk_q3['sfr'],edge_q3,face_q3,A=8,axis=[-1.6,1.5,0,1500],title='Star Formation Rates third quarter mass')
SFR_hist(disk_q4['sfr'],edge_q4,face_q4,A=8,axis=[-1.5,1.6,0,400],title='Star Formation Rates fourth quarter mass')

SFR_hist(disk_ei1['sfr'],edge_ei1,face_ei1,A=6,axis=[-2.7,1.5,0,1250],title='Star Formation Rates first eighth mass')
SFR_hist(disk_ei2['sfr'],edge_ei2,face_ei2,A=6,axis=[-2,1,0,1300],title='Star Formation Rates second eight mass')
SFR_hist(disk_ei3['sfr'],edge_ei3,face_ei3,A=6,axis=[-1.8,1.2,0,1000],title='Star Formation Rates third eighth mass')
SFR_hist(disk_ei4['sfr'],edge_ei4,face_ei4,A=6,axis=[-1.7,1.5,0,1000],title='Star Formation Rates fourth eighth mass')
SFR_hist(disk_ei5['sfr'],edge_ei5,face_ei5,A=6,axis=[-2,1.5,0,1000],title='Star Formation Rates fifth eighth mass')
SFR_hist(disk_ei6['sfr'],edge_ei6,face_ei6,A=6,axis=[-1.5,1.5,0,650],title='Star Formation Rates sixth eighth mass')
SFR_hist(disk_ei7['sfr'],edge_ei7,face_ei7,A=6,axis=[-1.4,1.5,0,350],title='Star Formation Rates seventh eighth mass')
SFR_hist(disk_ei8['sfr'],edge_ei8,face_ei8,A=6,axis=[-1.1,1.6,0,150],title='Star Formation Rates final eighth mass')

plt.scatter(disk['logM'], disk['z'], c='g', marker='.')
plt.ylabel('Redshift', fontsize=16, fontweight='bold')
plt.xlabel('Mass', fontsize=16, fontweight='bold')
plt.title('Redshift versus Mass', y=1.05,fontsize=17, fontweight='bold')
plt.axis([6,12.6,0,0.35],fontsize=30, fontweight='bold')
plt.show()

plt.scatter(disk['z'],disk['logM'], c='m',marker='.')
plt.ylabel('Mass', fontsize=16, fontweight='bold')
plt.xlabel('Redshift', fontsize=16, fontweight='bold')
plt.title('Mass versus Redshift', y=1.05,fontsize=17, fontweight='bold')
plt.axis([0,0.35,6,12.6],fontsize=30, fontweight='bold')
plt.show()
'''
cat=open_line_mass_catalog()
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['B/T_g']<0.35)
cut=np.logical_and(cut,cat['sfr_flag']==0)
disk=cat[cut]
#bin disk cat by mass
cut_mass1=np.logical_and(disk['logM']>9.5,disk['logM']<=10)
cut_mass2=np.logical_and(disk['logM']>10,disk['logM']<=10.5)
cut_mass3=np.logical_and(disk['logM']>10.5,disk['logM']<=11)
cut_mass4=np.logical_and(disk['logM']>11,disk['logM']<=11.5)
disk_mass1=disk[cut_mass1]
disk_mass2=disk[cut_mass2]
disk_mass3=disk[cut_mass3]
disk_mass4=disk[cut_mass4]
#bin cat mass by z
cut_mass1_z1=np.logical_and(disk_mass1['z']>0.05,disk_mass1['z']<=0.10)
cut_mass1_z2=np.logical_and(disk_mass1['z']>0.10,disk_mass1['z']<=0.15)
cut_mass1_z3=np.logical_and(disk_mass1['z']>0.15,disk_mass1['z']<=0.20)
cut_mass2_z1=np.logical_and(disk_mass2['z']>0.05,disk_mass2['z']<=0.10)
cut_mass2_z2=np.logical_and(disk_mass2['z']>0.10,disk_mass2['z']<=0.15)
cut_mass2_z3=np.logical_and(disk_mass2['z']>0.15,disk_mass2['z']<=0.20)
cut_mass3_z1=np.logical_and(disk_mass3['z']>0.05,disk_mass3['z']<=0.10)
cut_mass3_z2=np.logical_and(disk_mass3['z']>0.10,disk_mass3['z']<=0.15)
cut_mass3_z3=np.logical_and(disk_mass3['z']>0.15,disk_mass3['z']<=0.20)
cut_mass4_z1=np.logical_and(disk_mass4['z']>0.05,disk_mass4['z']<=0.10)
cut_mass4_z2=np.logical_and(disk_mass4['z']>0.10,disk_mass4['z']<=0.15)
cut_mass4_z3=np.logical_and(disk_mass4['z']>0.15,disk_mass4['z']<=0.20)
disk_mass1_z1=disk_mass1[cut_mass1_z1]
disk_mass1_z2=disk_mass1[cut_mass1_z2]
disk_mass1_z3=disk_mass1[cut_mass1_z3]
disk_mass2_z1=disk_mass2[cut_mass1_z1]
disk_mass2_z2=disk_mass2[cut_mass2_z2]
disk_mass2_z3=disk_mass2[cut_mass2_z3]
disk_mass3_z1=disk_mass3[cut_mass2_z1]
disk_mass3_z2=disk_mass3[cut_mass3_z2]
disk_mass3_z3=disk_mass3[cut_mass3_z3]
disk_mass4_z1=disk_mass4[cut_mass4_z1]
disk_mass4_z2=disk_mass4[cut_mass4_z2]
disk_mass4_z3=disk_mass4[cut_mass4_z3]
edge_mass1_z1=disk_mass1_z1['cosi']<0.30	
face_mass1_z1=disk_mass1_z1['cosi']>0.86
edge_mass4_z1=disk_mass4_z1['cosi']<0.30	
face_mass4_z1=disk_mass4_z1['cosi']>0.86
edge_mass4_z2=disk_mass4_z2['cosi']<0.30	
face_mass4_z2=disk_mass4_z2['cosi']>0.89
edge_mass4_z3=disk_mass4_z3['cosi']<0.30	
face_mass4_z3=disk_mass4_z3['cosi']>0.96

#bin disk by z
cut_z1=np.logical_and(disk['z']>0.05,disk['z']<=0.10)
cut_z2=np.logical_and(disk['z']>0.10,disk['z']<=0.15)
cut_z3=np.logical_and(disk['z']>0.15,disk['z']<=0.20)
disk_z1=disk[cut_z1]
disk_z2=disk[cut_z2]
disk_z3=disk[cut_z3]
cut_z1_mass1=np.logical_and(disk_z1['logM']>9.5,disk_z1['logM']<=10)
cut_z1_mass2=np.logical_and(disk_z1['logM']>10,disk_z1['logM']<=10.5)
cut_z1_mass3=np.logical_and(disk_z1['logM']>10.5,disk_z1['logM']<=11)
cut_z1_mass4=np.logical_and(disk_z1['logM']>11,disk_z1['logM']<=11.5)
cut_z2_mass1=np.logical_and(disk_z2['logM']>9.5,disk_z2['logM']<=10)
cut_z2_mass2=np.logical_and(disk_z2['logM']>10,disk_z2['logM']<=10.5)
cut_z2_mass3=np.logical_and(disk_z2['logM']>10.5,disk_z2['logM']<=11)
cut_z2_mass4=np.logical_and(disk_z2['logM']>11,disk_z2['logM']<=11.5)
cut_z3_mass1=np.logical_and(disk_z3['logM']>9.5,disk_z3['logM']<=10)
cut_z3_mass2=np.logical_and(disk_z3['logM']>10,disk_z3['logM']<=10.5)
cut_z3_mass3=np.logical_and(disk_z3['logM']>10.5,disk_z3['logM']<=11)
cut_z3_mass4=np.logical_and(disk_z3['logM']>11,disk_z3['logM']<=11.5)
disk_z1_mass1=disk_z1[cut_z1_mass1]
disk_z1_mass2=disk_z1[cut_z1_mass2]
disk_z1_mass3=disk_z1[cut_z1_mass3]
disk_z1_mass4=disk_z1[cut_z1_mass4]
disk_z2_mass1=disk_z2[cut_z2_mass1]
disk_z2_mass2=disk_z2[cut_z2_mass2]
disk_z2_mass3=disk_z2[cut_z2_mass3]
disk_z2_mass4=disk_z2[cut_z2_mass4]
disk_z3_mass1=disk_z3[cut_z3_mass1]
disk_z3_mass2=disk_z3[cut_z3_mass2]
disk_z3_mass3=disk_z3[cut_z3_mass3]
disk_z3_mass4=disk_z3[cut_z3_mass4]
edge_z1_mass1=disk_z1_mass1['cosi']<0.30	
face_z1_mass1=disk_z1_mass1['cosi']>0.86
edge_z1_mass2=disk_z1_mass2['cosi']<0.30	
face_z1_mass2=disk_z1_mass2['cosi']>0.86
edge_z1_mass3=disk_z1_mass3['cosi']<0.30	
face_z1_mass3=disk_z1_mass3['cosi']>0.86
edge_z1_mass4=disk_z1_mass4['cosi']<0.30	
face_z1_mass4=disk_z1_mass4['cosi']>0.86
'''
print len(disk_z1_mass1[edge_z1_mass1])
print len(disk_z1_mass1[face_z1_mass1])
print len(disk_z1_mass2[edge_z1_mass2])
print len(disk_z1_mass2[face_z1_mass2])
print len(disk_z1_mass3[edge_z1_mass3])
print len(disk_z1_mass3[face_z1_mass3])
print len(disk_z1_mass4[edge_z1_mass4])
print len(disk_z1_mass4[face_z1_mass4])
'''
def SFR_hist(data,edge,face, A=None, axis=None,title=None):
	plt.hist(data[edge],bins=A, color='r',  histtype='stepfilled', alpha=0.7, label='edge-on galaxies')
	plt.hist(data[face], bins=A, color='b',  histtype='stepfilled', alpha=0.7, label='face-on galaxies')
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean edge-on: %.3f'% np.mean(data[edge]))
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=3, label='mean face-on: %.3f'% np.mean(data[face]))
	#plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.legend(loc=2)
	plt.xlabel('log SFR (M$\odot$/yr)',fontsize=16, fontweight='bold')
	plt.ylabel('Number',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.09, .5)
	plt.show()

#SFR_hist(disk_mass1_z1['sfr'],edge_mass1_z1,face_mass1_z1,A=10,title='Star Formation Rates of Nearby Low Mass Galaxies') #axis=[-1.7,1.4,0,6000]
#SFR_hist(disk_mass4_z1['sfr'],edge_mass4_z1,face_mass4_z1,A=10,title='Star Formation Rates of Nearby High Mass Galaxies')
SFR_hist(disk_mass4_z3['sfr'],edge_mass4_z3,face_mass4_z3,axis=[-1.5,1.7,0,200],A=10,title='Star Formation Rates of Far High Mass Galaxies')
'''
print len(disk_z1)
print len(cut_z1_mass1)
print len(cut_z1_mass2)
print len(cut_z1_mass3)
print len(cut_z1_mass4)
print ''
print len(disk_z2)
print len(cut_z2_mass1)
print len(cut_z2_mass2)
print len(cut_z2_mass3)
print len(cut_z2_mass4)
print ''
print len(disk_z3)
print len(cut_z3_mass1)
print len(cut_z3_mass2)
print len(cut_z3_mass3)
print len(cut_z3_mass4)
print ''

print len(disk)
print len(disk_mass1)+len(disk_mass2)+len(disk_mass3)+len(disk_mass4)
print len(disk_z1)+len(disk_z2)+len(disk_z3)
print len(disk_mass1_z1)+len(disk_mass1_z2)+len(disk_mass1_z3)+len(disk_mass2_z1)+len(disk_mass2_z2)+len(disk_mass2_z3)+len(disk_mass3_z1)+len(disk_mass3_z2)+len(disk_mass3_z3)+len(disk_mass4_z1)+len(disk_mass4_z2)+len(disk_mass4_z3)
print len(disk_z1_mass1)+len(disk_z1_mass2)+len(disk_z1_mass3)+len(disk_z1_mass4)+len(disk_z2_mass1)+len(disk_z2_mass2)+len(disk_z2_mass3)+len(disk_z2_mass4)+len(disk_z3_mass1)+len(disk_z3_mass2)+len(disk_z3_mass3)+len(disk_z3_mass4)

#B/T
cat=open_line_mass_catalog()
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['B/T_g']<0.35)
cut=np.logical_and(cut,cat['sfr_flag']==0)
disk=cat[cut]
edge=disk['cosi']<0.30	
face=disk['cosi']>0.87
#print len(disk[edge])
#print len(disk[face])
cut_BTg1=np.logical_and(disk['B/T_g']>0,disk['B/T_g']<=0.07)
cut_BTg2=np.logical_and(disk['B/T_g']>0.07,disk['B/T_g']<=0.14)
cut_BTg3=np.logical_and(disk['B/T_g']>0.14,disk['B/T_g']<=0.21)
cut_BTg4=np.logical_and(disk['B/T_g']>0.21,disk['B/T_g']<=0.28)
cut_BTg5=np.logical_and(disk['B/T_g']>0.28,disk['B/T_g']<=0.35)
cut_BTg6=np.logical_and(disk['B/T_g']>0.05,disk['B/T_g']<=0.15)
cut_BTg7=np.logical_and(disk['B/T_g']>0.25,disk['B/T_g']<=0.35)
disk_BTg1=disk[cut_BTg1]
disk_BTg2=disk[cut_BTg2]
disk_BTg3=disk[cut_BTg3]
disk_BTg4=disk[cut_BTg4]
disk_BTg5=disk[cut_BTg5]
disk_BTg6=disk[cut_BTg6]
disk_BTg7=disk[cut_BTg7]

edge_BTg1=disk_BTg1['cosi']<0.28	
face_BTg1=disk_BTg1['cosi']>0.86
edge_BTg2=disk_BTg2['cosi']<0.28	
face_BTg2=disk_BTg2['cosi']>0.86
edge_BTg3=disk_BTg3['cosi']<0.27	
face_BTg3=disk_BTg3['cosi']>0.86
edge_BTg4=disk_BTg4['cosi']<0.28	
face_BTg4=disk_BTg4['cosi']>0.86
edge_BTg5=disk_BTg5['cosi']<0.30	
face_BTg5=disk_BTg5['cosi']>0.90
edge_BTg6=disk_BTg6['cosi']<0.28	
face_BTg6=disk_BTg6['cosi']>0.86
edge_BTg7=disk_BTg7['cosi']<0.30	
face_BTg7=disk_BTg7['cosi']>0.89

print 'B/T g 1 edge, face'
print len(disk_BTg1[edge_BTg1])	
print len(disk_BTg1[face_BTg1])
print ''
print 'B/T g 2 edge, face'
print len(disk_BTg2[edge_BTg2])	
print len(disk_BTg2[face_BTg2])
print ''
print 'B/T g 3 edge, face'
print len(disk_BTg3[edge_BTg3])	
print len(disk_BTg3[face_BTg3])
print ''
print 'B/T g 4 edge, face'
print len(disk_BTg4[edge_BTg4])	
print len(disk_BTg4[face_BTg4])
print ''
print 'B/T g 5 edge, face'
print len(disk_BTg5[edge_BTg5])	
print len(disk_BTg5[face_BTg5])
print ''
print 'B/T g 6 edge, face'
print len(disk_BTg6[edge_BTg6])	
print len(disk_BTg6[face_BTg6])
print ''
print 'B/T g 7 edge, face'
print len(disk_BTg7[edge_BTg7])	
print len(disk_BTg7[face_BTg7])
print ''

cut_BTr1=np.logical_and(disk['B/T_r']>0,disk['B/T_r']<=0.07)
cut_BTr2=np.logical_and(disk['B/T_r']>0.07,disk['B/T_r']<=0.14)
cut_BTr3=np.logical_and(disk['B/T_r']>0.14,disk['B/T_r']<=0.21)
cut_BTr4=np.logical_and(disk['B/T_r']>0.21,disk['B/T_r']<=0.28)
cut_BTr5=np.logical_and(disk['B/T_r']>0.28,disk['B/T_r']<=0.35)
cut_BTr6=np.logical_and(disk['B/T_r']>0.05,disk['B/T_r']<=0.15)
cut_BTr7=np.logical_and(disk['B/T_r']>0.25,disk['B/T_r']<=0.35)
disk_BTr1=disk[cut_BTr1]
disk_BTr2=disk[cut_BTr2]
disk_BTr3=disk[cut_BTr3]
disk_BTr4=disk[cut_BTr4]
disk_BTr5=disk[cut_BTr5]
disk_BTr6=disk[cut_BTr6]
disk_BTr7=disk[cut_BTr7]

edge_BTr1=disk_BTr1['cosi']<0.30	
face_BTr1=disk_BTr1['cosi']>0.89
edge_BTr2=disk_BTr2['cosi']<0.28	
face_BTr2=disk_BTr2['cosi']>0.86
edge_BTr3=disk_BTr3['cosi']<0.28	
face_BTr3=disk_BTr3['cosi']>0.86
edge_BTr4=disk_BTr4['cosi']<0.28	
face_BTr4=disk_BTr4['cosi']>0.86
edge_BTr5=disk_BTr5['cosi']<0.28	
face_BTr5=disk_BTr5['cosi']>0.86
edge_BTr6=disk_BTr6['cosi']<0.28	
face_BTr6=disk_BTr6['cosi']>0.86
edge_BTr7=disk_BTr7['cosi']<0.28	
face_BTr7=disk_BTr7['cosi']>0.86

print 'B/T r 1 edge, face'
print len(disk_BTr1[edge_BTr1])	
print len(disk_BTr1[face_BTr1])
print ''
print 'B/T r 2 edge, face'
print len(disk_BTr2[edge_BTr2])	
print len(disk_BTr2[face_BTr2])
print ''
print 'B/T r 3 edge, face'
print len(disk_BTr3[edge_BTr3])	
print len(disk_BTr3[face_BTr3])
print ''
print 'B/T r 4 edge, face'
print len(disk_BTr4[edge_BTr4])	
print len(disk_BTr4[face_BTr4])
print ''
print 'B/T r 5 edge, face'
print len(disk_BTr5[edge_BTr5])	
print len(disk_BTr5[face_BTr5])
print ''
print 'B/T r 6 edge, face'
print len(disk_BTr6[edge_BTr6])	
print len(disk_BTr6[face_BTr6])
print ''
print 'B/T r 7 edge, face'
print len(disk_BTr7[edge_BTr7])	
print len(disk_BTr7[face_BTr7])
print ''

def SFR_hist(data,edge,face, A=None, axis=None,title=None):
	plt.hist(data[edge],bins=A, color='r',  histtype='stepfilled', alpha=0.7, label='edge-on galaxies')
	plt.hist(data[face], bins=A, color='b',  histtype='stepfilled', alpha=0.7, label='face-on galaxies')
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean edge-on: %.3f'% np.mean(data[edge]))
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=3, label='mean face-on: %.3f'% np.mean(data[face]))
	#plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.legend(loc=1)
	plt.xlabel('log SFR (M$\odot$/yr)',fontsize=16, fontweight='bold')
	plt.ylabel('Number',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.09, .5)
	plt.show()

SFR_hist(disk_BTg1['sfr'],edge_BTg1,face_BTg1,A=10,title='Star Formation Rates of 0 < B/T <= 0.07 G Band')
SFR_hist(disk_BTg2['sfr'],edge_BTg2,face_BTg2,A=10,title='Star Formation Rates of 0.07 < B/T <= 0.14 G Band')
SFR_hist(disk_BTg3['sfr'],edge_BTg3,face_BTg3,A=10,title='Star Formation Rates of 0.14 < B/T <= 0.21 G Band')
SFR_hist(disk_BTg4['sfr'],edge_BTg4,face_BTg4,A=10,title='Star Formation Rates of 0.21 < B/T <= 0.28 G Band')
SFR_hist(disk_BTg5['sfr'],edge_BTg5,face_BTg5,A=10,title='Star Formation Rates of 0.28 < B/T <= 0.35 G Band')

SFR_hist(disk_BTr1['sfr'],edge_BTr1,face_BTr1,A=10,title='Star Formation Rates of 0 < B/T <= 0.07 R Band')
SFR_hist(disk_BTr2['sfr'],edge_BTr2,face_BTr2,A=10,title='Star Formation Rates of 0.07 < B/T <= 0.14 R Band')
SFR_hist(disk_BTr3['sfr'],edge_BTr3,face_BTr3,A=10,title='Star Formation Rates of 0.14 < B/T <= 0.21 R Band')
SFR_hist(disk_BTr4['sfr'],edge_BTr4,face_BTr4,A=10,title='Star Formation Rates of 0.21 < B/T <= 0.28 R Band')
SFR_hist(disk_BTr5['sfr'],edge_BTr5,face_BTr5,A=10,title='Star Formation Rates of 0.28 < B/T <= 0.35 R Band')

SFR_hist(disk_BTg6['sfr'],edge_BTg6,face_BTg6,A=10,axis=[-2,3,0,1700],title='G Band Star Formation Rates of Small Bulge Disk Galaxies')
SFR_hist(disk_BTg7['sfr'],edge_BTg7,face_BTg7,A=10,axis=[-2.5,2.5,0,510],title='G Band Star Formation Rates of Large Bulge Disk Galaxies')
SFR_hist(disk_BTr6['sfr'],edge_BTr6,face_BTr6,A=10,axis=[-2,3.5,0,1100],title='R Band Star Formation Rates of Small Bulge Disk Galaxies')
SFR_hist(disk_BTr7['sfr'],edge_BTr7,face_BTr7,A=10,axis=[-2,3,0,1500],title='R Band Star Formation Rates of Small Bulge Disk Galaxies')

#Rhl
cat=open_line_mass_catalog()
print 'Rhlg_S1 min, max, mean'
print np.amin(cat['Rhlg_S1'])
print np.amax(cat['Rhlg_S1'])
print np.mean(cat['Rhlg_S1'])
print ''
print 'Rhlg_S3 min, max, mean'
print np.amin(cat['Rhlg_S3'])
print np.amax(cat['Rhlg_S3'])
print np.mean(cat['Rhlg_S3'])
print ''
print 'Rhlr min, max, mean'
print np.amin(cat['Rhlr'])
print np.amax(cat['Rhlr'])
print np.mean(cat['Rhlr'])
print ''
cut=np.logical_and(cat['B/T_r']<0.35,cat['PpS']<0.32)
cut=np.logical_and(cut,cat['B/T_g']<0.35)
cut=np.logical_and(cut,cat['sfr_flag']==0)
disk=cat[cut]

edge=disk['cosi']<0.30	
face=disk['cosi']>0.87
print len(disk[edge])
print len(disk[face])

print 'Rhlg_S1 min, max, mean'
print np.amin(disk['Rhlg_S1'])
print np.amax(disk['Rhlg_S1'])
print np.mean(disk['Rhlg_S1'])
print ''
plt.hist(disk['Rhlg_S1'], bins=4, color='y', range=(3,11))
plt.title('Rhlg_S1 range 3-11', y=1.05,fontsize=17, fontweight='bold')
plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
plt.xlabel('Rhlg_S1',fontsize=16, fontweight='bold')
plt.show()
print 'Rhlg_S3 min, max, mean'
print np.amin(disk['Rhlg_S3'])
print np.amax(disk['Rhlg_S3'])
print np.mean(disk['Rhlg_S3'])
print ''
plt.hist(disk['Rhlg_S3'], bins=4, color='g', range=(3,11))
plt.title('Rhlg_S3 range 3-11', y=1.05,fontsize=17, fontweight='bold')
plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
plt.xlabel('Rhlg_S3',fontsize=16, fontweight='bold')
plt.show()
print 'Rhlr min, max, mean'
print np.amin(disk['Rhlr'])
print np.amax(disk['Rhlr'])
print np.mean(disk['Rhlr'])
print ''
plt.hist(disk['Rhlr'], bins=4, color='m', range=(3,11))
plt.title('Rhlr range 3-11', y=1.05,fontsize=17, fontweight='bold')
plt.ylabel('Frequency (counts)',fontsize=16, fontweight='bold')
plt.xlabel('Rhlr',fontsize=16, fontweight='bold')
plt.show()


cut_Rhlg_S1_1=np.logical_and(disk['Rhlg_S1']>3,disk['Rhlg_S1']<=5)
cut_Rhlg_S1_2=np.logical_and(disk['Rhlg_S1']>5,disk['Rhlg_S1']<=7)
cut_Rhlg_S1_3=np.logical_and(disk['Rhlg_S1']>7,disk['Rhlg_S1']<=9)
cut_Rhlg_S1_4=np.logical_and(disk['Rhlg_S1']>9,disk['Rhlg_S1']<=11)
disk_Rhlg_S1_1=disk[cut_Rhlg_S1_1]
disk_Rhlg_S1_2=disk[cut_Rhlg_S1_2]
disk_Rhlg_S1_3=disk[cut_Rhlg_S1_3]
disk_Rhlg_S1_4=disk[cut_Rhlg_S1_4]

edge_Rhlg_S1_1=disk_Rhlg_S1_1['cosi']<0.30	
face_Rhlg_S1_1=disk_Rhlg_S1_1['cosi']>0.88
edge_Rhlg_S1_2=disk_Rhlg_S1_2['cosi']<0.30	
face_Rhlg_S1_2=disk_Rhlg_S1_2['cosi']>0.87
edge_Rhlg_S1_3=disk_Rhlg_S1_3['cosi']<0.28	
face_Rhlg_S1_3=disk_Rhlg_S1_3['cosi']>0.86
edge_Rhlg_S1_4=disk_Rhlg_S1_4['cosi']<0.26	
face_Rhlg_S1_4=disk_Rhlg_S1_4['cosi']>0.86

print 'Rhlg_S1 1 edge, face'
print len(disk_Rhlg_S1_1[edge_Rhlg_S1_1])	
print len(disk_Rhlg_S1_1[face_Rhlg_S1_1])
print ''
print 'Rhlg_S1 2 edge, face'
print len(disk_Rhlg_S1_2[edge_Rhlg_S1_2])	
print len(disk_Rhlg_S1_2[face_Rhlg_S1_2])
print ''
print 'Rhlg_S1 3 edge, face'
print len(disk_Rhlg_S1_3[edge_Rhlg_S1_3])	
print len(disk_Rhlg_S1_3[face_Rhlg_S1_3])
print ''
print 'Rhlg_S1 4 edge, face'
print len(disk_Rhlg_S1_4[edge_Rhlg_S1_4])	
print len(disk_Rhlg_S1_4[face_Rhlg_S1_4])
print ''

cut_Rhlg_S3_1=np.logical_and(disk['Rhlg_S3']>3,disk['Rhlg_S3']<=5)
cut_Rhlg_S3_2=np.logical_and(disk['Rhlg_S3']>5,disk['Rhlg_S3']<=7)
cut_Rhlg_S3_3=np.logical_and(disk['Rhlg_S3']>7,disk['Rhlg_S3']<=9)
cut_Rhlg_S3_4=np.logical_and(disk['Rhlg_S3']>9,disk['Rhlg_S3']<=11)
disk_Rhlg_S3_1=disk[cut_Rhlg_S3_1]
disk_Rhlg_S3_2=disk[cut_Rhlg_S3_2]
disk_Rhlg_S3_3=disk[cut_Rhlg_S3_3]
disk_Rhlg_S3_4=disk[cut_Rhlg_S3_4]

edge_Rhlg_S3_1=disk_Rhlg_S3_1['cosi']<0.30	
face_Rhlg_S3_1=disk_Rhlg_S3_1['cosi']>0.91
edge_Rhlg_S3_2=disk_Rhlg_S3_2['cosi']<0.30	
face_Rhlg_S3_2=disk_Rhlg_S3_2['cosi']>0.89
edge_Rhlg_S3_3=disk_Rhlg_S3_3['cosi']<0.30	
face_Rhlg_S3_3=disk_Rhlg_S3_3['cosi']>0.87
edge_Rhlg_S3_4=disk_Rhlg_S3_4['cosi']<0.27	
face_Rhlg_S3_4=disk_Rhlg_S3_4['cosi']>0.87

print 'Rhlg_S3 1 edge, face'
print len(disk_Rhlg_S3_1[edge_Rhlg_S3_1])	
print len(disk_Rhlg_S3_1[face_Rhlg_S3_1])
print ''
print 'Rhlg_S3 2 edge, face'
print len(disk_Rhlg_S3_2[edge_Rhlg_S3_2])	
print len(disk_Rhlg_S3_2[face_Rhlg_S3_2])
print ''
print 'Rhlg_S3 3 edge, face'
print len(disk_Rhlg_S3_3[edge_Rhlg_S3_3])	
print len(disk_Rhlg_S3_3[face_Rhlg_S3_3])
print ''
print 'Rhlg_S3 4 edge, face'
print len(disk_Rhlg_S3_4[edge_Rhlg_S3_4])	
print len(disk_Rhlg_S3_4[face_Rhlg_S3_4])
print ''

cut_Rhlr_1=np.logical_and(disk['Rhlr']>3,disk['Rhlr']<=5)
cut_Rhlr_2=np.logical_and(disk['Rhlr']>5,disk['Rhlr']<=7)
cut_Rhlr_3=np.logical_and(disk['Rhlr']>7,disk['Rhlr']<=9)
cut_Rhlr_4=np.logical_and(disk['Rhlr']>9,disk['Rhlr']<=11)
disk_Rhlr_1=disk[cut_Rhlr_1]
disk_Rhlr_2=disk[cut_Rhlr_2]
disk_Rhlr_3=disk[cut_Rhlr_3]
disk_Rhlr_4=disk[cut_Rhlr_4]

edge_Rhlr_1=disk_Rhlr_1['cosi']<0.30	
face_Rhlr_1=disk_Rhlr_1['cosi']>0.88
edge_Rhlr_2=disk_Rhlr_2['cosi']<0.30	
face_Rhlr_2=disk_Rhlr_2['cosi']>0.86
edge_Rhlr_3=disk_Rhlr_3['cosi']<0.28	
face_Rhlr_3=disk_Rhlr_3['cosi']>0.86
edge_Rhlr_4=disk_Rhlr_4['cosi']<0.26	
face_Rhlr_4=disk_Rhlr_4['cosi']>0.86

print 'Rhlr 1 edge, face'
print len(disk_Rhlr_1[edge_Rhlr_1])	
print len(disk_Rhlr_1[face_Rhlr_1])
print ''
print 'Rhlr 2 edge, face'
print len(disk_Rhlr_2[edge_Rhlr_2])	
print len(disk_Rhlr_2[face_Rhlr_2])
print ''
print 'Rhlg_S3 3 edge, face'
print len(disk_Rhlr_3[edge_Rhlr_3])	
print len(disk_Rhlr_3[face_Rhlr_3])
print ''
print 'Rhlr 4 edge, face'
print len(disk_Rhlr_4[edge_Rhlr_4])	
print len(disk_Rhlr_4[face_Rhlr_4])
print ''


def SFR_hist(data,edge,face, A=None, axis=None,title=None):
	plt.hist(data[edge],bins=A, color='r',  histtype='stepfilled', alpha=0.7, label='edge-on galaxies')
	plt.hist(data[face], bins=A, color='b',  histtype='stepfilled', alpha=0.7, label='face-on galaxies')
	plt.axvline(data[edge].mean(), color='k', linestyle='--', linewidth=2, label='mean edge-on: %.3f'% np.mean(data[edge]))
	plt.axvline(data[face].mean(), color='k', linestyle='-', linewidth=3, label='mean face-on: %.3f'% np.mean(data[face]))
	#plt.title(title, y=1.05,fontsize=17, fontweight='bold')
	plt.axis(axis,fontsize=30, fontweight='bold')
	plt.legend(loc=1)
	plt.xlabel('log SFR (M$\odot$/yr)',fontsize=16, fontweight='bold')
	plt.ylabel('Number',fontsize=16, fontweight='bold')
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.09, .5)
	plt.show()

SFR_hist(disk_Rhlg_S1_1['sfr'],edge_Rhlg_S1_1,face_Rhlg_S1_1,A=10,title='Star Formation Rates of 3 < Rhlg_S1 <= 5')
SFR_hist(disk_Rhlg_S1_2['sfr'],edge_Rhlg_S1_2,face_Rhlg_S1_2,A=10,title='Star Formation Rates of 5 < Rhlg_S1 <= 7')
SFR_hist(disk_Rhlg_S1_3['sfr'],edge_Rhlg_S1_3,face_Rhlg_S1_3,A=10,title='Star Formation Rates of 7 < Rhlg_S1 <= 9')
SFR_hist(disk_Rhlg_S1_4['sfr'],edge_Rhlg_S1_4,face_Rhlg_S1_4,A=10,title='Star Formation Rates of 9 < Rhlg_S1 <= 11')

SFR_hist(disk_Rhlg_S3_1['sfr'],edge_Rhlg_S3_1,face_Rhlg_S3_1,A=10,title='Star Formation Rates of 3 < Rhlg_S3 <= 5')
SFR_hist(disk_Rhlg_S3_2['sfr'],edge_Rhlg_S3_2,face_Rhlg_S3_2,A=10,title='Star Formation Rates of 5 < Rhlg_S3 <= 7')
SFR_hist(disk_Rhlg_S3_3['sfr'],edge_Rhlg_S3_3,face_Rhlg_S3_3,A=10,title='Star Formation Rates of 7 < Rhlg_S3 <= 9')
SFR_hist(disk_Rhlg_S3_4['sfr'],edge_Rhlg_S3_4,face_Rhlg_S3_4,A=10,title='Star Formation Rates of 9 < Rhlg_S3 <= 11')

SFR_hist(disk_Rhlr_1['sfr'],edge_Rhlr_1,face_Rhlr_1,A=10,title='Star Formation Rates of 3 < Rhlr <= 5')
SFR_hist(disk_Rhlr_2['sfr'],edge_Rhlr_2,face_Rhlr_2,A=10,title='Star Formation Rates of 5 < Rhlr <= 7')
SFR_hist(disk_Rhlr_3['sfr'],edge_Rhlr_3,face_Rhlr_3,A=10,title='Star Formation Rates of 7 < Rhlr <= 9')
SFR_hist(disk_Rhlr_4['sfr'],edge_Rhlr_4,face_Rhlr_4,A=10,title='Star Formation Rates of 9 < Rhlr <= 11')


SFR_hist(disk_Rhlg_S1_1['sfr'],edge_Rhlg_S1_1,face_Rhlg_S1_1,A=10,axis=[-2,3.5,0,750],title='G Band Table 1 Star Formation Rates of Small Disk Galaxies')
SFR_hist(disk_Rhlg_S1_4['sfr'],edge_Rhlg_S1_4,face_Rhlg_S1_4,A=10,axis=[-1.5,3.5,0,550],title='G Band Table 1 Star Formation Rates of Large Disk Galaxies')
SFR_hist(disk_Rhlg_S3_1['sfr'],edge_Rhlg_S3_1,face_Rhlg_S3_1,A=10,axis=[-2,3.5,0,550],title='G Band Table 3 Star Formation Rates of Small Disk Galaxies')
SFR_hist(disk_Rhlg_S3_4['sfr'],edge_Rhlg_S3_4,face_Rhlg_S3_4,A=10,axis=[-1.5,3,0,600],title='G Band Table 3 Star Formation Rates of Large Disk Galaxies')
SFR_hist(disk_Rhlr_1['sfr'],edge_Rhlr_1,face_Rhlr_1,A=10,axis=[-2,3.5,0,700],title='R Band Star Formation Rates of Small Disk Galaxies')
SFR_hist(disk_Rhlr_4['sfr'],edge_Rhlr_4,face_Rhlr_4,A=10,axis=[-2,3.5,0,550],title='R Band Star Formation Rates of Large Disk Galaxies')

#Histogram of our catalog cosine inclinations of old and new catalogs together 
def angle_hist(data,cut):
	plt.hist(data, bins=8, range=(0,1), color='y',alpha=0.6,label='Full Catalog')
	plt.hist(data[cut], bins=8, range=(0,1), color='m',alpha=0.6,label='Disk Sample')	
	plt.axis([-0.1,1.1,0,120000],fontsize=30, fontweight='bold')
	plt.ylabel('Number',fontsize=16, fontweight='bold')
	plt.xlabel('Inclination Cosine Angle',fontsize=16, fontweight='bold')
	#plt.title('Inclination Cosine Angles', y=1.05,fontsize=17, fontweight='bold')
	plt.legend(loc=2)
	plt.gca().xaxis.set_label_coords(.5, -.09)
	plt.gca().yaxis.set_label_coords(-.12, .5)
	plt.show()

#angle_hist(cat['cosi'],cut)
'''