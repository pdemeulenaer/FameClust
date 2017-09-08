import os
import linecache
import re
import math
import subprocess
import shutil
import time
import random
from random import gauss
import numpy
import numexpr as ne
import scipy
from scipy.integrate import romb
from scipy.integrate import simps
from scipy.integrate import trapz
from scipy import integrate
from scipy.integrate import quad
import scipy.optimize
from pylab import *

from sys       import stdout
from time      import sleep
from random    import randint
from threading import Thread
import time
import timeit
import threading

import Lecture_module


#Launch example for one node: nohup taskset -c 0 python GenGridModels.py 800 400 n00 &


#The script is in Dropbox, but we want to go to the right directory
#os.chdir('/home/philippe/Desktop/Discrete_models_comparaison_jtao/SC_Parameters_20/Creation_of_cluster_models_HRS_GRIDs')

start_time = time.time()

#cut_isochrone = 1  #with isochrone cut (1==rapid generation), or without (0==slow generation)

# Part ONE : INPUT GENERAL PARAMETERS
# --------------------------------------------------

age_indice  = int(sys.argv[1])
mass_indice = int(sys.argv[2])
Z_indice    = str(sys.argv[3]) 
aa,mm,zz = Lecture_module.indices_to_aa_mm_zz(age_indice, mass_indice, Z_indice)
age, mass, Z, age_indice, mass_indice, Z_indice = Lecture_module.age_mass_Z_December2011(aa,mm,zz)
#print age, mass, Z, aa, mm, zz, age_indice, mass_indice, Z_indice
#raw_input()



number_models = 1000 #int(lines[8])

if mass < 1e5:
 number_models = 1000
if mass >= 1e5 and mass < 3e5:
 number_models = 300
if mass >= 3e5 and mass < 6e5:
 number_models = 100
if mass >= 6e5 and mass < 1e6:
 number_models = 30

common_age = age #8 #float(lines[16])
common_mass = mass #10000 #float(lines[20])
#Selection of metallicity
#Z = 0.01900 #Z for all clusters
Z_indice,zz_bidon = Lecture_module.Z_to_Zindex_and_zz(Z)
#Z_indice_list = ['p02','n00','m02','m04','m06','m08','m10','m12','m14','m16','m18','m20','m22']
logAge_desired  = np.zeros(number_models)
logMass_desired = np.zeros(number_models)
Z_desired       = np.zeros(number_models)
logAge_desired[:]  = np.around(common_age, decimals=2) #We round to decimals, to ensure that.
logMass_desired[:] = np.log10(common_mass)
Z_desired[:]       = Z



# Part TWO :Sampling the stellar mass distributions using FRS or RRS schemes
# --------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#For the building of magnitudes: loading of the photometry (isochrones)
choice_metallicity = 'z{0}'.format(str(Z)[2:])     #raw_input('Which metallicity do you select? (format: Z = 0.008 --> z008) : ')
logAge_iso,M_ini,photometry = Lecture_module.Read_isochrones(choice_metallicity)
number_of_filters = len(photometry[0,:])

#Preparing the output file containing the integrated magnitudes of all clusters
data_integrated_clusters = np.zeros((number_models,2+number_of_filters+1))
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Building the stellar masses: preparation
choice_FRS_or_RRS = 1 #Full Random Sampling or Reduced Random Sampling selected ??? (1/2) 
Expansion_factor = 6.
power_index=0.2
choice_pa_or_w = 1 #Pflamm-Altenburg & Kroupa (2007) or Weidner (2013) for the Optimal Sampling relation? (1/2) '))
#May be useful
m_low = 0.01
m_up = 150.
array_CDF, array_CDF_inverse         = Lecture_module.Array_Analytical_nCDF_and_inverse(m_low,m_up)
normalization_of_number_CDF          = Lecture_module.Analytical_nCDF(m_low,m_up)
array_massCDF, array_massCDF_inverse = Lecture_module.Array_Analytical_mCDF_and_inverse(m_low,m_up)
normalization_of_mass_CDF            = Lecture_module.Analytical_mCDF(m_low,m_up)
#print array_CDF, array_CDF_inverse
#print normalization_of_number_CDF 
#print 
#print array_massCDF, array_massCDF_inverse
#print normalization_of_mass_CDF
#raw_input()

#print CDF_imf_function_number(m_low, 150.)[0] / normalization_of_number_CDF
#os.chdir('/home/philippe/Desktop/Discrete_models_comparaison_jtao/SC_Parameters_20/Creation_of_cluster_models_HRS_GRIDs/'.format(choice_metallicity))
#path = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/SC_Parameters_20/Creation_of_cluster_models_HRS_GRIDs/Clusters_built/'
path = '/home/philippe/Desktop/Clusters_built/'
mass_stoch_half = np.zeros(number_models)
mass_stoch_full = np.zeros(number_models)


#BUILDING OF OS CLUSTER
#In that case we can build only one OS cluster for all RRS clusters!
#And thus, we can define the intervals of resampling for all RRS cluster!

Mecl = 10**logMass_desired[0]

#name_output =  path+'cluster_0.txt'
#if choice_pa_or_w == 1:	#We take the theoretical Mmax-Mecl relation of Kroupa 
# os.system('./optimal {0} "0.01(pow:-0.3)0.08(pow:-1.3)0.5(pow:-2.3)150"  150 > {1}'.format(Mecl, name_output))
#if choice_pa_or_w == 2:	#We take the Observational Mmax-Mecl relation of Weidner et al. (2013) 
# m_max_weidner = -0.66 + 1.08*np.log10(Mecl) - 0.15*np.log10(Mecl)**2 + 0.0084*np.log10(Mecl)**3
# m_max_weidner = 10**m_max_weidner
# os.system('./optimal_modified {0} "0.01(pow:-0.3)0.08(pow:-1.3)0.5(pow:-2.3)150"  150  {1} > {2}'.format(Mecl,m_max_weidner,name_output))

#mass_array = np.arange(2.,5.51,0.05)
#mass_array_indice = np.linspace(200,550,71).astype(np.int)

#if choice_pa_or_w==1: path = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/SC_Parameters_20/Creation_of_cluster_models_HRS_GRIDs/Cluster_OS_grid_PK07/'
if choice_pa_or_w==1: path = '/media/philippe/a36c9bac-04fd-4046-8af2-2038962a2127/philippe/Documents/PhD/Discrete_models_comparaison_jtao/SC_Parameters_20/Creation_of_cluster_models_HRS_GRIDs/Cluster_OS_grid_PK07/'
#if choice_pa_or_w==1: path = '/opt/Cluster_OS_grid_PK07/' #NOW IN SSD, 7 APRIL 2016
if choice_pa_or_w==2: path = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/SC_Parameters_20/Creation_of_cluster_models_HRS_GRIDs/Cluster_OS_grid_W13/'

#aaaa,mmmm,zzzz=Lecture_module.indices_to_aa_mm_zz(700, 400, 'm04')
#name_output =  path+'cluster_M{0}.txt'.format(mass_array_indice[mmmm-1]) #load from the txt files
#name_output =  path+'cluster_M{0}.npy'.format(mass_array_indice[mmmm-1])  #load from the npy files
#name_output =  path+'cluster_M{0}.txt'.format(mass_indice) #load from the txt files
name_output =  path+'cluster_M{0}.npy'.format(mass_indice)  #load from the npy files
#print name_output

#print 'Flag 1'

#data_cluster = np.genfromtxt(name_output) #load from the txt files
data_cluster = np.load(name_output)        #load from the npy files
star_optimal_index = data_cluster[:,0]
star_optimal_mass = data_cluster[:,1]
star_optimal_cumulative_mass = data_cluster[:,2]
star_half_stochastic_mass = np.zeros(len(star_optimal_index))
star_full_stochastic_mass = np.zeros(len(star_optimal_index))
m_low_used_array          = np.zeros(len(star_optimal_index)) 
m_up_used_array           = np.zeros(len(star_optimal_index)) 

#Determination of how much the bins should be multiplied, in order to be referenced to M300 OS bins case
Expansion_factor_used = Expansion_factor*(Mecl*0.001)**power_index

#Building the intervals of resampling for the stochastic sampling
if choice_FRS_or_RRS==1: #if 1, then FRS
 m_low_used_array[:] = 0.010000001 #FULL RANDOM SAMPLING
 m_up_used_array[:]  = 149.9999999 #FULL RANDOM SAMPLING
if choice_FRS_or_RRS==2: #if 2, then RRS
 m_low_used_array,m_up_used_array, average_mass_bin = Lecture_module.Intervals_of_resampling(star_optimal_index,star_optimal_mass,array_massCDF, array_massCDF_inverse, Expansion_factor_used, m_low,m_up)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#print 'Flag 2'


#########################################################################################################################
#IF CUT OF ISOCHRONES TO A THRESHOLD MASS (~1/5 of the stars ALIVE built, rest integrated)
#The idea is to take the 1/5 OS stars ALIVE, derive their flux, sum them, transform into mag, and build only the more massive masses


#PART PHOTOMETRY

#Selection of good piece of isochrone (for the desired age)
#It is important to ensure that the ages fall on the isochrone ages (from 6.6 to 10.13 by steps of 0.01)
index_photo = ne.evaluate("(logAge_iso>age-0.001) & (logAge_iso<age+0.001)") #We select only the isochrones with the desired age
masses_selected = M_ini[index_photo] 
photometry_selected = np.zeros((len(masses_selected),number_of_filters))
photometry_selected = photometry[index_photo]

#Interpolation of the magnitudes with masses of stars (for all bands)
photometry_interpolated_star = np.zeros((len(star_optimal_mass),number_of_filters))
for pp in range(0,number_of_filters):
 photometry_interpolated_star[:,pp] = np.interp(star_optimal_mass,masses_selected,photometry_selected[:,pp], left=50., right=50.)
U_interpolated_star = photometry_interpolated_star[:,2]
#I_interpolated_star = photometry_interpolated_star[:,6]

#Selection of the models with mag brighter than 50 (for all bands)
index_useful_mags = ne.evaluate("U_interpolated_star<50") #TRADITIONAL 
photometry_flux = np.zeros((len(U_interpolated_star[index_useful_mags]),number_of_filters))
photometry_interpolated_star = photometry_interpolated_star[index_useful_mags,:]

#Creating the integrated magnitudes (for all bands)
photometry_flux = 10**(-0.4*photometry_interpolated_star)

most_massive_alive_star = star_optimal_mass[index_useful_mags].max()
most_massive_alive_star_index = np.where(star_optimal_mass==most_massive_alive_star)[0][0]
number_to_generate = int(0.2* (len(star_optimal_mass)-most_massive_alive_star_index) )
index_threshold = most_massive_alive_star_index+number_to_generate
#this is the range of stars alive above the threshold: [most_massive_alive_star_index:index_threshold]
#print most_massive_alive_star_index,index_threshold

star_optimal_mass_alive_LowPart_sorted = np.sort(star_optimal_mass[index_threshold:],kind='quicksort')
present_day_cluster_mass_LowPart = star_optimal_mass_alive_LowPart_sorted.sum()
print( present_day_cluster_mass_LowPart )
all_cluster_mass_LowPart = star_optimal_mass[index_threshold:].sum()

#Putting data in the output array (for all bands)
photometry_flux_LowPart = np.zeros(number_of_filters) 
for pp in range(0,number_of_filters):
 photometry_flux_LowPart[pp] = np.sum(photometry_flux[index_threshold:,pp])
#print data_integrated_clusters_LowPart
#########################################################################################################################


print( star_optimal_mass[:index_threshold].sum(), star_optimal_mass[index_threshold:].sum() )
#print 'Flag 3'
input()



#Loop on all the clusters
for jj in range(0,number_models): # range(number_models):
 age = logAge_desired[jj]

 #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 #APPLYING THE STOCHASTIC SAMPLING!!!
 if choice_FRS_or_RRS==1: #if 2, then RRS
  m_low_used_array[:] = m_low_used_array[:]*0.  + 0.010000001 #FULL RANDOM SAMPLING
  m_up_used_array[:]  = m_up_used_array[:]*0.   + 149.9999999 #FULL RANDOM SAMPLING
 star_half_stochastic_mass = Lecture_module.Star_Cluster_stellar_masses_generator_N_star(m_low,m_up,m_low_used_array,m_up_used_array,array_CDF, array_CDF_inverse)
 #mass_stoch_half[jj]=star_half_stochastic_mass.sum()

 ##star_half_stochastic_mass_HighPart = np.sort(star_half_stochastic_mass,kind='mergesort')[::-1]
 ##star_half_stochastic_mass_HighPart = np.sort(star_half_stochastic_mass,kind='quicksort')[::-1]
 star_half_stochastic_mass_HighPart = star_half_stochastic_mass[::-1]  #7April2016
 star_half_stochastic_mass_HighPart = star_half_stochastic_mass_HighPart[:index_threshold]    #In absolute, mistake here!
 #Should be cut not at index_threshold, but at the mass where is the cut, because now new stars generated! But little influence? 
 mass_stoch_half[jj]=star_half_stochastic_mass_HighPart.sum() + all_cluster_mass_LowPart
 #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 '''
 #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 #PART PHOTOMETRY

 #Interpolation of the magnitudes with masses of stars (for all bands)
 photometry_interpolated_star = np.zeros((len(star_half_stochastic_mass_HighPart),number_of_filters))
 for pp in range(0,number_of_filters):
  photometry_interpolated_star[:,pp] = np.interp(star_half_stochastic_mass_HighPart,masses_selected,photometry_selected[:,pp], left=50., right=50.)
 U_interpolated_star = photometry_interpolated_star[:,2]
 #I_interpolated_star = photometry_interpolated_star[:,6]

 #Selection of the models with mag brighter than 50 (for all bands)
 index_useful_mags = ne.evaluate("U_interpolated_star<50") #TRADITIONAL 
 #index_useful_mags = ne.evaluate("(I_interpolated_star<50) & (I_interpolated_star>-3)") #CUTTING MASSIVE STARS WITH I BAND BRIGHTER THAN -3
 photometry_flux = np.zeros((len(U_interpolated_star[index_useful_mags]),number_of_filters))
 photometry_interpolated_star = photometry_interpolated_star[index_useful_mags,:]

 #Creating the integrated magnitudes (for all bands)
 photometry_flux = 10**(-0.4*photometry_interpolated_star)

 most_massive_alive_star = star_half_stochastic_mass_HighPart[index_useful_mags].max()
 select_stellar_mass_alive = ne.evaluate('star_half_stochastic_mass_HighPart<=most_massive_alive_star')
 star_half_stochastic_mass_alive = star_half_stochastic_mass_HighPart[select_stellar_mass_alive][::-1]
 present_day_cluster_mass = star_half_stochastic_mass_alive.sum() + present_day_cluster_mass_LowPart

 #Putting data in the output array (for all bands)
 data_integrated_clusters[jj,1] = np.log10(mass_stoch_half[jj] ) #the real (sampled) initial mass of the cluster, not just the desired one
 data_integrated_clusters[jj,-1]= np.log10(present_day_cluster_mass)        #the real (sampled) PRESENT-DAY mass of the cluster
 for pp in range(0,number_of_filters):
  data_integrated_clusters[jj,2+pp] = -2.5*np.log10( np.sum(photometry_flux[:,pp]) + photometry_flux_LowPart[pp] )

 #print jj, star_half_stochastic_mass.sum()#, age, data_integrated_clusters[jj,2],most_massive_alive_star, present_day_cluster_mass
 print jj, mass_stoch_half[jj]#, age, data_integrated_clusters[jj,2],most_massive_alive_star, present_day_cluster_mass
 #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 '''
 print( jj )
print( mass_stoch_half.mean() )

input()







# Part FOUR : Printing in final file
# ----------------------------------
data_integrated_clusters[:,0] = logAge_desired
data_integrated_clusters_Notextincted = np.zeros((number_models,5+number_of_filters+1))
for ii in range(0,number_models):
 data_integrated_clusters_Notextincted[ii,0] = ii+1        

for pp in range(0,number_of_filters):
 data_integrated_clusters_Notextincted[:,5+pp] = data_integrated_clusters[:,2+pp] 

data_integrated_clusters_Notextincted[:,0] = data_integrated_clusters_Notextincted[:,0] #ID 
data_integrated_clusters_Notextincted[:,1] = data_integrated_clusters[:,0]              #Age
data_integrated_clusters_Notextincted[:,2] = data_integrated_clusters[:,1]              #INITIAL Mass (real, sampled, not desired)
data_integrated_clusters_Notextincted[:,3] = 0.                                         #Extinction
data_integrated_clusters_Notextincted[:,4] = Z                                          #Metallicity
data_integrated_clusters_Notextincted[:,-1]= data_integrated_clusters[:,-1]             #PRESENT-DAY Mass

path_grid_out = '/home/philippe/Desktop/'
#path_grid_out = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grid_FRS_allZ_Kroupa_generated/GRID_NPZ/{0}/'.format(Z_indice)
#path_grid_out = '/mnt/storage/philippe/Grid_FRS_allZ_Kroupa_generated/GRID_NPZ/{0}/'.format(Z_indice)
outfile_NoEbv = open(path_grid_out+'Clusters_t{0}_M{1}_Z{2}'.format(age_indice,mass_indice,Z_indice),'w')  #If we need all the filters
#print >> outfile_NoEbv, '# ID        Age       Mass      Ext       Z         FUV       NUV       U         B         V         R         I         J         H         K         F275W     F336W     F475W     F814W     F110W     F160W     u_CFHT    g_CFHT    r_CFHT    i_CFHT    z_CFHT    u_SDSS    g_SDSS    r_SDSS    i_SDSS    z_SDSS    J_2MASS   H_2MASS   Ks_2MASS  a_BATC    b_BATC    c_BATC    d_BATC    e_BATC    f_BATC    g_BATC    h_BATC    i_BATC    j_BATC    k_BATC    m_BATC    n_BATC    o_BATC    p_BATC    t_BATC    IRAC36    IRAC45    IRAC58    IRAC80    mip24     mip70     mip160    Mass_present'
sss = '# ID        Age       Mass      Ext       Z         FUV       NUV       U         B         V         R         I         J         H         K         F275W     F336W     F475W     F814W     F110W     F160W     u_CFHT    g_CFHT    r_CFHT    i_CFHT    z_CFHT    u_SDSS    g_SDSS    r_SDSS    i_SDSS    z_SDSS    J_2MASS   H_2MASS   Ks_2MASS  a_BATC    b_BATC    c_BATC    d_BATC    e_BATC    f_BATC    g_BATC    h_BATC    i_BATC    j_BATC    k_BATC    m_BATC    n_BATC    o_BATC    p_BATC    t_BATC    IRAC36    IRAC45    IRAC58    IRAC80    mip24     mip70     mip160    Mass_present'
print(sss, end="", file=outfile_NoEbv)
np.savetxt(outfile_NoEbv,data_integrated_clusters_Notextincted, fmt='%9.5f') #
outfile_NoEbv.close()


print("--- %s seconds ---" % (time.time() - start_time))



















