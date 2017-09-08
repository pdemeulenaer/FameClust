# module de lecture des fichiers .sum de SimClust pour python
# Marie je te confie tout ca...

import os
import linecache
import re
import math
import subprocess
import shutil
import time
import random
from random import gauss
import math

import numpy as np
import numexpr as ne
import scipy
from scipy import integrate
from scipy.integrate import quad
from sklearn import mixture

#import numpy as np        #bizarre, pas de numpy ni scipy sur jt serveur!!


def Lecture_sumfile(file_name):

        infile1 = open(file_name,'r')

        UUU = []        #refresh important ###
        BBB = []
        VVV = []
        RRR = []
        III = []
        JJJ = []
        HHH = []
        KKK = []

        U = []        #refresh important ###
        B = []
        V = []
        R = []
        I = []
        J = []
        H = []
        K = []

        U_V = []        #refresh important ###
        B_V = []
        V_R = []
        V_I = []
        V_J = []
        V_H = []
        V_K = []

        kk=0
        for line in infile1:
         kk = kk+1
         ll = line.split()
         UUU.append(ll[4])
         BBB.append(ll[5])
         VVV.append(ll[6])
         RRR.append(ll[7])
         III.append(ll[8])
         JJJ.append(ll[9])
         HHH.append(ll[10])
         KKK.append(ll[11])
        infile1.close()

        for ii in range(1,kk) :
         U.append(float(UUU[ii]))
         B.append(float(BBB[ii]))
         V.append(float(VVV[ii]))
         R.append(float(RRR[ii]))
         I.append(float(III[ii]))
         J.append(float(JJJ[ii]))
         H.append(float(HHH[ii]))
         K.append(float(KKK[ii]))

        for ii in range(0,kk-1) :
         U_V.append(U[ii] - V[ii])
         B_V.append(B[ii] - V[ii])
         V_R.append(V[ii] - R[ii])
         V_I.append(V[ii] - I[ii])
         V_J.append(V[ii] - J[ii])
         V_H.append(V[ii] - H[ii])
         V_K.append(V[ii] - K[ii])


        print( 'Il y a', kk-1, 'clusters dans le fichier' )
        #print 'En voici ces magnitudes, normalisees par la bande V'
        #print str(U_V[0])+' 'str(B_V[0])+' 'str(0.)+' 'str(V_R[0])+' 'str(V_I[0])+' 'str(V_J[0])+' 'str(V_H[0])+' 'str(V_K[0])
        #for ii in range(0,kk) :
        #print U_V[0], B_V[0], 0., V_R[0], V_I[0], V_J[0], V_H[0], V_K[0]
        #print U_V[ii], B_V[ii], 0., V_R[ii], V_I[ii], V_J[ii], V_H[ii], V_K[ii]

        return U_V, B_V, V_R, V_I, V_J, V_H, V_K



def Lecture_sumfile_gCMD(file_name):

        infile1 = open(file_name,'r')

        UUU = []        #refresh important ###
        BBB = []
        VVV = []
        RRR = []
        III = []

        U = []        #refresh important ###
        B = []
        V = []
        R = []
        I = []

        U_V = []        #refresh important ###
        U_B = []
        B_V = []
        V_R = []
        V_I = []

        kk=0
        for line in infile1:
         kk = kk+1
         ll = line.split()
         UUU.append(ll[0])
         BBB.append(ll[1])
         VVV.append(ll[2])
         RRR.append(ll[3])
         III.append(ll[4])
        infile1.close()

        for ii in range(1,kk) :
         U.append(float(UUU[ii]))
         B.append(float(BBB[ii]))
         V.append(float(VVV[ii]))
         R.append(float(RRR[ii]))
         I.append(float(III[ii]))

        for ii in range(0,kk-1) :
         U_V.append(U[ii] - V[ii])
         U_B.append(U[ii] - B[ii])
         B_V.append(B[ii] - V[ii])
         V_R.append(V[ii] - R[ii])
         V_I.append(V[ii] - I[ii])


        print( 'Il y a', kk-1, 'clusters dans le fichier' )
        #print 'En voici ces magnitudes, normalisees par la bande V'
        #print str(U_V[0])+' 'str(B_V[0])+' 'str(0.)+' 'str(V_R[0])+' 'str(V_I[0])+' 'str(V_J[0])+' 'str(V_H[0])+' 'str(V_K[0])
        #for ii in range(0,kk) :
         #print U_V[0], B_V[0], 0., V_R[0], V_I[0], V_J[0], V_H[0], V_K[0]
         #print U_V[ii], B_V[ii], 0., V_R[ii], V_I[ii], V_J[ii], V_H[ii], V_K[ii]

        #return U_V, B_V, V_R, V_I
        return U, B, V, R, I, U_B, B_V




def Numberlines_sumfile(file_name):

        infile1 = open(file_name,'r')

        UUU = []        #refresh important ###
        BBB = []
        VVV = []
        RRR = []
        III = []
        JJJ = []
        HHH = []
        KKK = []

        kk=0
        for line in infile1:
         kk = kk+1
         ll = line.split()
         UUU.append(ll[4])
         BBB.append(ll[5])
         VVV.append(ll[6])
         RRR.append(ll[7])
         III.append(ll[8])
         JJJ.append(ll[9])
         HHH.append(ll[10])
         KKK.append(ll[11])
        infile1.close()

        #print 'Il y a', kk-1, 'clusters dans le fichier'

        return kk-1        #number of lines of file 


def age_mass_Z_May2011(a,m,zz):

        #Old Grid May 2011
        if   a ==1: age = 7.0                #10 Myr
        elif a ==2: age = 7.5                #30 Myr
        elif a ==3: age = 8.0                #100 Myr
        elif a ==4: age = 8.5                #300 Myr
        elif a ==5: age = 9.0                #1000 Myr
        elif a ==6: age = 10.0                #10000 Myr
        
        if   m==1: mass = 1000.                #Solar masses
        elif m==2: mass = 5000.                #Solar masses
        elif m==3: mass = 10000.        #Solar masses
        elif m==4: mass = 50000.        #Solar masses
        elif m==5: mass = 100000.        #Solar masses
        
        if   zz ==1: Z = 0.001                #10 Myr        
        elif zz ==2: Z = 0.004                #30 Myr        
        elif zz ==3: Z = 0.008                #100 Myr        
        elif zz ==4: Z = 0.01                #300 Myr        
        elif zz ==5: Z = 0.02                #1000 Myr
        elif zz ==6: Z = 0.03                #10000 Myr

        return age, mass, Z

def age_mass_Z_June2011(a,m,zz):

        #Old Age grid June 2011
        if   a ==1: age, age_indice = 7.00, '700'        
        elif a ==2: age, age_indice = 7.15, '715'        
        elif a ==3: age, age_indice = 7.30, '730'
        elif a ==4: age, age_indice = 7.45, '745'        
        elif a ==5: age, age_indice = 7.60, '760'        
        elif a ==6: age, age_indice = 7.75, '775'        
        elif a ==7: age, age_indice = 7.90, '790'        
        elif a ==8: age, age_indice = 8.05, '805'        
        elif a ==9: age, age_indice = 8.20, '820'        
        elif a ==10: age, age_indice = 8.35, '835'        
        elif a ==11: age, age_indice = 8.50, '850'        
        elif a ==12: age, age_indice = 8.65, '865'        
        elif a ==13: age, age_indice = 8.80, '880'        
        elif a ==14: age, age_indice = 8.95, '895'        
        elif a ==15: age, age_indice = 9.10, '910'        
        elif a ==16: age, age_indice = 9.25, '925'        
        elif a ==17: age, age_indice = 9.40, '940'        
        elif a ==18: age, age_indice = 9.55, '955'        
        elif a ==19: age, age_indice = 9.70, '970'        
        elif a ==20: age, age_indice = 9.85, '985'        
        elif a ==21: age, age_indice = 10.0, '1000'

        #Old Mass grid June 2011
        if   m ==1: mass, mass_indice = 10**3.00, '300' #Solar masses        
        elif m ==2: mass, mass_indice = 10**3.25, '325'        #Solar masses
        elif m ==3: mass, mass_indice = 10**3.50, '350'        #Solar masses                        
        elif m ==4: mass, mass_indice = 10**3.75, '375'        #Solar masses
        elif m ==5: mass, mass_indice = 10**4.00, '400'        #Solar masses
        elif m ==6: mass, mass_indice = 10**4.25, '425'        #Solar masses
        elif m ==7: mass, mass_indice = 10**4.50, '450'        #Solar masses
        elif m ==8: mass, mass_indice = 10**4.75, '475'        #Solar masses
        elif m ==9: mass, mass_indice = 10**5.00, '500'        #Solar masses

        #Old metallicity grid (June 2011), based on Zs=0.020
        # a tester : [-2.3, -1.9, -1.6, -1.3, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, +0.2]
        #if   zz ==1: Z, Z_indice = 0.00010, 'm23'        # [M/H] = -2.3
        #elif zz ==2: Z, Z_indice = 0.00025, 'm19'        # [M/H] = -1.9
        #elif zz ==3: Z, Z_indice = 0.00050, 'm16'        # [M/H] = -1.6
        #elif zz ==4: Z, Z_indice = 0.00100, 'm13'        # [M/H] = -1.3
        #elif zz ==5: Z, Z_indice = 0.00200, 'm10'        # [M/H] = -1.0
        #elif zz ==6: Z, Z_indice = 0.00317, 'm08'        # [M/H] = -0.8
        #elif zz ==7: Z, Z_indice = 0.00502, 'm06'        # [M/H] = -0.6
        ##elif zz ==8: Z, Z_indice = 0.00796, 'm04'        # [M/H] = -0.4                #BUG!!
        #elif zz ==8: Z, Z_indice = 0.00800, 'm04'        # [M/H] = -0.4
        #elif zz ==9: Z, Z_indice = 0.01262, 'm02'        # [M/H] = -0.2
        #elif zz ==10: Z, Z_indice = 0.02000, 'n00'        # [M/H] = +0.0
        #elif zz ==11: Z, Z_indice = 0.03000, 'p02'        # [M/H] = +0.2 ABOUT, not exactly

        #New metallicities, based on Zs=0.019, 11 metal
        #if   zz ==1:  Z, Z_indice = 0.00010, 'm23'        # [M/H] = -2.3
        #elif zz ==2:  Z, Z_indice = 0.00024, 'm19'        # [M/H] = -1.9
        #elif zz ==3:  Z, Z_indice = 0.00050, 'm16'        # [M/H] = -1.6
        #elif zz ==4:  Z, Z_indice = 0.00095, 'm13'        # [M/H] = -1.3
        #elif zz ==5:  Z, Z_indice = 0.00190, 'm10'        # [M/H] = -1.0
        #elif zz ==6:  Z, Z_indice = 0.00300, 'm08'        # [M/H] = -0.8
        #elif zz ==7:  Z, Z_indice = 0.00480, 'm06'        # [M/H] = -0.6
        #elif zz ==8:  Z, Z_indice = 0.00800, 'm04'        # [M/H] = -0.4                
        #elif zz ==9:  Z, Z_indice = 0.01200, 'm02'        # [M/H] = -0.2
        #elif zz ==10: Z, Z_indice = 0.01900, 'n00'        # [M/H] = +0.0
        #elif zz ==11: Z, Z_indice = 0.03000, 'p02'        # [M/H] = +0.2 
        
        
        #New metallicities, based on Zs=0.019, 13 metal
        if   zz ==1:  Z, Z_indice = 0.00012, 'm22'        # [M/H] = -2.2
        elif zz ==2:  Z, Z_indice = 0.00019, 'm20'        # [M/H] = -2.0
        elif zz ==3:  Z, Z_indice = 0.00030, 'm18'        # [M/H] = -1.8
        elif zz ==4:  Z, Z_indice = 0.00050, 'm16'        # [M/H] = -1.6
        elif zz ==5:  Z, Z_indice = 0.00075, 'm14'        # [M/H] = -1.4
        elif zz ==6:  Z, Z_indice = 0.00120, 'm12'        # [M/H] = -1.2
        elif zz ==7:  Z, Z_indice = 0.00190, 'm10'        # [M/H] = -1.0
        elif zz ==8:  Z, Z_indice = 0.00300, 'm08'        # [M/H] = -0.8
        elif zz ==9:  Z, Z_indice = 0.00480, 'm06'        # [M/H] = -0.6
        elif zz ==10: Z, Z_indice = 0.00800, 'm04'        # [M/H] = -0.4                
        elif zz ==11: Z, Z_indice = 0.01200, 'm02'        # [M/H] = -0.2
        elif zz ==12: Z, Z_indice = 0.01900, 'n00'        # [M/H] = +0.0
        elif zz ==13: Z, Z_indice = 0.03000, 'p02'        # [M/H] = +0.2 

        return age, mass, Z, age_indice, mass_indice, Z_indice


def age_mass_Z_September2011(a,m,zz):
        
        #Age grid September 2011
        if   a ==1: age, age_indice = 7.00, '700'        #10 Myr
        elif a ==2: age, age_indice = 7.10, '710'        
        elif a ==3: age, age_indice = 7.20, '720'
        elif a ==4: age, age_indice = 7.30, '730'        
        elif a ==5: age, age_indice = 7.40, '740'        
        elif a ==6: age, age_indice = 7.50, '750'        
        elif a ==7: age, age_indice = 7.60, '760'        
        elif a ==8: age, age_indice = 7.70, '770'        
        elif a ==9: age, age_indice = 7.80, '780'        
        elif a ==10: age, age_indice = 7.90, '790'        
        elif a ==11: age, age_indice = 8.00, '800'        
        elif a ==12: age, age_indice = 8.10, '810'        
        elif a ==13: age, age_indice = 8.20, '820'        
        elif a ==14: age, age_indice = 8.30, '830'        
        elif a ==15: age, age_indice = 8.40, '840'        
        elif a ==16: age, age_indice = 8.50, '850'        
        elif a ==17: age, age_indice = 8.60, '860'        
        elif a ==18: age, age_indice = 8.70, '870'        
        elif a ==19: age, age_indice = 8.80, '880'        
        elif a ==20: age, age_indice = 8.90, '890'        
        elif a ==21: age, age_indice = 9.00, '900'
        elif a ==22: age, age_indice = 9.10, '910'
        elif a ==23: age, age_indice = 9.20, '920'
        elif a ==24: age, age_indice = 9.30, '930'
        elif a ==25: age, age_indice = 9.40, '940'
        elif a ==26: age, age_indice = 9.50, '950'
        elif a ==27: age, age_indice = 9.60, '960'
        elif a ==28: age, age_indice = 9.70, '970'
        elif a ==29: age, age_indice = 9.80, '980'
        elif a ==30: age, age_indice = 9.90, '990'
        elif a ==31: age, age_indice = 10.0, '1000'
        
        #Mass grid September 2011
        if   m ==1 : mass, mass_indice = 10**3.00, '300'         #Solar masses        
        elif m ==2 : mass, mass_indice = 10**3.10, '310'        #Solar masses
        elif m ==3 : mass, mass_indice = 10**3.20, '320'        #Solar masses                        
        elif m ==4 : mass, mass_indice = 10**3.30, '330'        #Solar masses
        elif m ==5 : mass, mass_indice = 10**3.40, '340'        #Solar masses
        elif m ==6 : mass, mass_indice = 10**3.50, '350'        #Solar masses
        elif m ==7 : mass, mass_indice = 10**3.60, '360'        #Solar masses
        elif m ==8 : mass, mass_indice = 10**3.70, '370'        #Solar masses
        elif m ==9 : mass, mass_indice = 10**3.80, '380'        #Solar masses
        elif m ==10: mass, mass_indice = 10**3.90, '390'        #Solar masses
        elif m ==11: mass, mass_indice = 10**4.00, '400'        #Solar masses
        elif m ==12: mass, mass_indice = 10**4.10, '410'        #Solar masses
        elif m ==13: mass, mass_indice = 10**4.20, '420'        #Solar masses
        elif m ==14: mass, mass_indice = 10**4.30, '430'        #Solar masses
        elif m ==15: mass, mass_indice = 10**4.40, '440'        #Solar masses
        elif m ==16: mass, mass_indice = 10**4.50, '450'        #Solar masses
        elif m ==17: mass, mass_indice = 10**4.60, '460'        #Solar masses
        elif m ==18: mass, mass_indice = 10**4.70, '470'        #Solar masses
        elif m ==19: mass, mass_indice = 10**4.80, '480'        #Solar masses
        elif m ==20: mass, mass_indice = 10**4.90, '490'        #Solar masses
        elif m ==21: mass, mass_indice = 10**5.00, '500'        #Solar masses


        #September 2011, end
        #if   zz ==1:  Z, Z_indice = 0.00010    , 'm23'        # [M/H] = -2.3
        #elif zz ==2:  Z, Z_indice = 0.00012619 , 'm22'        # [M/H] = -2.2
        #elif zz ==3:  Z, Z_indice = 0.00015887 , 'm21'        # [M/H] = -2.1
        #elif zz ==4:  Z, Z_indice = 0.00020000 , 'm20'        # [M/H] = -2.0
        #elif zz ==5:  Z, Z_indice = 0.00025    , 'm19'        # [M/H] = -1.9
        #elif zz ==6:  Z, Z_indice = 0.00031698 , 'm18'        # [M/H] = -1.8
        #elif zz ==7:  Z, Z_indice = 0.00039905 , 'm17'        # [M/H] = -1.7
        #elif zz ==8:  Z, Z_indice = 0.00050    , 'm16'        # [M/H] = -1.6
        #elif zz ==9:  Z, Z_indice = 0.00063246 , 'm15'        # [M/H] = -1.5
        #elif zz ==10: Z, Z_indice = 0.00079621 , 'm14'        # [M/H] = -1.4
        #elif zz ==11: Z, Z_indice = 0.00100    , 'm13'        # [M/H] = -1.3
        #elif zz ==12: Z, Z_indice = 0.0012619  , 'm12'        # [M/H] = -1.2
        #elif zz ==13: Z, Z_indice = 0.0015887  , 'm11'        # [M/H] = -1.1
        #elif zz ==14: Z, Z_indice = 0.00200    , 'm10'        # [M/H] = -1.0
        #elif zz ==15: Z, Z_indice = 0.0025179  , 'm09'        # [M/H] = -0.9
        #elif zz ==16: Z, Z_indice = 0.00317    , 'm08'        # [M/H] = -0.8
        #elif zz ==17: Z, Z_indice = 0.0039905  , 'm07'        # [M/H] = -0.7
        #elif zz ==18: Z, Z_indice = 0.00502    , 'm06'        # [M/H] = -0.6
        #elif zz ==19: Z, Z_indice = 0.0063246  , 'm05'        # [M/H] = -0.5
        #elif zz ==20: Z, Z_indice = 0.00796    , 'm04'        # [M/H] = -0.4
        #elif zz ==21: Z, Z_indice = 0.010024   , 'm03'        # [M/H] = -0.3
        #elif zz ==22: Z, Z_indice = 0.01262    , 'm02'        # [M/H] = -0.2
        #elif zz ==23: Z, Z_indice = 0.015887   , 'm01'        # [M/H] = -0.1
        #elif zz ==24: Z, Z_indice = 0.02000    , 'n00'        # [M/H] = +0.0
        #elif zz ==25: Z, Z_indice = 0.025179   , 'p01'        # [M/H] = +0.1
        #elif zz ==26: Z, Z_indice = 0.03000    , 'p02'        # [M/H] = +0.2 ABOUT, not exactly




        #September 2011, end
        if   zz ==1:  Z, Z_indice = 0.00010    , 'm23'        # [M/H] = -2.3
        elif zz ==2:  Z, Z_indice = 0.00013    , 'm22'        # [M/H] = -2.2    #Should be 0.00012
        elif zz ==3:  Z, Z_indice = 0.00015    , 'm21'        # [M/H] = -2.1
        elif zz ==4:  Z, Z_indice = 0.00019    , 'm20'        # [M/H] = -2.0
        elif zz ==5:  Z, Z_indice = 0.00024    , 'm19'        # [M/H] = -1.9
        elif zz ==6:  Z, Z_indice = 0.00030    , 'm18'        # [M/H] = -1.8
        elif zz ==7:  Z, Z_indice = 0.00040    , 'm17'        # [M/H] = -1.7
        elif zz ==8:  Z, Z_indice = 0.00050    , 'm16'        # [M/H] = -1.6
        elif zz ==9:  Z, Z_indice = 0.00060    , 'm15'        # [M/H] = -1.5
        elif zz ==10: Z, Z_indice = 0.00080    , 'm14'        # [M/H] = -1.4    #Should be 0.00075
        elif zz ==11: Z, Z_indice = 0.00100    , 'm13'        # [M/H] = -1.3
        elif zz ==12: Z, Z_indice = 0.00130    , 'm12'        # [M/H] = -1.2    #Should be 0.00120
        elif zz ==13: Z, Z_indice = 0.00150    , 'm11'        # [M/H] = -1.1
        elif zz ==14: Z, Z_indice = 0.00190    , 'm10'        # [M/H] = -1.0
        elif zz ==15: Z, Z_indice = 0.00240    , 'm09'        # [M/H] = -0.9
        elif zz ==16: Z, Z_indice = 0.00300    , 'm08'        # [M/H] = -0.8
        elif zz ==17: Z, Z_indice = 0.00400    , 'm07'        # [M/H] = -0.7    #Should be 0.00380
        elif zz ==18: Z, Z_indice = 0.00500    , 'm06'        # [M/H] = -0.6    #Should be 0.00480
        elif zz ==19: Z, Z_indice = 0.00600    , 'm05'        # [M/H] = -0.5
        elif zz ==20: Z, Z_indice = 0.00800    , 'm04'        # [M/H] = -0.4
        elif zz ==21: Z, Z_indice = 0.01000    , 'm03'        # [M/H] = -0.3
        elif zz ==22: Z, Z_indice = 0.01300    , 'm02'        # [M/H] = -0.2    #Should be 0.01200
        elif zz ==23: Z, Z_indice = 0.01500    , 'm01'        # [M/H] = -0.1
        elif zz ==24: Z, Z_indice = 0.01900    , 'n00'        # [M/H] = +0.0
        elif zz ==25: Z, Z_indice = 0.02400    , 'p01'        # [M/H] = +0.1
        elif zz ==26: Z, Z_indice = 0.03000    , 'p02'        # [M/H] = +0.2 


        return age, mass, Z, age_indice, mass_indice, Z_indice


def age_mass_Z_November2011(a,m,zz):
        
        #Age grid Novtember 2011
        if   a ==1: age, age_indice = 7.00, '700'        #10 Myr
        elif a ==2: age, age_indice = 7.10, '710'        
        elif a ==3: age, age_indice = 7.20, '720'
        elif a ==4: age, age_indice = 7.30, '730'        
        elif a ==5: age, age_indice = 7.40, '740'        
        elif a ==6: age, age_indice = 7.50, '750'        
        elif a ==7: age, age_indice = 7.60, '760'        
        elif a ==8: age, age_indice = 7.70, '770'        
        elif a ==9: age, age_indice = 7.80, '780'        
        elif a ==10: age, age_indice = 7.90, '790'        
        elif a ==11: age, age_indice = 8.00, '800'        
        elif a ==12: age, age_indice = 8.10, '810'        
        elif a ==13: age, age_indice = 8.20, '820'        
        elif a ==14: age, age_indice = 8.30, '830'        
        elif a ==15: age, age_indice = 8.40, '840'        
        elif a ==16: age, age_indice = 8.50, '850'        
        elif a ==17: age, age_indice = 8.60, '860'        
        elif a ==18: age, age_indice = 8.70, '870'        
        elif a ==19: age, age_indice = 8.80, '880'        
        elif a ==20: age, age_indice = 8.90, '890'        
        elif a ==21: age, age_indice = 9.00, '900'
        elif a ==22: age, age_indice = 9.10, '910'
        elif a ==23: age, age_indice = 9.20, '920'
        elif a ==24: age, age_indice = 9.30, '930'
        elif a ==25: age, age_indice = 9.40, '940'
        elif a ==26: age, age_indice = 9.50, '950'
        elif a ==27: age, age_indice = 9.60, '960'
        elif a ==28: age, age_indice = 9.70, '970'
        elif a ==29: age, age_indice = 9.80, '980'
        elif a ==30: age, age_indice = 9.90, '990'
        elif a ==31: age, age_indice = 10.0, '1000'
        
        #Mass grid November 2011
        if   m ==1 : mass, mass_indice = 10**2.00, '200'         #Solar masses        
        elif m ==2 : mass, mass_indice = 10**2.10, '210'        #Solar masses
        elif m ==3 : mass, mass_indice = 10**2.20, '220'        #Solar masses                        
        elif m ==4 : mass, mass_indice = 10**2.30, '230'        #Solar masses
        elif m ==5 : mass, mass_indice = 10**2.40, '240'        #Solar masses
        elif m ==6 : mass, mass_indice = 10**2.50, '250'        #Solar masses
        elif m ==7 : mass, mass_indice = 10**2.60, '260'        #Solar masses
        elif m ==8 : mass, mass_indice = 10**2.70, '270'        #Solar masses
        elif m ==9 : mass, mass_indice = 10**2.80, '280'        #Solar masses
        elif m ==10: mass, mass_indice = 10**2.90, '290'        #Solar masses
        elif m ==11: mass, mass_indice = 10**3.00, '300'         #Solar masses        
        elif m ==12: mass, mass_indice = 10**3.10, '310'        #Solar masses
        elif m ==13: mass, mass_indice = 10**3.20, '320'        #Solar masses                        
        elif m ==14: mass, mass_indice = 10**3.30, '330'        #Solar masses
        elif m ==15: mass, mass_indice = 10**3.40, '340'        #Solar masses
        elif m ==16: mass, mass_indice = 10**3.50, '350'        #Solar masses
        elif m ==17: mass, mass_indice = 10**3.60, '360'        #Solar masses
        elif m ==18: mass, mass_indice = 10**3.70, '370'        #Solar masses
        elif m ==19: mass, mass_indice = 10**3.80, '380'        #Solar masses
        elif m ==20: mass, mass_indice = 10**3.90, '390'        #Solar masses
        elif m ==21: mass, mass_indice = 10**4.00, '400'        #Solar masses
        elif m ==22: mass, mass_indice = 10**4.10, '410'        #Solar masses
        elif m ==23: mass, mass_indice = 10**4.20, '420'        #Solar masses
        elif m ==24: mass, mass_indice = 10**4.30, '430'        #Solar masses
        elif m ==25: mass, mass_indice = 10**4.40, '440'        #Solar masses
        elif m ==26: mass, mass_indice = 10**4.50, '450'        #Solar masses
        elif m ==27: mass, mass_indice = 10**4.60, '460'        #Solar masses
        elif m ==28: mass, mass_indice = 10**4.70, '470'        #Solar masses
        elif m ==29: mass, mass_indice = 10**4.80, '480'        #Solar masses
        elif m ==30: mass, mass_indice = 10**4.90, '490'        #Solar masses
        elif m ==31: mass, mass_indice = 10**5.00, '500'        #Solar masses
        elif m ==32: mass, mass_indice = 10**5.10, '510'        #Solar masses
        elif m ==33: mass, mass_indice = 10**5.20, '520'        #Solar masses
        elif m ==34: mass, mass_indice = 10**5.30, '530'        #Solar masses
        elif m ==35: mass, mass_indice = 10**5.40, '540'        #Solar masses
        elif m ==36: mass, mass_indice = 10**5.50, '550'        #Solar masses
        elif m ==37: mass, mass_indice = 10**5.60, '560'        #Solar masses
        elif m ==38: mass, mass_indice = 10**5.70, '570'        #Solar masses
        elif m ==39: mass, mass_indice = 10**5.80, '580'        #Solar masses
        elif m ==40: mass, mass_indice = 10**5.90, '590'        #Solar masses
        elif m ==41: mass, mass_indice = 10**6.00, '600'        #Solar masses


        #November 2011, end
        if   zz ==1:  Z, Z_indice = 0.00010    , 'm23'        # [M/H] = -2.3
        elif zz ==2:  Z, Z_indice = 0.00013    , 'm22'        # [M/H] = -2.2    #Should be 0.00012
        elif zz ==3:  Z, Z_indice = 0.00015    , 'm21'        # [M/H] = -2.1
        elif zz ==4:  Z, Z_indice = 0.00019    , 'm20'        # [M/H] = -2.0
        elif zz ==5:  Z, Z_indice = 0.00024    , 'm19'        # [M/H] = -1.9
        elif zz ==6:  Z, Z_indice = 0.00030    , 'm18'        # [M/H] = -1.8
        elif zz ==7:  Z, Z_indice = 0.00040    , 'm17'        # [M/H] = -1.7
        elif zz ==8:  Z, Z_indice = 0.00050    , 'm16'        # [M/H] = -1.6
        elif zz ==9:  Z, Z_indice = 0.00060    , 'm15'        # [M/H] = -1.5
        elif zz ==10: Z, Z_indice = 0.00080    , 'm14'        # [M/H] = -1.4    #Should be 0.00075
        elif zz ==11: Z, Z_indice = 0.00100    , 'm13'        # [M/H] = -1.3
        elif zz ==12: Z, Z_indice = 0.00130    , 'm12'        # [M/H] = -1.2    #Should be 0.00120
        elif zz ==13: Z, Z_indice = 0.00150    , 'm11'        # [M/H] = -1.1
        elif zz ==14: Z, Z_indice = 0.00190    , 'm10'        # [M/H] = -1.0
        elif zz ==15: Z, Z_indice = 0.00240    , 'm09'        # [M/H] = -0.9
        elif zz ==16: Z, Z_indice = 0.00300    , 'm08'        # [M/H] = -0.8
        elif zz ==17: Z, Z_indice = 0.00400    , 'm07'        # [M/H] = -0.7    #Should be 0.00380
        elif zz ==18: Z, Z_indice = 0.00500    , 'm06'        # [M/H] = -0.6    #Should be 0.00480
        elif zz ==19: Z, Z_indice = 0.00600    , 'm05'        # [M/H] = -0.5
        elif zz ==20: Z, Z_indice = 0.00800    , 'm04'        # [M/H] = -0.4
        elif zz ==21: Z, Z_indice = 0.01000    , 'm03'        # [M/H] = -0.3
        elif zz ==22: Z, Z_indice = 0.01300    , 'm02'        # [M/H] = -0.2    #Should be 0.01200
        elif zz ==23: Z, Z_indice = 0.01500    , 'm01'        # [M/H] = -0.1
        elif zz ==24: Z, Z_indice = 0.01900    , 'n00'        # [M/H] = +0.0
        elif zz ==25: Z, Z_indice = 0.02400    , 'p01'        # [M/H] = +0.1
        elif zz ==26: Z, Z_indice = 0.03000    , 'p02'        # [M/H] = +0.2 


        return age, mass, Z, age_indice, mass_indice, Z_indice





def age_mass_Z_December2011(a,m,zz):
        
        #Age grid December 2011
        '''
        if   a ==1: age, age_indice = 7.00, '700'        #10 Myr
        elif a ==2: age, age_indice = 7.05, '705'        
        elif a ==3: age, age_indice = 7.10, '710'
        elif a ==4: age, age_indice = 7.15, '715'        
        elif a ==5: age, age_indice = 7.20, '720'
        elif a ==6: age, age_indice = 7.25, '725'
        elif a ==7: age, age_indice = 7.30, '730'
        elif a ==8: age, age_indice = 7.35, '735'        
        elif a ==9: age, age_indice = 7.40, '740'
        elif a ==10: age, age_indice = 7.45, '745'        
        elif a ==11: age, age_indice = 7.50, '750'        
        elif a ==12: age, age_indice = 7.55, '755'
        elif a ==13: age, age_indice = 7.60, '760'        
        elif a ==14: age, age_indice = 7.65, '765'
        elif a ==15: age, age_indice = 7.70, '770'
        elif a ==16: age, age_indice = 7.75, '775'        
        elif a ==17: age, age_indice = 7.80, '780'
        elif a ==18: age, age_indice = 7.85, '785'        
        elif a ==19: age, age_indice = 7.90, '790'
        elif a ==20: age, age_indice = 7.95, '795'        
        elif a ==21: age, age_indice = 8.00, '800'
        elif a ==22: age, age_indice = 8.05, '805'        
        elif a ==23: age, age_indice = 8.10, '810'
        elif a ==24: age, age_indice = 8.15, '815'        
        elif a ==25: age, age_indice = 8.20, '820'
        elif a ==26: age, age_indice = 8.25, '825'        
        elif a ==27: age, age_indice = 8.30, '830'
        elif a ==28: age, age_indice = 8.35, '835'        
        elif a ==29: age, age_indice = 8.40, '840'
        elif a ==30: age, age_indice = 8.45, '845'        
        elif a ==31: age, age_indice = 8.50, '850'
        elif a ==32: age, age_indice = 8.55, '855'        
        elif a ==33: age, age_indice = 8.60, '860'        
        elif a ==34: age, age_indice = 8.65, '865'
        elif a ==35: age, age_indice = 8.70, '870'
        elif a ==36: age, age_indice = 8.75, '875'        
        elif a ==37: age, age_indice = 8.80, '880'        
        elif a ==38: age, age_indice = 8.85, '885'
        elif a ==39: age, age_indice = 8.90, '890'
        elif a ==40: age, age_indice = 8.95, '895'        
        elif a ==41: age, age_indice = 9.00, '900'
        elif a ==42: age, age_indice = 9.05, '905'
        elif a ==43: age, age_indice = 9.10, '910'
        elif a ==44: age, age_indice = 9.15, '915'
        elif a ==45: age, age_indice = 9.20, '920'
        elif a ==46: age, age_indice = 9.25, '925'
        elif a ==47: age, age_indice = 9.30, '930'
        elif a ==48: age, age_indice = 9.35, '935'
        elif a ==49: age, age_indice = 9.40, '940'
        elif a ==50: age, age_indice = 9.45, '945'
        elif a ==51: age, age_indice = 9.50, '950'
        elif a ==52: age, age_indice = 9.55, '955'
        elif a ==53: age, age_indice = 9.60, '960'
        elif a ==54: age, age_indice = 9.65, '965'
        elif a ==55: age, age_indice = 9.70, '970'
        elif a ==56: age, age_indice = 9.75, '975'
        elif a ==57: age, age_indice = 9.80, '980'
        elif a ==58: age, age_indice = 9.85, '985'
        elif a ==59: age, age_indice = 9.90, '990'
        elif a ==60: age, age_indice = 9.95, '995'
        elif a ==61: age, age_indice = 10.0, '1000'
        elif a ==62: age, age_indice = 10.05, '1005'
        elif a ==63: age, age_indice = 10.10, '1010'
        '''

        #Age grid December 2011
        age  = 6.60
        #for ii in range(1,67):
        # if (a == ii) :
        age = age + (a-1)*0.05
        age_int = int(round(age*100))
        age_indice = str(age_int)
        
        #Mass grid December 2011
        mass  = 2.00
        #for ii in range(1,82):
        # if (m == ii) :
        mass = mass + (m-1)*0.05
        mass_int = int(round(mass*100.))
        mass = 10**mass
        mass_indice = str(mass_int)
        '''
        #Mass grid December 2011
        if   m ==1: mass, mass_indice = 10**2.00, '200'         #Solar masses        
        elif m ==2: mass, mass_indice = 10**2.05, '205'         #Solar masses        
        elif m ==3: mass, mass_indice = 10**2.10, '210'                #Solar masses
        elif m ==4: mass, mass_indice = 10**2.15, '215'         #Solar masses        
        elif m ==5: mass, mass_indice = 10**2.20, '220'                #Solar masses        
        elif m ==6: mass, mass_indice = 10**2.25, '225'         #Solar masses                        
        elif m ==7: mass, mass_indice = 10**2.30, '230'                #Solar masses
        elif m ==8: mass, mass_indice = 10**2.35, '235'         #Solar masses        
        elif m ==9: mass, mass_indice = 10**2.40, '240'                #Solar masses
        elif m ==10: mass, mass_indice = 10**2.45, '245'         #Solar masses        
        elif m ==11: mass, mass_indice = 10**2.50, '250'        #Solar masses
        elif m ==12: mass, mass_indice = 10**2.55, '255'        #Solar masses
        elif m ==13: mass, mass_indice = 10**2.60, '260'        #Solar masses
        elif m ==14: mass, mass_indice = 10**2.65, '265'        #Solar masses
        elif m ==15: mass, mass_indice = 10**2.70, '270'        #Solar masses
        elif m ==16: mass, mass_indice = 10**2.75, '275'        #Solar masses
        elif m ==17: mass, mass_indice = 10**2.80, '280'        #Solar masses
        elif m ==18: mass, mass_indice = 10**2.85, '285'        #Solar masses
        elif m ==19: mass, mass_indice = 10**2.90, '290'        #Solar masses
        elif m ==20: mass, mass_indice = 10**2.95, '295'        #Solar masses
        elif m ==21: mass, mass_indice = 10**3.00, '300'         #Solar masses        
        elif m ==22: mass, mass_indice = 10**3.05, '305'        #Solar masses
        elif m ==23: mass, mass_indice = 10**3.10, '310'        #Solar masses
        elif m ==24: mass, mass_indice = 10**3.15, '315'        #Solar masses
        elif m ==25: mass, mass_indice = 10**3.20, '320'        #Solar masses        
        elif m ==26: mass, mass_indice = 10**3.25, '325'        #Solar masses                
        elif m ==27: mass, mass_indice = 10**3.30, '330'        #Solar masses
        elif m ==28: mass, mass_indice = 10**3.35, '335'        #Solar masses
        elif m ==29: mass, mass_indice = 10**3.40, '340'        #Solar masses
        elif m ==30: mass, mass_indice = 10**3.45, '345'        #Solar masses
        elif m ==31: mass, mass_indice = 10**3.50, '350'        #Solar masses
        elif m ==32: mass, mass_indice = 10**3.55, '355'        #Solar masses
        elif m ==33: mass, mass_indice = 10**3.60, '360'        #Solar masses
        elif m ==34: mass, mass_indice = 10**3.65, '365'        #Solar masses
        elif m ==35: mass, mass_indice = 10**3.70, '370'        #Solar masses
        elif m ==36: mass, mass_indice = 10**3.75, '375'        #Solar masses
        elif m ==37: mass, mass_indice = 10**3.80, '380'        #Solar masses
        elif m ==38: mass, mass_indice = 10**3.85, '385'        #Solar masses
        elif m ==39: mass, mass_indice = 10**3.90, '390'        #Solar masses
        elif m ==40: mass, mass_indice = 10**3.95, '395'        #Solar masses
        elif m ==41: mass, mass_indice = 10**4.00, '400'        #Solar masses
        elif m ==42: mass, mass_indice = 10**4.05, '405'        #Solar masses
        elif m ==43: mass, mass_indice = 10**4.10, '410'        #Solar masses
        elif m ==44: mass, mass_indice = 10**4.15, '415'        #Solar masses
        elif m ==45: mass, mass_indice = 10**4.20, '420'        #Solar masses
        elif m ==46: mass, mass_indice = 10**4.25, '425'        #Solar masses
        elif m ==47: mass, mass_indice = 10**4.30, '430'        #Solar masses
        elif m ==48: mass, mass_indice = 10**4.35, '435'        #Solar masses
        elif m ==49: mass, mass_indice = 10**4.40, '440'        #Solar masses
        elif m ==50: mass, mass_indice = 10**4.45, '445'        #Solar masses
        elif m ==51: mass, mass_indice = 10**4.50, '450'        #Solar masses
        elif m ==52: mass, mass_indice = 10**4.55, '455'        #Solar masses
        elif m ==53: mass, mass_indice = 10**4.60, '460'        #Solar masses
        elif m ==54: mass, mass_indice = 10**4.65, '465'        #Solar masses
        elif m ==55: mass, mass_indice = 10**4.70, '470'        #Solar masses
        elif m ==56: mass, mass_indice = 10**4.75, '475'        #Solar masses
        elif m ==57: mass, mass_indice = 10**4.80, '480'        #Solar masses
        elif m ==58: mass, mass_indice = 10**4.85, '485'        #Solar masses
        elif m ==59: mass, mass_indice = 10**4.90, '490'        #Solar masses
        elif m ==60: mass, mass_indice = 10**4.95, '495'        #Solar masses
        elif m ==61: mass, mass_indice = 10**5.00, '500'        #Solar masses
        elif m ==62: mass, mass_indice = 10**5.05, '505'        #Solar masses
        elif m ==63: mass, mass_indice = 10**5.10, '510'        #Solar masses
        elif m ==64: mass, mass_indice = 10**5.15, '515'        #Solar masses
        elif m ==65: mass, mass_indice = 10**5.20, '520'        #Solar masses
        elif m ==66: mass, mass_indice = 10**5.25, '525'        #Solar masses
        elif m ==67: mass, mass_indice = 10**5.30, '530'        #Solar masses
        elif m ==68: mass, mass_indice = 10**5.35, '535'        #Solar masses
        elif m ==69: mass, mass_indice = 10**5.40, '540'        #Solar masses
        elif m ==70: mass, mass_indice = 10**5.45, '545'        #Solar masses
        elif m ==71: mass, mass_indice = 10**5.50, '550'        #Solar masses
        elif m ==72: mass, mass_indice = 10**5.55, '555'        #Solar masses
        elif m ==73: mass, mass_indice = 10**5.60, '560'        #Solar masses
        elif m ==74: mass, mass_indice = 10**5.65, '565'        #Solar masses
        elif m ==75: mass, mass_indice = 10**5.70, '570'        #Solar masses
        elif m ==76: mass, mass_indice = 10**5.75, '575'        #Solar masses
        elif m ==77: mass, mass_indice = 10**5.80, '580'        #Solar masses
        elif m ==78: mass, mass_indice = 10**5.85, '585'        #Solar masses
        elif m ==79: mass, mass_indice = 10**5.90, '590'        #Solar masses
        elif m ==80: mass, mass_indice = 10**5.95, '595'        #Solar masses
        elif m ==81: mass, mass_indice = 10**6.00, '600'        #Solar masses
        '''

        #November 2011, end
        #'''
        if   zz ==1:  Z, Z_indice = 0.00010    , 'm23'        # [M/H] = -2.3
        elif zz ==2:  Z, Z_indice = 0.00013    , 'm22'        # [M/H] = -2.2    #Should be 0.00012
        elif zz ==3:  Z, Z_indice = 0.00015    , 'm21'        # [M/H] = -2.1
        elif zz ==4:  Z, Z_indice = 0.00019    , 'm20'        # [M/H] = -2.0
        elif zz ==5:  Z, Z_indice = 0.00024    , 'm19'        # [M/H] = -1.9
        elif zz ==6:  Z, Z_indice = 0.00030    , 'm18'        # [M/H] = -1.8
        elif zz ==7:  Z, Z_indice = 0.00040    , 'm17'        # [M/H] = -1.7
        elif zz ==8:  Z, Z_indice = 0.00050    , 'm16'        # [M/H] = -1.6
        elif zz ==9:  Z, Z_indice = 0.00060    , 'm15'        # [M/H] = -1.5
        elif zz ==10: Z, Z_indice = 0.00080    , 'm14'        # [M/H] = -1.4    #Should be 0.00075
        elif zz ==11: Z, Z_indice = 0.00100    , 'm13'        # [M/H] = -1.3
        elif zz ==12: Z, Z_indice = 0.00130    , 'm12'        # [M/H] = -1.2    #Should be 0.00120
        elif zz ==13: Z, Z_indice = 0.00150    , 'm11'        # [M/H] = -1.1
        elif zz ==14: Z, Z_indice = 0.00190    , 'm10'        # [M/H] = -1.0
        elif zz ==15: Z, Z_indice = 0.00240    , 'm09'        # [M/H] = -0.9
        elif zz ==16: Z, Z_indice = 0.00300    , 'm08'        # [M/H] = -0.8
        elif zz ==17: Z, Z_indice = 0.00400    , 'm07'        # [M/H] = -0.7    #Should be 0.00380
        elif zz ==18: Z, Z_indice = 0.00500    , 'm06'        # [M/H] = -0.6    #Should be 0.00480
        elif zz ==19: Z, Z_indice = 0.00600    , 'm05'        # [M/H] = -0.5
        elif zz ==20: Z, Z_indice = 0.00800    , 'm04'        # [M/H] = -0.4
        elif zz ==21: Z, Z_indice = 0.01000    , 'm03'        # [M/H] = -0.3
        elif zz ==22: Z, Z_indice = 0.01300    , 'm02'        # [M/H] = -0.2    #Should be 0.01200
        elif zz ==23: Z, Z_indice = 0.01500    , 'm01'        # [M/H] = -0.1
        elif zz ==24: Z, Z_indice = 0.01900    , 'n00'        # [M/H] = +0.0
        elif zz ==25: Z, Z_indice = 0.02400    , 'p01'        # [M/H] = +0.1
        elif zz ==26: Z, Z_indice = 0.03000    , 'p02'        # [M/H] = +0.2 
        elif zz ==27: Z, Z_indice = 0.04000    , 'p03'        # [M/H] = +0.3
        elif zz ==28: Z, Z_indice = 0.05000    , 'p04'        # [M/H] = +0.4 

        return age, mass, Z, age_indice, mass_indice, Z_indice

def mass_to_mass_indice(mass): #Mass should be entered in LOG! 4 in place of 10000Ms 
        mass_int = int(round(mass*100.))
        mass_indice = str(mass_int)
        return mass_indice


def age_mass_Z_December2011_inverse(age, mass, Z, age_indice, mass_indice, Z_indice):

        #Age grid December 2011
        age = (age - 6.60)/0.05
        a = int(round(age)) + 1
        age_indice = (float(age_indice)/100. - 6.60)/0.05
        a_indice = int(round(age_indice)) + 1
        #if a != a_indice:
        # print a, a_indice
        # print 'error'
        # raw_input()

        #Mass grid December 2011
        mass = (math.log10(mass) - 2.0)/0.05
        m = int(round(mass)) + 1
        mass_indice = (float(mass_indice)/100. - 2.0)/0.05
        m_indice = int(round(mass_indice)) + 1
        #if m != m_indice:
        # print m, m_indice
        # print 'error'
        # raw_input()

        #November 2011, end
        if    Z == 0.00010 or Z_indice == 'm23' : zz = 1        # [M/H] = -2.3
        elif  Z == 0.00013 or Z_indice == 'm22' : zz = 2        # [M/H] = -2.2    #Should be 0.00012
        elif  Z == 0.00015 or Z_indice == 'm21' : zz = 3        # [M/H] = -2.1
        elif  Z == 0.00019 or Z_indice == 'm20' : zz = 4        # [M/H] = -2.0
        elif  Z == 0.00024 or Z_indice == 'm19' : zz = 5        # [M/H] = -1.9
        elif  Z == 0.00030 or Z_indice == 'm18' : zz = 6        # [M/H] = -1.8
        elif  Z == 0.00040 or Z_indice == 'm17' : zz = 7        # [M/H] = -1.7
        elif  Z == 0.00050 or Z_indice == 'm16' : zz = 8        # [M/H] = -1.6
        elif  Z == 0.00060 or Z_indice == 'm15' : zz = 9        # [M/H] = -1.5
        elif  Z == 0.00080 or Z_indice == 'm14' : zz = 10        # [M/H] = -1.4    #Should be 0.00075
        elif  Z == 0.00100 or Z_indice == 'm13' : zz = 11        # [M/H] = -1.3
        elif  Z == 0.00130 or Z_indice == 'm12' : zz = 12        # [M/H] = -1.2    #Should be 0.00120
        elif  Z == 0.00150 or Z_indice == 'm11' : zz = 13        # [M/H] = -1.1
        elif  Z == 0.00190 or Z_indice == 'm10' : zz = 14        # [M/H] = -1.0
        elif  Z == 0.00240 or Z_indice == 'm09' : zz = 15        # [M/H] = -0.9
        elif  Z == 0.00300 or Z_indice == 'm08' : zz = 16        # [M/H] = -0.8
        elif  Z == 0.00400 or Z_indice == 'm07' : zz = 17        # [M/H] = -0.7    #Should be 0.00380
        elif  Z == 0.00500 or Z_indice == 'm06' : zz = 18        # [M/H] = -0.6    #Should be 0.00480
        elif  Z == 0.00600 or Z_indice == 'm05' : zz = 19        # [M/H] = -0.5
        elif  Z == 0.00800 or Z_indice == 'm04' : zz = 20        # [M/H] = -0.4
        elif  Z == 0.01000 or Z_indice == 'm03' : zz = 21        # [M/H] = -0.3
        elif  Z == 0.01300 or Z_indice == 'm02' : zz = 22        # [M/H] = -0.2    #Should be 0.01200
        elif  Z == 0.01500 or Z_indice == 'm01' : zz = 23        # [M/H] = -0.1
        elif  Z == 0.01900 or Z_indice == 'n00' : zz = 24        # [M/H] = +0.0
        elif  Z == 0.02400 or Z_indice == 'p01' : zz = 25        # [M/H] = +0.1
        elif  Z == 0.03000 or Z_indice == 'p02' : zz = 26        # [M/H] = +0.2 
        elif  Z == 0.04000 or Z_indice == 'p03' : zz = 27        # [M/H] = +0.3
        elif  Z == 0.05000 or Z_indice == 'p04' : zz = 28        # [M/H] = +0.4 

        #return age, mass, Z, age_indice, mass_indice, Z_indice
        return a, m, zz


def indices_to_aa_mm_zz(age_indice, mass_indice, Z_indice):

        #Age grid December 2011
        #age = (age - 6.60)/0.05
        #a = int(round(age)) + 1
        age_indice = (float(age_indice)/100. - 6.60)/0.05
        a = int(round(age_indice)) + 1
        #if a != a_indice:
        # print a, a_indice
        # print 'error'
        # raw_input()

        #Mass grid December 2011
        #mass = (math.log10(mass) - 2.0)/0.05
        #m = int(round(mass)) + 1
        mass_indice = (float(mass_indice)/100. - 2.0)/0.05
        m = int(round(mass_indice)) + 1
        #if m != m_indice:
        # print m, m_indice
        # print 'error'
        # raw_input()

        #November 2011, end
        if    Z_indice == 'm23' : zz = 1        # [M/H] = -2.3
        elif  Z_indice == 'm22' : zz = 2        # [M/H] = -2.2    #Should be 0.00012
        elif  Z_indice == 'm21' : zz = 3        # [M/H] = -2.1
        elif  Z_indice == 'm20' : zz = 4        # [M/H] = -2.0
        elif  Z_indice == 'm19' : zz = 5        # [M/H] = -1.9
        elif  Z_indice == 'm18' : zz = 6        # [M/H] = -1.8
        elif  Z_indice == 'm17' : zz = 7        # [M/H] = -1.7
        elif  Z_indice == 'm16' : zz = 8        # [M/H] = -1.6
        elif  Z_indice == 'm15' : zz = 9        # [M/H] = -1.5
        elif  Z_indice == 'm14' : zz = 10        # [M/H] = -1.4    #Should be 0.00075
        elif  Z_indice == 'm13' : zz = 11        # [M/H] = -1.3
        elif  Z_indice == 'm12' : zz = 12        # [M/H] = -1.2    #Should be 0.00120
        elif  Z_indice == 'm11' : zz = 13        # [M/H] = -1.1
        elif  Z_indice == 'm10' : zz = 14        # [M/H] = -1.0
        elif  Z_indice == 'm09' : zz = 15        # [M/H] = -0.9
        elif  Z_indice == 'm08' : zz = 16        # [M/H] = -0.8
        elif  Z_indice == 'm07' : zz = 17        # [M/H] = -0.7    #Should be 0.00380
        elif  Z_indice == 'm06' : zz = 18        # [M/H] = -0.6    #Should be 0.00480
        elif  Z_indice == 'm05' : zz = 19        # [M/H] = -0.5
        elif  Z_indice == 'm04' : zz = 20        # [M/H] = -0.4
        elif  Z_indice == 'm03' : zz = 21        # [M/H] = -0.3
        elif  Z_indice == 'm02' : zz = 22        # [M/H] = -0.2    #Should be 0.01200
        elif  Z_indice == 'm01' : zz = 23        # [M/H] = -0.1
        elif  Z_indice == 'n00' : zz = 24        # [M/H] = +0.0
        elif  Z_indice == 'p01' : zz = 25        # [M/H] = +0.1
        elif  Z_indice == 'p02' : zz = 26        # [M/H] = +0.2 
        elif  Z_indice == 'p03' : zz = 27        # [M/H] = +0.3
        elif  Z_indice == 'p04' : zz = 28        # [M/H] = +0.4 

        #return age, mass, Z, age_indice, mass_indice, Z_indice
        return a, m, zz




def Zindex_to_Z_and_zz(Z_indice):

        if    Z_indice == 'm23' : Z,zz = 0.00010, 1        # [M/H] = -2.3
        elif  Z_indice == 'm22' : Z,zz = 0.00013, 2        # [M/H] = -2.2    #Should be 0.00012
        elif  Z_indice == 'm21' : Z,zz = 0.00015, 3        # [M/H] = -2.1
        elif  Z_indice == 'm20' : Z,zz = 0.00019, 4        # [M/H] = -2.0
        elif  Z_indice == 'm19' : Z,zz = 0.00024, 5        # [M/H] = -1.9
        elif  Z_indice == 'm18' : Z,zz = 0.00030, 6        # [M/H] = -1.8
        elif  Z_indice == 'm17' : Z,zz = 0.00040, 7        # [M/H] = -1.7
        elif  Z_indice == 'm16' : Z,zz = 0.00050, 8        # [M/H] = -1.6
        elif  Z_indice == 'm15' : Z,zz = 0.00060, 9        # [M/H] = -1.5
        elif  Z_indice == 'm14' : Z,zz = 0.00080, 10        # [M/H] = -1.4    #Should be 0.00075
        elif  Z_indice == 'm13' : Z,zz = 0.00100, 11        # [M/H] = -1.3
        elif  Z_indice == 'm12' : Z,zz = 0.00130, 12        # [M/H] = -1.2    #Should be 0.00120
        elif  Z_indice == 'm11' : Z,zz = 0.00150, 13        # [M/H] = -1.1
        elif  Z_indice == 'm10' : Z,zz = 0.00190, 14        # [M/H] = -1.0
        elif  Z_indice == 'm09' : Z,zz = 0.00240, 15        # [M/H] = -0.9
        elif  Z_indice == 'm08' : Z,zz = 0.00300, 16        # [M/H] = -0.8
        elif  Z_indice == 'm07' : Z,zz = 0.00400, 17        # [M/H] = -0.7    #Should be 0.00380
        elif  Z_indice == 'm06' : Z,zz = 0.00500, 18        # [M/H] = -0.6    #Should be 0.00480
        elif  Z_indice == 'm05' : Z,zz = 0.00600, 19        # [M/H] = -0.5
        elif  Z_indice == 'm04' : Z,zz = 0.00800, 20        # [M/H] = -0.4
        elif  Z_indice == 'm03' : Z,zz = 0.01000, 21        # [M/H] = -0.3
        elif  Z_indice == 'm02' : Z,zz = 0.01300, 22        # [M/H] = -0.2    #Should be 0.01200
        elif  Z_indice == 'm01' : Z,zz = 0.01500, 23        # [M/H] = -0.1
        elif  Z_indice == 'n00' : Z,zz = 0.01900, 24        # [M/H] = +0.0
        elif  Z_indice == 'p01' : Z,zz = 0.02400, 25        # [M/H] = +0.1
        elif  Z_indice == 'p02' : Z,zz = 0.03000, 26        # [M/H] = +0.2 
        elif  Z_indice == 'p03' : Z,zz = 0.04000, 27        # [M/H] = +0.3
        elif  Z_indice == 'p04' : Z,zz = 0.05000, 28        # [M/H] = +0.4 

        return Z,zz



def Z_to_Zindex_and_zz(Z):

        if    Z==0.00010  : Z_indice ,zz = 'm23', 1        # [M/H] = -2.3
        elif  Z==0.00013  : Z_indice ,zz = 'm22', 2        # [M/H] = -2.2    #Should be 0.00012
        elif  Z==0.00015  : Z_indice ,zz = 'm21', 3        # [M/H] = -2.1
        elif  Z==0.00019  : Z_indice ,zz = 'm20', 4        # [M/H] = -2.0
        elif  Z==0.00024  : Z_indice ,zz = 'm19', 5        # [M/H] = -1.9
        elif  Z==0.00030  : Z_indice ,zz = 'm18', 6        # [M/H] = -1.8
        elif  Z==0.00040  : Z_indice ,zz = 'm17', 7        # [M/H] = -1.7
        elif  Z==0.00050  : Z_indice ,zz = 'm16', 8        # [M/H] = -1.6
        elif  Z==0.00060  : Z_indice ,zz = 'm15', 9        # [M/H] = -1.5
        elif  Z==0.00080  : Z_indice ,zz = 'm14', 10        # [M/H] = -1.4    #Should be 0.00075
        elif  Z==0.00100  : Z_indice ,zz = 'm13', 11        # [M/H] = -1.3
        elif  Z==0.00130  : Z_indice ,zz = 'm12', 12        # [M/H] = -1.2    #Should be 0.00120
        elif  Z==0.00150  : Z_indice ,zz = 'm11', 13        # [M/H] = -1.1
        elif  Z==0.00190  : Z_indice ,zz = 'm10', 14        # [M/H] = -1.0
        elif  Z==0.00240  : Z_indice ,zz = 'm09', 15        # [M/H] = -0.9
        elif  Z==0.00300  : Z_indice ,zz = 'm08', 16        # [M/H] = -0.8
        elif  Z==0.00400  : Z_indice ,zz = 'm07', 17        # [M/H] = -0.7    #Should be 0.00380
        elif  Z==0.00500  : Z_indice ,zz = 'm06', 18        # [M/H] = -0.6    #Should be 0.00480
        elif  Z==0.00600  : Z_indice ,zz = 'm05', 19        # [M/H] = -0.5
        elif  Z==0.00800  : Z_indice ,zz = 'm04', 20        # [M/H] = -0.4
        elif  Z==0.01000  : Z_indice ,zz = 'm03', 21        # [M/H] = -0.3
        elif  Z==0.01300  : Z_indice ,zz = 'm02', 22        # [M/H] = -0.2    #Should be 0.01200
        elif  Z==0.01500  : Z_indice ,zz = 'm01', 23        # [M/H] = -0.1
        elif  Z==0.01900  : Z_indice ,zz = 'n00', 24        # [M/H] = +0.0
        elif  Z==0.02400  : Z_indice ,zz = 'p01', 25        # [M/H] = +0.1
        elif  Z==0.03000  : Z_indice ,zz = 'p02', 26        # [M/H] = +0.2 
        elif  Z==0.04000  : Z_indice ,zz = 'p03', 27        # [M/H] = +0.3
        elif  Z==0.05000  : Z_indice ,zz = 'p04', 28        # [M/H] = +0.4 

        return Z_indice,zz


def Z_to_Zindex_and_zz_and_jj(Z):

        if    Z==0.00010  : Z_indice ,zz,jj = 'm23', 1 , 99        # [M/H] = -2.3
        elif  Z==0.00013  : Z_indice ,zz,jj = 'm22', 2 , 12        # [M/H] = -2.2    #Should be 0.00012
        elif  Z==0.00015  : Z_indice ,zz,jj = 'm21', 3 , 99        # [M/H] = -2.1
        elif  Z==0.00019  : Z_indice ,zz,jj = 'm20', 4 , 11        # [M/H] = -2.0
        elif  Z==0.00024  : Z_indice ,zz,jj = 'm19', 5 , 99        # [M/H] = -1.9
        elif  Z==0.00030  : Z_indice ,zz,jj = 'm18', 6 , 10        # [M/H] = -1.8
        elif  Z==0.00040  : Z_indice ,zz,jj = 'm17', 7 , 99        # [M/H] = -1.7
        elif  Z==0.00050  : Z_indice ,zz,jj = 'm16', 8 , 9        # [M/H] = -1.6
        elif  Z==0.00060  : Z_indice ,zz,jj = 'm15', 9 , 99        # [M/H] = -1.5
        elif  Z==0.00080  : Z_indice ,zz,jj = 'm14', 10, 8        # [M/H] = -1.4    #Should be 0.00075
        elif  Z==0.00100  : Z_indice ,zz,jj = 'm13', 11, 99        # [M/H] = -1.3
        elif  Z==0.00130  : Z_indice ,zz,jj = 'm12', 12, 7        # [M/H] = -1.2    #Should be 0.00120
        elif  Z==0.00150  : Z_indice ,zz,jj = 'm11', 13, 99        # [M/H] = -1.1
        elif  Z==0.00190  : Z_indice ,zz,jj = 'm10', 14, 6        # [M/H] = -1.0
        elif  Z==0.00240  : Z_indice ,zz,jj = 'm09', 15, 99        # [M/H] = -0.9
        elif  Z==0.00300  : Z_indice ,zz,jj = 'm08', 16, 5        # [M/H] = -0.8
        elif  Z==0.00400  : Z_indice ,zz,jj = 'm07', 17, 99        # [M/H] = -0.7    #Should be 0.00380
        elif  Z==0.00500  : Z_indice ,zz,jj = 'm06', 18, 4        # [M/H] = -0.6    #Should be 0.00480
        elif  Z==0.00600  : Z_indice ,zz,jj = 'm05', 19, 99        # [M/H] = -0.5
        elif  Z==0.00800  : Z_indice ,zz,jj = 'm04', 20, 3        # [M/H] = -0.4
        elif  Z==0.01000  : Z_indice ,zz,jj = 'm03', 21, 99        # [M/H] = -0.3
        elif  Z==0.01300  : Z_indice ,zz,jj = 'm02', 22, 2        # [M/H] = -0.2    #Should be 0.01200
        elif  Z==0.01500  : Z_indice ,zz,jj = 'm01', 23, 99        # [M/H] = -0.1
        elif  Z==0.01900  : Z_indice ,zz,jj = 'n00', 24, 1        # [M/H] = +0.0
        elif  Z==0.02400  : Z_indice ,zz,jj = 'p01', 25, 99        # [M/H] = +0.1
        elif  Z==0.03000  : Z_indice ,zz,jj = 'p02', 26, 0        # [M/H] = +0.2 
        elif  Z==0.04000  : Z_indice ,zz,jj = 'p03', 27, 99        # [M/H] = +0.3
        elif  Z==0.05000  : Z_indice ,zz,jj = 'p04', 28, 99        # [M/H] = +0.4 

        return Z_indice,zz,jj


def zz_to_Z_and_Zindex(zz):

        if   zz ==1:  Z, Z_indice = 0.00010    , 'm23'        # [M/H] = -2.3
        elif zz ==2:  Z, Z_indice = 0.00013    , 'm22'        # [M/H] = -2.2    #Should be 0.00012
        elif zz ==3:  Z, Z_indice = 0.00015    , 'm21'        # [M/H] = -2.1
        elif zz ==4:  Z, Z_indice = 0.00019    , 'm20'        # [M/H] = -2.0
        elif zz ==5:  Z, Z_indice = 0.00024    , 'm19'        # [M/H] = -1.9
        elif zz ==6:  Z, Z_indice = 0.00030    , 'm18'        # [M/H] = -1.8
        elif zz ==7:  Z, Z_indice = 0.00040    , 'm17'        # [M/H] = -1.7
        elif zz ==8:  Z, Z_indice = 0.00050    , 'm16'        # [M/H] = -1.6
        elif zz ==9:  Z, Z_indice = 0.00060    , 'm15'        # [M/H] = -1.5
        elif zz ==10: Z, Z_indice = 0.00080    , 'm14'        # [M/H] = -1.4    #Should be 0.00075
        elif zz ==11: Z, Z_indice = 0.00100    , 'm13'        # [M/H] = -1.3
        elif zz ==12: Z, Z_indice = 0.00130    , 'm12'        # [M/H] = -1.2    #Should be 0.00120
        elif zz ==13: Z, Z_indice = 0.00150    , 'm11'        # [M/H] = -1.1
        elif zz ==14: Z, Z_indice = 0.00190    , 'm10'        # [M/H] = -1.0
        elif zz ==15: Z, Z_indice = 0.00240    , 'm09'        # [M/H] = -0.9
        elif zz ==16: Z, Z_indice = 0.00300    , 'm08'        # [M/H] = -0.8
        elif zz ==17: Z, Z_indice = 0.00400    , 'm07'        # [M/H] = -0.7    #Should be 0.00380
        elif zz ==18: Z, Z_indice = 0.00500    , 'm06'        # [M/H] = -0.6    #Should be 0.00480
        elif zz ==19: Z, Z_indice = 0.00600    , 'm05'        # [M/H] = -0.5
        elif zz ==20: Z, Z_indice = 0.00800    , 'm04'        # [M/H] = -0.4
        elif zz ==21: Z, Z_indice = 0.01000    , 'm03'        # [M/H] = -0.3
        elif zz ==22: Z, Z_indice = 0.01300    , 'm02'        # [M/H] = -0.2    #Should be 0.01200
        elif zz ==23: Z, Z_indice = 0.01500    , 'm01'        # [M/H] = -0.1
        elif zz ==24: Z, Z_indice = 0.01900    , 'n00'        # [M/H] = +0.0
        elif zz ==25: Z, Z_indice = 0.02400    , 'p01'        # [M/H] = +0.1
        elif zz ==26: Z, Z_indice = 0.03000    , 'p02'        # [M/H] = +0.2 
        elif zz ==27: Z, Z_indice = 0.04000    , 'p03'        # [M/H] = +0.3
        elif zz ==28: Z, Z_indice = 0.05000    , 'p04'        # [M/H] = +0.4 

        return Z,Z_indice


def node_from_NPZ(aa,mm,zz):  #Function to extract a node from a NPZ grid
    U,B,V = np.zeros(1000), np.zeros(1000), np.zeros(1000) 
    age, mass, Z, age_indice, mass_indice, Z_indice = Lecture_module.age_mass_Z_December2011(aa,mm,zz)
    #Open the npz node file   
    path_in = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grid_FRS_allZ_Kroupa_NPZ/NPZ_files/{0}/'.format(Z_indice)
    node = np.load(path_in+'Clusters_t{0}_M{1}_Z{2}.npz'.format(age_indice,mass_indice,Z_indice))

    #Extract the means, covariances and weights
    means       = node['means']
    covariances = node['covariances']
    weights     = node['weights']
    components_number = len(weights)

    #Build the node using means, covariances and weights
    #Here small block to ensure that we build 1000 models per bin
    int_weights = np.zeros(10, dtype=int)
    for ii in range(0,components_number-1):
     int_weights[ii] = int(round(weights[ii]*1000,0))   
    if np.sum(int_weights) != 1000: #If more or less than 1000models, remove or add the exceeding/missing models randomly
     if np.sum(int_weights) > 1000: sign = -1  #remove some models
     if np.sum(int_weights) < 1000: sign = 1   #add some models
     excess = abs(1000-np.sum(int_weights))
     for ex in range(0,excess):
      random_integer = np.random.randint(0,10)
      int_weights[random_integer] = int_weights[random_integer] + sign*1

    #Building of nodes
    if components_number ==1:
     data = np.random.multivariate_normal(means[0],covariances[0],int_weights)   
    if components_number >1:
     for ii in range(0,components_number):
      if ii==0:
       data = np.random.multivariate_normal(means[0,:],covariances[0],int_weights[0])             
      if ii>0:
       data_new = np.random.multivariate_normal(means[ii,:],covariances[ii],int_weights[ii])             
       data = np.vstack((data,data_new)) 
       del data_new 
    
    #Here i show how to extract the UBV from the data of the node                    
    #U[:] = data[:,2]
    #B[:] = data[:,3]
    #V[:] = data[:,4] 
    return data #U,B,V #temporary, i can extract all the photometry if needed


def age_mass_Z_March2012(a,m,zz):


        #Age grid March 2012
        age  = 6.60
        age = age + (a-1)*0.01
        age_int = int(round(age*100))
        age_indice = str(age_int)
        
        #Mass grid March 2012
        mass  = 2.00
        mass = mass + (m-1)*0.01
        mass_int = int(round(mass*100.))
        mass = 10**mass
        mass_indice = str(mass_int)

        #November 2011, end

        if   zz ==1:  Z, Z_indice = 0.00010    , 'm23'        # [M/H] = -2.3
        elif zz ==2:  Z, Z_indice = 0.00013    , 'm22'        # [M/H] = -2.2    #Should be 0.00012
        elif zz ==3:  Z, Z_indice = 0.00015    , 'm21'        # [M/H] = -2.1
        elif zz ==4:  Z, Z_indice = 0.00019    , 'm20'        # [M/H] = -2.0
        elif zz ==5:  Z, Z_indice = 0.00024    , 'm19'        # [M/H] = -1.9
        elif zz ==6:  Z, Z_indice = 0.00030    , 'm18'        # [M/H] = -1.8
        elif zz ==7:  Z, Z_indice = 0.00040    , 'm17'        # [M/H] = -1.7
        elif zz ==8:  Z, Z_indice = 0.00050    , 'm16'        # [M/H] = -1.6
        elif zz ==9:  Z, Z_indice = 0.00060    , 'm15'        # [M/H] = -1.5
        elif zz ==10: Z, Z_indice = 0.00080    , 'm14'        # [M/H] = -1.4    #Should be 0.00075
        elif zz ==11: Z, Z_indice = 0.00100    , 'm13'        # [M/H] = -1.3
        elif zz ==12: Z, Z_indice = 0.00130    , 'm12'        # [M/H] = -1.2    #Should be 0.00120
        elif zz ==13: Z, Z_indice = 0.00150    , 'm11'        # [M/H] = -1.1
        elif zz ==14: Z, Z_indice = 0.00190    , 'm10'        # [M/H] = -1.0
        elif zz ==15: Z, Z_indice = 0.00240    , 'm09'        # [M/H] = -0.9
        elif zz ==16: Z, Z_indice = 0.00300    , 'm08'        # [M/H] = -0.8
        elif zz ==17: Z, Z_indice = 0.00400    , 'm07'        # [M/H] = -0.7    #Should be 0.00380
        elif zz ==18: Z, Z_indice = 0.00500    , 'm06'        # [M/H] = -0.6    #Should be 0.00480
        elif zz ==19: Z, Z_indice = 0.00600    , 'm05'        # [M/H] = -0.5
        elif zz ==20: Z, Z_indice = 0.00800    , 'm04'        # [M/H] = -0.4
        elif zz ==21: Z, Z_indice = 0.01000    , 'm03'        # [M/H] = -0.3
        elif zz ==22: Z, Z_indice = 0.01300    , 'm02'        # [M/H] = -0.2    #Should be 0.01200
        elif zz ==23: Z, Z_indice = 0.01500    , 'm01'        # [M/H] = -0.1
        elif zz ==24: Z, Z_indice = 0.01900    , 'n00'        # [M/H] = +0.0
        elif zz ==25: Z, Z_indice = 0.02400    , 'p01'        # [M/H] = +0.1
        elif zz ==26: Z, Z_indice = 0.03000    , 'p02'        # [M/H] = +0.2 
        elif zz ==27: Z, Z_indice = 0.04000    , 'p03'        # [M/H] = +0.3
        elif zz ==28: Z, Z_indice = 0.05000    , 'p04'        # [M/H] = +0.4 


        return age, mass, Z, age_indice, mass_indice, Z_indice


def age_mass_Z_March2012_inverse(age, mass, Z, age_indice, mass_indice, Z_indice):

        #Age
        age = (age - 6.60)/0.01
        a = int(round(age)) + 1
        age_indice = (float(age_indice)/100. - 6.60)/0.01
        a_indice = int(round(age_indice)) + 1
        #if a != a_indice:
        # print a, a_indice
        # print 'error'
        # raw_input()

        #Mass
        mass = (math.log10(mass) - 2.0)/0.01
        m = int(round(mass)) + 1
        mass_indice = (float(mass_indice)/100. - 2.0)/0.01
        m_indice = int(round(mass_indice)) + 1
        #if m != m_indice:
        # print m, m_indice
        # print 'error'
        # raw_input()

        #November 2011, end
        if    Z == 0.00010 or Z_indice == 'm23' : zz = 1        # [M/H] = -2.3
        elif  Z == 0.00013 or Z_indice == 'm22' : zz = 2        # [M/H] = -2.2    #Should be 0.00012
        elif  Z == 0.00015 or Z_indice == 'm21' : zz = 3        # [M/H] = -2.1
        elif  Z == 0.00019 or Z_indice == 'm20' : zz = 4        # [M/H] = -2.0
        elif  Z == 0.00024 or Z_indice == 'm19' : zz = 5        # [M/H] = -1.9
        elif  Z == 0.00030 or Z_indice == 'm18' : zz = 6        # [M/H] = -1.8
        elif  Z == 0.00040 or Z_indice == 'm17' : zz = 7        # [M/H] = -1.7
        elif  Z == 0.00050 or Z_indice == 'm16' : zz = 8        # [M/H] = -1.6
        elif  Z == 0.00060 or Z_indice == 'm15' : zz = 9        # [M/H] = -1.5
        elif  Z == 0.00080 or Z_indice == 'm14' : zz = 10        # [M/H] = -1.4    #Should be 0.00075
        elif  Z == 0.00100 or Z_indice == 'm13' : zz = 11        # [M/H] = -1.3
        elif  Z == 0.00130 or Z_indice == 'm12' : zz = 12        # [M/H] = -1.2    #Should be 0.00120
        elif  Z == 0.00150 or Z_indice == 'm11' : zz = 13        # [M/H] = -1.1
        elif  Z == 0.00190 or Z_indice == 'm10' : zz = 14        # [M/H] = -1.0
        elif  Z == 0.00240 or Z_indice == 'm09' : zz = 15        # [M/H] = -0.9
        elif  Z == 0.00300 or Z_indice == 'm08' : zz = 16        # [M/H] = -0.8
        elif  Z == 0.00400 or Z_indice == 'm07' : zz = 17        # [M/H] = -0.7    #Should be 0.00380
        elif  Z == 0.00500 or Z_indice == 'm06' : zz = 18        # [M/H] = -0.6    #Should be 0.00480
        elif  Z == 0.00600 or Z_indice == 'm05' : zz = 19        # [M/H] = -0.5
        elif  Z == 0.00800 or Z_indice == 'm04' : zz = 20        # [M/H] = -0.4
        elif  Z == 0.01000 or Z_indice == 'm03' : zz = 21        # [M/H] = -0.3
        elif  Z == 0.01300 or Z_indice == 'm02' : zz = 22        # [M/H] = -0.2    #Should be 0.01200
        elif  Z == 0.01500 or Z_indice == 'm01' : zz = 23        # [M/H] = -0.1
        elif  Z == 0.01900 or Z_indice == 'n00' : zz = 24        # [M/H] = +0.0
        elif  Z == 0.02400 or Z_indice == 'p01' : zz = 25        # [M/H] = +0.1
        elif  Z == 0.03000 or Z_indice == 'p02' : zz = 26        # [M/H] = +0.2 
        elif  Z == 0.04000 or Z_indice == 'p03' : zz = 27        # [M/H] = +0.3
        elif  Z == 0.05000 or Z_indice == 'p04' : zz = 28        # [M/H] = +0.4 
        return a, m, zz


def Filters_Index_Name(Filter_index):
        if Filter_index == 1:   name,e_name = 'FUV','e_FUV'                         # Galex
        if Filter_index == 2:   name,e_name = 'NUV','e_NUV'                         # Galex
        if Filter_index == 3:   name,e_name = 'U','e_U'                        # Johnson  
        if Filter_index == 4:   name,e_name = 'B','e_B'                        # Johnson  
        if Filter_index == 5:   name,e_name = 'V','e_V'                        # Johnson  
        if Filter_index == 6:   name,e_name = 'R','e_R'                        # Cousins  
        if Filter_index == 7:   name,e_name = 'I','e_I'                        # Cousins  
        if Filter_index == 8:   name,e_name = 'J','e_J'                        # Bessel 
        if Filter_index == 9:   name,e_name = 'H','e_H'                        # Bessel 
        if Filter_index == 10:  name,e_name = 'K','e_K'                        # Bessel 
        if Filter_index == 11:  name,e_name = 'F275W','e_275'                # WFC3
        if Filter_index == 12:  name,e_name = 'F336W','e_336'                # WFC3
        if Filter_index == 13:  name,e_name = 'F475W','e_475'                # WFC3
        if Filter_index == 14:  name,e_name = 'F814W','e_814'                # WFC3
        if Filter_index == 15:  name,e_name = 'F110W','e_110'                # WFC3
        if Filter_index == 16:  name,e_name = 'F160W','e_160'                # WFC3
        if Filter_index == 17:  name,e_name = 'u_CFHT','e_uCFHT'                # CFHT     
        if Filter_index == 18:  name,e_name = 'g_CFHT','e_gCFHT'                # CFHT     
        if Filter_index == 19:  name,e_name = 'r_CFHT','e_rCFHT'                # CFHT     
        if Filter_index == 20:  name,e_name = 'i_CFHT','e_iCFHT'                # CFHT     
        if Filter_index == 21:  name,e_name = 'z_CFHT','e_zCFHT'                # CFHT     
        if Filter_index == 22:  name,e_name = 'u_SDSS','e_uSDSS'                # SDSS   
        if Filter_index == 23:  name,e_name = 'g_SDSS','e_gSDSS'                  # SDSS   
        if Filter_index == 24:  name,e_name = 'r_SDSS','e_rSDSS'                # SDSS   
        if Filter_index == 25:  name,e_name = 'i_SDSS','e_iSDSS'                # SDSS   
        if Filter_index == 26:  name,e_name = 'z_SDSS','e_zSDSS'                # SDSS   
        if Filter_index == 27:  name,e_name = 'J_2MASS','e_J2MASS'                 # 2MASS
        if Filter_index == 28:  name,e_name = 'H_2MASS','e_H2MASS'                 # 2MASS
        if Filter_index == 29:  name,e_name = 'Ks_2MASS','e_Ks2MASS'                # 2MASS
        if Filter_index == 30:  name,e_name = 'a_BATC','e_aBATC'                   # BATC
        if Filter_index == 31:  name,e_name = 'b_BATC','e_bBATC'                 # BATC
        if Filter_index == 32:  name,e_name = 'c_BATC','e_cBATC'                 # BATC
        if Filter_index == 33:  name,e_name = 'd_BATC','e_dBATC'                 # BATC
        if Filter_index == 34:  name,e_name = 'e_BATC','e_eBATC'                   # BATC
        if Filter_index == 35:  name,e_name = 'f_BATC','e_fBATC'                   # BATC
        if Filter_index == 36:  name,e_name = 'g_BATC','e_gBATC'                   # BATC
        if Filter_index == 37:  name,e_name = 'h_BATC','e_hBATC'                   # BATC
        if Filter_index == 38:  name,e_name = 'i_BATC','e_iBATC'                   # BATC
        if Filter_index == 39:  name,e_name = 'j_BATC','e_jBATC'                   # BATC
        if Filter_index == 40:  name,e_name = 'k_BATC','e_kBATC'                   # BATC
        if Filter_index == 41:  name,e_name = 'm_BATC','e_mBATC'                   # BATC
        if Filter_index == 42:  name,e_name = 'n_BATC','e_nBATC'                   # BATC
        if Filter_index == 43:  name,e_name = 'o_BATC','e_oBATC'                   # BATC
        if Filter_index == 44:  name,e_name = 'p_BATC','e_pBATC'                   # BATC
        if Filter_index == 45:  name,e_name = 't_BATC','e_tBATC'                   # BATC
        if Filter_index == 46:  name,e_name = 'IRAC36','e_36'                         # Spitzer
        if Filter_index == 47:  name,e_name = 'IRAC45','e_45'                         # Spitzer
        if Filter_index == 48:  name,e_name = 'IRAC58','e_58'                         # Spitzer
        if Filter_index == 49:  name,e_name = 'IRAC80','e_80'                         # Spitzer
        if Filter_index == 50:  name,e_name = 'mips24','e_m24'                         # Spitzer
        if Filter_index == 51:  name,e_name = 'mips70','e_m70'                         # Spitzer
        if Filter_index == 52:  name,e_name = 'mips160','e_m10'                 # Spitzer
        return name,e_name

def Filters_Name_Index(name):
        if name == 'FUV'          : Filter_index = 1  # Galex
        if name == 'NUV'          : Filter_index = 2  # Galex
        if name == 'U'         : Filter_index = 3  # Johnson  
        if name == 'B'         : Filter_index = 4  # Johnson  
        if name == 'V'         : Filter_index = 5  # Johnson  
        if name == 'R'         : Filter_index = 6  # Cousins  
        if name == 'I'         : Filter_index = 7  # Cousins  
        if name == 'J'         : Filter_index = 8  # Bessel 
        if name == 'H'         : Filter_index = 9  # Bessel 
        if name == 'K'         : Filter_index = 10 # Bessel 
        if name == 'F275W'           : Filter_index = 11 # WFC3
        if name == 'F336W'           : Filter_index = 12 # WFC3
        if name == 'F475W'           : Filter_index = 13 # WFC3
        if name == 'F814W'           : Filter_index = 14 # WFC3
        if name == 'F110W'           : Filter_index = 15 # WFC3
        if name == 'F160W'           : Filter_index = 16 # WFC3
        if name == 'u_CFHT'         : Filter_index = 17 # CFHT     
        if name == 'g_CFHT'         : Filter_index = 18 # CFHT     
        if name == 'r_CFHT'         : Filter_index = 19 # CFHT     
        if name == 'i_CFHT'         : Filter_index = 20 # CFHT     
        if name == 'z_CFHT'         : Filter_index = 21 # CFHT     
        if name == 'u_SDSS'         : Filter_index = 22 # SDSS   
        if name == 'g_SDSS'           : Filter_index = 23 # SDSS   
        if name == 'r_SDSS'         : Filter_index = 24 # SDSS   
        if name == 'i_SDSS'         : Filter_index = 25 # SDSS   
        if name == 'z_SDSS'         : Filter_index = 26 # SDSS   
        if name == 'J_2MASS'         : Filter_index = 27 # 2MASS
        if name == 'H_2MASS'          : Filter_index = 28 # 2MASS
        if name == 'Ks_2MASS'         : Filter_index = 29 # 2MASS
        if name == 'a_BATC'         : Filter_index = 30 # BATC
        if name == 'b_BATC'         : Filter_index = 31 # BATC
        if name == 'c_BATC'         : Filter_index = 32 # BATC
        if name == 'd_BATC'         : Filter_index = 33 # BATC
        if name == 'e_BATC'           : Filter_index = 34 # BATC
        if name == 'f_BATC'           : Filter_index = 35 # BATC
        if name == 'g_BATC'           : Filter_index = 36 # BATC
        if name == 'h_BATC'           : Filter_index = 37 # BATC
        if name == 'i_BATC'           : Filter_index = 38 # BATC
        if name == 'j_BATC'           : Filter_index = 39 # BATC
        if name == 'k_BATC'           : Filter_index = 40 # BATC
        if name == 'm_BATC'           : Filter_index = 41 # BATC
        if name == 'n_BATC'           : Filter_index = 42 # BATC
        if name == 'o_BATC'           : Filter_index = 43 # BATC
        if name == 'p_BATC'           : Filter_index = 44 # BATC
        if name == 't_BATC'           : Filter_index = 45 # BATC
        if name == 'IRAC36'         : Filter_index = 46 # Spitzer
        if name == 'IRAC45'         : Filter_index = 47 # Spitzer
        if name == 'IRAC58'         : Filter_index = 48 # Spitzer
        if name == 'IRAC80'         : Filter_index = 49 # Spitzer
        if name == 'mips24'         : Filter_index = 50 # Spitzer
        if name == 'mips70'         : Filter_index = 51 # Spitzer
        if name == 'mips160'         : Filter_index = 52 # Spitzer
        return Filter_index









#=============================================================================================================

def imf_function(m):
 #IMF in number unit

 #Popescu_IMF(m):
 #k = 11.836901316 #1./0.0844815694  
 #if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)
 #if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.3) * m**(-1.3)
 #if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.35) * m**(-2.35)
 #if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.35) * m**(-2.35)

 #myGrid_IMF(m): #Kroupa01CB
 #k = 26.578777009 #1./0.03762400353
 #if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)                                 # imf = 12.5558*m**(-0.3)        
 #if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.8) * m**(-1.8)                          # imf = 0.284104*m**(-1.8)        
 #if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.8) * (1./0.5)**(-2.7) * m**(-2.7)         # imf = 0.152248*m**(-2.7)    
 #if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.8) * (1./0.5)**(-2.7) * m**(-2.3)         # imf = 0.152248*m**(-2.3) 

 #Kroupa01(m): #Canonical IMF, or Kr01NCB
 #k = 1 #Not normalized
 #k = 4.2391287 #1./0.09053131901837103  #Number normalization, so that CDF_imf_function_number(0.01,150.) = 1
 k = 11.045901141 #1./0.09053131901837103   #Mass normalization, so that CDF_imf_function_mass(0.01,150.) = 1    
 if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)                                 # imf = 5.177595858*m**(-0.3)        
 if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.3) * m**(-1.3)                                  # imf = 0.414207669*m**(-1.3)        
 if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3) * m**(-2.3)         # imf = 0.207103834*m**(-2.3)    
 if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3) * m**(-2.3)         # imf = 0.207103834*m**(-2.3)  

 #Kroupa02CB(m): #Kroupa 2002 corrected for binaries
 #k = 14.724042932 #1./0.06791612905804065
 #if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)                                 # imf = 6.901668114*m**(-0.3)        
 #if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.3) * m**(-1.3)                          # imf = 0.552133449*m**(-1.3)        
 #if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3) * m**(-2.3)         # imf = 0.276066725*m**(-2.3)    
 #if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3) * m**(-2.7)         # imf = 0.276066725*m**(-2.7) 

 #Kroupa93(m):  #Kroupa 1993 
 #k = 14.229724031 #1./0.07027543175430642
 #if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)                                 # imf = 6.66996375*m**(-0.3)        
 #if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.3) * m**(-1.3)                          # imf = 0.5335971*m**(-1.3)        
 #if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.2) * m**(-2.2)         # imf = 0.285947606*m**(-2.2)    
 #if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.2) * m**(-2.7)         # imf = 0.285947606*m**(-2.7) 

 #Strange: before, for Kr93, i used this:                where is the mistake??? 
 #if(x >= 0.01 and x < 0.08):     imf = 16.98*x**(-0.3)        
 #if(x >= 0.08 and x < 0.5):      imf = 0.5518*x**(-1.3)        
 #if(x >= 0.5 and x < 1.):        imf = 0.2956*x**(-2.2)    
 #if(x >= 1.):                    imf = 0.2956*x**(-2.7)  

 #return m*imf #To check that integral of IMF is really 1 (for a mass normalization of the IMF)
 return imf


def imf_function_m(m):  
        #IMF in mass unit
        return m*imf_function(m)


def CDF_imf_function_mass(m_low, m):                #CDF of the MASS IMF
        args = ()
        Integral_IMF = integrate.quad(imf_function_m, m_low, m, args)
        return Integral_IMF
        #print CDF_imf_function(150.) #The integral of the IMF over the whole mass range is the total mass 


def CDF_imf_function_number(m_low, m):         #CDF of the NUMBER IMF!!
        args = ()
        Integral_IMF = integrate.quad(imf_function, m_low, m, args)
        return Integral_IMF
        #print CDF_imf_function(150.) #The integral of the IMF over the whole mass range is 1, of course. 


    
def Analytical_nCDF(m_inf,m_sup):    
        k = 11.045901141 #1./0.09053131901837103   #Mass normalization, so that CDF_imf_function_mass(0.01,150.) = 1    
        C0 = k * (1./0.08)**(-0.3)
        C1 = k * (1./0.08)**(-1.3)
        C2 = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3)
        C3 = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3)
        alpha0 = 0.3
        alpha1 = 1.3
        alpha2 = 2.3
        alpha3 = 2.3
        m0 = 0.01
        m1 = 0.08
        m2 = 0.5
        m3 = 1.

        if(m_inf >= 0.01 and m_inf < 0.08 ):   nCDF_inf = C0/(1-alpha0) * (m_inf**(1-alpha0) - m0**(1-alpha0))
        if(m_inf >= 0.08 and m_inf < 0.5  ):   nCDF_inf = C0/(1-alpha0) * (m1**(1-alpha0) - m0**(1-alpha0))   +   C1/(1-alpha1) * (m_inf**(1-alpha1) - m1**(1-alpha1))        
        if(m_inf >= 0.5  and m_inf < 1.   ):   nCDF_inf = C0/(1-alpha0) * (m1**(1-alpha0) - m0**(1-alpha0))   +   C1/(1-alpha1) * (m2**(1-alpha1) - m1**(1-alpha1))   +   C2/(1-alpha2) * (m_inf**(1-alpha2) - m2**(1-alpha2))  
        if(m_inf >= 1.   and m_inf <= 150.):   nCDF_inf = C0/(1-alpha0) * (m1**(1-alpha0) - m0**(1-alpha0))   +   C1/(1-alpha1) * (m2**(1-alpha1) - m1**(1-alpha1))   +   C2/(1-alpha2) * (m3**(1-alpha2) - m2**(1-alpha2))   +   C3/(1-alpha3) * (m_inf**(1-alpha3) - m3**(1-alpha3))    

        if(m_sup >= 0.01 and m_sup < 0.08 ):   nCDF_sup = C0/(1-alpha0) * (m_sup**(1-alpha0) - m0**(1-alpha0))
        if(m_sup >= 0.08 and m_sup < 0.5  ):   nCDF_sup = C0/(1-alpha0) * (m1**(1-alpha0) - m0**(1-alpha0))   +   C1/(1-alpha1) * (m_sup**(1-alpha1) - m1**(1-alpha1))        
        if(m_sup >= 0.5  and m_sup < 1.   ):   nCDF_sup = C0/(1-alpha0) * (m1**(1-alpha0) - m0**(1-alpha0))   +   C1/(1-alpha1) * (m2**(1-alpha1) - m1**(1-alpha1))   +   C2/(1-alpha2) * (m_sup**(1-alpha2) - m2**(1-alpha2))  
        if(m_sup >= 1.   and m_sup <= 150.):   nCDF_sup = C0/(1-alpha0) * (m1**(1-alpha0) - m0**(1-alpha0))   +   C1/(1-alpha1) * (m2**(1-alpha1) - m1**(1-alpha1))   +   C2/(1-alpha2) * (m3**(1-alpha2) - m2**(1-alpha2))   +   C3/(1-alpha3) * (m_sup**(1-alpha3) - m3**(1-alpha3))                                       
        return nCDF_sup-nCDF_inf
                
                
def Analytical_mCDF(m_inf,m_sup):    
        k = 11.045901141 #1./0.09053131901837103   #Mass normalization, so that CDF_imf_function_mass(0.01,150.) = 1    
        C0 = k * (1./0.08)**(-0.3)
        C1 = k * (1./0.08)**(-1.3)
        C2 = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3)
        C3 = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3)
        alpha0 = 0.3
        alpha1 = 1.3
        alpha2 = 2.3
        alpha3 = 2.3
        m0 = 0.01
        m1 = 0.08
        m2 = 0.5
        m3 = 1.

        if(m_inf >= 0.01 and m_inf < 0.08 ):   mCDF_inf = C0/(2-alpha0) * (m_inf**(2-alpha0) - m0**(2-alpha0))
        if(m_inf >= 0.08 and m_inf < 0.5  ):   mCDF_inf = C0/(2-alpha0) * (m1**(2-alpha0) - m0**(2-alpha0))   +   C1/(2-alpha1) * (m_inf**(2-alpha1) - m1**(2-alpha1))        
        if(m_inf >= 0.5  and m_inf < 1.   ):   mCDF_inf = C0/(2-alpha0) * (m1**(2-alpha0) - m0**(2-alpha0))   +   C1/(2-alpha1) * (m2**(2-alpha1) - m1**(2-alpha1))   +   C2/(2-alpha2) * (m_inf**(2-alpha2) - m2**(2-alpha2))  
        if(m_inf >= 1.   and m_inf <= 150.):   mCDF_inf = C0/(2-alpha0) * (m1**(2-alpha0) - m0**(2-alpha0))   +   C1/(2-alpha1) * (m2**(2-alpha1) - m1**(2-alpha1))   +   C2/(2-alpha2) * (m3**(2-alpha2) - m2**(2-alpha2))   +   C3/(2-alpha3) * (m_inf**(2-alpha3) - m3**(2-alpha3))    

        if(m_sup >= 0.01 and m_sup < 0.08 ):   mCDF_sup = C0/(2-alpha0) * (m_sup**(2-alpha0) - m0**(2-alpha0))
        if(m_sup >= 0.08 and m_sup < 0.5  ):   mCDF_sup = C0/(2-alpha0) * (m1**(2-alpha0) - m0**(2-alpha0))   +   C1/(2-alpha1) * (m_sup**(2-alpha1) - m1**(2-alpha1))        
        if(m_sup >= 0.5  and m_sup < 1.   ):   mCDF_sup = C0/(2-alpha0) * (m1**(2-alpha0) - m0**(2-alpha0))   +   C1/(2-alpha1) * (m2**(2-alpha1) - m1**(2-alpha1))   +   C2/(2-alpha2) * (m_sup**(2-alpha2) - m2**(2-alpha2))  
        if(m_sup >= 1.   and m_sup <= 150.):   mCDF_sup = C0/(2-alpha0) * (m1**(2-alpha0) - m0**(2-alpha0))   +   C1/(2-alpha1) * (m2**(2-alpha1) - m1**(2-alpha1))   +   C2/(2-alpha2) * (m3**(2-alpha2) - m2**(2-alpha2))   +   C3/(2-alpha3) * (m_sup**(2-alpha3) - m3**(2-alpha3))      
        return mCDF_sup-mCDF_inf                
                


def Array_Analytical_nCDF_and_inverse(m_low,m_up):
        #This function returns only the array of the masses (array_CDF), and the array of the CDF (array_CDF_inverse)
        normalization_of_Analytical_nCDF = Analytical_nCDF(m_low,m_up)
        array_Analytical_nCDF = np.linspace(-2.,2.176091259,10000) #in log10  #(0.15, m_up., num=1000)  in linear  120Ms = 2.079181246 in logm
        array_Analytical_nCDF = 10**array_Analytical_nCDF
        array_Analytical_nCDF_inverse = np.zeros(10000)
        for ii in range(0, len(array_Analytical_nCDF)):
         array_Analytical_nCDF_inverse[ii] = Analytical_nCDF(m_low,array_Analytical_nCDF[ii]) / normalization_of_Analytical_nCDF
        return array_Analytical_nCDF, array_Analytical_nCDF_inverse    

def Array_Analytical_mCDF_and_inverse(m_low,m_up):
        #This function returns only the array of the masses (array_CDF), and the array of the CDF (array_CDF_inverse)
        normalization_of_Analytical_mCDF = Analytical_mCDF(m_low,m_up)
        array_Analytical_mCDF = np.linspace(-2.,2.176091259,10000) #in log10  #(0.15, m_up., num=1000)  in linear  120Ms = 2.079181246 in logm
        array_Analytical_mCDF = 10**array_Analytical_mCDF
        array_Analytical_mCDF_inverse = np.zeros(10000)
        for ii in range(0, len(array_Analytical_mCDF)):
         array_Analytical_mCDF_inverse[ii] = Analytical_mCDF(m_low,array_Analytical_mCDF[ii]) / normalization_of_Analytical_mCDF
        return array_Analytical_mCDF, array_Analytical_mCDF_inverse   
    

'''
#UNTIL 7 APRIL 2016
def Star_Cluster_stellar_masses_generator_N_star(m_low,m_up,m_low_used_array,m_up_used_array,array_CDF, array_CDF_inverse):
        #Now we can use it: Generate a random number between 0 and 1 and then interpolate in the array_CDF_inverse to get the mass. 
        N = len(m_low_used_array)
        cdf_low = np.zeros(N)
        cdf_up = np.zeros(N)
        cdf_low[:] = np.interp(m_low_used_array[:], array_CDF, array_CDF_inverse)
        cdf_up[:]  = np.interp(m_up_used_array[:] , array_CDF, array_CDF_inverse)
        #normalization_of_Analytical_nCDF = Analytical_nCDF(m_low,m_up)
        #for ii in range(0,len(cdf_low)):
        # cdf_low[ii] = Analytical_nCDF(m_low,m_low_used_array[ii])
        # cdf_up[ii]  = Analytical_nCDF(m_low,m_up_used_array[ii])
        #cdf_low = cdf_low/normalization_of_Analytical_nCDF
        #cdf_up  = cdf_up/normalization_of_Analytical_nCDF
        stars_random_number_generated = np.random.rand(N) #*cdf_range
        stars_random_number_generated[:] =  stars_random_number_generated[:]*(cdf_up[:]-cdf_low[:])
        stars_random_number_generated = stars_random_number_generated + cdf_low
        #stars_random_number_generated_sorted = np.sort(stars_random_number_generated, axis=-1, kind='mergesort')
        stars_random_number_generated_sorted = np.sort(stars_random_number_generated, axis=-1, kind='quicksort')
        stars_mass_generated_sorted = np.interp(stars_random_number_generated_sorted, array_CDF_inverse, array_CDF)
        sort_choice = 0 #not sorted masses!
        if sort_choice == 1:
         return stars_mass_generated_sorted #If we want the stars sorted by mass, from lowest to highest
        if sort_choice == 0:
         stars_mass_generated_unsorted = np.zeros(N)
         stars_mass_generated_unsorted[np.argsort(stars_random_number_generated)] = stars_mass_generated_sorted
         return stars_mass_generated_unsorted #If we want the stars UNsorted by mass!
'''

#7 APRIL 2016
def Star_Cluster_stellar_masses_generator_N_star(m_low,m_up,m_low_used_array,m_up_used_array,array_CDF, array_CDF_inverse):
        #Now we can use it: Generate a random number between 0 and 1 and then interpolate in the array_CDF_inverse to get the mass. 
        N = len(m_low_used_array)
        cdf_low = np.zeros(N)
        cdf_up = np.zeros(N)
        #cdf_low[:] = np.interp(m_low_used_array[:], array_CDF, array_CDF_inverse)
        #cdf_up[:]  = np.interp(m_up_used_array[:] , array_CDF, array_CDF_inverse)
        cdf_low[:] = np.interp(m_low_used_array[0], array_CDF, array_CDF_inverse) #FULL RANDOM SAMPLING 
        cdf_up[:]  = np.interp(m_up_used_array[0] , array_CDF, array_CDF_inverse) #FULL RANDOM SAMPLING 
        #cdf_low[:] = 7.90935267438e-09  #FULL RANDOM SAMPLING
        #cdf_up[:]  = 1.                 #FULL RANDOM SAMPLING 

        #stars_random_number_generated        = np.random.rand(N) #*cdf_range
        #stars_random_number_generated[:]     = stars_random_number_generated[:]*(cdf_up[:]-cdf_low[:])
        #stars_random_number_generated        = stars_random_number_generated + cdf_low

        stars_random_number_generated        = ( np.random.rand(N) * (cdf_up[0]-cdf_low[0]) ) + cdf_low[0]

        #stars_random_number_generated_sorted = np.sort(stars_random_number_generated, axis=-1, kind='mergesort')
        stars_random_number_generated_sorted = np.sort(stars_random_number_generated, axis=-1, kind='quicksort')
        stars_mass_generated_sorted          = np.interp(stars_random_number_generated_sorted, array_CDF_inverse, array_CDF)

        return stars_mass_generated_sorted

def extrap(x, xp, yp):
        """np.interp function with linear extrapolation"""
        y = np.interp(x, xp, yp)
        y = np.where(x<xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
        y = np.where(x>xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]), y)
        return y


def Intervals_of_resampling(star_optimal_index,star_optimal_mass,array_massCDF, array_massCDF_inverse, Expansion_factor,m_low,m_up):
        #I will need to use the limits given by the optimal sampling
        #So i need to take the positions of the stars in OS (their masses in fact)
        #Then i compute a vector from which the elements are the middle masses between stars. 
        #The first element is m_low, and the last is has to be derived!
        m_low_used_array          = np.zeros(len(star_optimal_index))
        m_up_used_array           = np.zeros(len(star_optimal_index))
        m_low_used_array_expanded = np.zeros(len(star_optimal_index))
        m_up_used_array_expanded  = np.zeros(len(star_optimal_index))
        interval_size             = np.zeros(len(star_optimal_index))
        ratio_high_on_low_parts   = np.zeros(len(star_optimal_index))
        rule_50pc_left                  = np.zeros(len(star_optimal_index))
        rule_50pc_right                  = np.zeros(len(star_optimal_index))
        CDF_star_optimal_mass          = np.zeros(len(star_optimal_index))
        CDF_star_optimal_mass_low = np.zeros(len(star_optimal_index)) 
        CDF_star_optimal_mass_up  = np.zeros(len(star_optimal_index))
        CDF_average_mass_bin          = np.zeros(len(star_optimal_index))
        average_mass_bin          = np.zeros(len(star_optimal_index))
 
         
        #Creation of this vector of half-masses elements, bins limits 
        #Method "Limits of intervals given by optimal sampling"  
        for ii in range(0,len(star_optimal_index)):
         if ii==0:                                 #First bin 
          m_low_used = 0.5*(star_optimal_mass[ii]+star_optimal_mass[ii+1])                        #In linear!!!
          m_up_used  = 150. #???? Could not be true!!! Modified in the following :)
         if ii>0 and ii<len(star_optimal_index)-1:                                #Common bins
          m_low_used = 0.5*(star_optimal_mass[ii]+star_optimal_mass[ii+1])                        #In linear!!!
          m_up_used  = 0.5*(star_optimal_mass[ii-1]+star_optimal_mass[ii])                        #In linear!!!
         if ii==len(star_optimal_index)-1:        #Last bin 
          m_low_used = m_low
          m_up_used  = 0.5*(star_optimal_mass[ii-1]+star_optimal_mass[ii])                        #In linear!!!
         m_low_used_array[ii] = max(m_low_used,m_low)
         m_up_used_array[ii] = min(m_up_used,m_up)


        #Realizing the extrapolation for the most massive star superior limit
        index_regime_star_1 = 1 #index of the second most massive star of the last IMF regime
        #index_regime_star_2 = int(len(star_optimal_index)*0.03) #index of a star close to the inferior limit of the last IMF regime
         
        #index of a star close to the inferior limit of the relevant IMF regime. 
        #This is to garantee that we take the good regime for the most massive star interval!
        if star_optimal_mass[0]>1:
         itemindex=np.where(star_optimal_mass>1)
         index_regime_star_2 = itemindex[0][-1]
        if star_optimal_mass[0]<=1 and star_optimal_mass[0]>0.5:
         itemindex=np.where(star_optimal_mass>0.5)
         index_regime_star_2 = itemindex[0][-1]
        if star_optimal_mass[0]<=0.5 and star_optimal_mass[0]>0.08:
         itemindex=np.where(star_optimal_mass>0.08)
         index_regime_star_2 = itemindex[0][-1]
        if len(itemindex[0])<=2: index_regime_star_2=index_regime_star_2+2

        mass_1 = star_optimal_mass[index_regime_star_1] 
        mass_2 = star_optimal_mass[index_regime_star_2]
        mass_regime = star_optimal_mass[index_regime_star_1:index_regime_star_2+1]
        mass_regime = np.log10(mass_regime)
        mass_regime = mass_regime[::-1] #Reverse the order of the array, because need to be croissant numbers, for interpolation
        m_up_used_array_regime = m_up_used_array[index_regime_star_1:index_regime_star_2+1] - star_optimal_mass[index_regime_star_1:index_regime_star_2+1]
        m_up_used_array_regime = np.log10(m_up_used_array_regime)
        m_up_used_array_regime = m_up_used_array_regime[::-1]
        mass_extrapolated = extrap(np.log10(star_optimal_mass[0]), mass_regime, m_up_used_array_regime) 
        mass_extrapolated = (10**mass_extrapolated) + star_optimal_mass[0]
        m_up_used_array[0] = mass_extrapolated
        interval_size[:] = m_up_used_array[:] - m_low_used_array[:]
        #ratio_high_on_low_parts[:] = (m_up_used_array[:]-star_optimal_mass[:]) / (star_optimal_mass[:]-m_low_used_array[:])


        normalization_of_number_CDF = Analytical_nCDF(m_low,m_up)
        normalization_of_mass_CDF   = Analytical_mCDF(m_low,m_up)
        #Expansion_factor = 6
        for ii in range(0,len(star_optimal_index)):
         #Linear growth of the intervals size
         m_low_used_array_expanded[ii] = max(star_optimal_mass[ii] - Expansion_factor*(star_optimal_mass[ii]-m_low_used_array[ii]), m_low)
         m_up_used_array_expanded[ii]  = min(star_optimal_mass[ii] + Expansion_factor*(m_up_used_array[ii]-star_optimal_mass[ii]),  m_up)
        m_low_used_array = m_low_used_array_expanded
        m_up_used_array = m_up_used_array_expanded

        #'''
        m_low_used_array = m_low_used_array + 0.000000001 #Security! If 0.01, it crashes.
        saturated = np.zeros(len(star_optimal_mass))
        criteria = np.zeros(len(star_optimal_mass))
        delta_m = np.zeros(len(star_optimal_mass)) 
        bins_low_mask  = ne.evaluate('star_optimal_mass <  0.3837739123138731') 
        bins_high_mask = ne.evaluate('star_optimal_mass >= 0.3837739123138731') 
        star_effective_mass     = np.zeros(len(star_optimal_mass))
        saturated[:] = 0 #By default, all stars are unsaturated

        #Classifying saturated and unsaturated bins. Computation of missing mass
        for ii in range(0,len(star_optimal_mass)):
           
         if star_optimal_mass[ii] >= 0.3837739123138731: #Case where bins are above the mirror
          criteria[ii] = Analytical_mCDF(m_low_used_array[ii], 150.)/Analytical_nCDF(m_low_used_array[ii], 150.)
          if criteria[ii] < star_optimal_mass[ii]:  #Case saturated
           saturated[ii] = 1
           delta_m[ii] = star_optimal_mass[ii] - criteria[ii]  #Mass lost in that bin. Shown as positive because it will be added later
           star_effective_mass[ii] = criteria[ii]
           m_up_used_array[ii] = 150.
          if criteria[ii] >= star_optimal_mass[ii]: #Case valid
           saturated[ii] = 0
           m_up_used_array[ii] = scipy.optimize.brentq(mCDF_on_nCDF_up, 0.01, 150., args=(star_optimal_mass[ii], m_low_used_array[ii]), xtol=1e-12,disp=True)
           star_effective_mass[ii] = star_optimal_mass[ii]

         if star_optimal_mass[ii] < 0.3837739123138731: #Case where bins are below the mirror 
          criteria[ii] = Analytical_mCDF(0.01,m_up_used_array[ii])/Analytical_nCDF(0.01,m_up_used_array[ii])
          if criteria[ii] > star_optimal_mass[ii]:  #Case saturated
           saturated[ii] = -1
           delta_m[ii] = star_optimal_mass[ii] - criteria[ii]        #Mass gained in that bin. Shown as negative because it will be removed later
           star_effective_mass[ii] = criteria[ii]
           m_low_used_array[ii] = 0.010000001
          if criteria[ii] <= star_optimal_mass[ii]: #Case valid
           saturated[ii] = 0
           m_low_used_array[ii] = scipy.optimize.brentq(mCDF_on_nCDF_low, 0.01, 150., args=(star_optimal_mass[ii], m_up_used_array[ii]), xtol=1e-12,disp=True)
           star_effective_mass[ii] = star_optimal_mass[ii]

        #print 'saturated sum, (security) ', saturated.sum()
        mask_delta = ne.evaluate('delta_m<0')
        #print 'Delta negative sum ', delta_m[mask_delta].sum()
        #print 'Delta positive sum ', delta_m[~mask_delta].sum()
        Delta_all = delta_m.sum()
        #print 'Delta_all = ', Delta_all
        #print 'star_optimal_mass.sum()   ', star_optimal_mass.sum()  #Sum of the effective mass contained in each bin
        #print 'star_effective_mass.sum() ', star_effective_mass.sum()  #Sum of the effective mass contained in each bin
        #print
         
        
        #mask_saturated = ne.evaluate('saturated != 0')
        mask_saturated_right = ne.evaluate('saturated > 0')
        if len(star_optimal_mass[mask_saturated_right]) > 0:
         mean_saturated_right = np.mean(star_optimal_mass[mask_saturated_right])
         m_up_used_array[mask_saturated_right] = 150.
         #print 'len(star_optimal_mass[mask_saturated_right]) ', len(star_optimal_mass[mask_saturated_right])
         #print 'star_optimal_mass[mask_saturated_right] ', star_optimal_mass[mask_saturated_right]
         #print 'mean_saturated_right ', mean_saturated_right
         #print
         m_low_used_array[mask_saturated_right] = scipy.optimize.brentq(mCDF_on_nCDF_low, 0.01, 150., args=(mean_saturated_right, 149.9999), xtol=1e-12,disp=True)

        mask_saturated_left = ne.evaluate('saturated < 0')
        if len(star_optimal_mass[mask_saturated_left]) > 0:
         mean_saturated_left = np.mean(star_optimal_mass[mask_saturated_left])
         m_low_used_array[mask_saturated_left] = 0.010000001

         #print 'len(star_optimal_mass[mask_saturated_left]) ', len(star_optimal_mass[mask_saturated_left])
         #print 'star_optimal_mass[mask_saturated_left] ', star_optimal_mass[mask_saturated_left]
         #print 'mean_saturated_left ', mean_saturated_left
         #print 
         m_up_used_array[mask_saturated_left] = scipy.optimize.brentq(mCDF_on_nCDF_up, 0.01, 150., args=(mean_saturated_left, 0.0100000001), xtol=1e-12,disp=True)

        #mask = ne.evaluate('m_up_used_array>=149.9')
        #fig = plt.figure(1)
        #plt.plot(star_optimal_mass,-delta_m,'bo-')
        #plt.plot(star_optimal_mass[mask],-delta_m[mask],'ro')
        #plt.axhline(y=0.,color='k',ls='dashed')
        #plt.title('Mean mass per bin - OS mass per bin ; EF = 5000 ; {0} clusters'.format(choice_number_SC))
        #ylabel('Mean mass per bin - OS mass per bin')
        #xlabel('OS mass')
        #plt.show()
        #'''
        return m_low_used_array, m_up_used_array, average_mass_bin





def Read_isochrones(choice_metallicity):
        logAge,M_ini,M_act,logL,logTe,logG,mbol=[],[],[],[],[],[],[]
        FUV,NUV,U,B,V,R,I,J,H,K=[],[],[],[],[],[],[],[],[],[]
        u_CFHT,g_CFHT,r_CFHT,i_CFHT,z_CFHT,IRAC36,IRAC45,IRAC58,IRAC80=[],[],[],[],[],[],[],[],[]
        u_SDSS,g_SDSS,r_SDSS,i_SDSS,z_SDSS,J_2MASS,H_2MASS,K_2MASS=[],[],[],[],[],[],[],[]
        mip24,mip70,mip160,CO,M_hec,period,pmode,logMdot,int_IMF=[],[],[],[],[],[],[],[],[]
        a_BATC,b_BATC,c_BATC,d_BATC,e_BATC,f_BATC,g_BATC,h_BATC,i_BATC,j_BATC,k_BATC,m_BATC,n_BATC,o_BATC,p_BATC,t_BATC=[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]
        F275W,F336W,F475W,F814W,F110W,F160W=[],[],[],[],[],[]
        #isochrone_name = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Girardi2010_isochrones_preparation/{0}/{0}_AllBands_HST_CB_step_logt01_Girardi2010_CMD25L.dat'.format(choice_metallicity) 
        #isochrone_name = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Girardi2010_isochrones_preparation/{0}/{0}_AllBands_HST_CB_step_logt01_Girardi2010_CMD25L_withACS.dat'.format(choice_metallicity) 
        isochrone_name = '/media/philippe/a36c9bac-04fd-4046-8af2-2038962a2127/philippe/Documents/PhD/Discrete_models_comparaison_jtao/SC_Parameters_20/Girardi2010_isochrones_preparation/{0}/{0}_AllBands_HST_CB_step_logt01_Girardi2010_CMD25L_withACS.dat'.format(choice_metallicity)
        infile_isochrone = open(isochrone_name,'r')
        for line in infile_isochrone:
          ll = line.split()
          if ll[0]!='#' and ll[0]!='#log(age/yr)':
            logAge.append(float(ll[0]))
            M_ini.append(float(ll[1]))
            M_act.append(float(ll[2]))
            logL.append(float(ll[3]))
            logTe.append(float(ll[4]))
            logG.append(float(ll[5]))
            mbol.append(float(ll[6]))
            FUV.append(float(ll[7]))
            NUV.append(float(ll[8]))
            U.append(float(ll[9]))
            B.append(float(ll[10]))
            V.append(float(ll[11]))
            R.append(float(ll[12]))
            I.append(float(ll[13]))
            J.append(float(ll[14]))
            H.append(float(ll[15]))
            K.append(float(ll[16]))
            F275W.append(float(ll[17]))
            F336W.append(float(ll[18]))
            F475W.append(float(ll[19]))
            F814W.append(float(ll[20]))
            F110W.append(float(ll[21]))
            F160W.append(float(ll[22]))
            u_CFHT.append(float(ll[23]))
            g_CFHT.append(float(ll[24]))
            r_CFHT.append(float(ll[25]))
            i_CFHT.append(float(ll[26]))
            z_CFHT.append(float(ll[27]))
            u_SDSS.append(float(ll[28]))
            g_SDSS.append(float(ll[29]))
            r_SDSS.append(float(ll[30]))
            i_SDSS.append(float(ll[31]))
            z_SDSS.append(float(ll[32]))
            J_2MASS.append(float(ll[33]))
            H_2MASS.append(float(ll[34]))
            K_2MASS.append(float(ll[35]))
            a_BATC.append(float(ll[36]))
            b_BATC.append(float(ll[37]))
            c_BATC.append(float(ll[38]))
            d_BATC.append(float(ll[39]))
            e_BATC.append(float(ll[40]))
            f_BATC.append(float(ll[41]))
            g_BATC.append(float(ll[42]))
            h_BATC.append(float(ll[43]))
            i_BATC.append(float(ll[44]))
            j_BATC.append(float(ll[45]))
            k_BATC.append(float(ll[46]))
            m_BATC.append(float(ll[47]))
            n_BATC.append(float(ll[48]))
            o_BATC.append(float(ll[49]))
            p_BATC.append(float(ll[50]))
            t_BATC.append(float(ll[51]))
            IRAC36.append(float(ll[52]))
            IRAC45.append(float(ll[53]))
            IRAC58.append(float(ll[54]))
            IRAC80.append(float(ll[55]))
            mip24.append(float(ll[56]))
            mip70.append(float(ll[57]))
            mip160.append(float(ll[58]))
        infile_isochrone.close()
        logAge = np.array(logAge)
        M_ini = np.array(M_ini)
        FUV = np.array(FUV)
        NUV = np.array(NUV)
        U = np.array(U)
        B = np.array(B)
        V = np.array(V)
        R = np.array(R)
        I = np.array(I)
        J = np.array(J)
        H = np.array(H)
        K = np.array(K)
        F275W = np.array(F275W)
        F336W = np.array(F336W)
        F475W = np.array(F475W)
        F814W = np.array(F814W)
        F110W = np.array(F110W)
        F160W = np.array(F160W)
        u_CFHT = np.array(u_CFHT)
        g_CFHT = np.array(g_CFHT)
        r_CFHT = np.array(r_CFHT)
        i_CFHT = np.array(i_CFHT)
        z_CFHT = np.array(z_CFHT)
        u_SDSS = np.array(u_SDSS)
        g_SDSS = np.array(g_SDSS)
        r_SDSS = np.array(r_SDSS)
        i_SDSS = np.array(i_SDSS)
        z_SDSS = np.array(z_SDSS)
        J_2MASS = np.array(J_2MASS)
        H_2MASS = np.array(H_2MASS)
        K_2MASS = np.array(K_2MASS)
        a_BATC = np.array(a_BATC)
        b_BATC = np.array(b_BATC)
        c_BATC = np.array(c_BATC)
        d_BATC = np.array(d_BATC)
        e_BATC = np.array(e_BATC)
        f_BATC = np.array(f_BATC)
        g_BATC = np.array(g_BATC)
        h_BATC = np.array(h_BATC)
        i_BATC = np.array(i_BATC)
        j_BATC = np.array(j_BATC)
        k_BATC = np.array(k_BATC)
        m_BATC = np.array(m_BATC)
        n_BATC = np.array(n_BATC)
        o_BATC = np.array(o_BATC)
        p_BATC = np.array(p_BATC)
        t_BATC = np.array(t_BATC)
        IRAC36 = np.array(IRAC36)
        IRAC45 = np.array(IRAC45)
        IRAC58 = np.array(IRAC58)
        IRAC80 = np.array(IRAC80)
        mip24 = np.array(mip24)
        mip70 = np.array(mip70)
        mip160 = np.array(mip160)

        photometry = np.zeros((len(U),52))
        photometry[:,0] = FUV
        photometry[:,1] = NUV
        photometry[:,2] = U
        photometry[:,3] = B
        photometry[:,4] = V
        photometry[:,5] = R
        photometry[:,6] = I
        photometry[:,7] = J
        photometry[:,8] = H
        photometry[:,9] = K
        photometry[:,10] = F275W
        photometry[:,11] = F336W
        photometry[:,12] = F475W
        photometry[:,13] = F814W
        photometry[:,14] = F110W
        photometry[:,15] = F160W
        photometry[:,16] = u_CFHT
        photometry[:,17] = g_CFHT
        photometry[:,18] = r_CFHT
        photometry[:,19] = i_CFHT
        photometry[:,20] = z_CFHT
        photometry[:,21] = u_SDSS
        photometry[:,22] = g_SDSS
        photometry[:,23] = r_SDSS
        photometry[:,24] = i_SDSS
        photometry[:,25] = z_SDSS
        photometry[:,26] = J_2MASS
        photometry[:,27] = H_2MASS
        photometry[:,28] = K_2MASS
        photometry[:,29] = a_BATC
        photometry[:,30] = b_BATC
        photometry[:,31] = c_BATC
        photometry[:,32] = d_BATC
        photometry[:,33] = e_BATC
        photometry[:,34] = f_BATC
        photometry[:,35] = g_BATC
        photometry[:,36] = h_BATC
        photometry[:,37] = i_BATC
        photometry[:,38] = j_BATC
        photometry[:,39] = k_BATC
        photometry[:,40] = m_BATC
        photometry[:,41] = n_BATC
        photometry[:,42] = o_BATC
        photometry[:,43] = p_BATC
        photometry[:,44] = t_BATC
        photometry[:,45] = IRAC36
        photometry[:,46] = IRAC45
        photometry[:,47] = IRAC58
        photometry[:,48] = IRAC80
        photometry[:,49] = mip24
        photometry[:,50] = mip70
        photometry[:,51] = mip160
        return logAge,M_ini,photometry



def Read_isochrones_April2016(choice_metallicity):
        #isochrone_name = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Girardi2010_isochrones_preparation/{0}/{0}_AllBands_HST_CB_step_logt01_Girardi2010_CMD25L.dat'.format(choice_metallicity) 
        #isochrone_name = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Girardi2010_isochrones_preparation/{0}/{0}_AllBands_HST_CB_step_logt01_Girardi2010_CMD25L_withACS.dat'.format(choice_metallicity) 
        isochrone_name = '/media/philippe/a36c9bac-04fd-4046-8af2-2038962a2127/philippe/Documents/PhD/Discrete_models_comparaison_jtao/SC_Parameters_20/Girardi2010_isochrones_preparation/{0}/{0}_AllBands_HST_CB_step_logt01_Girardi2010_CMD25L_withACS.dat'.format(choice_metallicity) 

        photometry = np.zeros((len(U),52))
        iso_data = np.genfromtxt(isochrone_name,comments='#')
        #photometry = np.zeros((size_all_photometry,52))
        #logAge     = np.zeros(size_all_photometry)
        #M_ini      = np.zeros(size_all_photometry)
        logAge       = iso_data[:,0]
        M_ini        = iso_data[:,1]
        photometry   = iso_data[:,7:59] #float(ll[7:59])

        return logAge,M_ini,photometry





def FilterName_to_FilterIndex(filter_name):
    '''
    This function converts the name of a filter to its index number    
    USE: 
    filter_name = 'V'
    filter_name_index = FilterName_to_FilterIndex(filter_name)
    print filter_name_index
    >> 5
    '''
    
    filters_names_list = ['FUV','NUV',
                          'U','B','V','R','I','J','H','K',
                          'F275W','F336W','F475W','F814W','F110W','F160W',
                          'u_CFHT','g_CFHT','r_CFHT','i_CFHT','z_CFHT',
                          'u_SDSS','g_SDSS','r_SDSS','i_SDSS','z_SDSS',
                          'J_2MASS','H_2MASS','Ks_2MASS',
                          'a_BATC','b_BATC','c_BATC','d_BATC','e_BATC','f_BATC','g_BATC','h_BATC','i_BATC','j_BATC','k_BATC','m_BATC','n_BATC','o_BATC','p_BATC','t_BATC',
                          'IRAC36','IRAC45','IRAC58','IRAC80',
                          'mip24','mip70','mip160']
    
    ii=1
    for name_list in filters_names_list:
        if filter_name == name_list:
            filter_name_index = ii
        ii=ii+1
    
    return filter_name_index

def FilterIndex_to_FilterName(filter_index):    
    '''
    This function converts the index number of a filter to its name
    USE: 
    filter_index = 5
    filter_name = FilterIndex_to_FilterName(filter_index)
    print filter_name
    >>> V
    '''
    
    filters_names_list = ['FUV','NUV',
                          'U','B','V','R','I','J','H','K',
                          'F275W','F336W','F475W','F814W','F110W','F160W',
                          'u_CFHT','g_CFHT','r_CFHT','i_CFHT','z_CFHT',
                          'u_SDSS','g_SDSS','r_SDSS','i_SDSS','z_SDSS',
                          'J_2MASS','H_2MASS','Ks_2MASS',
                          'a_BATC','b_BATC','c_BATC','d_BATC','e_BATC','f_BATC','g_BATC','h_BATC','i_BATC','j_BATC','k_BATC','m_BATC','n_BATC','o_BATC','p_BATC','t_BATC',
                          'IRAC36','IRAC45','IRAC58','IRAC80',
                          'mip24','mip70','mip160']
    
    filter_name = filters_names_list[filter_index-1]
    return filter_name




def mCDF_on_nCDF_up(x, OS, a):
        return OS - Analytical_mCDF(a,x)/Analytical_nCDF(a,x)

def mCDF_on_nCDF_low(x, OS, b):
        return OS - Analytical_mCDF(x,b)/Analytical_nCDF(x,b)

#=============================================================================================================


def straight_line_array_from_two_points(array_x,x1,x2,y1,y2):
    #example of use: (creates an array_x and an array_y)
    #array_x = np.linspace(6.6,10.1,71)
    #array_y = straight_line_array_from_two_points(array_x,x1,x2,y1,y2)
    c = (y2-y1)/(x2-x1)       #ordonnee a l'origine and slope of first line    
    y0  = y2-c*x2   #ordonnee a l'origine of first line (derived using y2-y0 = c*(x2-x0) where x0=0 by definition of abscissa at origine)
    return y0 + array_x*c 


#===============================================================================================================
#This function is designed to build three arrays: one of the mass, one of the NUMBER IMF and one of the MASS IMF
#Just for purpose of plotting.
#It requires the imf_function (which is defined from a Kroupa multi power-law)
def IMF_function_for_plot(Mass_min,Mass_max):
    #Here I built the IMF function because we will plot it
    IMF_abs        = np.zeros(100)
    numberIMF_oord = np.zeros(100)
    massIMF_oord   = np.zeros(100)
    #Mass_min = 0.01 #Ms
    #Mass_max = 150. #Ms    
    step = (np.log10(Mass_max) - np.log10(Mass_min)) / 100.
    for ii in range(0,len(IMF_abs)):
     m = 10**(np.log10(Mass_min)+ii*step)
     IMF_abs[ii]        = m
     numberIMF_oord[ii] = imf_function(m)
     massIMF_oord[ii]   = imf_function_m(m) #which is exactly same as imf_function(m)*m
    return IMF_abs, numberIMF_oord, massIMF_oord

#USAGE:
# IMF_abs, numberIMF_oord, massIMF_oord = IMF_function(0.01,150.)
# 0.01,150. are the minimum and maximum mass of the IMF, in linear mass
#===============================================================================================================



#===============================================================================================================
def Saving_MASS_of_All_Stars_AliveDead_of_EachCluster(logMass_desired_jj,jj,star_half_stochastic_mass,choice_mass_distribution,mass_indice,Z_indice):
 '''Saving masses of cluster ALL stars (also dead and very low masses) in separated files'''
 
 # logMass_desired_jj : logMass_desired[jj], the mass of the cluster, in log10
 # jj : the number of the cluster

 #if (logMass_desired[jj] <= 3. and jj < 1000) or (logMass_desired[jj] <= 4. and jj < 100) or (logMass_desired[jj] <= 5. and jj < 10) : 
 if (logMass_desired_jj <= 3. and jj < 1000) or (logMass_desired_jj <= 4. and jj < 100) or (logMass_desired_jj <= 5. and jj < 10) : 
  #only for a reduced numbar of clusters, because a lot of data! (number depends on mass of clusters)
  all_stars_created = np.sort(star_half_stochastic_mass, axis=-1, kind='mergesort')
  all_stars_created = all_stars_created[::-1]
  
  if choice_mass_distribution == 1:
   file_output_CMD_all_stars = open('/home/philippe/Desktop/results/CMD_all_stars_M{1}_Z{2}/cluster_{0}_CMD_all_stars'.format(jj,mass_indice,Z_indice) , 'w' )
  elif choice_mass_distribution > 1:
   file_output_CMD_all_stars = open('/home/philippe/Desktop/results/CMD_all_stars/cluster_{0}_CMD_all_stars'.format(jj) , 'w' )

  print >> file_output_CMD_all_stars, '# Mass'
  np.savetxt(file_output_CMD_all_stars,all_stars_created)
  file_output_CMD_all_stars.close()
  del all_stars_created
#===============================================================================================================



#===============================================================================================================
def Saving_MASS_and_PHOTOMETRY_of_All_Stars_Alive_of_EachCluster(logMass_desired_jj,jj,star_half_stochastic_mass,choice_mass_distribution,age,mass_indice,Z_indice,photometry_flux,index_useful_mags,number_of_filters):
 '''Saving photometry of cluster alive stars in separated files. From this we can build resolved CMDs of each cluster model''' 

 # logMass_desired_jj : logMass_desired[jj], the mass of the cluster, in log10
 # jj : the number of the cluster

 #if (logMass_desired[jj] <= 3. and jj < 1000) or (logMass_desired[jj] <= 4. and jj < 100) or (logMass_desired[jj] <= 5. and jj < 10) : 
 if (logMass_desired_jj <= 3. and jj < 1000) or (logMass_desired_jj <= 4. and jj < 100) or (logMass_desired_jj <= 5. and jj < 10) :
  #only for a reduced number of clusters, because a lot of data! (number depends on mass of clusters)
  # It would be good to sort the stars from the most massive to the less massive...
  data_print_CMD = np.zeros((len(star_half_stochastic_mass[index_useful_mags]),number_of_filters+1+1)) #lines: number of stars ; columns: id, age, mass, photometry
  data_print_CMD[:,0]  = age
  data_print_CMD[:,1]  = star_half_stochastic_mass[index_useful_mags]
  data_print_CMD[:,2:] = -2.5*np.log10(photometry_flux[:,:])

  if choice_mass_distribution == 1:
   file_output_CMD = open('/home/philippe/Desktop/results/CMD_M{1}_Z{2}/cluster_{0}_CMD'.format(jj,mass_indice,Z_indice) , 'w' )
  elif choice_mass_distribution > 1:
   file_output_CMD = open('/home/philippe/Desktop/results/CMD/cluster_{0}_CMD'.format(jj) , 'w' )

  print >> file_output_CMD, '# Age     Mass    FUV     NUV     U       B       V       R       I       J       H       K       F275W   F336W   F475W   F814W   F110W   F160W   u_CFHT  g_CFHT  r_CFHT  i_CFHT  z_CFHT  u_SDSS  g_SDSS  r_SDSS  i_SDSS  z_SDSS  J_2Mass H_2Mass K_2Mass a_BATC  b_BATC  c_BATC  d_BATC  e_BATC  f_BATC  g_BATC  h_BATC  i_BATC  j_BATC  k_BATC  m_BATC  n_BATC  o_BATC  p_BATC  t_BATC  IRAC36  IRAC45  IRAC58  IRAC80   mip24   mip70  mip160'
  np.savetxt(file_output_CMD,data_print_CMD, fmt='%.5f')
  file_output_CMD.close()
#===============================================================================================================








def Reading_of_InputFile(InputFile_Name):
        '''This function is the reading function of input files for FameClust program'''
        InputFile = open('/home/philippe/Desktop/Discrete_models_comparaison_jtao/SC_Parameters_20/'+InputFile_Name).readlines()
        my_list = []
        for line in InputFile:
            item = str.split(line)
            print( item )
            if item[0][0] != '#':
             #my_list.append(item[0])
             for item_of_item in item:
              my_list.append(item_of_item)

        number_filters_GLOBAL = int(my_list[0])
        print( 'Number of filters selected:                 ', number_filters_GLOBAL )

        filters_selected_index = np.arange(number_filters_GLOBAL)         #integer array containing the indexes of the selected filters 
        for ff in range(0,number_filters_GLOBAL):
         filters_selected_index[ff] = int(my_list[ff+1])
        print( 'Indexes of filters selected:               ', filters_selected_index )

        Distance_modulus_host_galaxy = my_list[number_filters_GLOBAL+1]
        print( 'Distance modulus of the host galaxy:        ', Distance_modulus_host_galaxy )
        app_or_abs = int(my_list[number_filters_GLOBAL+2])
        print( 'Apparent mags [1], Absolute mags [2]:       ', app_or_abs )
        file_observed_clusters = my_list[number_filters_GLOBAL+3]
        print( 'Input file of the observed clusters:        ' )
        print( '    ',file_observed_clusters )
        number_cluster_observed = int(my_list[number_filters_GLOBAL+4])
        print( 'Number of observed clusters:                ',  number_cluster_observed )        #[obsolete?] used with MagLim!
        choice_extinction = int(my_list[number_filters_GLOBAL+5])
        print( 'Cluster(s) extincted [1], not extincted [2]:', choice_extinction )
        choice_extinction_law = int(my_list[number_filters_GLOBAL+6])
        print( 'Extinction law of MW [1], of LMC [2]:       ', choice_extinction_law )
        path_file_out_cluster = my_list[number_filters_GLOBAL+7]
        print( 'Path of output files for derived parameters:' )
        print( '    ',path_file_out_cluster )
        Grid_path1 = my_list[number_filters_GLOBAL+8]
        Grid_path2 = my_list[number_filters_GLOBAL+9]       
        print( 'Path and name of the grid of models:' )
        print( Grid_path1 )
        print( Grid_path2 )
        number_nodes_age,number_nodes_mass = int(my_list[number_filters_GLOBAL+10]),int(my_list[number_filters_GLOBAL+11])
        print( 'Number of age nodes and of mass nodes:', number_nodes_age,number_nodes_mass )
        number_models_per_node = int(my_list[number_filters_GLOBAL+12])
        print( 'Number of models per node:', number_models_per_node )
        number_of_filters_in_grid = int(my_list[number_filters_GLOBAL+13])
        print( 'Number of filters available in the grid of models:', number_of_filters_in_grid )
        File_information_filters_and_ExtCurve = my_list[number_filters_GLOBAL+14]
        print( 'Name of the file containing infos of filters and extinction curve:', File_information_filters_and_ExtCurve )
        Min_extinction, Max_extinction = float(my_list[number_filters_GLOBAL+15]),float(my_list[number_filters_GLOBAL+16])
        print( 'Minimum and maximum extinction:',Min_extinction, Max_extinction )
        choice_RealClusters_or_ArtificialTest = int(my_list[number_filters_GLOBAL+17])
        print( 'Artificial test or real clusters studied (1/2):', choice_RealClusters_or_ArtificialTest )

        return number_filters_GLOBAL,        \
               filters_selected_index,       \
               Distance_modulus_host_galaxy, \
               app_or_abs,                   \
               file_observed_clusters,       \
               number_cluster_observed,      \
               choice_extinction,            \
               choice_extinction_law,        \
               path_file_out_cluster,        \
               Grid_path1,                   \
               Grid_path2,                   \
               number_nodes_age,             \
               number_nodes_mass,            \
               number_models_per_node,       \
               number_of_filters_in_grid,    \
               File_information_filters_and_ExtCurve, \
               Min_extinction,               \
               Max_extinction,               \
               choice_RealClusters_or_ArtificialTest




#======================================================================================================================
#  CARDELLI, CLAYTON & MATHIS (1989) EXTINCTION LAW
#======================================================================================================================
#Here i will create a Cardelli (1989) law function between 0.1 and 3.33 microns

def a(x):
 #x is the wavelength inverse, in (micron)**-1     
 if x >=0.3 and x<=1.1: #optical range
  return 0.574*x**1.61     
 if x >=1.1 and x<=3.3: #optical/NIR range 
  y=x-1.82   
  return 1. + 0.17699*y - 0.50447*y**2 - 0.02427*y**3 + 0.72085*y**4 + 0.01979*y**5 - 0.77530*y**6 + 0.32999*y**7       
 if x >=3.3 and x<5.9:   #UV 
  return 1.752 - 0.316*x - 0.104/( (x-4.67)**2 +0.341 )          
 if x >=5.9 and x<=8:    #FUV
  return 1.752 - 0.316*x - 0.104/( (x-4.67)**2 +0.341 )  - 0.04473*(x-5.9)**2 - 0.009779*(x-5.9)**3
 if x >=8 and x<=10:    #very FUV
  return -1.073 - 0.628*(x-8) + 0.137*(x-8)**2 - 0.070*(x-8)**3      
        
    
def b(x):
 #x is the wavelength inverse, in (micron)**-1     
 if x >=0.3 and x<=1.1: #optical range
  return -0.527*x**1.61  
 if x >=1.1 and x<=3.3: #optical/NIR range 
  y=x-1.82  
  return 1.41338*y + 2.28305*y**2 + 1.07233*y**3 - 5.38434*y**4 - 0.62251*y**5 + 5.30260*y**6 - 2.09002*y**7        
 if x >=3.3 and x<5.9:   #UV 
  return -3.090 + 1.825*x + 1.206/( (x-4.62)**2 +0.263 )   
 if x >=5.9 and x<=8:    #FUV 
  return -3.090 + 1.825*x + 1.206/( (x-4.62)**2 +0.263 )  + 0.2130*(x-5.9)**2 + 0.1207*(x-5.9)**3  
 if x >=8 and x<=10:    #very FUV    
  return 13.670 + 4.257*(x-8) - 0.420*(x-8)**2 + 0.374*(x-8)**3  
    
def Cardelli(x):  
 '''
 Cardelli, Clayton & Mathis (1989) extinction law.
 Rv is fixed to 3.1 (standard value for MW, diffuse MW) 
 x is the wavelength inverse, in (micron)**-1 
 '''
 Rv=3.1    
 return a(x) + b(x)/Rv

def Cardelli_freeRv(x,Rv):  
 '''
 Cardelli Clayton & Mathis (1989) extinction law.
 Rv is free. Examples:
 Rv = 2.5 : diffuse (~ M31)
 Rv = 3.1 : standard value for MW, diffuse MW) 
 Rv = 5   : dense medium
 x is the wavelength inverse, in (micron)**-1 
 '''   
 return a(x) + b(x)/Rv

#EXAMPLE OF USE
'''
#Creation of the lambda and lambda inverse
lamb = np.arange(0.1,3.33,0.01) #lambda in micron
x = 1./lamb[:]                  #in micron**-1
#print lamb[:]
#print x, len(x)

A_lambda_on_Av = np.zeros(len(x)) #Cardelli(x) 
for ii in range(0,len(x)):
 A_lambda_on_Av[ii] = Cardelli(x[ii])


A_lambda_on_Av_Rv25 = np.zeros(len(x))
for ii in range(0,len(x)):
 A_lambda_on_Av_Rv25[ii] = Cardelli_freeRv(x[ii],2.5)
'''
#======================================================================================================================


#======================================================================================================================
# CARDELLI (1989) + O'DONNEL (1994) EXTINCTION LAW
#======================================================================================================================
def Cardelli_COD(x, **kwargs):    
    """
    NAME:
     CCM_UNRED
    PURPOSE:
     Deredden a flux vector using the CCM 1989 parameterization
    EXPLANATION:
     The reddening curve is that of Cardelli, Clayton, and Mathis (1989 ApJ.
     345, 245), including the update for the near-UV given by O'Donnell 
     (1994, ApJ, 422, 158).   Parameterization is valid from the IR to the 
     far-UV (3.5 microns to 0.1 microns).    

     Users might wish to consider using the alternate procedure FM_UNRED
     which uses the extinction curve of Fitzpatrick (1999).
    
    CALLING SEQUENCE:
     ccm_unred(wave, flux, ebv [, R_V = ])      
     
    INPUT:
     WAVE - wavelength vector (Angstroms)
     FLUX - calibrated flux vector, same number of elements as WAVE
             If only 3 parameters are supplied, then this vector will
             updated on output to contain the dereddened flux.
     EBV  - color excess E(B-V), scalar.  If a negative EBV is supplied,
             then fluxes will be reddened rather than deredenned.

    OUTPUT:
     FUNRED - unreddened flux vector, same units and number of elements
             as FLUX

    OPTIONAL INPUT KEYWORD
     R_V - scalar specifying the ratio of total selective extinction
             R(V) = A(V) / E(B - V). If not specified, then R_V = 3.1
             Extreme values of R(V) range from 2.75 to 5.3

    EXAMPLE:
     Determine how a flat spectrum (in wavelength) between 1200 A and 3200 A
     is altered by a reddening of E(B-V) = 0.1.  Assume an "average"
     reddening for the diffuse interstellar medium (R(V) = 3.1)

       >>> w = 1200 + arange(40)*50       #Create a wavelength vector
       >>> f = w*0 + 1                    #Create a "flat" flux vector
       >>> fnew = ccm_unred(w, f, -0.1)   #Redden (negative E(B-V)) flux vector
       >>> plot(w,fnew)                   

    NOTES:
     (1) The CCM curve shows good agreement with the Savage & Mathis (1979)
             ultraviolet curve shortward of 1400 A, but is probably
             preferable between 1200 and 1400 A.
     (2)  Many sightlines with peculiar ultraviolet interstellar extinction 
             can be represented with a CCM curve, if the proper value of 
             R(V) is supplied.
     (3)  Curve is extrapolated between 912 and 1000 A as suggested by
             Longo et al. (1989, ApJ, 339,474)
     (4)  Use the 4 parameter calling sequence if you wish to save the 
             original flux vector.
     (5)  Valencic et al. (2004, ApJ, 616, 912) revise the ultraviolet CCM
             curve (3.3 -- 8.0 um-1).    But since their revised curve does
             not connect smoothly with longer and shorter wavelengths, it is
             not included here.

    REQUIRED MODULES:
     scipy, numpy
    REVISION HISTORY:
       Written   W. Landsman        Hughes/STX   January, 1992
       Extrapolate curve for wavelengths between 900 and 1000 A   Dec. 1993
       Use updated coefficients for near-UV from O'Donnell   Feb 1994
       Allow 3 parameter calling sequence      April 1998
       Converted to IDLV5.0                    April 1998
       Ported to Python        C. Theissen    August 2012
    """
    
    # Import modules   
    import numpy as n

    # Set defaults   
    R_V = 3.1

    for key in kwargs:
        if key.lower() == 'r_v':
            R_V = kwargs[key]
    '''                
    if isinstance(wave, int) or isinstance(wave, float):
        x = 10000. / n.array([wave])              # Convert to inverse microns
    else:
        x = 10000. / n.array(wave)                # Convert to inverse microns 
    '''
    npts = len( x )

    a = n.zeros((npts))  
    b = n.zeros((npts))
    
    ###############################

    good = n.where( (x > 0.3) & (x < 1.1) )     # Infrared
    Ngood = len(x[good])


    if Ngood > 0:
        a[good] =  0.574 * x[good]**(1.61)
        b[good] = -0.527 * x[good]**(1.61)

    ###############################

    good = n.where( (x >= 1.1) & (x < 3.3) )     # Optical/NIR
    Ngood = len(good[0])

    if Ngood > 0:                               # Use new constants from O'Donnell (1994)
        y = x[good] - 1.82
        #c1 = n.array([ 0.32999, -0.77530, 0.01979, 0.72085,        # Original
        #               -0.02427,  -0.50447, 0.17699, 1. ])         # coefficients              
        #c2 = n.array([ -2.09002, 5.30260, -0.62251, -5.38434,       # from CCM89
        #               1.07233, 2.28305, 1.41338, 0. ]) 
        c1 = n.array([ -0.505 , 1.647, -0.827, -1.718,              # New coefficients
                       1.137, 0.701, -0.609, 0.104, 1. ])           # from O'Donnell
        c2 = n.array([ 3.347,  -10.805, 5.491, 11.102,              # (1994)
                       -7.985, -3.989, 2.908, 1.952, 0. ])

        a[good] = n.polyval(c1, y)
        b[good] = n.polyval(c2, y)

    ###############################

    good = n.where( (x >= 3.3) & (x < 8) )                # Mid-UV
    Ngood = len(x[good])

    if Ngood > 0:
        y = x[good]
        F_a = n.zeros((Ngood))
        F_b = n.zeros((Ngood))
        good1 = n.where( (y > 5.9) )
        Ngood1 = len(y[good1])
        
        if Ngood1 > 0:
            y1 = y[good1] - 5.9
            F_a[good1] = -0.04473 * y1**2 - 0.009779 * y1**3
            F_b[good1] = 0.2130 * y1**2  +  0.1207 * y1**3
    
        a[good] = 1.752 - 0.316*y - (0.104 / ( (y-4.67)**2 + 0.341 )) + F_a
        b[good] = -3.090 + 1.825*y + (1.206 / ( (y-4.62)**2 + 0.263 )) + F_b

    ###############################

    good = n.where( (x >= 8) & (x < 11) )              #Far-UV
    Ngood = len(x[good])

    if Ngood > 0:
        y = x[good] - 8.
        c1 = [ -0.070, 0.137, -0.628, -1.073 ]
        c2 = [ 0.374, -0.420, 4.257, 13.670 ]
        
        a[good] = n.polyval(c1, y)
        b[good] = n.polyval(c2, y)

    ###############################

    # Now apply extinction correction to input flux vector

    A_V = 1. #R_V * ebv
    A_lambda = A_V * (a + b / R_V)
    return A_lambda
    #return flux * 10.**(0.4 * A_lambda)       # Derive unreddened flux


#EXAMPLE OF USE
#A_lambda_on_Av_COD = Cardelli_COD(x) 
#======================================================================================================================

































##### OLD! ######


'''
def imf_function(m):
 #IMF in number unit

 #Popescu_IMF(m):
 #k = 11.836901316 #1./0.0844815694  
 #if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)
 #if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.3) * m**(-1.3)
 #if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.35) * m**(-2.35)
 #if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.35) * m**(-2.35)

 #myGrid_IMF(m): #Kroupa01CB
 #k = 26.578777009 #1./0.03762400353
 #if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)                                 # imf = 12.5558*m**(-0.3)        
 #if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.8) * m**(-1.8)                          # imf = 0.284104*m**(-1.8)        
 #if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.8) * (1./0.5)**(-2.7) * m**(-2.7)         # imf = 0.152248*m**(-2.7)    
 #if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.8) * (1./0.5)**(-2.7) * m**(-2.3)         # imf = 0.152248*m**(-2.3) 

 #Kroupa01(m): #Canonical IMF, or Kr01NCB
 k = 11.045901141 #1./0.09053131901837103
 if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)                                 # imf = 5.177595858*m**(-0.3)        
 if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.3) * m**(-1.3)                                  # imf = 0.414207669*m**(-1.3)        
 if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3) * m**(-2.3)         # imf = 0.207103834*m**(-2.3)    
 if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3) * m**(-2.3)         # imf = 0.207103834*m**(-2.3) 

 #Kroupa02CB(m): #Kroupa 2002 corrected for binaries
 #k = 14.724042932 #1./0.06791612905804065
 #if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)                                 # imf = 6.901668114*m**(-0.3)        
 #if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.3) * m**(-1.3)                          # imf = 0.552133449*m**(-1.3)        
 #if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3) * m**(-2.3)         # imf = 0.276066725*m**(-2.3)    
 #if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.3) * m**(-2.7)         # imf = 0.276066725*m**(-2.7) 

 #Kroupa93(m):  #Kroupa 1993 
 #k = 14.229724031 #1./0.07027543175430642
 #if(m >= 0.01 and m < 0.08):   imf = k * (1./0.08)**(-0.3) * m**(-0.3)                                 # imf = 6.66996375*m**(-0.3)        
 #if(m >= 0.08 and m < 0.5):    imf = k * (1./0.08)**(-1.3) * m**(-1.3)                          # imf = 0.5335971*m**(-1.3)        
 #if(m >= 0.5 and m < 1.):      imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.2) * m**(-2.2)         # imf = 0.285947606*m**(-2.2)    
 #if(m >= 1. and m < 150.):     imf = k * (0.5/0.08)**(-1.3) * (1./0.5)**(-2.2) * m**(-2.7)         # imf = 0.285947606*m**(-2.7) 

 #Strange: before, for Kr93, i used this:                where is the mistake??? 
 #if(x >= 0.01 and x < 0.08):     imf = 16.98*x**(-0.3)        
 #if(x >= 0.08 and x < 0.5):      imf = 0.5518*x**(-1.3)        
 #if(x >= 0.5 and x < 1.):        imf = 0.2956*x**(-2.2)    
 #if(x >= 1.):                    imf = 0.2956*x**(-2.7)  

 #return m*imf #To check that integral of IMF is really 1. 
 return imf


def imf_function_m(m):  
        #IMF in mass unit
        return m*imf_function(m)


def CDF_imf_function_mass(m):                #CDF of the MASS IMF
        args = ()
        Integral_IMF = integrate.quad(imf_function_m, 0.01, m, args)
        return Integral_IMF
        #print CDF_imf_function(150.) #The integral of the IMF over the whole mass range is 1, of course. 


def CDF_imf_function_number(m_low, m):         #CDF of the NUMBER IMF!!
        args = ()
        Integral_IMF = integrate.quad(imf_function, m_low, m, args)
        #Integral_IMF = integrate.quad(imf_function, 0.01, m, args)
        return Integral_IMF
        #print CDF_imf_function(150.) #The integral of the IMF over the whole mass range is 1, of course. 


def Star_Cluster_stellar_masses_generator_N_star(m_low,m_up,m_low_used_array,m_up_used_array,array_CDF, array_CDF_inverse):
        #Now we can use it: Generate a random number between 0 and 1 and then interpolate in the array_CDF_inverse to get the mass. 
        N = len(m_low_used_array)
        cdf_low = np.zeros(N)
        cdf_up = np.zeros(N)
        cdf_low[:] = np.interp(m_low_used_array[:], array_CDF, array_CDF_inverse)
        cdf_up[:]  = np.interp(m_up_used_array[:] , array_CDF, array_CDF_inverse)
        stars_random_number_generated = np.random.rand(N) #*cdf_range
        stars_random_number_generated[:] =  stars_random_number_generated[:]*(cdf_up[:]-cdf_low[:])
        stars_random_number_generated = stars_random_number_generated + cdf_low
        stars_random_number_generated_sorted = np.sort(stars_random_number_generated, axis=-1, kind='mergesort')
        stars_mass_generated_sorted = np.interp(stars_random_number_generated_sorted, array_CDF_inverse, array_CDF)
        sort_choice = 0 #not sorted masses!
        if sort_choice == 1:
         return stars_mass_generated_sorted #If we want the stars sorted by mass, from lowest to highest
        if sort_choice == 0:
         stars_mass_generated_unsorted = np.zeros(N)
         stars_mass_generated_unsorted[np.argsort(stars_random_number_generated)] = stars_mass_generated_sorted
         return stars_mass_generated_unsorted #If we want the stars UNsorted by mass!


def Array_CDF_and_inverse(m_low,m_up):
        #This function returns only the array of the masses (array_CDF), and the array of the CDF (array_CDF_inverse)
        #mav = 1/CDF_imf_function_number(m_low,m_up)[0]
        normalization_of_number_CDF = CDF_imf_function_number(m_low,m_up)[0]
        array_CDF = np.linspace(-2.,2.176091259,1000) #in log10  #(0.15, m_up., num=1000)  in linear  120Ms = 2.079181246 in logm
        array_CDF = 10**array_CDF
        array_CDF_inverse = np.zeros(1000)
        for ii in range(0, len(array_CDF)):
         array_CDF_inverse[ii] = CDF_imf_function_number(m_low,array_CDF[ii])[0] / normalization_of_number_CDF
        return array_CDF, array_CDF_inverse


def extrap(x, xp, yp):
        """np.interp function with linear extrapolation"""
        y = np.interp(x, xp, yp)
        y = np.where(x<xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
        y = np.where(x>xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]), y)
        return y


def Intervals_of_resampling(m_low,m_up,star_optimal_index,star_optimal_mass):
        array_CDF, array_CDF_inverse = Array_CDF_and_inverse(m_low,m_up)  #(arrays needed: there are x and y of CDF)
        #I will need to use the limits given by the optimal sampling
        #So i need to take the positions of the stars in OS (their masses in fact)
        #Then i compute a vector from which the elements are the middle masses between stars. 
        #The first element is m_low, and the last is has to be derived!
        m_low_used_array          = np.zeros(len(star_optimal_index))
        m_up_used_array           = np.zeros(len(star_optimal_index))
        m_low_used_array_expanded = np.zeros(len(star_optimal_index))
        m_up_used_array_expanded  = np.zeros(len(star_optimal_index))
        interval_size             = np.zeros(len(star_optimal_index))
        ratio_high_on_low_parts   = np.zeros(len(star_optimal_index))
        rule_50pc_left                  = np.zeros(len(star_optimal_index))
        rule_50pc_right                  = np.zeros(len(star_optimal_index))
        CDF_star_optimal_mass          = np.zeros(len(star_optimal_index))
        CDF_star_optimal_mass_low = np.zeros(len(star_optimal_index)) 
        CDF_star_optimal_mass_up  = np.zeros(len(star_optimal_index)) 
 
         
        #Creation of this vector of half-masses elements, bins limits 
        #Method "Limits of intervals given by optimal sampling"  
        for ii in range(0,len(star_optimal_index)):
         if ii==0:                                 #First bin 
          m_low_used = 0.5*(star_optimal_mass[ii]+star_optimal_mass[ii+1])                        #In linear!!!
          m_up_used  = 150. #???? Could not be true!!! 
         if ii>0 and ii<len(star_optimal_index)-1:                                #Common bins
          m_low_used = 0.5*(star_optimal_mass[ii]+star_optimal_mass[ii+1])                        #In linear!!!
          m_up_used  = 0.5*(star_optimal_mass[ii-1]+star_optimal_mass[ii])                        #In linear!!!
         if ii==len(star_optimal_index)-1:        #Last bin 
          m_low_used = m_low
          m_up_used  = 0.5*(star_optimal_mass[ii-1]+star_optimal_mass[ii])                        #In linear!!!
         m_low_used_array[ii] = max(m_low_used,m_low)
         m_up_used_array[ii] = min(m_up_used,m_up)


        #Realizing the extrapolation for the most massive star superior limit
        index_regime_star_1 = 1 #index of the second most massive star of the last IMF regime
        #index_regime_star_2 = int(len(star_optimal_index)*0.03) #index of a star close to the inferior limit of the last IMF regime
         
        #index of a star close to the inferior limit of the relevant IMF regime. 
        #This is to garantee that we take the good regime for the most massive star interval!
        if star_optimal_mass[0]>1:
         itemindex=np.where(star_optimal_mass>1)
         index_regime_star_2 = itemindex[0][-1]
        if star_optimal_mass[0]<=1 and star_optimal_mass[0]>0.5:
         itemindex=np.where(star_optimal_mass>0.5)
         index_regime_star_2 = itemindex[0][-1]
        if star_optimal_mass[0]<=0.5 and star_optimal_mass[0]>0.08:
         itemindex=np.where(star_optimal_mass>0.08)
         index_regime_star_2 = itemindex[0][-1]
        if len(itemindex[0])<=2: index_regime_star_2=index_regime_star_2+2

        mass_1 = star_optimal_mass[index_regime_star_1] 
        mass_2 = star_optimal_mass[index_regime_star_2]
        mass_regime = star_optimal_mass[index_regime_star_1:index_regime_star_2+1]
        mass_regime = np.log10(mass_regime)
        mass_regime = mass_regime[::-1] #Reverse the order of the array, because need to be croissant numbers, for interpolation
        m_up_used_array_regime = m_up_used_array[index_regime_star_1:index_regime_star_2+1] - star_optimal_mass[index_regime_star_1:index_regime_star_2+1]
        m_up_used_array_regime = np.log10(m_up_used_array_regime)
        m_up_used_array_regime = m_up_used_array_regime[::-1]
        mass_extrapolated = extrap(np.log10(star_optimal_mass[0]), mass_regime, m_up_used_array_regime) 
        mass_extrapolated = (10**mass_extrapolated) + star_optimal_mass[0]
        m_up_used_array[0] = mass_extrapolated
        interval_size[:] = m_up_used_array[:] - m_low_used_array[:]
        ratio_high_on_low_parts[:] = (m_up_used_array[:]-star_optimal_mass[:]) / (star_optimal_mass[:]-m_low_used_array[:])


        normalization_of_number_CDF = CDF_imf_function_number(m_low,m_up)[0]
        Expansion_factor = 6.
        for ii in range(0,len(star_optimal_index)):
         #Linear growth of the intervals size
         m_low_used_array_expanded[ii] = max(star_optimal_mass[ii] - Expansion_factor*(star_optimal_mass[ii]-m_low_used_array[ii]), m_low)
         m_up_used_array_expanded[ii]  = min(star_optimal_mass[ii] + Expansion_factor*(m_up_used_array[ii]-star_optimal_mass[ii]),  m_up)
         
         rule_50pc_left[ii]  = CDF_imf_function_number(m_low_used_array_expanded[ii], star_optimal_mass[ii])[0] / normalization_of_number_CDF
         rule_50pc_right[ii] = CDF_imf_function_number(star_optimal_mass[ii],  m_up_used_array_expanded[ii])[0] / normalization_of_number_CDF

         CDF_star_optimal_mass[ii] = CDF_imf_function_number(m_low, star_optimal_mass[ii])[0] / normalization_of_number_CDF
         CDF_star_optimal_mass_low[ii] = CDF_star_optimal_mass[ii] - rule_50pc_right[ii]
         CDF_star_optimal_mass_up[ii]  = CDF_star_optimal_mass[ii] + rule_50pc_right[ii]
         CDF_star_optimal_mass_low[ii] = max(CDF_star_optimal_mass_low[ii],0.)  #Saturation to low masses
         CDF_star_optimal_mass_up[ii]  = min(CDF_star_optimal_mass_up[ii] ,1.)  #Saturation to high masses
         #Now i need to retransform in mass the CDF_star_optimal_mass_low and CDF_star_optimal_mass_up...
         m_low_used_array_expanded[ii] = np.interp(CDF_star_optimal_mass_low[ii], array_CDF_inverse, array_CDF)

         m_low_used_array_expanded[ii] = max(m_low_used_array_expanded[ii],m_low)
         m_up_used_array_expanded[ii] = min(m_up_used_array_expanded[ii],m_up)
         rule_50pc_left[ii]  = CDF_imf_function_number(m_low_used_array_expanded[ii], star_optimal_mass[ii])[0] / normalization_of_number_CDF
         rule_50pc_right[ii] = CDF_imf_function_number(star_optimal_mass[ii],  m_up_used_array_expanded[ii])[0] / normalization_of_number_CDF


        ##Print for every cluster these data:  (could be useful for analyze of the goodness of the code)
        #array  = np.zeros((len(star_optimal_index),10))
        #array[:,0] = m_low_used_array
        #array[:,1] = m_up_used_array
        #array[:,2] = star_optimal_mass
        #array[:,3] = interval_size
        #array[:,4] = ratio_high_on_low_parts
        #array[:,5] = m_low_used_array_expanded
        #array[:,6] = m_up_used_array_expanded
        #array[:,7] = rule_50pc_left
        #array[:,8] = rule_50pc_right
        #array[:,9] = CDF_star_optimal_mass
        #np.savetxt('test.out', array,fmt='%.6f')
        #raw_input() 


        m_low_used_array = m_low_used_array_expanded
        m_up_used_array = m_up_used_array_expanded

        return m_low_used_array, m_up_used_array
'''




















