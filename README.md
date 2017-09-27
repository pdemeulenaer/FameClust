# FameClust

==============================================================================

         FameClust: Finding Age, Mass, Extinction of star CLUSTers

==============================================================================

Philippe de Meulenaer
Vilnius University Observatory, FTMC
PhD in Astrophysics


Codes are running under Ubuntu 16.04,
Python 3.5.3 |Anaconda custom (64-bit)| (default, Mar  6 2017, 11:58:13) 
[GCC 4.4.7 20120313 (Red Hat 4.4.7-1)] on linux


==========================
1. CLUSTER GENERATION CODE
==========================

Launch: python GenGridModels.py 800 400 n00

Time: log(t) = 8 -> 800 [ranges from 6.60 to 10.10, step of 0.05]

Mass: log(m) = 4 -> 400 [ranges from 2.00 to 7.00, step of 0.05]

Metallicity: solar -> n00 [-2.2=m22,...,-0.2=m02, then solar=n00, then +0.2=p02, step of 0.2]



==========================
2. FAMEClust CODE
==========================

Compilation: gfortran -O3 -march=native -ffast-math -funroll-loops -mcmodel=medium Module_lecture.f90 FameClust_V101_observation.f90 -o FameClust_V101_observation.exe

Execution: time ./FameClust_V101_observation.exe InputFameClustNEW_WFC3_PHAT_ALL_new 1 10 n00

where:
- InputFameClustNEW_WFC3_PHAT_ALL_new is the input file
- 1 and 10 are the first and last clusters analysed from the list (sometimes no need to analyse them all)
- n00 is the metallicity of the grid used, here solar [-2.2=m22,...,-0.2=m02, then solar=n00, then +0.2=p02, step of 0.2]

The InputFameClustNEW_WFC3_PHAT_ALL_new contains the path of the cluster data to analyse, the model grid, and others. It is important to adapt these paths to your system to make the code work. The grids of models for different metallicities can be given on request (write to me, pdemeulenaer@gmail.com)


