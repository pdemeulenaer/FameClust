# FameClust

==============================================================================

         FameClust: Finding Age, Mass, Extinction of star CLUSTers (astrophysical package written in python and Fortran)

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


Articles based on this method:

-------------------------------------

* de Meulenaer, P.; Narbutis, D.; Mineikis, T.; Vansevičius, V.	
  Deriving physical parameters of unresolved star clusters. I. Age, mass, and extinction degeneracies (http://adsabs.harvard.edu/abs/2013A%26A...550A..20D . PDF here: https://www.aanda.org/articles/aa/pdf/2013/02/aa20674-12.pdf)

* de Meulenaer, P.; Narbutis, D.; Mineikis, T.; Vansevičius, V.	
  Deriving physical parameters of unresolved star clusters. II. The degeneracies of age, mass, extinction, and metallicity (http://adsabs.harvard.edu/abs/2014A%26A...569A...4D . PDF here: https://www.aanda.org/articles/aa/pdf/2014/09/aa23988-14.pdf)

* de Meulenaer, P.; Narbutis, D.; Mineikis, T.; Vansevičius, V.	
  Deriving physical parameters of unresolved star clusters. III. Application to M 31 PHAT clusters (http://adsabs.harvard.edu/abs/2015A%26A...574A..66D . PDF here: https://www.aanda.org/articles/aa/pdf/2015/02/aa25121-14.pdf)

* de Meulenaer, P.; Narbutis, D.; Mineikis, T.; Vansevičius, V.	
  Deriving physical parameters of unresolved star clusters. IV. The M 33 star cluster system (http://adsabs.harvard.edu/abs/2015A%26A...581A.111D . PDF here: https://www.aanda.org/articles/aa/pdf/2015/09/aa26544-15.pdf)

* de Meulenaer, P.; Stonkutė, R.; Vansevičius, V.	
  Deriving physical parameters of unresolved star clusters. V. M 31 PHAT star clusters (http://adsabs.harvard.edu/abs/2017A%26A...602A.112D . PDF here: https://ui.adsabs.harvard.edu/link_gateway/2017A%26A...602A.112D/PUB_PDF)


