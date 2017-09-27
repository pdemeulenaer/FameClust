 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !											    !
 !											    !
 !		FameClust 10.1 and Chi Square methods 		 			    !	
 !		(Finding of Age Mass and Extinction of star Clusters)			    !
 !		Philippe de Meulenaer, PhD Student in Astrophysics (year three)		    !
 !		Astronomical Observatory, Vilnius University				    !
 !											    !
 !											    !
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
 ! Marie je te confie cela...
 !
 ! Date of conception: 28 november 2013
 ! Last update : 28 november 2013
 ! OLD Compilation : gfortran -O2 -march=native -mcmodel=medium Module_lecture.f90 FameClust_V101_observation.f90 -o FameClust_V101_observation.exe
 ! FASTER: gfortran -O3 -march=native -ffast-math -funroll-loops -mcmodel=medium Module_lecture.f90 FameClust_V101_observation.f90 -o FameClust_V101_observation.exe
 ! Compilation VHD : -mcmodel=large 
 ! Example: time ./FameClust_V101_observation.exe InputFameClustNEW_UBVRI_Z01900_M400 1 10 n00
 ! 
 !	   (No pain, no gain!)
 !         (Let's DO it!)
 !        O 
 !    ___o
 !   (*,*)
 !   (   )
 !---"--"----
 !
 !
 !!!!!!!!!!!!!!!!!!!!!!!!!!
 !IMPORTATIONS
 !!!!!!!!!!!!!!!!!!!!!!!!!!

 Program main
 use Lecture_module_fortran
 implicit none

 !!!!!!!!!!!!!!!!
 !INITIALISATIONS
 !!!!!!!!!!!!!!!!
 integer, parameter :: N_bin_all_nodes = 7171000 !5041000 !4331000 !5751000 !4331000 !106253000
 real(4), allocatable :: magnitude_GRID(:,:)
 real(4), allocatable :: age_list(:), mass_list(:), Ebv_list(:), sigma_obs_vector(:)
 real(4), allocatable :: age_max_chi2(:),mass_max_chi2(:),Ebv_max_chi2(:), chi2(:)
 real(4), allocatable :: A_lambda_filters_selected(:),M1_input(:,:),M1(:),M2(:),gaps_in_filter(:,:)
 real(4), allocatable :: A_lambda_filters_selected_foreground(:) !For foreground extinction!!!
 real(4), allocatable :: M1_input_initial(:,:)
 real(4) :: chi2_min, chi2_min_old
 real(4) :: age, mass, Z, Ebv, age_maximal, mass_maximal, Z_maximal 
 real(4) :: distance_modulus
 real(4) :: node_solution(1:71,1:101,1:301,1:9) !node_solution(1:71,1:101,1:121,1:9)
 real(4) :: lambda, lambda_f_MW, lambda_f_LMC,lambda_f_M31
 real(4) :: Probability_node,Probability_node_2,Probability_node_3,Probability_node_4
 real(4) :: GRID_read(1:1000,1:52),data_input(1:52),data_sigma_input(1:52),xx
 real(4) :: Z_selected, Rv, index_real
 real(4), allocatable :: Grid_completed(:,:), M0(:,:), M1_M0(:,:),M2_M1(:),M0_prime(:,:), M1_M0_prime(:,:)
 real(4), allocatable :: d_direct_square(:),d_parallel_square(:)
 real(4), allocatable :: d_perpendicular_square(:),d_perpendicular(:)
 real(4), allocatable :: Ebv_vector(:), Grid_completed_selection(:,:)
 real(4), allocatable :: age_max_CG(:),mass_max_CG(:),Ebv_max_CG(:),Proba_max_CG(:)
 real(4), allocatable :: age_max_d2(:),mass_max_d2(:),Ebv_max_d2(:),Proba_max_d2(:)
 real(4), allocatable :: k(:), sigma_filter(:,:), sigma_filter_inverse(:,:) !General case
 real(4), allocatable :: sum_sigma_inverse_square(:), sum_sigma_inverse(:),proba_model_vector(:)
 real(4) :: M2_M1_square, M2_M1_square_inverse
 real(4) :: sigma_magnitude,counting_inverse, proba_model, sigma_filter_automatic
 real(4) :: k_lower, k_higher, a_sigma, b_sigma, c_sigma, sigma_exp
 real(4) :: histo_age(1:401), histo_mass(1:401), histo_Ebv(1:401), histo_Z(1:401)
 real(4) :: Rv_foreground, Ebv_foreground
 integer, allocatable :: sigma1_total(:), sigma2_total(:), sigma3_total(:)
 integer, allocatable ::  filters_selected(:)
 integer, allocatable :: Cluster_ID(:),counting(:),OB_size(:), ID_models(:),nodata_all(:)
 integer(8) :: compteur
 integer :: sigma1_number, sigma2_number, sigma3_number
 integer :: Ext_limit1,Ext_limit2,filter_ID
 integer :: ff,hh,ii,jj,kk,ll,a,m,zz,Ext, number_cluster,n_lines,node
 integer :: list, number_cluster_observed
 integer :: multi_1,multi_2, method, choice, choice_extinction, weight
 integer :: number_filters, choice_filters, filters(1:52), choice_extinction_law
 integer :: idum, choice_noise, app_or_abs, switch, age_int, mass_int, aa,mm,zz_bis
 integer ::  number_begin, number_end, in_box, sigma_factor, choice_sigma
 integer :: max_proba_position(1:3)
 character(len=50) :: file_name,arg
 character(len=200) :: InputFile_Name
 character(len=50) :: age_indice,mass_indice,Z_indice,Ebv_indice,noise_flag,extinction_flag
 character(len=50) :: extinction_file_name
 character(len=50) :: multi_char,jj_char, file_out_cluster2
 character(len=300) :: file_out_cluster, file_out_cluster_f90, file_name_grid,file_Test_1000_random_clusters
 character(len=300) :: file_observed_clusters,file_name_nodes_bin,forma, file_out_cluster_models
 character(len=300) :: file_out_cluster_node,file_out_cluster_histo
 CHARACTER(len=3) :: Z_indice_selected
 write(*,*)
 CALL system('date')

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !GETTING THE ARGUMENTS OF THE COMMAND-LINE 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 call getarg(1, InputFile_Name)
 !InputFile_Name = '/home/philippe/Desktop/Discrete_models_comparaison_jtao&
 !	&/SC_Parameters_20/'//trim(adjustl(adjustr(InputFile_Name)))
 InputFile_Name = '/home/philippe/Desktop/FameClust/'//trim(adjustl(adjustr(InputFile_Name))) 
 call getarg(2, arg) 
 read (arg,'(I10)') number_begin
 call getarg(3, arg)
 read (arg,'(I10)') number_end
 call getarg(4, Z_indice_selected) 
 call age_mass_Z_December_Inverse(aa,mm,zz,Ext,age,mass,Z,Ebv, Z_indice_selected)
 call age_mass_Z_December(a,m,zz,Ext, age, mass, Z_selected, Ebv, age_indice, mass_indice, Z_indice, Ebv_indice)
 write(*,*)Z_indice,zz,Z_selected




 !# -------------------------
 !# Loading of the input file
 !# -------------------------
 filters(:) = 0
 open(unit=10,file=InputFile_Name)
 READ(10,*)
 READ(10,*)
 READ(10,*)number_filters
 READ(10,*)
 allocate(filters_selected(1:number_filters),A_lambda_filters_selected(1:number_filters))
 allocate(A_lambda_filters_selected_foreground(1:number_filters)) !For foreground extinction!!!
 DO ii = 1,number_filters
  READ(10,*) choice_filters
  filters_selected(ii)=choice_filters
  filters(choice_filters) = 1
 ENDDO
 READ(10,*)
 READ(10,*) distance_modulus  		!M33: 24.54 (McConnachie2004;2005) !M31: 24.47 (Narbutis2008)
 READ(10,*)
 READ(10,*) app_or_abs			!apparent or absolute magnitude? (1/2)
 READ(10,*)
 READ(10,*) file_observed_clusters	!name of the file containing the observed clusters
 READ(10,*)
 READ(10,*) number_cluster_observed	!How many clusters are there in your file
 READ(10,*)
 READ(10,*) choice_extinction		!cluster(s) studied extincted or not ? (1/2)
 READ(10,*)
 READ(10,*) choice_extinction_law
 READ(10,*)
 !READ(10,*) choice_sigma		!Automatic sigma (1) input file sigma (2)
 !READ(10,*)
 READ(10,'(a)') file_out_cluster  	!This format '(a)' because of the slashes in the name of directories
 write(*,*) file_out_cluster		!Path where to store the output files
 close(10)


 !# -----------------------------------
 !# Allocations of the different tables
 !# -----------------------------------
 ALLOCATE(Cluster_ID(1:number_cluster_observed))
 ALLOCATE(Grid_completed(1:N_bin_all_nodes,1:4+number_filters+2))
 ALLOCATE(M0(1:N_bin_all_nodes,1:number_filters),M0_prime(1:N_bin_all_nodes,1:number_filters))
 ALLOCATE(M1_input(1:number_cluster_observed,1:number_filters))
 ALLOCATE(M1_input_initial(1:number_cluster_observed,1:number_filters))
 ALLOCATE(M1(1:number_filters),M2(1:number_filters))
 ALLOCATE(M1_M0(1:N_bin_all_nodes,1:number_filters),M1_M0_prime(1:N_bin_all_nodes,1:number_filters))
 ALLOCATE(M2_M1(1:number_filters))
 ALLOCATE(d_direct_square(1:N_bin_all_nodes))
 ALLOCATE(d_parallel_square(1:N_bin_all_nodes))
 ALLOCATE(d_perpendicular_square(1:N_bin_all_nodes))
 ALLOCATE(d_perpendicular(1:N_bin_all_nodes))
 ALLOCATE(Ebv_vector(1:N_bin_all_nodes))
 ALLOCATE(Grid_completed_selection(1:7171000,1:4+number_filters+2+number_filters))
 ALLOCATE(age_max_CG(1:number_cluster_observed),mass_max_CG(1:number_cluster_observed))
 ALLOCATE(Ebv_max_CG(1:number_cluster_observed),Proba_max_CG(1:number_cluster_observed))
 ALLOCATE(age_max_d2(1:number_cluster_observed),mass_max_d2(1:number_cluster_observed))
 ALLOCATE(Ebv_max_d2(1:number_cluster_observed),Proba_max_d2(1:number_cluster_observed))
 ALLOCATE(k(1:N_bin_all_nodes),counting(1:number_cluster_observed))
 ALLOCATE(OB_size(1:number_cluster_observed))
 ALLOCATE(sigma_filter(1:number_cluster_observed,1:number_filters))	!General case
 ALLOCATE(sigma_filter_inverse(1:number_cluster_observed,1:number_filters))	!General case
 ALLOCATE(sum_sigma_inverse_square(1:number_cluster_observed))
 ALLOCATE(sum_sigma_inverse(1:number_cluster_observed))
 ALLOCATE(ID_models(1:N_bin_all_nodes))
 ALLOCATE(gaps_in_filter(1:number_cluster_observed,1:number_filters))
 ALLOCATE(nodata_all(1:number_cluster_observed),proba_model_vector(1:number_filters))


 !# -----------------------------------------------------------------------------------------
 !#Loading of the observations
 !# -----------------------------------------------------------------------------------------
 !choice_sigma=1	!Automatic sigmas
 !choice_sigma=2	!maximum between data sigmas and power-law sigmas
 choice_sigma=3	! data sigmas and power-law sigmas
 CALL chdir(file_out_cluster)
 OPEN(unit=11,file=file_observed_clusters)
 READ(11,*)
 M1_input(:,:)=0.
 DO list = 1, number_cluster_observed
  !If we want input sigma, for each filter, for each cluster observed
  if (choice_sigma==2 .or. choice_sigma==3) then 
   READ(11,*) Cluster_ID(list), xx,xx,xx,xx, data_input(1:52), xx,xx,xx,xx, data_sigma_input(1:52)
   do ff = 1,number_filters   
    sigma_filter(list,ff) = data_sigma_input(filters_selected(ff)) 
   enddo
  !If choice_sigma=1, automated sigma are the ones loaded from Filters_information.dat file
  else				
   READ(11,*) index_real, xx,xx,xx,xx, data_input(1:52)
   Cluster_ID(list) = nint(index_real)
  endif

  !write(*,*)number_cluster_observed, Cluster_ID(list), data_input(11:16)
  !read(*,*)

  do ff = 1,number_filters   
   M1_input(list,ff) = data_input(filters_selected(ff))
   M1_input_initial(list,ff) = data_input(filters_selected(ff))
  enddo
  IF (app_or_abs == 1) then
   do ii = 1, number_filters
    !M1_input(ii,list) = M1_input(ii,list) - distance_modulus
    M1_input(list,ii) = M1_input(list,ii) - distance_modulus
   enddo
  ENDIF
 ENDDO
 CLOSE(11)
 !write(*,*) M1_input(1,1:2)

   

 !FOR GALEX: AB -> vega system!!! (GALEX in my grid is given in vega system)
 !DO ff = 1,number_filters
 ! if (filters_selected(ff)==1) M1_input_initial(:,ff) = M1_input_initial(:,ff) + 2.128
 ! if (filters_selected(ff)==2) M1_input_initial(:,ff) = M1_input_initial(:,ff) + 1.662 
 !ENDDO

 !FOR GALEX: ARTIFICIAL BRIGHTENING OF FUV AND NUV (test)
 !DO ff = 1,number_filters
 ! if (filters_selected(ff)==1) M1_input_initial(:,ff) = M1_input_initial(:,ff) -1.
 ! if (filters_selected(ff)==2) M1_input_initial(:,ff) = M1_input_initial(:,ff) -1.
 !ENDDO

 !FOR CFHT of SR10: CALIBRATION TO FAN14
 !if (filters_selected(6)==17) write(*,*)'CALIBRATION CFHT ON'
 !DO ff = 1,number_filters
 ! if (filters_selected(ff)==17) M1_input(:,ff) = M1_input(:,ff) - 0.342
 ! if (filters_selected(ff)==18) M1_input(:,ff) = M1_input(:,ff) - 0.132 
 ! if (filters_selected(ff)==19) M1_input(:,ff) = M1_input(:,ff) - 0.152 
 ! if (filters_selected(ff)==20) M1_input(:,ff) = M1_input(:,ff) - 0.080 
 ! if (filters_selected(ff)==21) M1_input(:,ff) = M1_input(:,ff) - 0.080 
 !ENDDO
 !FOR UBVRI of Ma12: CALIBRATION TO FAN14
 !if (filters_selected(1)==3) write(*,*)'CALIBRATION UBVRI ON'
 !DO ff = 1,number_filters
 ! if (filters_selected(ff)==3) M1_input(:,ff) = M1_input(:,ff) - 0.235
 ! if (filters_selected(ff)==4) M1_input(:,ff) = M1_input(:,ff) - 0.170 
 ! if (filters_selected(ff)==5) M1_input(:,ff) = M1_input(:,ff) - 0.102 
 ! if (filters_selected(ff)==6) M1_input(:,ff) = M1_input(:,ff) - 0.136 
 ! if (filters_selected(ff)==7) M1_input(:,ff) = M1_input(:,ff) - 0.182 
 !ENDDO

 !read(*,*)


 !Indication of the possible gaps in input data of observed clusters
 gaps_in_filter(:,:) = 1
 DO list = 1, number_cluster_observed
  do ff=1,number_filters
   if (M1_input(list,ff) >= 40.) then
    gaps_in_filter(list,ff) = 0
   endif
  enddo
 Enddo



 !# -----------------------------------------------------------------------------------------
 !#Loading of the A_lambda (extinction parameters) for the filters selected in the input file
 !# -----------------------------------------------------------------------------------------
 !CALL chdir('/home/philippe/Desktop/Discrete_models_comparaison_jtao/SC_Parameters_20/Source/')
 CALL chdir('/home/philippe/Desktop/FameClust/Source/')
 !open(unit=18, file = 'Filters_information.dat')
 !open(unit=18, file = 'Filters_information_observations.dat')
 open(unit=18, file = 'Filters_information_observations_WFC3_from_STSCI.dat')
 !open(unit=18, file = 'Filters_information_observations_WFC3_from_STSCI_test.dat')
 read(18,*)
 jj=0
 do ii = 1,52 !46
  read(18,*)lambda, lambda_f_MW, lambda_f_LMC, lambda_f_M31, filter_ID, sigma_filter_automatic, xx, a_sigma, b_sigma, c_sigma
  if (filters(ii) == 1) then
   jj=jj+1
   if (choice_extinction_law == 1) then		!Case MW
    Rv = 3.1
    A_lambda_filters_selected(jj)=lambda_f_MW
   elseif (choice_extinction_law == 2) then	!Case LMC average (Gordon 2003)
    Rv = 3.4
    A_lambda_filters_selected(jj)=lambda_f_LMC
   elseif (choice_extinction_law == 3) then	!Case M31
    Rv = 5.2 !2.5 for M31
    A_lambda_filters_selected(jj)=lambda_f_M31
   endif

   A_lambda_filters_selected_foreground(jj)=lambda_f_MW !For the foreground MW extinction!!!

   if (choice_sigma==1) then	  !Automatic sigmas, same for all observed clusters  
    sigma_filter(:,jj) = sigma_filter_automatic

   elseif (choice_sigma==3) then  !sigmas taken from data
    !If sigma lower than 0.05mag, --> 0.05mag.
    do list = 1, number_cluster_observed
     if (gaps_in_filter(list,jj) == 1) then  !only for the ones with data
      !sigma_filter(list,jj) = sigma_filter(list,jj)*2 	!Criteria decrease importance of bad data
      !if (sigma_filter(list,jj)<0.01) sigma_filter(list,jj) = 0.01  	!Minimum criteria  
      if (sigma_filter(list,jj)<0.03) sigma_filter(list,jj) = 0.03  	!Very low criteria  
      !if (sigma_filter(list,jj)<0.05) sigma_filter(list,jj) = 0.05  	!Classical criteria  
      !if (sigma_filter(list,jj)<0.07) sigma_filter(list,jj) = 0.07   	!large criteria 
      !sigma_filter(list,jj) = sigma_filter(list,jj) + 0.05 	!additional uncertainty
     elseif (gaps_in_filter(list,jj) == 0) then  !prevent gaps
      sigma_filter(list,jj) = 0.
     endif
    enddo
   endif 

  endif
 enddo
 close(18)


 !DO list = number_begin,number_end
  !if (gaps_in_filter(list,1)==1) sigma_filter(list,1) = sigma_filter(list,1) + 0.03
  !if (gaps_in_filter(list,2)==1) sigma_filter(list,2) = sigma_filter(list,2) + 0.03
  !if (gaps_in_filter(list,3)==1) sigma_filter(list,3) = sigma_filter(list,3) + 0.03
  !if (gaps_in_filter(list,4)==1) sigma_filter(list,4) = sigma_filter(list,4) + 0.05
  !if (gaps_in_filter(list,5)==1) sigma_filter(list,5) = sigma_filter(list,5) + 0.30
  !if (gaps_in_filter(list,6)==1) sigma_filter(list,6) = sigma_filter(list,6) + 0.50
 !ENDDO



 sum_sigma_inverse(:)=0.
 sum_sigma_inverse_square(:)=0.
 do list=1,10
  write(*,*)
  do ff=1,number_filters
   write(*,*)ff,A_lambda_filters_selected(ff),sigma_filter(list,ff), M1_input_initial(list,ff)
   sigma_filter_inverse(:,ff) = 1/sigma_filter(:,ff)
  enddo
 enddo

 DO list = number_begin,number_end
  do ff=1,number_filters
   if (gaps_in_filter(list,ff)==1) then
    sum_sigma_inverse_square(list) = sum_sigma_inverse_square(list) + sigma_filter_inverse(list,ff)
   endif
  enddo
 ENDDO

 sum_sigma_inverse(:) = sum_sigma_inverse_square(:)**0.5
 !read(*,*)




 !# ----------------------------------------------
 !# Loading of the grid of models (from .bin file)
 !# ----------------------------------------------
 !CALL chdir('/home/philippe/Desktop/Discrete_models_comparaison_jtao/SC_Parameters_20/Source/')
 CALL chdir('/home/philippe/Desktop/FameClust/Source/')
 Grid_completed(:,:) = 0.
 CALL system('date')
 WRITE(*,*) '' 
 WRITE(*,*) ' Loading of the grid in the program'
 ii=1
 Do aa = 1,71  		!Loop on the age. 
  Do mm = 1,101 ![2.7-5.0] !81  !Loop on the mass.  
    Ext = 1
    call age_mass_Z_December(aa,mm,zz,Ext,age,mass,Z,Ebv,age_indice,mass_indice,Z_indice,Ebv_indice)

    !gCMD Grid!
    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grid_gCMD014_71_81_Z'
    !file_name_nodes_bin=trim(file_name_nodes_bin)//trim(Z_indice)//'_AllB_HST_20pc_untruncated_iso_binary/'
    !file_name_nodes_bin=trim(file_name_nodes_bin)//trim(Z_indice)//'_AllB_HST_20pc_untruncated_iso_Kr01NCB_binary/'

    !FRS Grid!
    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grids_with_ACS/Grid_FRS_Z'   !HDD
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_Kroupa_1000models_per_node_with_ACS_binary/'
    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grids_with_ACS/'
    file_name_nodes_bin = '/media/philippe/a36c9bac-04fd-4046-8af2-2038962a2127/philippe/Documents/'&
      &//'PhD/Discrete_models_comparaison_jtao/Grids_with_ACS/'
    file_name_nodes_bin = trim(file_name_nodes_bin)//'/Grid_FRS_1000models_per_node_OldGenerationCode/'
    !file_name_nodes_bin = trim(file_name_nodes_bin)//'/Grid_FRS_1000models_per_node_NewGenerationCode/'
    file_name_nodes_bin = trim(file_name_nodes_bin)//'GRID_files/Grid_FRS_generated_from_npz/BINARY_files/' // trim(Z_indice) //'/' !Same grid as above, just different place

    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grids_with_ACS/PEGASE_FRS_grids/Grid_FRS_Z'   !HDD
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_Kroupa_1000models_per_node_with_ACS_PEGASE_binary/'

    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grid_Z'              !TRADITIONAL
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_test2_binary/'

    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grid_Z'
    !file_name_nodes_bin = trim(file_name_nodes_bin) //adjustl(trim(Z_indice))//'_from_npzGMM_binary/'

    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grids_with_ACS/'
    !file_name_nodes_bin = trim(file_name_nodes_bin) //'Grid_FRS_generated_not_from_npz/Grid_FRS_Z'
    !file_name_nodes_bin = trim(file_name_nodes_bin) //adjustl(trim(Z_indice))//'_Kroupa_1000models_per_node_with_ACS_binary/'



    !HRS Grid!
    !file_name_nodes_bin = '/mnt/storage/philippe/Grid_HRS_Z'
    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grid_HRS_Z'
    !file_name_nodes_bin = '/opt/Grid_HRS_Z'
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_ExpFactor6_Kroupa_binary/'
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_ExpFactor6_Weidner_corrected_binary/'
    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grid_HRS_Z'   !HDD
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_ExpFactor6_Weidner_correctedTEST_binary/'
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_ExpFactor6_Weidner_corrected_interpolated_binary/'
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_ExpFactor6_Weidner_corrected_with_ACS_binary/'
    !file_name_nodes_bin = '/opt/Grid_HRS_Z'   !SDD
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_ExpFactor6_Weidner_corrected_binary/'


    !SSP grid!!
    !file_name_nodes_bin = '/home/philippe/Desktop/Discrete_models_comparaison_jtao/Grid_SSP_CMD25/Grid_SSP_Z'
    !file_name_nodes_bin = trim(file_name_nodes_bin) // trim(Z_indice) // '_nonoise_binary/'

    file_name_nodes_bin = trim(file_name_nodes_bin) // 't' // trim(age_indice) // '_M' //&
	& trim(mass_indice) // '_Z' // trim(Z_indice)
    call Lecture_UBVRI_ugriz_GRID_Jan2013(file_name_nodes_bin,GRID_READ) 
    do ff=1,number_filters 
     Grid_completed(ii:ii+999,4+ff) = GRID_READ(:,filters_selected(ff))
    enddo
    Grid_completed(ii:ii+999,1) = age
    Grid_completed(ii:ii+999,2) = log10(mass)
    Grid_completed(ii:ii+999,4) = Z
    ii=ii+1000
  Enddo
  write(*,*)age
 Enddo
 do ff=1,number_filters 
  M0(:,ff) = Grid_completed(:,4+ff) 
 enddo
 
 
 WRITE(*,*) ' The grid has been loaded correctly, from: ', file_name_nodes_bin 
 WRITE(*,*) '' 
 CALL system('date')




 !# -------------------------------------------------------
 !# Preparation of the extinction (foreground used or not?)
 !# -------------------------------------------------------

 Ebv_foreground = 0. !0.06 !mag. If 0, no foreground extinction adopted
 Rv_foreground  = 3.1  !As the standard value for the MW (using CCM 1989)

 do list = number_begin,number_end
  !De-reddening of the observation by the foreground extinction
  M1_input(list,:) = M1_input(list,:) - A_lambda_filters_selected_foreground * Rv_foreground * Ebv_foreground
 enddo

 !FOR M31
 !k_lower  = 0.02 	!Ebv_foreground=0.04 == k_higher = 0.02	!M31
 !k_lower  = 0.03 	!Ebv_foreground=0.06 == k_higher = 0.03	!M31  !as in Fouesneau 2014, using Schlegel value, E(B-V)=0.062
 k_lower  = 0.03   !IN CASE OF FOREGROUND EXTINCTION ENABLED!
 k_higher = 0.55 !1.5 !0.55 !0.5 	!Ebv_limit = 1.0 == k_higher = 0.5	!M31   if 1. -> Ebv_limit = 2.

 !FOR M33
 !k_lower  = 0.02 	!Ebv_foreground=0.04 == k_higher = 0.02	!M31
 !k_lower  = 0.03 	!Ebv_foreground=0.06 == k_higher = 0.03	!M31  !as in Fouesneau 2014, using Schlegel value, E(B-V)=0.062
 !k_higher = 0.50 !0.15 !for Ebvlimit=0.30 !0.5 	!Ebv_limit = 1.0 == k_higher = 0.5	!M31   if 1. -> Ebv_limit = 2.

 !k_lower  = 0.0 	!Ebv_foreground=0.06 == k_higher = 0.03	!M31  !as in Fouesneau 2014, using Schlegel value, E(B-V)=0.062
 !k_higher = 0.50 !0.15 !for Ebvlimit=0.30 !0.5 	!Ebv_limit = 1.0 == k_higher = 0.5	!M31   if 1. -> Ebv_limit = 2.


 !just to check the foreground-dereddened photometry
 do list=1,10
  write(*,*)
  do ff=1,number_filters
   write(*,*)ff,M1_input_initial(list,ff) - distance_modulus, M1_input(list,ff)
  enddo
 enddo
 !read(*,*)



 !# --------------------------------------------------------
 !# THE MAIN LOOP OF PARAMETER DERIVATION (for all clusters)
 !# --------------------------------------------------------

 DO list = number_begin,number_end
	
	!# ------------------------------------------------------------------------
	!# Reddening of the observation, to obtain the line going through M1 and M2
	!# ------------------------------------------------------------------------
	M1(:) = M1_input(list,:)                      !We take the observation from the list of observed clusters:
	M2 = M1 - A_lambda_filters_selected * Rv*2    !We deredden observation by a quantity E(B-V)=2, to have the line (M2,M1)


	!# ----------------------------------------------------------
	!# Derivation of the distances d_perpendicular and d_parallel
	!# ----------------------------------------------------------
        do ff=1,number_filters
	 M1_M0(:,ff) = ( M1(ff) - M0(:,ff) ) * sigma_filter_inverse(list,ff)/sum_sigma_inverse(list)
         if (gaps_in_filter(list,ff)==0) then
	  M1_M0(:,ff) = 0.
	 endif
        enddo

	M2_M1 = (M2(:)-M1(:)) * sigma_filter_inverse(list,:)/sum_sigma_inverse(list)
        do ff=1,number_filters
         if (gaps_in_filter(list,ff)==0) then
	  M2_M1(ff) = 0.
	 endif
	enddo

	d_direct_square(:) = 0.
	do ff=1,number_filters
	 d_direct_square(:) = d_direct_square(:) + M1_M0(:,ff)*M1_M0(:,ff)
	enddo

	M2_M1_square = 0.
	do ff=1,number_filters
	 M2_M1_square = M2_M1_square + M2_M1(ff)*M2_M1(ff)
	enddo
	M2_M1_square_inverse = M2_M1_square**(-1)

	!#IF WE WANT THE k FACTOR!!!
	k(:)=0.
	do ff=1,number_filters
	 k(:) = k(:) - M1_M0(:,ff)*M2_M1(ff)
	enddo
	k=k*M2_M1_square_inverse
	d_parallel_square = k*k*M2_M1_square

	d_perpendicular_square = d_direct_square - d_parallel_square
	Ebv_vector = 2*(d_parallel_square*M2_M1_square_inverse)**0.5	
	 
        do ff=1,number_filters  !Now we want to compute the reddening of the M0 points, to get the M0' ones.
	 M0_prime(:,ff) = M0(:,ff) + A_lambda_filters_selected(ff) * Rv*Ebv_vector(:)
        enddo
	
	!Now we want to get the distances between M0' and the M1 point (i.e. the coordinates of M0' in the system of M1).  
        do ff=1,number_filters
	 M1_M0_prime(:,ff) = abs(M0_prime(:,ff) - M1(ff))
         if (gaps_in_filter(list,ff)==0) then
	  M1_M0_prime(:,ff) = 0.
	 endif
        enddo



	!# -----------------------------
	!# Selection of models in the OB
	!# -----------------------------
	!#Traditional OB
	counting(list) = 0
	OB_size(list) = 3
	sigma_factor = 3	!Indicates the size of the OB, in sigma units
	ID_models(:) = 0

	!If we want to explore all the models of the grid
	!counting(list) = N_bin_all_nodes
	!OB_size(list) = 7 !Means all nodes
	!If we want to explore an OB in place of all the models of the grid
	do 
	 do ii = 1,N_bin_all_nodes
	  in_box=0 !1
	 
	  !HERE PUT THE K FACTOR SELECTION,TO EXCLUDE MODELS OUT OF E(B-V)=[0,1], out of [0,k_max]
	  if (k(ii) >= k_lower .and. k(ii)<k_higher) then		!k is in unit of E(B-V) = 2! 
	   !do ff = 1,number_filters
	   ! if (gaps_in_filter(list,ff)==1) then 
	   !  if (M1_M0_prime(ii,ff)<sigma_factor*sigma_filter(list,ff)) then	!General case
           !  !if (abs(M1_M0(ii,ff))<sigma_factor*sigma_filter(list,ff)) then	!IF WE DO NOT WANT EXTINCTION!!!!!!!
	   !   in_box=1
	   !  else
	   !   in_box=0
	   !   exit
 	   !  endif
	   ! else
           !  in_box=1
	   ! endif
	   !enddo
    in_box=1  !!!!!!!!!!!!INFINITE OB SIZE!!!!!!!!!!!!!!!! (means that we remove the OB, take all models of all nodes!!!)
	  endif
	
	  if (k_lower == 0.) then
	   if (k(ii) < 0) then	!In that case we look if M0 ITSELF is in OB
	    do ff = 1,number_filters
	     if (gaps_in_filter(list,ff)==1) then
	      if (abs(M1_M0(ii,ff))<sigma_factor*sigma_filter(list,ff)) then		!General case
	       in_box=1
	      else
	       in_box=0
	       exit
 	      endif
	     else
              in_box=1
	     endif
	    enddo
	   endif
	  endif

	
	  if (in_box==1) then
	     counting(list) = counting(list) + 1
	     ID_models(counting(list)) = ii
	     Grid_completed_selection(counting(list),:) = Grid_completed(ii,:)
	     if (k(ii) >= 0 .and. k(ii)<k_higher) then	!Classical case, where k>=0.
	     !if (k(ii) >= 0) then	!Case for inspection of the whole grid 
	      !Grid_completed_selection(counting(list),3) = Ebv_vector(ii) ! 2*(d_parallel_square(ii)*M2_M1_square_inverse)**0.5
	      Grid_completed_selection(counting(list),3) = Ebv_vector(ii) + Ebv_foreground !ADDING FOREGROUND EXTINCTION!!!

	     else                                       !Case where k<0. --> M0 is in OB, no extinction.
	      !Grid_completed_selection(counting(list),3) = 0. 
	      Grid_completed_selection(counting(list),3) = 0. + Ebv_foreground !ADDING FOREGROUND EXTINCTION!!!
	     endif
	     Grid_completed_selection(counting(list),4+number_filters+1) = d_perpendicular_square(ii) !**0.5
	     Grid_completed_selection(counting(list),4+number_filters+2) = d_perpendicular_square(ii)**(-1)
	     Grid_completed_selection(counting(list),4+number_filters+2+1:4+number_filters+2+number_filters) = &
	& + M1_M0_prime(ii,1:number_filters) !Here add the M1-M0_prime table in Grid_completed_selection                   !IS THIS CORRECT??????????????????!



	  endif
	 enddo
	 if (counting(list) > 1000) then
	  exit
	 elseif (sigma_factor >= 6) then	!We do not explore further than 6 sigmas from observation
	  exit
	 else
	  sigma_factor = sigma_factor + 1
	  OB_size(list) = OB_size(list) +1
	 endif
	enddo


	counting_inverse = (counting(list)*1.)**(-1)
	write(*,*) 'Cluster  ', list, 'ID =', Cluster_ID(list), ' #in OB ',counting(list), ' OB_size = ', OB_size(list)
	!Here i could write a file containing all the models in the OB


	!# ------------------------------------
	!# Building the probabilities 
	!# ------------------------------------
	!The best thing to do now would be to create nodes files. I could launch separately the code on different Z. 
	!Then, with another code, i will read the results of all the Z (or only one, if desired) and produce final results
	!So, all what we have to do here is simply to put the results in a table age/mass/Ebv with their associated results:
	!sigma1, sigma2, sigma3(?), proba SG(?), proba CG and proba d2 of the whole node!

	!Writing of final solution 
	node_solution = 0.
	Proba_max_CG(list) = 0.
	Proba_max_d2(list) = 0.
        !histo_age(:) = 0.	
        !histo_mass(:)= 0.	
        !histo_Ebv(:) = 0.	
        !histo_Z(:)   = 0.
	!zz_bis = 14-nint(0.5*zz)



	do ii = 1,counting(list)
         !age_int = nint(Grid_completed_selection(ii,1)*100)
         !mass_int = nint(Grid_completed_selection(ii,2)*100)
	 aa  = nint(20*(Grid_completed_selection(ii,1)-6.60)+1)
	 mm  = nint(20*(Grid_completed_selection(ii,2)-2.00)+1)
	 !Ext = nint(50*Grid_completed_selection(ii,3)) + 1
	 Ext = nint(100*Grid_completed_selection(ii,3)) + 1      !HERE BE CAREFULL!! IS IT 100 OR 50?????



!	 proba_model=0.
!	 do ff = 1,number_filters
!	  if (gaps_in_filter(list,ff)==1) then
!	   proba_model = proba_model &
!	& + ( Grid_completed_selection(ii,4+number_filters+2+ff) * sigma_filter_inverse(list,ff) )**2	!General case
!         endif
!	 enddo
!	 node_solution(aa,mm,Ext,1) = node_solution(aa,mm,Ext,1) &
!	& + exp(-0.5*proba_model) * counting_inverse	 

	 !ALTERNATIVE CASE WITH PART 1/sqrt(2pi)sigma before proba (sqrt(2pi) not useful)  In fact it gives (almost) the same results as here above (not the same probabilities, but same parameters)
	 proba_model_vector(:)=0.
	 proba_model=1.
	 do ff = 1,number_filters
	  if (gaps_in_filter(list,ff)==1) then
	   proba_model_vector(ff) = &
	& + ( Grid_completed_selection(ii,4+number_filters+2+ff) * sigma_filter_inverse(list,ff) )**2	!General case
	   proba_model_vector(ff) = sigma_filter_inverse(list,ff) * exp(-0.5*proba_model_vector(ff))
	   proba_model=proba_model*proba_model_vector(ff)
          endif
 	 enddo
	 node_solution(aa,mm,Ext,1) = node_solution(aa,mm,Ext,1) &
	& + proba_model * counting_inverse     				!no Prior
	!& + proba_model * Grid_completed_selection(ii,1)**(-1) * &	!Prior
	!& Grid_completed_selection(ii,2)**(-1) * counting_inverse	!Prior   

	 node_solution(aa,mm,Ext,2) = node_solution(aa,mm,Ext,2) &
	& + Grid_completed_selection(ii,4+number_filters+2) * counting_inverse			!d2

	 if (node_solution(aa,mm,Ext,1) > Proba_max_CG(list) ) then !.and. Ext < 52) then
	  age_max_CG(list) = Grid_completed_selection(ii,1)
	  mass_max_CG(list) = Grid_completed_selection(ii,2)
	  Ebv_max_CG(list) = Grid_completed_selection(ii,3)
	  Proba_max_CG(list) = node_solution(aa,mm,Ext,1)
	 endif
	 !if (node_solution(age_int,mass_int,Ext,2) > Proba_max_d2(list) .and. Ext < 52) then
	 ! age_max_d2(list) = Grid_completed_selection(ii,1)
	 ! mass_max_d2(list) = Grid_completed_selection(ii,2)
	 ! Ebv_max_d2(list) = Grid_completed_selection(ii,3)
	 ! Proba_max_d2(list) = node_solution(age_int,mass_int,Ext,2)
	 !endif


	 !1D histograms!
	 !histo_age(aa)   = histo_age(aa)   + exp(-0.5*proba_model)
	 !histo_mass(mm)  = histo_mass(mm)  + exp(-0.5*proba_model) 
	 !histo_Ebv(Ext)  = histo_Ebv(Ext)  + exp(-0.5*proba_model) 
	 !histo_Z(zz_bis) = histo_Z(zz_bis) + exp(-0.5*proba_model) 

	enddo

	!max_proba_position = MAXLOC(node_solution(:,:,:,1)) !It gives the position (aa,mm,Ext) of the max 
	!write(*,*) max_proba_position
	
	!Option: output the ID and parameters of models located in the OB. Slows down the code... 
	WRITE(jj_char, '(i10)' ) Cluster_ID(list)
	!file_out_cluster_models=trim(file_out_cluster)//'models_files/'
	!OPEN(unit = 41,file=trim(file_out_cluster_models)//'Cluster_'//trim(adjustl(jj_char))//'_Z'//Z_indice_selected//'_V101_models.dat')
	!do ii = 1,counting(list)
        ! if (nint(Grid_completed_selection(ii,1)*100) ==  nint(age_max_CG(list)*100) ) then
        !  if (nint(Grid_completed_selection(ii,2)*100) == nint(mass_max_CG(list)*100) ) then
        !   if (nint(50*Grid_completed_selection(ii,3))  == nint(50*Ebv_max_CG(list)) ) then
	!    write(41,*)ID_models(ii), mod (ID_models(ii), 1000), Grid_completed_selection(ii,1), &
	!& Grid_completed_selection(ii,2), Grid_completed_selection(ii,3)
	!   endif
	!  endif
	! endif
	!enddo
 	!close(41)

	!# ----------------------------------------------
	!Writing of node file for each individual cluster  (takes a bit more time)
	!# ----------------------------------------------
	file_out_cluster_node=trim(file_out_cluster)//'node_files/'
	OPEN(unit = 40,file=trim(file_out_cluster_node)//'Cluster_'//trim(adjustl(jj_char))//'_node_Z'//Z_indice_selected//'.dat')
	write(40,*)'#aa  mm  Ext Proba'
	!OPEN(unit = 40,file=trim(file_out_cluster)//'_'//trim(adjustl(jj_char))//'_node_Z'//Z_indice_selected//'.bin', form='UNFORMATTED')
 	do aa = 1,71 !650, 1200
 	 do mm = 1,101 !200, 600
 	  do Ext = 1,121 !Ext_limit1,Ext_limit2
	   if (node_solution(aa,mm,Ext,1) > 0.) then
 	    write(40,'(i4,i4,i4,E12.3,E12.3,i4)') &
 	    !write(40) &
	& aa,mm, Ext, &
	& node_solution(aa,mm,Ext,1) !, &
	!& node_solution(age_int,mass_int,Ext,2), &
	!& nint(node_solution(age_int,mass_int,Ext,3))
	   endif
 	  enddo
 	 enddo
 	enddo
 	close(40)

	!# --------------------------------------------------
	!Writing of 1D histo file for each individual cluster 
	!# --------------------------------------------------
	!file_out_cluster_histo=trim(file_out_cluster)//'histo_files/'
	!OPEN(unit = 42,file=trim(file_out_cluster_histo)//'Cluster_'//trim(adjustl(jj_char))//'_histo1D_Z'//Z_indice_selected//'_V101.dat')
 	!do Ext = 1,401
	! write(42,'(E12.4,E12.4,E12.4,E12.4)') histo_age(Ext), histo_mass(Ext), histo_Ebv(Ext), histo_Z(Ext)
 	!enddo
 	!close(42)

 ENDDO


 !# ----------------------------------------------
 !Writing of solutions for all clusters (for 1Z)
 !# ----------------------------------------------
 file_out_cluster_f90 = trim(file_out_cluster)
 file_out_cluster = trim(file_out_cluster)//'file_out_cluster'
 file_out_cluster_f90 = trim(file_out_cluster_f90)//'All_clusters_parameters_results_f90_V100'
 file_out_cluster_f90 = trim(file_out_cluster_f90)//'_Z'//Z_indice_selected
 open(unit=41,file=file_out_cluster_f90)
 !write(41,*)'#   age3 mas3 Ebv3  Proba3       age4 mas4 Ebv4  Proba4         count    OB'
 write(41,*)'# ID  age3 mas3 Ebv3  Proba3         count    OB'
 do list = number_begin,number_end
  !write(41,'(i5,F7.2,F5.2,F5.2,E12.4,F7.2,F5.2,F5.2,E12.4,i10,i5)') &
  write(41,'(i5,F7.2,F5.2,F5.2,E12.4,i10,i5)') &
	& Cluster_ID(list), &
	& age_max_CG(list),mass_max_CG(list),Ebv_max_CG(list),Proba_max_CG(list), &
	!& age_max_d2(list),mass_max_d2(list),Ebv_max_d2(list),Proba_max_d2(list), &
	& counting(list), OB_size(list)
 enddo
 close(41)


 CALL system('date')
 WRITE(*,*)' Computation completed' 
 stop

 deALLOCATE(Cluster_ID)
 deALLOCATE(Grid_completed)
 deALLOCATE(M0)
 deALLOCATE(M1_input)
 deALLOCATE(M1,M2)
 deALLOCATE(M1_M0,M1_M0_prime)
 deALLOCATE(M2_M1)
 deALLOCATE(d_direct_square)
 deALLOCATE(d_parallel_square)
 deALLOCATE(d_perpendicular_square)
 deALLOCATE(d_perpendicular)
 deALLOCATE(Ebv_vector)
 deALLOCATE(Grid_completed_selection)
 deALLOCATE(age_max_CG,mass_max_CG)
 deALLOCATE(Ebv_max_CG,Proba_max_CG)
 deALLOCATE(age_max_d2,mass_max_d2)
 deALLOCATE(Ebv_max_d2,Proba_max_d2)
 deALLOCATE(k,counting)
 deALLOCATE(OB_size)
 deALLOCATE(sigma_filter)




 end program
