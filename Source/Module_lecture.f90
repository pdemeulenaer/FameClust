! Module python : module de lecture des fichiers .sum de SimClust 

Module Lecture_module_fortran
	
 contains 


!--------------------------------------------------------------------------------------
!SUBROUTINE Lecture_sumfile_colours(file_name, U_V, B_V, V_R, V_I, V_J, V_H, V_K)	!Grid SimClust
SUBROUTINE Lecture_sumfile_colours(file_name, U_V, B_V, V_R, V_I)			!Grid gCMD
!--------------------------------------------------------------------------------------

 implicit none

 real(4), dimension(1:10000) :: U, B, V, R, I !, J, H, K                  		!Only if SimClust
 real(4), dimension(0:10000) :: UUU, BBB, VVV, RRR, III !, JJJ, HHH, KKK 		!Only if SimClust
 real(4), dimension(1:10000), intent(out) :: U_V, B_V, V_R, V_I !, V_J, V_H, V_K  	!Only if SimClust
 character(len=50), intent(in) ::  file_name
 character(len=50) :: clsNr,clsMlum_ini,clsMlum_act,clsMbol,clsMu,clsMb,clsMv,clsMr,clsMi,clsMj,clsMh,clsMk 
 character(len=50) :: aaaa,bbbb,cccc,dddd
 integer :: ii,kk


 open(unit=19,file=file_name)
 !read(19,*)clsNr,clsMlum_ini,clsMlum_act,clsMbol,clsMu,clsMb,clsMv,clsMr,clsMi,clsMj,clsMh,clsMk  	!Grid SimClust
 do ii = 1,10000
	!Grid SimClust
	!read(19,*)aaaa,bbbb,cccc,dddd,UUU(ii),BBB(ii),VVV(ii),RRR(ii),III(ii),JJJ(ii),HHH(ii),KKK(ii)
	!Grid gCMD
	read(19,*)UUU(ii),BBB(ii),VVV(ii),RRR(ii),III(ii)
 enddo
 close(19)

 U(1:10000) = UUU(1:10000)
 B(1:10000) = BBB(1:10000)
 V(1:10000) = VVV(1:10000)
 R(1:10000) = RRR(1:10000)
 I(1:10000) = III(1:10000)
 !J(1:10000) = JJJ(1:10000)  	!Only if SimClust
 !H(1:10000) = HHH(1:10000)  	!Only if SimClust
 !K(1:10000) = KKK(1:10000)  	!Only if SimClust

 U_V(:) = U(:) - V(:)
 B_V(:) = B(:) - V(:)
 V_R(:) = V(:) - R(:)
 V_I(:) = V(:) - I(:)
 !V_J(:) = V(:) - J(:)	  	!Only if SimClust
 !V_H(:) = V(:) - H(:)  	!Only if SimClust
 !V_K(:) = V(:) - K(:)  	!Only if SimClust

 return 

 END SUBROUTINE


!--------------------------------------------------------------------------------------
!SUBROUTINE Lecture_sumfile_magnitudes(file_name, U, B, V, R, I, J, H, K)	!Grid SimClust
SUBROUTINE Lecture_sumfile_magnitudes(file_name, U, B, V, R, I)			!Grid gCMD
!Lecture_sumfile_magnitudes_2.0(file_name,U_read_index,B_read_index,V_read_index,R_read_index,I_read_index,height_bin, N_bin)
!--------------------------------------------------------------------------------------

 implicit none

 real(4), dimension(1:1000), intent(out) :: U, B, V, R, I !, J, H, K    !Only if SimClust
 real(4), dimension(0:1000) :: UUU, BBB, VVV, RRR, III !, JJJ, HHH, KKK !Only if SimClust
 real(4), dimension(1:1000) :: U_V, B_V, V_R, V_I !, V_J, V_H, V_K      !Only if SimClust
 character(len=50), intent(in) ::  file_name
 !character(len=50) :: clsNr,clsMlum_ini,clsMlum_act,clsMbol,clsMu,clsMb,clsMv,clsMr,clsMi,clsMj,clsMh,clsMk 
 character(len=50) :: aaaa,bbbb,cccc,dddd
 integer :: ii,kk

 !If ASCII files
 !open(unit=19,file=file_name)
 !read(19,*)clsNr,clsMlum_ini,clsMlum_act,clsMbol,clsMu,clsMb,clsMv,clsMr,clsMi,clsMj,clsMh,clsMk  
 !do ii = 1,10000
 !	read(19,*)aaaa,bbbb,cccc,dddd,U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii)
 !	!write(*,*)aaaa,bbbb,cccc,dddd,UUU(ii),BBB(ii),VVV(ii),RRR(ii),III(ii),JJJ(ii),HHH(ii),KKK(ii)
 !	!read(19,*)U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii)    	!Grid SimClust
 !	read(19,*)U(ii),B(ii),V(ii),R(ii),I(ii) !,J(ii),H(ii),K(ii)   	!Grid gCMD
 !enddo
 !close(19)

 
 !If binary file!!! (only magnitudes written)
 open(unit=8,file=trim(file_name)//'.bin', form='UNFORMATTED')
 do ii = 1,1000
 	!read(8)U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii)    	!Grid SimClust
  	read(8)U(ii),B(ii),V(ii),R(ii),I(ii) !,J(ii),H(ii),K(ii)   	!Grid gCMD
	!write(*,*)U(ii),B(ii),V(ii),R(ii),I(ii)
 enddo
 close(8)

 return 

 END SUBROUTINE



!--------------------------------------------------------------------------------------
!SUBROUTINE Lecture_sumfile_magnitudes(file_name, U, B, V, R, I, J, H, K)	!Grid SimClust
!SUBROUTINE Lecture_sumfile_magnitudes(file_name, U, B, V, R, I)		!Grid gCMDV 1.0
SUBROUTINE Lecture_sumfile_magnitudes_20(file_name,U,B,V,R,I,height_bin,N_bin) !Grid gCMDV 2.0
!--------------------------------------------------------------------------------------

 implicit none

 !integer, dimension(1:10000), intent(out) :: U, B, V, R, I !, J, H, K    !Only if SimClust
 !integer, dimension(0:10000) :: UUU, BBB, VVV, RRR, III !, JJJ, HHH, KKK !Only if SimClust
 !integer, dimension(1:10000) :: U_V, B_V, V_R, V_I !, V_J, V_H, V_K      !Only if SimClust
 real(4), dimension(1:10000), intent(out) :: U, B, V, R, I !, J, H, K    !Only if SimClust
 real(4), dimension(0:10000) :: UUU, BBB, VVV, RRR, III !, JJJ, HHH, KKK !Only if SimClust
 real(4), dimension(1:10000) :: U_V, B_V, V_R, V_I !, V_J, V_H, V_K      !Only if SimClust
 integer, intent (out) :: N_bin
 integer, dimension(1:10000) :: height_bin(:)
 real(4) :: minimum, BinSize
 character(len=50), intent(in) ::  file_name
 character(len=50) :: clsNr,clsMlum_ini,clsMlum_act,clsMbol,clsMu,clsMb,clsMv,clsMr,clsMi,clsMj,clsMh,clsMk 
 character(len=50) :: aaaa,bbbb,cccc,dddd
 integer :: ii,kk


 !If ASCII files
 !open(unit=19,file=file_name)
 !read(19,*)minimum, BinSize, N_bin  
 !do ii = 1,N_bin   !10000
 ! read(19,*)U(ii),B(ii),V(ii),R(ii),I(ii), height_bin(ii)  	!Grid gCMD
 !enddo
 !close(19)

 !If ASCII files
 open(unit=19,file=file_name, form='UNFORMATTED')
 read(19)minimum, BinSize, N_bin  
 do ii = 1,N_bin   !10000
  read(19)U(ii),B(ii),V(ii),R(ii),I(ii), height_bin(ii)  	!Grid gCMD
 enddo
 close(19)

 return 

 END SUBROUTINE




!--------------------------------------------------------------------------------------
SUBROUTINE Lecture_sumfile_magnitudes_AllBands(file_name,FUV,NUV,U,B,V,R,I,J,H,K,u_SDSS,g_SDSS,r_SDSS, &
	& i_SDSS,z_SDSS,IRAC36,IRAC45,IRAC58,IRAC80,mips24,mips70,mips160)
!--------------------------------------------------------------------------------------

 implicit none
 real(4), dimension(1:1000), intent(out) :: FUV,NUV,U,B,V,R,I,J,H,K,u_SDSS,g_SDSS,r_SDSS, &
	& i_SDSS,z_SDSS,IRAC36,IRAC45,IRAC58,IRAC80,mips24,mips70,mips160
 character(len=50), intent(in) ::  file_name
 integer :: ii,kk

 !If ASCII files
 !open(unit=19,file=file_name)
 !read(19,*)clsNr,clsMlum_ini,clsMlum_act,clsMbol,clsMu,clsMb,clsMv,clsMr,clsMi,clsMj,clsMh,clsMk  
 !do ii = 1,10000
 !	!read(19,*)U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii)    	!Grid SimClust
 !	read(19,*)U(ii),B(ii),V(ii),R(ii),I(ii) !,J(ii),H(ii),K(ii)   	!Grid gCMD
 !enddo
 !close(19)

 
 !If binary file!!! (only magnitudes written)
 open(unit=8,file=trim(file_name)//'.bin', form='UNFORMATTED')
 do ii = 1,1000
  	read(8)FUV(ii),NUV(ii),U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii), &
	& u_SDSS(ii),g_SDSS(ii),r_SDSS(ii),i_SDSS(ii),z_SDSS(ii),IRAC36(ii),IRAC45(ii),IRAC58(ii), &
	& IRAC80(ii),mips24(ii),mips70(ii),mips160(ii)
	!write(*,*)U(ii)
	!pause
 enddo
 close(8)

 return 

 END SUBROUTINE




! NEW VERSION, READS GRID CONTAINING ONLY FUV NUV UBVRIJHK WFC3 CFHT SDSS 2MASS
! (but in binary files, i leave 0 for other filters, to keep the same structure for FameClust)
!--------------------------------------------------------------------------------------
SUBROUTINE Lecture_sumfile_magnitudes_AllBands_ASCII_HRS_with_ACS(file_name,FUV,NUV,U,B,V,R,I,J,H,K,&
	& F275W,F336W,F475W,F814W,F110W,F160W,&
	& u_CFHT,g_CFHT,r_CFHT,i_CFHT,z_CFHT,&
	& u_SDSS,g_SDSS,r_SDSS,i_SDSS,z_SDSS,J_2MASS,H_2MASS,Ks_2MASS,number_models_per_node)
!--------------------------------------------------------------------------------------

 implicit none
 integer, intent(in) :: number_models_per_node
 real(4), dimension(1:number_models_per_node), intent(out) :: FUV,NUV,U,B,V,R,I,J,H,K,&
	& F275W,F336W,F475W,F814W,F110W,F160W,&
	& u_CFHT,g_CFHT,r_CFHT,i_CFHT,z_CFHT,&
	& u_SDSS,g_SDSS,r_SDSS,i_SDSS,z_SDSS,J_2MASS,H_2MASS,Ks_2MASS
 real(4), dimension(1:number_models_per_node) :: age, mass
 character(len=50), intent(in) ::  file_name
 integer :: ii,kk

 !If ASCII files
 open(unit=19,file=file_name)
 read(19,*) !To pass the header line
 do ii = 1,number_models_per_node !1000
  	read(19,*)age(ii),mass(ii),FUV(ii),NUV(ii),U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii), &
	& F275W(ii),F336W(ii),F475W(ii),F814W(ii),F110W(ii),F160W(ii), &
	& u_CFHT(ii),g_CFHT(ii),r_CFHT(ii),i_CFHT(ii),z_CFHT(ii), &
	& u_SDSS(ii),g_SDSS(ii),r_SDSS(ii),i_SDSS(ii),z_SDSS(ii),J_2MASS(ii),H_2MASS(ii),Ks_2MASS(ii)

	!write(*,*)FUV(ii),NUV(ii),U(ii),B(ii)
 enddo
 close(19)

 return 
 END SUBROUTINE



! [DEPRECATED]
!--------------------------------------------------------------------------------------
SUBROUTINE Lecture_UBVRI_ugriz_GRID_OLD(file_name,U,B,V,R,I,u_SDSS,g_SDSS,r_SDSS, &
	& i_SDSS,z_SDSS)
!--------------------------------------------------------------------------------------

 implicit none
 real(4), dimension(1:1000), intent(out) :: U,B,V,R,I,u_SDSS,g_SDSS,r_SDSS, &
	& i_SDSS,z_SDSS !,IRAC36,IRAC45,IRAC58,IRAC80,mips24,mips70,mips160
 real(4) :: FUV,NUV,J,H,K
 character(len=300), intent(in) :: file_name
 integer :: ii,kk

 
 !If binary file!!! (only magnitudes written)
 open(unit=8,file=trim(file_name)//'.bin', form='UNFORMATTED')
 do ii = 1,1000
  	read(8)FUV,NUV,U(ii),B(ii),V(ii),R(ii),I(ii),J,H,K, &
	& u_SDSS(ii),g_SDSS(ii),r_SDSS(ii),i_SDSS(ii),z_SDSS(ii)!,IRAC36(ii),IRAC45(ii),IRAC58(ii), &
	!& IRAC80(ii),mips24(ii),mips70(ii),mips160(ii)

  	!read(8)FUV(ii),NUV(ii),U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii), &
	!& u_SDSS(ii),g_SDSS(ii),r_SDSS(ii),i_SDSS(ii),z_SDSS(ii),IRAC36(ii),IRAC45(ii),IRAC58(ii), &
	!& IRAC80(ii),mips24(ii),mips70(ii),mips160(ii)
	!write(*,*)U(ii)
	!pause
 enddo
 close(8)

 return 

 END SUBROUTINE


! [DEPRECATED]
!--------------------------------------------------------------------------------------
!SUBROUTINE Lecture_UBVRI_ugriz_GRID(file_name,U,B,V,R,I,u_SDSS,g_SDSS,r_SDSS, &
!	& i_SDSS,z_SDSS)
!SUBROUTINE Lecture_UBVRI_ugriz_GRID(file_name,FUV,NUV,U,B,V,R,I,J,H,K,u_CFHT,g_CFHT,r_CFHT, &
!	& i_CFHT,z_CFHT,u_SDSS,g_SDSS,r_SDSS,i_SDSS,z_SDSS,J_2MASS,H_2MASS,Ks_2MASS,a_BATC,b_BATC, &
!	& c_BATC,d_BATC,e_BATC,f_BATC,g_BATC,h_BATC,i_BATC,j_BATC,k_BATC,m_BATC,n_BATC,o_BATC,p_BATC, &
!	& t_BATC,IRAC36,IRAC45,IRAC58,IRAC80,mips24,mips70,mips160)
SUBROUTINE Lecture_UBVRI_ugriz_GRID(file_name,GRID_read)
!--------------------------------------------------------------------------------------

 implicit none
! real(4), dimension(1:1000), intent(out) :: U,B,V,R,I,u_SDSS,g_SDSS,r_SDSS, &
!	& i_SDSS,z_SDSS !,IRAC36,IRAC45,IRAC58,IRAC80,mips24,mips70,mips160
! real(4), dimension(1:1000), intent(out) :: FUV,NUV,U,B,V,R,I,J,H,K,u_CFHT,g_CFHT,r_CFHT, &
! 	& i_CFHT,z_CFHT,u_SDSS,g_SDSS,r_SDSS,i_SDSS,z_SDSS,J_2MASS,H_2MASS,Ks_2MASS,a_BATC,b_BATC, &
! 	& c_BATC,d_BATC,e_BATC,f_BATC,g_BATC,h_BATC,i_BATC,j_BATC,k_BATC,m_BATC,n_BATC,o_BATC,p_BATC, &
! 	& t_BATC,IRAC36,IRAC45,IRAC58,IRAC80,mips24,mips70,mips160
 real(4), dimension(1:46,1:1000), intent(out) :: GRID_read
 !real(4) :: FUV,NUV,J,H,K
 character(len=300), intent(in) :: file_name
 integer :: ii,kk

 
 !If binary file!!! (only magnitudes written)
 open(unit=8,file=trim(file_name)//'.bin', form='UNFORMATTED')
 do ii = 1,1000
!  	read(8)FUV,NUV,U(ii),B(ii),V(ii),R(ii),I(ii),J,H,K, &
!	& u_SDSS(ii),g_SDSS(ii),r_SDSS(ii),i_SDSS(ii),z_SDSS(ii)!,IRAC36(ii),IRAC45(ii),IRAC58(ii), &
!	!& IRAC80(ii),mips24(ii),mips70(ii),mips160(ii)

  	!read(8)FUV(ii),NUV(ii),U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii), &
	!& u_CFHT(ii),g_CFHT(ii),r_CFHT(ii),i_CFHT(ii),z_CFHT(ii), &
	!& u_SDSS(ii),g_SDSS(ii),r_SDSS(ii),i_SDSS(ii),z_SDSS(ii),J_2MASS(ii),H_2MASS(ii),Ks_2MASS(ii), &
	!& a_BATC(ii),b_BATC(ii),c_BATC(ii),d_BATC(ii),e_BATC(ii),f_BATC(ii),g_BATC(ii),h_BATC(ii), &
	!& i_BATC(ii),j_BATC(ii),k_BATC(ii),m_BATC(ii),n_BATC(ii),o_BATC(ii),p_BATC(ii),t_BATC(ii), &
	!& IRAC36(ii),IRAC45(ii),IRAC58(ii),IRAC80(ii),mips24(ii),mips70(ii),mips160(ii)
  	read(8)GRID_read(1:46,ii)!,GRID_read(2,ii),GRID_read(3,ii),GRID_read(4,ii),GRID_read(5,ii), &
	!& GRID_read(6,ii),GRID_read(7,ii),GRID_read(1,ii),GRID_read(1,ii), &
	!& GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii), &
	!& GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii), &
	!& GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii), &
	!& GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii), &
	!& GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii), &
	!& GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii), &
	!& GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii),GRID_read(1,ii), &
	!& GRID_read(1,ii)

  	!read(8)FUV(ii),NUV(ii),U(ii),B(ii),V(ii),R(ii),I(ii),J(ii),H(ii),K(ii), &
	!& u_SDSS(ii),g_SDSS(ii),r_SDSS(ii),i_SDSS(ii),z_SDSS(ii),IRAC36(ii),IRAC45(ii),IRAC58(ii), &
	!& IRAC80(ii),mips24(ii),mips70(ii),mips160(ii)
	!write(*,*)U(ii)
	!pause
 enddo
 close(8)

 return 

 END SUBROUTINE



! [DEPRECATED]
!--------------------------------------------------------------------------------------
SUBROUTINE Lecture_UBVRI_ugriz_GRID_Jan2013(file_name,GRID_read)
!--------------------------------------------------------------------------------------
 implicit none
 real(4), dimension(1:1000,1:52), intent(out) :: GRID_read
 character(len=300), intent(in) :: file_name
 integer :: ii,kk
 open(unit=8,file=trim(file_name)//'.bin', form='UNFORMATTED')
 do ii = 1,1000
  read(8)GRID_read(ii,1:52)
 enddo
 close(8)
 return 
 END SUBROUTINE


!--------------------------------------------------------------------------------------
SUBROUTINE Lecture_GRID_March2016(file_name_nodes_bin,GRID_READ,number_models_per_node) 
!--------------------------------------------------------------------------------------
 implicit none
 integer, intent(in) :: number_models_per_node
 real(4), dimension(1:number_models_per_node,1:52), intent(out) :: GRID_read
 character(len=300), intent(in) :: file_name_nodes_bin
 integer :: ii,kk
 open(unit=8,file=trim(file_name_nodes_bin)//'.bin', form='UNFORMATTED')
 do ii = 1,number_models_per_node
  read(8)GRID_read(ii,1:52)
 enddo
 close(8)
 return 
 END SUBROUTINE





!--------------------------------------------------------------------------------------
SUBROUTINE Lecture_GMM_ASCII_GRID_March2016(path_GMM_ascii_file,means,covariances,weights,&
	& age_indice,mass_indice,Z_indice,components_number) 
!--------------------------------------------------------------------------------------
 implicit none
 integer, intent(in) :: components_number
 real(4), dimension(1:components_number,1:29), intent(out) :: means
 real(4), dimension(1:components_number,1:29,1:29), intent(out) :: covariances
 real(4), dimension(1:components_number), intent(out) :: weights
 character(len=300), intent(in) :: path_GMM_ascii_file
 character(len=50) , intent(in) :: age_indice,mass_indice,Z_indice
 character(len=300) :: file_name_means, file_name_covariances, file_name_weights
 integer :: ii,jj

 file_name_means = trim(path_GMM_ascii_file)// 'means/Clusters_t' //&
        & trim(age_indice) // &
        & '_M' // trim(mass_indice) //&
        & '_Z' // trim(Z_indice) // '_means'

 file_name_covariances = trim(path_GMM_ascii_file)// 'covariances/Clusters_t' //&
        & trim(age_indice) // &
        & '_M' // trim(mass_indice) //&
        & '_Z' // trim(Z_indice) // '_covariances'

 file_name_weights = trim(path_GMM_ascii_file)// 'weights/Clusters_t' //&
        & trim(age_indice) // &
        & '_M' // trim(mass_indice) //&
        & '_Z' // trim(Z_indice) // '_weights'

 !write(*,*)file_name_means
 !write(*,*)file_name_covariances
 !write(*,*)file_name_weights


 !Reading MEANS files
 open(unit=20,file=trim(file_name_means))
 do ii = 1,components_number
  read(20,*)means(ii,1:29)
 enddo
 close(20)


 !Reading COVARIANCES files
 open(unit=21,file=trim(file_name_covariances))
 do ii = 1,components_number
  do jj = 1,29
   read(21,*)covariances(ii,jj,1:29)
  enddo
  if (ii<components_number) read(21,*)
 enddo
 close(21)


 !Reading WEIGHTS files
 open(unit=22,file=trim(file_name_weights))
 do ii = 1,components_number
  read(22,*)weights(ii)
 enddo
 close(22)


 return 
 END SUBROUTINE








!--------------------------------------------------------------------------------------
SUBROUTINE Lecture_GMM_BIN_GRID_March2016(path_GMM_bin_file,means,covariances,weights,&
	& age_indice,mass_indice,Z_indice,components_number) 
!--------------------------------------------------------------------------------------
 implicit none
 integer, intent(in) :: components_number
 real(4), dimension(1:components_number,1:29), intent(out) :: means
 real(4), dimension(1:components_number,1:29,1:29), intent(out) :: covariances
 real(4), dimension(1:components_number), intent(out) :: weights
 character(len=300), intent(in) :: path_GMM_bin_file
 character(len=50) , intent(in) :: age_indice,mass_indice,Z_indice
 character(len=300) :: file_name_means, file_name_covariances, file_name_weights
 integer :: ii,jj

 file_name_means = trim(path_GMM_bin_file)// 'means/Clusters_t' //&
        & trim(age_indice) // &
        & '_M' // trim(mass_indice) //&
        & '_Z' // trim(Z_indice) // '_means.bin'

 file_name_covariances = trim(path_GMM_bin_file)// 'covariances/Clusters_t' //&
        & trim(age_indice) // &
        & '_M' // trim(mass_indice) //&
        & '_Z' // trim(Z_indice) // '_covariances.bin'

 file_name_weights = trim(path_GMM_bin_file)// 'weights/Clusters_t' //&
        & trim(age_indice) // &
        & '_M' // trim(mass_indice) //&
        & '_Z' // trim(Z_indice) // '_weights.bin'

 !write(*,*)file_name_means
 !write(*,*)file_name_covariances
 !write(*,*)file_name_weights


 !Reading MEANS files
 open(unit=8,file=trim(file_name_means), form='UNFORMATTED')
 do ii = 1,components_number
  read(8)means(ii,1:29)
 enddo
 close(8)


 !Reading COVARIANCES files
 open(unit=9,file=trim(file_name_covariances), form='UNFORMATTED')
 do ii = 1,components_number
  do jj = 1,29
   read(9)covariances(ii,jj,1:29)
  enddo
  if (ii<components_number) read(9)
 enddo
 close(9)


 !Reading WEIGHTS files
 open(unit=10,file=trim(file_name_weights), form='UNFORMATTED')
 do ii = 1,components_number
  read(10)weights(ii)
 enddo
 close(10)


 return 
 END SUBROUTINE








! [DEPRECATED]
!--------------------------------------------------------------------------------------
SUBROUTINE age_mass_Z(a,m,zz, age, mass, Z)
!--------------------------------------------------------------------------------------

 implicit none

 real(4), intent(out) :: age, mass, Z
 integer, intent(in) :: a, m, zz
 integer :: ii,kk


	if (a == 1) age = 7.0	!10 Myr
	if (a == 2) age = 7.5	!30 Myr
	if (a == 3) age = 8.0	!100 Myr
	if (a == 4) age = 8.5	!300 Myr
	if (a == 5) age = 9.0	!1000 Myr
	if (a == 6) age = 10.0	!10000 Myr
	
	if (m == 1) mass = 1000.	!Solar masses
	if (m == 2) mass = 5000.	!Solar masses
	if (m == 3) mass = 10000.	!Solar masses
	if (m == 4) mass = 50000.	!Solar masses
	if (m == 5) mass = 100000.	!Solar masses
	
	if (zz == 1) Z = 0.001	!10 Myr	
	if (zz == 2) Z = 0.004	!30 Myr	
	if (zz == 3) Z = 0.008	!100 Myr	
	if (zz == 4) Z = 0.01	!300 Myr	
	if (zz == 5) Z = 0.02	!1000 Myr
	if (zz == 6) Z = 0.03	!10000 Myr

 return 

 END SUBROUTINE



! [DEPRECATED]
!--------------------------------------------------------------------------------------
SUBROUTINE age_mass_Z_NEW_GRID(a,m,zz,Ext, age, mass, Z, Ebv, age_indice, mass_indice, Z_indice, Ebv_indice)
!--------------------------------------------------------------------------------------

 implicit none

 real(4), intent(out) :: age, mass, Z, Ebv
 integer, intent(in) :: a, m, zz, Ext
 integer :: ii,kk
 character(len=50), intent(out) :: age_indice, mass_indice, Z_indice, Ebv_indice


 ! Old grid, June 2011
!	if (a == 1) then
!	 age  = 7.00
!	 age_indice = '700'	!10 Myr
!	elseif (a == 2) then
!	 age = 7.15
!	 age_indice = '715'	
!	elseif (a == 3) then
!	 age = 7.30
!	 age_indice = '730'
!	elseif (a == 4) then
!	 age = 7.45
!	 age_indice = '745'	
!	elseif (a == 5) then
!	 age = 7.60
!	 age_indice = '760'	
!	elseif (a == 6) then
!	 age = 7.75
!	 age_indice = '775'	
!	elseif (a == 7) then
!	 age = 7.90
!	 age_indice = '790'	
!	elseif (a == 8) then
!	 age = 8.05
!	 age_indice = '805'	
!	elseif (a == 9) then
!	 age = 8.20
!	 age_indice = '820'	
!	elseif (a == 10) then
!	 age = 8.35
!	 age_indice = '835'	
!	elseif (a == 11) then
!	 age = 8.50
!	 age_indice = '850'	
!	elseif (a == 12) then
!	 age = 8.65
!	 age_indice = '865'	
!	elseif (a == 13) then
!	 age = 8.80
!	 age_indice = '880'	
!	elseif (a == 14) then
!	 age = 8.95
!	 age_indice = '895'	
!	elseif (a == 15) then
!	 age = 9.10
!	 age_indice = '910'	
!	elseif (a == 16) then
!	 age = 9.25
!	 age_indice = '925'	
!	elseif (a == 17) then
!	 age = 9.40
!	 age_indice = '940'	
!	elseif (a == 18) then
!	 age = 9.55
!	 age_indice = '955'	
!	elseif (a == 19) then
!	 age = 9.70
!	 age_indice = '970'	
!	elseif (a == 20) then
!	 age = 9.85
!	 age_indice = '985'	
!	elseif (a == 21) then
!	 age = 10.0
!	 age_indice = '1000'	
!	endif

 !New grid: September 2011

	if (a == 1) then
	 age  = 7.00
	 age_indice = '700'	!10 Myr
	elseif (a == 2) then
	 age = 7.10
	 age_indice = '710'	
	elseif (a == 3) then
	 age = 7.20
	 age_indice = '720'
	elseif (a == 4) then
	 age = 7.30
	 age_indice = '730'	
	elseif (a == 5) then
	 age = 7.40
	 age_indice = '740'	
	elseif (a == 6) then
	 age = 7.50
	 age_indice = '750'	
	elseif (a == 7) then
	 age = 7.60
	 age_indice = '760'	
	elseif (a == 8) then
	 age = 7.70
	 age_indice = '770'	
	elseif (a == 9) then
	 age = 7.80
	 age_indice = '780'	
	elseif (a == 10) then
	 age = 7.90
	 age_indice = '790'	
	elseif (a == 11) then
	 age = 8.00
	 age_indice = '800'	
	elseif (a == 12) then
	 age = 8.10
	 age_indice = '810'	
	elseif (a == 13) then
	 age = 8.20
	 age_indice = '820'	
	elseif (a == 14) then
	 age = 8.30
	 age_indice = '830'	
	elseif (a == 15) then
	 age = 8.40
	 age_indice = '840'	
	elseif (a == 16) then
	 age = 8.50
	 age_indice = '850'	
	elseif (a == 17) then
	 age = 8.60
	 age_indice = '860'	
	elseif (a == 18) then
	 age = 8.70
	 age_indice = '870'	
	elseif (a == 19) then
	 age = 8.80
	 age_indice = '880'	
	elseif (a == 20) then
	 age = 8.90
	 age_indice = '890'	
	elseif (a == 21) then
	 age = 9.00
	 age_indice = '900'	
	elseif (a == 22) then
	 age = 9.10
	 age_indice = '910'	
	elseif (a == 23) then
	 age = 9.20
	 age_indice = '920'	
	elseif (a == 24) then
	 age = 9.30
	 age_indice = '930'	
	elseif (a == 25) then
	 age = 9.40
	 age_indice = '940'	
	elseif (a == 26) then
	 age = 9.50
	 age_indice = '950'	
	elseif (a == 27) then
	 age = 9.60
	 age_indice = '960'	
	elseif (a == 28) then
	 age = 9.70
	 age_indice = '970'	
	elseif (a == 29) then
	 age = 9.80
	 age_indice = '980'	
	elseif (a == 30) then
	 age = 9.90
	 age_indice = '990'	
	elseif (a == 31) then
	 age = 10.0
	 age_indice = '1000'
	endif

 !Old grid: June 2011

	!if (m == 1) then 
	! mass = 10**3.00
	! mass_indice = '300'    !Solar masses	
	!elseif (m == 2) then
	! mass = 10**3.25
	! mass_indice = '325'	!Solar masses
	!elseif (m == 3) then
	! mass = 10**3.50
	! mass_indice = '350'	!Solar masses			
	!elseif (m == 4) then
	! mass = 10**3.75
	! mass_indice = '375'	!Solar masses
	!elseif (m == 5) then
	! mass = 10**4.00
	! mass_indice = '400'	!Solar masses
	!elseif (m == 6) then
	! mass = 10**4.25
	! mass_indice = '425'	!Solar masses
	!elseif (m == 7) then
	! mass = 10**4.50
	! mass_indice = '450'	!Solar masses
	!elseif (m == 8) then
	! mass = 10**4.75
	! mass_indice = '475'	!Solar masses
	!elseif (m == 9) then
	! mass = 10**5.00
	! mass_indice = '500'	!Solar masses
	!endif

 !New grid: September 2011

	if (m == 1) then 
	 mass = 10**3.00
	 mass_indice = '300'    !Solar masses	
	elseif (m == 2) then
	 mass = 10**3.10
	 mass_indice = '310'	!Solar masses
	elseif (m == 3) then
	 mass = 10**3.20
	 mass_indice = '320'	!Solar masses			
	elseif (m == 4) then
	 mass = 10**3.30
	 mass_indice = '330'	!Solar masses
	elseif (m == 5) then
	 mass = 10**3.40
	 mass_indice = '340'	!Solar masses
	elseif (m == 6) then
	 mass = 10**3.50
	 mass_indice = '350'	!Solar masses
	elseif (m == 7) then
	 mass = 10**3.60
	 mass_indice = '360'	!Solar masses
	elseif (m == 8) then
	 mass = 10**3.70
	 mass_indice = '370'	!Solar masses
	elseif (m == 9) then
	 mass = 10**3.80
	 mass_indice = '380'	!Solar masses
	elseif (m == 10) then
	 mass = 10**3.90
	 mass_indice = '390'	!Solar masses
	elseif (m == 11) then
	 mass = 10**4.00
	 mass_indice = '400'	!Solar masses
	elseif (m == 12) then
	 mass = 10**4.10
	 mass_indice = '410'	!Solar masses
	elseif (m == 13) then
	 mass = 10**4.20
	 mass_indice = '420'	!Solar masses
	elseif (m == 14) then
	 mass = 10**4.30
	 mass_indice = '430'	!Solar masses
	elseif (m == 15) then
	 mass = 10**4.40
	 mass_indice = '440'	!Solar masses
	elseif (m == 16) then
	 mass = 10**4.50
	 mass_indice = '450'	!Solar masses
	elseif (m == 17) then
	 mass = 10**4.60
	 mass_indice = '460'	!Solar masses
	elseif (m == 18) then
	 mass = 10**4.70
	 mass_indice = '470'	!Solar masses
	elseif (m == 19) then
	 mass = 10**4.80
	 mass_indice = '480'	!Solar masses
	elseif (m == 20) then
	 mass = 10**4.90
	 mass_indice = '490'	!Solar masses
	elseif (m == 21) then
	 mass = 10**5.00
	 mass_indice = '500'	!Solar masses
	endif


	! a tester : [-2.3, -1.9, -1.6, -1.3, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, +0.2]
	!if (zz == 1) then	
	! Z = 0.00010
	! Z_indice = 'm23'	! [M/H] = -2.3
	!elseif (zz == 2) then
	! Z = 0.00025
	! Z_indice = 'm19'	! [M/H] = -1.9
	!elseif (zz == 3) then
	! Z = 0.00050
	! Z_indice = 'm16'	! [M/H] = -1.6
	!elseif (zz == 4) then
	! Z = 0.00100
	! Z_indice = 'm13'	! [M/H] = -1.3
	!elseif (zz == 5) then
	! Z = 0.00200
	! Z_indice = 'm10'	! [M/H] = -1.0
	!elseif (zz == 6) then
	! Z = 0.00317
	! Z_indice = 'm08'	! [M/H] = -0.8
	!elseif (zz == 7) then
	! Z = 0.00502
	! Z_indice = 'm06'	! [M/H] = -0.6
	!elseif (zz == 8) then
	! Z = 0.00796
	! Z_indice = 'm04'	! [M/H] = -0.4
	!elseif (zz == 9) then
	! Z = 0.01262
	! Z_indice = 'm02'	! [M/H] = -0.2
	!elseif (zz == 10) then
	! Z = 0.02000
	! Z_indice = 'n00'	! [M/H] = +0.0
	!elseif (zz == 11) then
	! Z = 0.03000
	! Z_indice = 'p02'	! [M/H] = +0.2 ABOUT, not exactly
	!endif

	!September 2011
	! a tester :[-2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, +0.2]
	!if (zz == 1) then	
	! Z =  0.00012619
	! Z_indice = 'm22'	! [M/H] = -2.2
	!elseif (zz == 2) then
	! Z = 0.00020000
	! Z_indice = 'm20'	! [M/H] = -2.0
	!elseif (zz == 3) then
	! Z = 0.00031698
	! Z_indice = 'm18'	! [M/H] = -1.8
	!elseif (zz == 4) then
	! Z = 0.00050
	! Z_indice = 'm16'	! [M/H] = -1.6
	!elseif (zz == 5) then
	! Z = 0.00079621
	! Z_indice = 'm14'	! [M/H] = -1.4
	!elseif (zz == 6) then
	! Z = 0.0012619
	! Z_indice = 'm12'	! [M/H] = -1.2
	!elseif (zz == 7) then
	! Z = 0.00200
	! Z_indice = 'm10'	! [M/H] = -1.0
	!elseif (zz == 8) then
	! Z = 0.00317
	! Z_indice = 'm08'	! [M/H] = -0.8
	!elseif (zz == 9) then
	! Z = 0.00502
	! Z_indice = 'm06'	! [M/H] = -0.6
	!elseif (zz == 10) then
	! Z = 0.00796
	! Z_indice = 'm04'	! [M/H] = -0.4
	!elseif (zz == 11) then
	! Z = 0.01262
	! Z_indice = 'm02'	! [M/H] = -0.2
	!elseif (zz == 12) then
	! Z = 0.02000
	! Z_indice = 'n00'	! [M/H] = +0.0
	!elseif (zz == 13) then
	! Z = 0.03000
	! Z_indice = 'p02'	! [M/H] = +0.2 ABOUT, not exactly
	!endif


	!October 2011
	if (zz == 1) then	
	 Z = 0.00010
	 Z_indice = 'm23'	! [M/H] = -2.3
	elseif (zz == 2) then
	 Z =  0.00012619
	 Z_indice = 'm22'	! [M/H] = -2.2
        elseif (zz == 3) then
         Z = 0.00015887
         Z_indice = 'm21'	! [M/H] = -2.1
	elseif (zz == 4) then
	 Z = 0.00020000
	 Z_indice = 'm20'	! [M/H] = -2.0
	elseif (zz == 5) then
	 Z = 0.00025
	 Z_indice = 'm19'	! [M/H] = -1.9
	elseif (zz == 6) then
	 Z = 0.00031698
	 Z_indice = 'm18'	! [M/H] = -1.8
        elseif (zz == 7) then
         Z = 0.00039905
         Z_indice = 'm17'	! [M/H] = -1.7
	elseif (zz == 8) then
	 Z = 0.00050
	 Z_indice = 'm16'	! [M/H] = -1.6
        elseif (zz == 9) then 
         Z = 0.00063246
         Z_indice = 'm15'	! [M/H] = -1.5
	elseif (zz == 10) then
	 Z = 0.00079621
	 Z_indice = 'm14'	! [M/H] = -1.4
	elseif (zz == 11) then
	 Z = 0.00100
	 Z_indice = 'm13'	! [M/H] = -1.3
	elseif (zz == 12) then
	 Z = 0.0012619
	 Z_indice = 'm12'	! [M/H] = -1.2
        elseif (zz == 13) then
         Z = 0.0015887
         Z_indice = 'm11'	! [M/H] = -1.1
	elseif (zz == 14) then
	 Z = 0.00200
	 Z_indice = 'm10'	! [M/H] = -1.0
        elseif (zz == 15) then
         Z = 0.0025179 
         Z_indice = 'm09'	! [M/H] = -0.9
	elseif (zz == 16) then
	 Z = 0.00317
	 Z_indice = 'm08'	! [M/H] = -0.8
        elseif (zz == 17) then
         Z = 0.0039905
         Z_indice = 'm07'	! [M/H] = -0.7
	elseif (zz == 18) then
	 Z = 0.00502
	 Z_indice = 'm06'	! [M/H] = -0.6
        elseif (zz == 19) then
         Z = 0.0063246
         Z_indice = 'm05'	! [M/H] = -0.5
	elseif (zz == 20) then
	 Z = 0.00796
	 Z_indice = 'm04'	! [M/H] = -0.4
        elseif (zz == 21) then
         Z = 0.010024
         Z_indice = 'm03'	! [M/H] = -0.3
	elseif (zz == 22) then
	 Z = 0.01262
	 Z_indice = 'm02'	! [M/H] = -0.2
        elseif (zz == 23) then 
         Z = 0.015887
         Z_indice = 'm01'	! [M/H] = -0.1
	elseif (zz == 24) then
	 Z = 0.02000
	 Z_indice = 'n00'	! [M/H] = +0.0
        elseif (zz == 25) then
         Z = 0.025179
         Z_indice = 'p01'	! [M/H] = +0.1
	elseif (zz == 26) then
	 Z = 0.03000
	 Z_indice = 'p02'	! [M/H] = +0.2 ABOUT, not exactly
	endif




 	!Extinction 
        if  (Ext == 1) then
         Ebv = 0.0
         Ebv_indice = '00'
        elseif  (Ext == 2) then
         Ebv = 0.1
         Ebv_indice = '01'
        elseif  (Ext == 3) then
         Ebv = 0.2
         Ebv_indice = '02'
        elseif  (Ext == 4) then
         Ebv = 0.3
         Ebv_indice = '03'
        elseif  (Ext == 5) then
         Ebv = 0.4
         Ebv_indice = '04'
        elseif  (Ext == 6) then
         Ebv = 0.5
         Ebv_indice = '05'
        elseif  (Ext == 7) then
         Ebv = 0.6
         Ebv_indice = '06'
        elseif  (Ext == 8) then
         Ebv = 0.7
         Ebv_indice = '07'
        elseif  (Ext == 9) then
         Ebv = 0.8
         Ebv_indice = '08'
        elseif  (Ext == 10) then
         Ebv = 0.9
         Ebv_indice = '09'
        elseif  (Ext == 11) then
         Ebv = 1.0
         Ebv_indice = '10'
        elseif  (Ext == 12) then
         Ebv = 1.1
         Ebv_indice = '11'
        elseif  (Ext == 13) then
         Ebv = 1.2
         Ebv_indice = '12'
        elseif  (Ext == 14) then
         Ebv = 1.3
         Ebv_indice = '13'
        elseif  (Ext == 15) then
         Ebv = 1.4
         Ebv_indice = '14'
        elseif  (Ext == 16) then
         Ebv = 1.5
         Ebv_indice = '15'
        elseif  (Ext == 17) then
         Ebv = 1.6
         Ebv_indice = '16'
        elseif  (Ext == 18) then
         Ebv = 1.7
         Ebv_indice = '17'
        elseif  (Ext == 19) then
         Ebv = 1.8
         Ebv_indice = '18'
        elseif  (Ext == 20) then
         Ebv = 1.9
         Ebv_indice = '19'
        elseif  (Ext == 21) then
         Ebv = 2.0
         Ebv_indice = '20'
        endif
 return


 END SUBROUTINE


! [DEPRECATED]
!--------------------------------------------------------------------------------------
SUBROUTINE age_mass_Z_November(a,m,zz,Ext, age, mass, Z, Ebv, age_indice, mass_indice, Z_indice, Ebv_indice)
!--------------------------------------------------------------------------------------

 implicit none

 real(4), intent(out) :: age, mass, Z, Ebv
 integer, intent(in) :: a, m, zz, Ext
 integer :: ii,kk
 character(len=50), intent(out) :: age_indice, mass_indice, Z_indice, Ebv_indice


	if (a == 1) then
	 age  = 7.00
	 age_indice = '700'	!10 Myr
	elseif (a == 2) then
	 age = 7.10
	 age_indice = '710'	
	elseif (a == 3) then
	 age = 7.20
	 age_indice = '720'
	elseif (a == 4) then
	 age = 7.30
	 age_indice = '730'	
	elseif (a == 5) then
	 age = 7.40
	 age_indice = '740'	
	elseif (a == 6) then
	 age = 7.50
	 age_indice = '750'	
	elseif (a == 7) then
	 age = 7.60
	 age_indice = '760'	
	elseif (a == 8) then
	 age = 7.70
	 age_indice = '770'	
	elseif (a == 9) then
	 age = 7.80
	 age_indice = '780'	
	elseif (a == 10) then
	 age = 7.90
	 age_indice = '790'	
	elseif (a == 11) then
	 age = 8.00
	 age_indice = '800'	
	elseif (a == 12) then
	 age = 8.10
	 age_indice = '810'	
	elseif (a == 13) then
	 age = 8.20
	 age_indice = '820'	
	elseif (a == 14) then
	 age = 8.30
	 age_indice = '830'	
	elseif (a == 15) then
	 age = 8.40
	 age_indice = '840'	
	elseif (a == 16) then
	 age = 8.50
	 age_indice = '850'	
	elseif (a == 17) then
	 age = 8.60
	 age_indice = '860'	
	elseif (a == 18) then
	 age = 8.70
	 age_indice = '870'	
	elseif (a == 19) then
	 age = 8.80
	 age_indice = '880'	
	elseif (a == 20) then
	 age = 8.90
	 age_indice = '890'	
	elseif (a == 21) then
	 age = 9.00
	 age_indice = '900'	
	elseif (a == 22) then
	 age = 9.10
	 age_indice = '910'	
	elseif (a == 23) then
	 age = 9.20
	 age_indice = '920'	
	elseif (a == 24) then
	 age = 9.30
	 age_indice = '930'	
	elseif (a == 25) then
	 age = 9.40
	 age_indice = '940'	
	elseif (a == 26) then
	 age = 9.50
	 age_indice = '950'	
	elseif (a == 27) then
	 age = 9.60
	 age_indice = '960'	
	elseif (a == 28) then
	 age = 9.70
	 age_indice = '970'	
	elseif (a == 29) then
	 age = 9.80
	 age_indice = '980'	
	elseif (a == 30) then
	 age = 9.90
	 age_indice = '990'	
	elseif (a == 31) then
	 age = 10.0
	 age_indice = '1000'
	endif


	if (m == 1) then 
	 mass = 10**2.00
	 mass_indice = '200'    !Solar masses	
	elseif (m == 2) then
	 mass = 10**2.10
	 mass_indice = '210'	!Solar masses
	elseif (m == 3) then
	 mass = 10**2.20
	 mass_indice = '220'	!Solar masses			
	elseif (m == 4) then
	 mass = 10**2.30
	 mass_indice = '230'	!Solar masses
	elseif (m == 5) then
	 mass = 10**2.40
	 mass_indice = '240'	!Solar masses
	elseif (m == 6) then
	 mass = 10**2.50
	 mass_indice = '250'	!Solar masses
	elseif (m == 7) then
	 mass = 10**2.60
	 mass_indice = '260'	!Solar masses
	elseif (m == 8) then
	 mass = 10**2.70
	 mass_indice = '270'	!Solar masses
	elseif (m == 9) then
	 mass = 10**2.80
	 mass_indice = '280'	!Solar masses
	elseif (m == 10) then
	 mass = 10**2.90
	 mass_indice = '290'	!Solar masses
	elseif (m == 11) then 
	 mass = 10**3.00
	 mass_indice = '300'    !Solar masses	
	elseif (m == 12) then
	 mass = 10**3.10
	 mass_indice = '310'	!Solar masses
	elseif (m == 13) then
	 mass = 10**3.20
	 mass_indice = '320'	!Solar masses			
	elseif (m == 14) then
	 mass = 10**3.30
	 mass_indice = '330'	!Solar masses
	elseif (m == 15) then
	 mass = 10**3.40
	 mass_indice = '340'	!Solar masses
	elseif (m == 16) then
	 mass = 10**3.50
	 mass_indice = '350'	!Solar masses
	elseif (m == 17) then
	 mass = 10**3.60
	 mass_indice = '360'	!Solar masses
	elseif (m == 18) then
	 mass = 10**3.70
	 mass_indice = '370'	!Solar masses
	elseif (m == 19) then
	 mass = 10**3.80
	 mass_indice = '380'	!Solar masses
	elseif (m == 20) then
	 mass = 10**3.90
	 mass_indice = '390'	!Solar masses
	elseif (m == 21) then
	 mass = 10**4.00
	 mass_indice = '400'	!Solar masses
	elseif (m == 22) then
	 mass = 10**4.10
	 mass_indice = '410'	!Solar masses
	elseif (m == 23) then
	 mass = 10**4.20
	 mass_indice = '420'	!Solar masses
	elseif (m == 24) then
	 mass = 10**4.30
	 mass_indice = '430'	!Solar masses
	elseif (m == 25) then
	 mass = 10**4.40
	 mass_indice = '440'	!Solar masses
	elseif (m == 26) then
	 mass = 10**4.50
	 mass_indice = '450'	!Solar masses
	elseif (m == 27) then
	 mass = 10**4.60
	 mass_indice = '460'	!Solar masses
	elseif (m == 28) then
	 mass = 10**4.70
	 mass_indice = '470'	!Solar masses
	elseif (m == 29) then
	 mass = 10**4.80
	 mass_indice = '480'	!Solar masses
	elseif (m == 30) then
	 mass = 10**4.90
	 mass_indice = '490'	!Solar masses
	elseif (m == 31) then
	 mass = 10**5.00
	 mass_indice = '500'	!Solar masses
	elseif (m == 32) then
	 mass = 10**5.10
	 mass_indice = '510'	!Solar masses
	elseif (m == 33) then
	 mass = 10**5.20
	 mass_indice = '520'	!Solar masses
	elseif (m == 34) then
	 mass = 10**5.30
	 mass_indice = '530'	!Solar masses
	elseif (m == 35) then
	 mass = 10**5.40
	 mass_indice = '540'	!Solar masses
	elseif (m == 36) then
	 mass = 10**5.50
	 mass_indice = '550'	!Solar masses
	elseif (m == 37) then
	 mass = 10**5.60
	 mass_indice = '560'	!Solar masses
	elseif (m == 38) then
	 mass = 10**5.70
	 mass_indice = '570'	!Solar masses
	elseif (m == 39) then
	 mass = 10**5.80
	 mass_indice = '580'	!Solar masses
	elseif (m == 40) then
	 mass = 10**5.90
	 mass_indice = '590'	!Solar masses
	elseif (m == 41) then
	 mass = 10**6.00
	 mass_indice = '600'	!Solar masses
	endif


	if (zz == 1) then	
	 Z = 0.00010
	 Z_indice = 'm23'	! [M/H] = -2.3
	elseif (zz == 2) then
	 Z =  0.00012619
	 Z_indice = 'm22'	! [M/H] = -2.2
        elseif (zz == 3) then
         Z = 0.00015887
         Z_indice = 'm21'	! [M/H] = -2.1
	elseif (zz == 4) then
	 Z = 0.00020000
	 Z_indice = 'm20'	! [M/H] = -2.0
	elseif (zz == 5) then
	 Z = 0.00025
	 Z_indice = 'm19'	! [M/H] = -1.9
	elseif (zz == 6) then
	 Z = 0.00031698
	 Z_indice = 'm18'	! [M/H] = -1.8
        elseif (zz == 7) then
         Z = 0.00039905
         Z_indice = 'm17'	! [M/H] = -1.7
	elseif (zz == 8) then
	 Z = 0.00050
	 Z_indice = 'm16'	! [M/H] = -1.6
        elseif (zz == 9) then 
         Z = 0.00063246
         Z_indice = 'm15'	! [M/H] = -1.5
	elseif (zz == 10) then
	 Z = 0.00079621
	 Z_indice = 'm14'	! [M/H] = -1.4
	elseif (zz == 11) then
	 Z = 0.00100
	 Z_indice = 'm13'	! [M/H] = -1.3
	elseif (zz == 12) then
	 Z = 0.0012619
	 Z_indice = 'm12'	! [M/H] = -1.2
        elseif (zz == 13) then
         Z = 0.0015887
         Z_indice = 'm11'	! [M/H] = -1.1
	elseif (zz == 14) then
	 Z = 0.00200
	 Z_indice = 'm10'	! [M/H] = -1.0
        elseif (zz == 15) then
         Z = 0.0025179 
         Z_indice = 'm09'	! [M/H] = -0.9
	elseif (zz == 16) then
	 Z = 0.00317
	 Z_indice = 'm08'	! [M/H] = -0.8
        elseif (zz == 17) then
         Z = 0.0039905
         Z_indice = 'm07'	! [M/H] = -0.7
	elseif (zz == 18) then
	 Z = 0.00502
	 Z_indice = 'm06'	! [M/H] = -0.6
        elseif (zz == 19) then
         Z = 0.0063246
         Z_indice = 'm05'	! [M/H] = -0.5
	elseif (zz == 20) then
	 Z = 0.00796
	 Z_indice = 'm04'	! [M/H] = -0.4
        elseif (zz == 21) then
         Z = 0.010024
         Z_indice = 'm03'	! [M/H] = -0.3
	elseif (zz == 22) then
	 Z = 0.01262
	 Z_indice = 'm02'	! [M/H] = -0.2
        elseif (zz == 23) then 
         Z = 0.015887
         Z_indice = 'm01'	! [M/H] = -0.1
	elseif (zz == 24) then
	 Z = 0.02000
	 Z_indice = 'n00'	! [M/H] = +0.0
        elseif (zz == 25) then
         Z = 0.025179
         Z_indice = 'p01'	! [M/H] = +0.1
	elseif (zz == 26) then
	 Z = 0.03000
	 Z_indice = 'p02'	! [M/H] = +0.2 ABOUT, not exactly
	endif




 	!Extinction 
        if  (Ext == 1) then
         Ebv = 0.0
         Ebv_indice = '00'
        elseif  (Ext == 2) then
         Ebv = 0.1
         Ebv_indice = '01'
        elseif  (Ext == 3) then
         Ebv = 0.2
         Ebv_indice = '02'
        elseif  (Ext == 4) then
         Ebv = 0.3
         Ebv_indice = '03'
        elseif  (Ext == 5) then
         Ebv = 0.4
         Ebv_indice = '04'
        elseif  (Ext == 6) then
         Ebv = 0.5
         Ebv_indice = '05'
        elseif  (Ext == 7) then
         Ebv = 0.6
         Ebv_indice = '06'
        elseif  (Ext == 8) then
         Ebv = 0.7
         Ebv_indice = '07'
        elseif  (Ext == 9) then
         Ebv = 0.8
         Ebv_indice = '08'
        elseif  (Ext == 10) then
         Ebv = 0.9
         Ebv_indice = '09'
        elseif  (Ext == 11) then
         Ebv = 1.0
         Ebv_indice = '10'
        elseif  (Ext == 12) then
         Ebv = 1.1
         Ebv_indice = '11'
        elseif  (Ext == 13) then
         Ebv = 1.2
         Ebv_indice = '12'
        elseif  (Ext == 14) then
         Ebv = 1.3
         Ebv_indice = '13'
        elseif  (Ext == 15) then
         Ebv = 1.4
         Ebv_indice = '14'
        elseif  (Ext == 16) then
         Ebv = 1.5
         Ebv_indice = '15'
        elseif  (Ext == 17) then
         Ebv = 1.6
         Ebv_indice = '16'
        elseif  (Ext == 18) then
         Ebv = 1.7
         Ebv_indice = '17'
        elseif  (Ext == 19) then
         Ebv = 1.8
         Ebv_indice = '18'
        elseif  (Ext == 20) then
         Ebv = 1.9
         Ebv_indice = '19'
        elseif  (Ext == 21) then
         Ebv = 2.0
         Ebv_indice = '20'
        endif
 return


 END SUBROUTINE



!--------------------------------------------------------------------------------------
SUBROUTINE age_mass_Z_December(a,m,zz,Ext, age, mass, Z, Ebv, age_indice, mass_indice, Z_indice, Ebv_indice)
!--------------------------------------------------------------------------------------

 implicit none

 real(4), intent(out) :: age, mass, Z, Ebv
 integer, intent(in) :: a, m, zz, Ext
 integer :: ii,kk,age_int,mass_int
 character(len=50), intent(out) :: age_indice, mass_indice, Z_indice, Ebv_indice

        age  = 6.60
	!Do ii = 1,63
	! if (a == ii) then
	age = age + (a-1)*0.05
	age_int = nint(age*100)
        write(age_indice,'(I4)') age_int
        age_indice = adjustl(adjustr(age_indice))
	! endif
	!Enddo

        mass  = 2.00
	!Do ii = 1,81
	! if (m == ii) then
	mass = mass + (m-1)*0.05
	mass_int = nint(mass*100)
        mass = 10**mass
        write(mass_indice,'(I4)') mass_int
        mass_indice = adjustl(adjustr(mass_indice))
	! endif
	!Enddo

	if (zz == 1) then	
	 Z = 0.00010
	 Z_indice = 'm23'	! [M/H] = -2.3
	elseif (zz == 2) then
	 Z =  0.00012619
	 Z_indice = 'm22'	! [M/H] = -2.2
        elseif (zz == 3) then
         Z = 0.00015887
         Z_indice = 'm21'	! [M/H] = -2.1
	elseif (zz == 4) then
	 Z = 0.00020000
	 Z_indice = 'm20'	! [M/H] = -2.0
	elseif (zz == 5) then
	 Z = 0.00025
	 Z_indice = 'm19'	! [M/H] = -1.9
	elseif (zz == 6) then
	 Z = 0.00031698
	 Z_indice = 'm18'	! [M/H] = -1.8
        elseif (zz == 7) then
         Z = 0.00039905
         Z_indice = 'm17'	! [M/H] = -1.7
	elseif (zz == 8) then
	 Z = 0.00050
	 Z_indice = 'm16'	! [M/H] = -1.6
        elseif (zz == 9) then 
         Z = 0.00063246
         Z_indice = 'm15'	! [M/H] = -1.5
	elseif (zz == 10) then
	 Z = 0.00079621
	 Z_indice = 'm14'	! [M/H] = -1.4
	elseif (zz == 11) then
	 Z = 0.00100
	 Z_indice = 'm13'	! [M/H] = -1.3
	elseif (zz == 12) then
	 Z = 0.0012619
	 Z_indice = 'm12'	! [M/H] = -1.2
        elseif (zz == 13) then
         Z = 0.0015887
         Z_indice = 'm11'	! [M/H] = -1.1
	elseif (zz == 14) then
	 Z = 0.00200
	 Z_indice = 'm10'	! [M/H] = -1.0
        elseif (zz == 15) then
         Z = 0.0025179 
         Z_indice = 'm09'	! [M/H] = -0.9
	elseif (zz == 16) then
	 Z = 0.00317
	 Z_indice = 'm08'	! [M/H] = -0.8
        elseif (zz == 17) then
         Z = 0.0039905
         Z_indice = 'm07'	! [M/H] = -0.7
	elseif (zz == 18) then
	 Z = 0.00502
	 Z_indice = 'm06'	! [M/H] = -0.6
        elseif (zz == 19) then
         Z = 0.0063246
         Z_indice = 'm05'	! [M/H] = -0.5
	elseif (zz == 20) then
	 Z = 0.00796
	 Z_indice = 'm04'	! [M/H] = -0.4
        elseif (zz == 21) then
         Z = 0.010024
         Z_indice = 'm03'	! [M/H] = -0.3
	elseif (zz == 22) then
	 Z = 0.01262
	 Z_indice = 'm02'	! [M/H] = -0.2
        elseif (zz == 23) then 
         Z = 0.015887
         Z_indice = 'm01'	! [M/H] = -0.1
	elseif (zz == 24) then
	 Z = 0.02000
	 Z_indice = 'n00'	! [M/H] = +0.0
        elseif (zz == 25) then
         Z = 0.025179
         Z_indice = 'p01'	! [M/H] = +0.1
	elseif (zz == 26) then
	 Z = 0.03000
	 Z_indice = 'p02'	! [M/H] = +0.2 ABOUT, not exactly
	endif


 	!Extinction 

	Ebv = 0.
	!Do ii = 1,101
	! if (Ext == ii) then
	Ebv = Ebv + (Ext-1)*0.02
        Ebv_indice = adjustl(adjustr(Ebv_indice))
	! endif
	!Enddo

        
 return


 END SUBROUTINE



!--------------------------------------------------------------------------------------
SUBROUTINE age_mass_Z_December_Inverse(a,m,zz,Ext,age,mass,Z,Ebv, Z_indice)
!--------------------------------------------------------------------------------------

 implicit none

 real(4), intent(in) :: age, mass, Z, Ebv
 integer, intent(out) :: a, m, zz, Ext
 integer :: ii,kk,age_int,mass_int
 !character(len=50), intent(out) :: age_indice, mass_indice, Z_indice, Ebv_indice
 character(len=3), intent(in) :: Z_indice

	!Age
	a = 20*(age-6.60)+1
	!Mass
        m = 20*(mass-2.00)+1
	!Metallicity
	if (Z >= 0.00009 .and. Z <= 0.00011) zz = 1		!Z = 0.00010  ! [M/H] = -2.3
	if (Z_indice == 'm23') zz = 1
	if (Z >= 0.00012 .and. Z <= 0.00014) zz = 2		!Z = 0.00013  ! [M/H] = -2.2
	if (Z_indice == 'm22') zz = 2
	if (Z >= 0.00015 .and. Z <= 0.00017) zz = 3		!Z = 0.00016  ! [M/H] = -2.1
	if (Z_indice == 'm21') zz = 3
	if (Z >= 0.00018 .and. Z <= 0.00021) zz = 4		!Z = 0.00020  ! [M/H] = -2.0
	if (Z_indice == 'm20') zz = 4
	if (Z >= 0.00023 .and. Z <= 0.00027) zz = 5		!Z = 0.00025  ! [M/H] = -1.9
	if (Z_indice == 'm19') zz = 5
	if (Z >= 0.00028 .and. Z <= 0.00032) zz = 6		!Z = 0.00030  ! [M/H] = -1.8
	if (Z_indice == 'm18') zz = 6
	if (Z >= 0.00035 .and. Z <= 0.00044) zz = 7		!Z = 0.00040  ! [M/H] = -1.7
	if (Z_indice == 'm17') zz = 7
	if (Z >= 0.00045 .and. Z <= 0.00055) zz = 8		!Z = 0.00050  ! [M/H] = -1.6
	if (Z_indice == 'm16') zz = 8
	if (Z >= 0.00060 .and. Z <= 0.00070) zz = 9		!Z = 0.00063  ! [M/H] = -1.5
	if (Z_indice == 'm15') zz = 9
	if (Z >= 0.00074 .and. Z <= 0.00082) zz = 10		!Z = 0.00080  ! [M/H] = -1.4
	if (Z_indice == 'm14') zz = 10
	if (Z >= 0.00090 .and. Z <= 0.00105) zz = 11		!Z = 0.00100  ! [M/H] = -1.3
	if (Z_indice == 'm13') zz = 11
	if (Z >= 0.00110 .and. Z <= 0.00140) zz = 12		!Z = 0.00126  ! [M/H] = -1.2
	if (Z_indice == 'm12') zz = 12
	if (Z >= 0.00145 .and. Z <= 0.00180) zz = 13		!Z = 0.00160  ! [M/H] = -1.1
	if (Z_indice == 'm11') zz = 13
	if (Z >= 0.00185 .and. Z <= 0.00205) zz = 14		!Z = 0.00190  ! [M/H] = -1.0
	if (Z_indice == 'm10') zz = 14
	if (Z >= 0.00230 .and. Z <= 0.00280) zz = 15		!Z = 0.00250  ! [M/H] = -0.9
	if (Z_indice == 'm09') zz = 15
	if (Z >= 0.00300 .and. Z <= 0.00330) zz = 16		!Z = 0.00317  ! [M/H] = -0.8
	if (Z_indice == 'm08') zz = 16
	if (Z >= 0.00350 .and. Z <= 0.00450) zz = 17		!Z = 0.00400  ! [M/H] = -0.7
	if (Z_indice == 'm07') zz = 17
	if (Z >= 0.00470 .and. Z <= 0.00530) zz = 18		!Z = 0.00500  ! [M/H] = -0.6
	if (Z_indice == 'm06') zz = 18
	if (Z >= 0.00550 .and. Z <= 0.00700) zz = 19		!Z = 0.00632  ! [M/H] = -0.5
	if (Z_indice == 'm05') zz = 19
	if (Z >= 0.00750 .and. Z <= 0.00810) zz = 20		!Z = 0.00796  ! [M/H] = -0.4
	if (Z_indice == 'm04') zz = 20
	if (Z >= 0.00900 .and. Z <= 0.01100) zz = 21		!Z = 0.01000  ! [M/H] = -0.3
	if (Z_indice == 'm03') zz = 21
	if (Z >= 0.01200 .and. Z <= 0.01350) zz = 22		!Z = 0.01300  ! [M/H] = -0.2
	if (Z_indice == 'm02') zz = 22
	if (Z >= 0.01400 .and. Z <= 0.01700) zz = 23		!Z = 0.01600  ! [M/H] = -0.1
	if (Z_indice == 'm01') zz = 23
	if (Z >= 0.01800 .and. Z <= 0.02100) zz = 24		!Z = 0.02000  ! [M/H] = +0.0
	if (Z_indice == 'n00') zz = 24
	if (Z >= 0.02300 .and. Z <= 0.02700) zz = 25		!Z = 0.02500  ! [M/H] = +0.1  
	if (Z_indice == 'p01') zz = 25
	if (Z >= 0.02800 .and. Z <= 0.03500) zz = 26		!Z = 0.03000  ! [M/H] = +0.2 ABOUT, not exactly
	if (Z_indice == 'p02') zz = 26
 	!Extinction 
        Ext = 50*Ebv + 1       
 return
 END SUBROUTINE









!--------------------------------------------------------------------------------------
SUBROUTINE age_mass_Z_March2012(a,m,zz,Ext, age, mass, Z, Ebv, age_indice, mass_indice, Z_indice, Ebv_indice)
!--------------------------------------------------------------------------------------

 implicit none

 real(4), intent(out) :: age, mass, Z, Ebv
 integer, intent(in) :: a, m, zz, Ext
 integer :: ii,kk,age_int,mass_int
 character(len=50), intent(out) :: age_indice, mass_indice, Z_indice, Ebv_indice

        age  = 6.60
	!Do ii = 1,353
	! if (a == ii) then
	age = age + (a-1)*0.01
	age_int = nint(age*100)
        write(age_indice,'(I4)') age_int
        age_indice = adjustl(adjustr(age_indice))
	! endif
	!Enddo

        mass  = 2.00
	!Do ii = 1,301
	! if (m == ii) then
	mass = mass + (m-1)*0.01
	mass_int = nint(mass*100)
        mass = 10**mass
        write(mass_indice,'(I4)') mass_int
        mass_indice = adjustl(adjustr(mass_indice))
	! endif
	!Enddo


	if (zz == 1) then	
	 Z = 0.00010
	 Z_indice = 'm23'	! [M/H] = -2.3
	elseif (zz == 2) then
	 Z =  0.00012619
	 Z_indice = 'm22'	! [M/H] = -2.2
        elseif (zz == 3) then
         Z = 0.00015887
         Z_indice = 'm21'	! [M/H] = -2.1
	elseif (zz == 4) then
	 Z = 0.00020000
	 Z_indice = 'm20'	! [M/H] = -2.0
	elseif (zz == 5) then
	 Z = 0.00025
	 Z_indice = 'm19'	! [M/H] = -1.9
	elseif (zz == 6) then
	 Z = 0.00031698
	 Z_indice = 'm18'	! [M/H] = -1.8
        elseif (zz == 7) then
         Z = 0.00039905
         Z_indice = 'm17'	! [M/H] = -1.7
	elseif (zz == 8) then
	 Z = 0.00050
	 Z_indice = 'm16'	! [M/H] = -1.6
        elseif (zz == 9) then 
         Z = 0.00063246
         Z_indice = 'm15'	! [M/H] = -1.5
	elseif (zz == 10) then
	 Z = 0.00079621
	 Z_indice = 'm14'	! [M/H] = -1.4
	elseif (zz == 11) then
	 Z = 0.00100
	 Z_indice = 'm13'	! [M/H] = -1.3
	elseif (zz == 12) then
	 Z = 0.0012619
	 Z_indice = 'm12'	! [M/H] = -1.2
        elseif (zz == 13) then
         Z = 0.0015887
         Z_indice = 'm11'	! [M/H] = -1.1
	elseif (zz == 14) then
	 Z = 0.00200
	 Z_indice = 'm10'	! [M/H] = -1.0
        elseif (zz == 15) then
         Z = 0.0025179 
         Z_indice = 'm09'	! [M/H] = -0.9
	elseif (zz == 16) then
	 Z = 0.00317
	 Z_indice = 'm08'	! [M/H] = -0.8
        elseif (zz == 17) then
         Z = 0.0039905
         Z_indice = 'm07'	! [M/H] = -0.7
	elseif (zz == 18) then
	 Z = 0.00502
	 Z_indice = 'm06'	! [M/H] = -0.6
        elseif (zz == 19) then
         Z = 0.0063246
         Z_indice = 'm05'	! [M/H] = -0.5
	elseif (zz == 20) then
	 Z = 0.00796
	 Z_indice = 'm04'	! [M/H] = -0.4
        elseif (zz == 21) then
         Z = 0.010024
         Z_indice = 'm03'	! [M/H] = -0.3
	elseif (zz == 22) then
	 Z = 0.01262
	 Z_indice = 'm02'	! [M/H] = -0.2
        elseif (zz == 23) then 
         Z = 0.015887
         Z_indice = 'm01'	! [M/H] = -0.1
	elseif (zz == 24) then
	 Z = 0.02000
	 Z_indice = 'n00'	! [M/H] = +0.0
        elseif (zz == 25) then
         Z = 0.025179
         Z_indice = 'p01'	! [M/H] = +0.1
	elseif (zz == 26) then
	 Z = 0.03000
	 Z_indice = 'p02'	! [M/H] = +0.2 ABOUT, not exactly
	endif


 	!Extinction 

	Ebv = 0.
	!Do ii = 1,101
	! if (Ext == ii) then
	Ebv = Ebv + (Ext-1)*0.02
        Ebv_indice = adjustl(adjustr(Ebv_indice))
	! endif
	!Enddo
 return
 END SUBROUTINE



!--------------------------------------------------------------------------------------
SUBROUTINE age_mass_Z_March2012_Inverse(a,m,zz,Ext,age,mass,Z,Ebv)
!--------------------------------------------------------------------------------------

 implicit none

 real(4), intent(in) :: age, mass, Z, Ebv
 integer, intent(out) :: a, m, zz, Ext
 integer :: ii,kk,age_int,mass_int
 !character(len=50), intent(out) :: age_indice, mass_indice, Z_indice, Ebv_indice

	!Age
	a = 100*(age-6.60)+1
	!Mass
        m = 100*(mass-2.00)+1
	!Metallicity
	if (Z >= 0.00009 .and. Z <= 0.00011) zz = 1		!Z = 0.00010  ! [M/H] = -2.3
	if (Z >= 0.00012 .and. Z <= 0.00014) zz = 2		!Z = 0.00013  ! [M/H] = -2.2
	if (Z >= 0.00015 .and. Z <= 0.00017) zz = 3		!Z = 0.00016  ! [M/H] = -2.1
	if (Z >= 0.00018 .and. Z <= 0.00021) zz = 4		!Z = 0.00020  ! [M/H] = -2.0
	if (Z >= 0.00023 .and. Z <= 0.00027) zz = 5		!Z = 0.00025  ! [M/H] = -1.9
	if (Z >= 0.00028 .and. Z <= 0.00032) zz = 6		!Z = 0.00030  ! [M/H] = -1.8
	if (Z >= 0.00035 .and. Z <= 0.00044) zz = 7		!Z = 0.00040  ! [M/H] = -1.7
	if (Z >= 0.00045 .and. Z <= 0.00055) zz = 8		!Z = 0.00050  ! [M/H] = -1.6
	if (Z >= 0.00060 .and. Z <= 0.00070) zz = 9		!Z = 0.00063  ! [M/H] = -1.5
	if (Z >= 0.00074 .and. Z <= 0.00082) zz = 10		!Z = 0.00080  ! [M/H] = -1.4
	if (Z >= 0.00090 .and. Z <= 0.00105) zz = 11		!Z = 0.00100  ! [M/H] = -1.3
	if (Z >= 0.00110 .and. Z <= 0.00140) zz = 12		!Z = 0.00126  ! [M/H] = -1.2
	if (Z >= 0.00145 .and. Z <= 0.00180) zz = 13		!Z = 0.00160  ! [M/H] = -1.1
	if (Z >= 0.00185 .and. Z <= 0.00205) zz = 14		!Z = 0.00190  ! [M/H] = -1.0
	if (Z >= 0.00230 .and. Z <= 0.00280) zz = 15		!Z = 0.00250  ! [M/H] = -0.9
	if (Z >= 0.00300 .and. Z <= 0.00330) zz = 16		!Z = 0.00317  ! [M/H] = -0.8
	if (Z >= 0.00350 .and. Z <= 0.00450) zz = 17		!Z = 0.00400  ! [M/H] = -0.7
	if (Z >= 0.00470 .and. Z <= 0.00530) zz = 18		!Z = 0.00500  ! [M/H] = -0.6
	if (Z >= 0.00550 .and. Z <= 0.00700) zz = 19		!Z = 0.00632  ! [M/H] = -0.5
	if (Z >= 0.00750 .and. Z <= 0.00810) zz = 20		!Z = 0.00796  ! [M/H] = -0.4
	if (Z >= 0.00900 .and. Z <= 0.01100) zz = 21		!Z = 0.01000  ! [M/H] = -0.3
	if (Z >= 0.01200 .and. Z <= 0.01350) zz = 22		!Z = 0.01300  ! [M/H] = -0.2
	if (Z >= 0.01400 .and. Z <= 0.01700) zz = 23		!Z = 0.01600  ! [M/H] = -0.1
	if (Z >= 0.01800 .and. Z <= 0.02100) zz = 24		!Z = 0.02000  ! [M/H] = +0.0
	if (Z >= 0.02300 .and. Z <= 0.02700) zz = 25		!Z = 0.02500  ! [M/H] = +0.1  
	if (Z >= 0.02800 .and. Z <= 0.03500) zz = 26		!Z = 0.03000  ! [M/H] = +0.2 ABOUT, not exactly
 	!Extinction 
        Ext = 50*Ebv + 1       
 return
 END SUBROUTINE


END MODULE Lecture_module_fortran

