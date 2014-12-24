PROGRAM Fisher_Distance
  USE cosmo
  USE growth
  USE linearpk
  USE angular_distance
  ! Ref: Seo & Eisenstein, ApJ, 665, 14 (2007) .. Eq. (25)
  IMPLICIT none
  integer :: npara=2 ! # of parameters
  double precision, dimension(2)   :: der
  double precision, dimension(2,2) :: cov,fis,fistot
  double precision, dimension(6,6) :: fis3x3, fis6x6
  integer, dimension(7):: list
  integer, dimension(2)   :: work
  double precision :: sigma0=12.4d0*0.817d0/0.9d0 ! non-linearity scale in units of h Mpc^-1, rescaled by WMAP5's sigma8 value (0.817)
  double precision :: BAO_AMP=0.5817d0 ! A_0 in Seo&Eisenstein
  double precision :: kmax_ov_h,k0_ov_h, kmin_ov_h, k0_new, k_new! h Mpc^-1
  double precision :: linear_pk,pk,p02,p01,Pb, fistest
  integer:: status
  double precision:: k_full, fbao, pk_ref, dfbao, k_ov_h
  double precision :: wdamp_perp, wdamp_para, wdamp_silk, wdamp_factor, wdamp_sinterm
  double precision :: sigma_para,sigma_perp,sigma_silk=8.38d0 ! h^-1 Mpc
  double precision :: dlnk,mu,mu2,dmu,factor,wkmu,dummy,wdamp,Rmu
  double precision :: z,zin,beta,bias,g,dgdlna,sigma_z,delta_v,fz, w,dvdz, Dz
  double precision :: zmin,zmax,area,Vsurvey,Vtot,Ntot,Ngl, ngal,dndz,da, one_over_h, Vsurvey2, dvdz2
  character(len=128) :: filename
  integer :: n,i,j,ibin,nbins,ifile, n1, n2
  external linear_pk,dgdlna,g,dvdz,da, one_over_h
  ! Specify three cosmological parameters for computing the growth factor and the angular diameter distance. 
  ! The data type has been defined in MODULE cosmo.
  !############################################
  
  !******************* PLANCK ********************************
  H_0 = 67.d0
  c   = 299792.458d0
  om0 = 0.267d0
  ode0= 0.686d0
  ok0 = 0d0
  ob0 = 0.049
  ok0 = 0.d0
  w0  =-1.d0
  w_a = 0.d0
  CALL setup_growth ! tabulate the growth factor
  CALL setup_da       ! tabulate the angular diameter distance
  ! ask for survey specific parameters
  !==============Enter the File name ================================== 
  !filename = '0n_Corrected_nz_Srms_2.txt'
  filename = 'MID_MK_B2_Real_Corrected_nz_Srms_2.txt'
  !filename = 'SKA1MID_B2_Opt_Corrected_nz_Srms_2.txt'
  !filename = 'SKA1MID_B2_Pess_Corrected_nz_Srms_2.txt'
  !filename = 'SKA2_Opt_Corrected_nz_Srms_2.txt'
  !filename = 'SKA2_Real_Corrected_nz_Srms_2.txt'
  !filename =  'SKA2_Pess_Corrected_nz_Srms_2.txt'!
  !filename= '/home/ubuntu/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_bao_SKA_Euclid/inputs/number_EuclidmJy_ref.txt'
  open(2,file=filename,status='unknown')
  
  !=============  Read in linear P(k) ===================================
  !filename='linear_matterpower_1.dat' 
  filename= 'deriv_mikrom_dfbao_dk.dat'
  !open(22,file='deriv_mikrom_dfbao_dk.dat',status='unknown')
  !filename= 'wmap5baosn_max_likelihood_matterpower_at_z=30.dat'
  n  =129! no. of lines in the file
  zin=1d0 ! redshift of the input power spectrum
  CALL open_linearpk(filename,n) ! This fuction just to interpolate the values in the filename	
  
  !================================================================  
  ! loop over redshift bins
  fistot=0d0   ! total Fisher matrix integrated over all the redshift bins
  fis6x6=0d0 ! total Fisher matrix for w, Omega_k, and Omega_m
  Vtot  =0d0   ! total survey volume
  Ntot  =0d0
  nbins =0    ! initialize the number of bins
  !=========== Loop over redshift bins =======================================
  do ibin=1,100
     !read(2,*,end=10)z ,dndz,bias,kmax_ov_h,kmin_ov_h,delta_v, Vsurvey , dvdz2! uncomment for Euclid ref
     read(2,*,end=10) z ,dndz,bias, kmax_ov_h, Vsurvey, dvdz2! Uncomment for SKA
     nbins=nbins+1
     !--------------------------------------------------------------------
     beta  = (1d0+dgdlna(z)/g(z))/bias
     w     = w0 + w_a*(z/(1d0+z))
     !ngal  = dndz  		  		! Uncomment for  Euclid
     dvdz2 = dvdz2*(H_0/100d0)**3 		! convert the comoving volume to Mpc^3 h^-3 Uncomment for SKA
     ngal = dndz/((3.14159d0/180d0)**2*dvdz2) 	! Uncomment for SKA or Euclid new
     Ngl   = (dndz * Vsurvey) ! for Euclid
     Ntot  = Ntot + Ngl 
     Vtot  = Vtot+Vsurvey ! h^-3 Mpc^3.
     !==================Test the growth =================================
     !open(101,file='test_all_terms_Euclid.txt' ,status='unknown')
     !write(101,*) z, ngal, Vsurvey, g(z) , beta*bias, bias, beta!, D_plus!/(1d0+z), (1d0+dgdlna(z)/g(z))
     !write(101,*) z, Vsurvey, dvdz2, dndz, ngal, bias, beta, beta*bias, g(z)
     !=================  computing the Fisher matrix... =================================
     open(11,file=filename,status='old')
     read(11,*)k0_ov_h,dummy, dummy, dummy
     k_ov_h=k0_ov_h
     fis=0d0
     fistest= 0d0
     !================ loop over k ===================================
     do while (k_ov_h <=kmax_ov_h)
        read(11,*)k_ov_h, fbao, dfbao,pk_ref
        !dlnk=dlog(k_ov_h)-dlog(k0_ov_h) ! uncomment for widder bins
        !print*, 'dlnk', dlnk 
        !===========================================================
         k0_new=6.8492999999999996E-005 
         k_new =9.5589999999999998E-005
         dlnk=dlog(k_new)-dlog(k0_new) !uncomment to use with Plank's parameters 
         !factor=(k_ov_h)**3d0*dlnk/(2d0 * 2d0 *3.1415926535d0**2d0) ! h^3 Mpc^-3
         factor=(k_ov_h)**3d0*dlnk/(8.d0 *3.1415926535d0**2d0) ! h^3 Mpc^-3
         !print*, k_ov_h
         !=========== loop over mu..========================
         mu=0d0
         dmu=1d-3
         !print*, g(z)/g(zin)
           do while (mu<=1d0)
             mu2= mu*mu
	     
             !P02 is P(k) at k=0.2 h Mpc^-1 in unitys of h^-3 Mpc^3
             p02 = 848.536902895
             p02 =p02*(g(z)/g(zin)*(1d0+zin)/(1d0+z))**2d0 &
              *bias**2d0*(1d0+beta*mu2)**2d0
             
	     ! P(k) the galaxy power spectrum  in units of h^-3 Mpc^3  
             pk=pk_ref *(1d0 + fbao)*((g(z)/g(zin))*((1d0+zin)/(1d0+z)))**2d0 &
               *bias**2d0*(1d0+beta*mu2)**2d0      
             !print*, pk
             !==========================================     
             wkmu=Vsurvey*(ngal/(1d0+ngal*pk))**2d0
             !print*, k_ov_h ,  wkmu
             wdamp = ((dfbao/(1d0 + fbao))*pk)**2d0
             der(1)=(mu2-1d0)*k_ov_h ! dPkdlogDa
             der(2)=mu2*k_ov_h ! dPkdlogH
             !print*,  der(1), der(2)
             do i=1,npara
               fistest=factor*Vsurvey*(ngal*pk/(1d0+ngal*pk))**2d0
               fis(i,:)= fis(i,:) +factor* wkmu*wdamp*der(i)*der(:)*dmu 
           enddo
           mu=mu+dmu
        enddo
        !================= save files ====================================== 
        !open(5,file='Pk_Euclid_testPlank.txt' ,status='unknown')
        !write(5,*) k_ov_h, pk! , 1d0/ngal
      
        !open(15,file='gz.txt' ,status='unknown')
        !write(15,*)  z, g(z)
       
    	!open(16,file='cosmic_v_k_Euclid_fis_test_7.txt' ,status='unknown')
    	!write(16,*) k_ov_h,  dsqrt(1d0/fistest)!
    	
        open(100,file='cosmic_variance_5.txt' ,status='unknown')
        write(100,*) z , ((1d0/P02)*dvdz2 ), (dndz*((3.14159d0/180d0)**2*dvdz2))
       !======================================================================
    k0_ov_h=k_ov_h
    enddo
    fistot=fistot+fis
    print*,'=== redshift bin#=(',ibin,') ==='
    print*, ' z = ', z
    print*,'Vsur =', Vsurvey,' h^-3 Mpc^3'
    print*,'ngal =',ngal,'h^3 Mpc^-3'
    print*,'1/ngal =',1d0/ngal,'h^-3 Mpc^3'
    print*,'Vsur * n(z)=', ngal* Vsurvey!,' h^-3 Mpc^3'
    print'(1A7,1F8.5)','Bias =',bias
    print'(1A7,1F8.5)','Beta =',beta
    print'(1A7,1F8.5)','g(z) =',g(z)
    print'(1A7,1F8.5)','gin a =',g(zin)/(1d0 +zin)
    print'(1A7,1F8.5)','f(z) =',1d0+dgdlna(z)/g(z)
    print'(1A7,1F8.5,1A9)','sigz =',sigma_z,' h^-1 Mpc'
    print'(1A7,1F8.5,1A9)','kmin =',kmin_ov_h,' h Mpc^-1'
    print'(1A7,1F8.5,1A9)','kmax =',kmax_ov_h,' h Mpc^-1'
    CALL report_result(z,bias,npara,fis)
    print*,''
    close(11)
    CALL transform_fisher(z,fis,fis6x6)
 enddo
10 close(2)
CALL close_linearpk
 if(nbins>1)then
     write(*,*) ,'=== combined ==='
     write(*,*) ,'Vsur =',Vtot,' h^-3 Mpc^3'
     write(*,*) , 'N(z)=',Ntot
     write(12,*),'=== combined ==='
     write(12,*),'Vsur =',Vtot,' h^-3 Mpc^3'
     write(12,*), 'N(z)=',Ntot
     CALL report_result(z, bias,npara,fistot)
  endif
  print*,''
 CALL report_result3x3(fis6x6)
 print*,''
END PROGRAM Fisher_Distance

!====================================================================
!===================! SUBROUTINES  ===================================

SUBROUTINE report_result(z,bias,npara,fis)
  IMPLICIT none
  integer, intent(IN) :: npara
  double precision, intent(IN) :: fis(npara,npara), z,bias
  double precision, allocatable, dimension(:,:) :: cov
  double precision, allocatable, dimension(:)   :: work
  double precision :: r12,err_lnda,err_lnh,err_lnR,beta, linear_pk,dgdlna,g
  integer :: i,j
  external linear_pk,dgdlna,g
  ALLOCATE(cov(2,2),work(2))
  cov=fis
  CALL DVERT(cov,2,2,work)
  beta=(1d0+dgdlna(z)/g(z))/bias
  r12=cov(1,2)/sqrt(cov(1,1)*cov(2,2))
  err_lnda=sqrt(cov(1,1))
  err_lnh=sqrt(cov(2,2))
  err_lnR=err_lnda*sqrt((1d0-r12**2d0) &
       /(1d0+2d0*r12*err_lnda/err_lnh+(err_lnda/err_lnh)**2d0))
  !==============save desired fisher files =======================
  !open(12,file='Fisher_bao_Real_150_Srms_corrected_nz_test.txt' ,status='unknown')
  open(12,file='Fisher_bao_0n_Srms_s30k.txt' ,status='unknown')
  !open(12,file='Fisher_bao_Euclid_s15k_3.txt' ,status='unknown')
  !open(12,file='Fisher_bao_Euclid_new_2.txt' ,status='unknown') 
  !write(12,*) 'z',  'fis(1,1)', 'fis(1,2)', 'fis(2,1)', 'fis(2,2)'
  !write(12,*) z,  fis(1,1), fis(1,2), fis(2,1), fis(2,2)
  write(12,'(4F18.5)') err_lnda,err_lnh,err_lnR
  !open(13,file='output_Fisher_bao_Pess_23_Srms_corrected_nz_2.txt', status='unknown')
  !open(13,file='output_Fisher_bao_0n_rms_s30k.txt' ,status='unknown')
  !open(13,file='output_Fisher_bao_Euclid_s15k_3.txt', status='unknown')
  !open(13,file='output_Fisher_bao_Euclid_new_2.txt', status='unknown')
  !write(13,'(6F18.5)') z, err_lnda*1d2,err_lnh*1d2,err_lnR*1d2, beta
  !========================================================
  print'(1A15,1F9.5)','Err[lnDa](%) =',err_lnda*1d2
  print'(1A15,1F9.5)','Err[lnH](%)  =',err_lnh*1d2
  print'(1A15,1F9.5)','r(lnDa,lnH)  =',r12
  print'(1A15,1F9.5)','Err[lnR](%)  =',err_lnR*1d2
  DEALLOCATE(cov,work)
  return
END SUBROUTINE report_result

!===========================================================
SUBROUTINE report_result3x3(fis)
  USE cosmo
  IMPLICIT none
  double precision, intent(IN) :: fis(6,6)
  double precision :: work(6),cov(6,6), A(6,6), M55DET, DET5x5,DET2x2,A05(6,6)
  integer :: i, j 
  cov=fis
 !==============================Calculate DET 2x2 ===============
  CALL DVERT(cov,6,6,work)
  A = fis
  DET2x2 =  A(1,1)*A(2,2) - A(1,2)*A(2,1)  
  write(12,*)'*** Error on w, with  Om0, ob0, ok0 and  H0 marginalized over ***'
  write(12,*) 'The Fisher Matrix = ', '[', 0.0, 0.0 ,    0.,  0., 0. , 0.,  0. ,  ';'&
    ,0., A(1,1) ,    A(1,2),  A(1,3) , A(1,4),  A(1,5), A(1,6),  ';'&
   ,0.,  A(2,1) ,  A(2,2),  A(2,3), A(2,4),  A(2,5), A(2,6), ';' &
   ,0.,  A(3,1),  A(3,2),  A(3,3), A(3,4),  A(3,5),  A(3,6),  ';'  &
   ,0.,  A(4,1),  A(4,2),  A(4,3), A(4,4),  A(4,5),  A(4,6), ';'    &
   , 0., A(5,1),  A(5,2),  A(5,3), A(5,4),  A(5,5),  A(5,6), ';'  &
   , 0., A(6,1),  A(6,2),  A(6,3), A(6,4),  A(6,5),  A(6,6), ']'     
   
 !====================write the resutls to file 12 ========================
  write(12,*),'Figure of Merit  =' ,  1d0/dsqrt((cov(1,1) * cov(2,2) - cov(1,2)* cov(2,1)))    
  write(12,'(1A20,1F9.5)')'Err[w] =',dsqrt(cov(1,1))
  write(12,'(1A20,1F9.5)')'Err[wa] =',dsqrt(cov(2,2))
  print*,'Figure of Merit  =' , (1d0/ dsqrt((cov(1,1) * cov(2,2) - cov(1,2)* cov(2,1))) )
  print'(1A20,1F9.5)', 'Err[w] =',dsqrt(cov(1,1))
  print'(1A20,1F9.5)', 'Err[wa] =',dsqrt(cov(2,2))
  write(12,'(1A20,1F9.5)')'Err[wa_w] =',(cov(1,2))
  write(12,'(1A20,1F9.5)')'Err[ Omega_m]=',(dsqrt(cov(5,5)))
  write(12,* )'Err[ Omega_b]=',(dsqrt(cov(6,6))* ob0)
  write(12, '(1A20,1F9.5)')'Err[Omega_m_w]=',(cov(1,5))
  write(12,'(1A20,1F9.5)')'Err[Omega_k]=',dsqrt(cov(3,3))
  write(12,'(1A20,1F9.5)')'Err[Omega_k_w]=',(cov(3,1))
 write(12,'(1A20,1F9.5)') 'Err[H_0] =',dsqrt(cov(4,4))
  return
END SUBROUTINE report_result3x3
!================================================================
SUBROUTINE transform_fisher(z,fisDH,fis6x6)
  USE cosmo
  IMPLICIT none
  integer :: a,b,i,j
  double precision, intent(IN)    :: z,fisDH(2,2)
  double precision, intent(INOUT) :: fis6x6(6,6)
  double precision :: dpdq(6,2)
  double precision :: chi,h2,func0,func1,func2,func3,func4,rombint,fz
  external h2,func0,func1,func2,func3,func4,rombint
  chi=rombint(func0,0d0,z,1d-7)
   fz =(1d0 + z)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(z/(1d0+ z)))
  dpdq(1,1)=-1.5d0*ode0*rombint(func1,0d0,z,1d-7)/chi         ! dlnDa/dw
  dpdq(1,2)= 1.5d0*ode0*dlog(1d0+z)*fz/h2(z)                       !dlnH/dw
  dpdq(2,1)=-1.5d0*ode0*rombint(func4,0d0,z,1d-7)/chi        !dlnDa/dwa
  dpdq(2,2)= 1.5d0*ode0*(dlog(1d0+z) - z/(1d0+ z))*fz/h2(z) !dlnH/dwa
  dpdq(3,1)=-0.5d0*ob0*rombint(func3,0d0,z,1d-7)/chi               !dlnDa/d(Omega_b)
  dpdq(3,2)= 0.5d0*ob0*((1d0+z)**3d0 - fz)/h2(z)                       !dlnH/dOmega_b
  dpdq(4,1)=-0.5d0*rombint(func2,0d0,z,1d-7)/chi+chi**2d0/6d0 !dlnDa/dOmega_k
  dpdq(4,2)= 0.5d0*((1d0+z)**2d0 - fz)/h2(z)                         !dlnH/dOmega_k
  dpdq(5,1)=-0.5d0*rombint(func3,0d0,z,1d-7)/chi                !dlnDa/d(Omega_m)
  dpdq(5,2)= 0.5d0*((1d0+z)**3d0 - fz)/h2(z)                        !dlnH/dOmega_m
  dpdq(6,1)= -100d0/H_0		                                      !dlnDa/dH_0
  dpdq(6,2)= 100d0/H_0		                                              !dlnH/dH_0

  do a=1,6
     do b=1,6
        do i=1,2
           do j=1,2
              ! transform and accumulate fis6x6
              fis6x6(a,b)=fis6x6(a,b)+dpdq(a,i)*dpdq(b,j)*fisDH(i,j) 
           enddo
        enddo
     enddo
  enddo
  return
END SUBROUTINE transform_fisher
!==============! Functions ==================================
DOUBLE PRECISION FUNCTION h2(redshift)
  USE cosmo
  ! h2(z) = Omega_matter(1+z)^3+Omega_lambda
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: fz
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  h2 =  ( (om0+ob0)*(1d0+redshift)**3d0+ ok0 * (1d0+redshift)**2d0+ode0* fz)
  return
END FUNCTION h2
DOUBLE PRECISION FUNCTION func0(redshift)
  USE cosmo
  ! func0(z) = 1/[h2(z)]^0.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2
  external :: h2
  func0 = 1d0/dsqrt(h2(redshift))
  return
END FUNCTION func0
DOUBLE PRECISION FUNCTION func1(redshift)
  USE cosmo
  ! func1(z) = ln(1+z)/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func1 = dlog(1d0+redshift) * fz /h2(redshift)**1.5d0
  return
END FUNCTION func1
DOUBLE PRECISION FUNCTION func2(redshift)
  USE cosmo
  ! func2(z) = (1+z)^2/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2,fz
  external :: h2
   fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func2 =( (1d0+redshift)**2d0 - fz)/h2(redshift)**1.5d0
  return
END FUNCTION func2
DOUBLE PRECISION FUNCTION func3(redshift)
  USE cosmo
  ! func3(z) = (1+z)^3/[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
   fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func3 =( (1d0+redshift)**3d0 - fz) /h2(redshift)**1.5d0
  return
END FUNCTION func3
DOUBLE PRECISION FUNCTION func4(redshift)
  USE cosmo
  ! func4(z) = ln(1+z) - z/1+z /[h2(z)]^1.5
  IMPLICIT none
  DOUBLE PRECISION, intent(IN) :: redshift
  DOUBLE PRECISION :: h2, fz
  external :: h2
  fz =(1d0 + redshift)**(3d0*(1d0 + w0+ w_a)) * exp(-3d0 * w_a *(redshift/(1d0+ redshift)))
  func4 = fz * (dlog(1d0+redshift) - (redshift/(1d0+ redshift)))/h2(redshift)**1.5d0
  return
END FUNCTION func4
