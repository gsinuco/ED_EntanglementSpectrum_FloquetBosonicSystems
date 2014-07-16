! En este programa calculamos de manera exacta el spectro de Floquet
! en un tigh-binding Hamiltonian PRL...
! Single particle!


SUBROUTINE FLOQUET_BOSE_HUBBARD_SUB(q_,U_T,U_FTY, EIGENVALUES_H,INFO)

  USE physical_constants  
  USE LATTICE
  USE INTERACTION
  USE TIME_DEPENDENT
  USE INTERFACE

  IMPLICIT NONE

  COMPLEX*16, DIMENSION(:,:),INTENT(OUT):: U_T
  COMPLEX*16, DIMENSION(:), INTENT(OUT) :: U_FTY
  DOUBLE PRECISION, DIMENSION(:),INTENT(OUT):: EIGENVALUES_H
  INTEGER, INTENT(IN) :: q_
  INTEGER, INTENT(INOUT) :: INFO

  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: U_X, U_Y,U_AUX
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: U,EIGENVALUES
  !DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: EIGENVALUES_H

  INTEGER BASIS_DIM,i,i_,k_,q
  DOUBLE PRECISION D_H
  COMPLEX*16 :: z_Y,z_X

  !OPEN(UNIT=1,FILE='SPECTRUMb.dat',ACTION='WRITE')
  !OPEN(UNIT=2,FILE='EIGENVECTORSb.dat',ACTION='WRITE')

  BASIS_DIM =  INT(D_H(N_SITES,N_PARTICLES))
  !GAMMA(1.0*N_PARTICLES+1.0*N_SITES_X*N_SITES_Y-1.0+1.0)/&
  !& (GAMMA(1.0*N_PARTICLES+1)*GAMMA(1.0*N_SITES_X*N_SITES_Y-1.0+1.0))

  ALLOCATE(U_X(N_SITES_X,N_SITES_X))
  ALLOCATE(U_Y(N_SITES_X,N_SITES_X))
  !ALLOCATE(U_T(N_SITES_X,N_SITES_X))
  ALLOCATE(U_AUX(N_SITES_X,N_SITES_X))
  !ALLOCATE(U_FTY(N_SITES_Y))
  ALLOCATE(EIGENVALUES(N_SITES_X))
  !ALLOCATE(EIGENVALUES_H(N_SITES_X))
  
  INFO  = 0

  z_X       = DCMPLX(0.0,-1.0)*tau
  z_Y       = DCMPLX(0.0,-1.0)*(T-tau)
  U_X       = DCMPLX(0.0,0.0)
  U_AUX     = DCMPLX(0.0,0.0)

  !The Hamiltonian H_X and its diagonalization
  DO i_=1,N_SITES_X-1        
     U_AUX(i_,i_+1)  = -J_X
  END DO
  !U_AUX(1,N_SITES_X) = -J_X  ! periodic Boundary conditions
  U_AUX = U_AUX + transpose(conjg(U_AUX))
  CALL LAPACK_FULLEIGENVALUES(U_AUX,N_SITES_X,EIGENVALUES_H,INFO)

!  CALL MANY_BODY_SPECTRUM(U_AUX,EIGENVALUES_H,U_MB,E_MB,INFO)

  !Evolution operator U(0,tau) =  exp(-ihbar tau H_X)
  DO i_=1,N_SITES_X	
     U_X(i_,i_) = exp(z_X*EIGENVALUES_H(i_))
  END DO

  
  DO q=q_,q_!-N_SITES_Y/2,N_SITES_Y/2-1

     U_T       = DCMPLX(0.0,0.0)
     U_Y       = DCMPLX(0.0,0.0)
       
     DO i_=1,N_SITES_X
        U_Y(i_,i_) = exp(z_Y*2.0*cos(2.0*pi*i_*alpha + 2*pi*q/N_SITES_Y))
     END DO
     U_Y = MATMUL(TRANSPOSE(CONJG(U_AUX)),MATMUL(U_Y,U_AUX))
     U_T = MATMUL(U_Y,U_X)      

     CALL LAPACK_ZGEEV(U_T,N_SITES_X,EIGENVALUES,INFO)
     !WRITE(2,2064) U_T ! EIGENVECTORS STORED AS COLUMNS OF U_T; U_T(:,1) IN THE MOMENTUM BASIS |k,q>

     EIGENVALUES_H = 0.0
     CALL EV_QUASIENERGY(EIGENVALUES,EIGENVALUES_H,SIZE(EIGENVALUES_H),INFO)
     !WRITE(1,30128) EIGENVALUES_H

      U_T   = MATMUL(U_AUX,U_T)!   see notes!   
      DO i_=1,N_SITES_Y
        U_FTY(i_) = exp(-DCMPLX(0.0,1.0)*2*pi*q*i_/N_SITES_Y)/SQRT(1.0*N_SITES_Y) ! phase factor for the Fourier Transform along y direction
      END DO

  END DO
201 format(2e15.6)
208 format(8e15.6,8e15.6)
2064 format(32e15.6,32e15.6)
30128 format(128e15.6)  
301 format(1e15.6)  
END SUBROUTINE FLOQUET_BOSE_HUBBARD_SUB

! SUBROUTINE LAPACK_FULLEIGENVALUES(H,N,W_SPACE,INFO)
! 
!   IMPLICIT NONE
!   INTEGER,                        INTENT(IN) :: N
!   COMPLEX*16, DIMENSION(:,:),     INTENT(INOUT) :: H
!   COMPLEX*16, DIMENSION(:), INTENT(OUT) :: W_SPACE
!   INTEGER,                        INTENT(OUT):: INFO
!     
! 
!   !---SETTING  LAPACK VARIABLES: START ---------!
!   CHARACTER         JOBZ, UPLO
!   INTEGER           LWORK,LDA
!   INTEGER,          DIMENSION(:),   ALLOCATABLE :: IWORK,ISUPPZ
!   DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: RWORK
!   COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: WORK
!   COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_AUX  
!   JOBZ = 'V'  ! NOTE THAT COLUMNS IN RESULTING H ARE THE EIGENVECTORS H(:.1)
!   UPLO = 'L'  
! 
!   ALLOCATE(H_AUX(N,N))
!   H_AUX = 0.0
!   LDA   =  N
!   LWORK = 2*N
!   ALLOCATE(WORK(2*N))
!   ALLOCATE(RWORK(3*N-2))
!   !---- use zheev to get the optimun value of LWORK
!   WORK = -1
!   CALL ZHEEV(JOBZ,UPLO,N,H_AUX,LDA,W_SPACE,WORK,LWORK,RWORK,INFO)
!   LWORK = INT(WORK(1))
!   DEALLOCATE(WORK)
!   ALLOCATE(WORK(LWORK))
!   IF(INFO /= 0) WRITE(*,*) "DIAG FAIL"
!   DEALLOCATE(H_AUX)  
!   CALL ZHEEV(JOBZ,UPLO,N,H,LDA,W_SPACE, WORK, LWORK, RWORK,INFO)       
! 	
!   DEALLOCATE(WORK)
!   DEALLOCATE(RWORK) 
!   
! END SUBROUTINE LAPACK_FULLEIGENVALUES
! 
! 
! SUBROUTINE LAPACK_ZGEEV(H,N,W_SPACE,INFO)
!  
!    IMPLICIT NONE
!    INTEGER,                        INTENT(IN)    :: N
!    COMPLEX*16, DIMENSION(:,:),     INTENT(INOUT) :: H
!    DOUBLE PRECISION, DIMENSION(:), INTENT(OUT)   :: W_SPACE
!    INTEGER,                        INTENT(INOUT)   :: INFO
! 
! ! 
! !   !---SETTING  LAPACK VARIABLES: START ---------!
! !.. Scalar Arguments ..
!       CHARACTER          JOBVL, JOBVR
!       INTEGER            LDA, LDVL, LDVR, LWORK,i_
! !*     ..
! !*     .. Array Arguments ..
!       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
!       COMPLEX*16, DIMENSION(:,:),ALLOCATABLE      :: VL,VR
!       COMPLEX*16, DIMENSION(:),ALLOCATABLE        :: WORK
! 
! 
!    JOBVL = 'N'
!    JOBVR = 'V'
! 
!    LDA   =  N
!    LDVL  =  1
!    LDVR  =  N
!    ALLOCATE(VL(LDVL,LDVL))
!    ALLOCATE(VR(LDVR,LDVR))
!    ALLOCATE(WORK(4*N))
!    ALLOCATE(RWORK(2*N))
!    !---- use zgeev to get the optimun value of LWORK
!    LWORK = -1
!    CALL ZGEEV( JOBVL, JOBVR, N, H, LDA, W_SPACE, VL, LDVL, VR, LDVR,&
!      &                  WORK, LWORK, RWORK, INFO )
!    IF(INFO /= 0) WRITE(*,*) "LWORK REQUEST FAILED"
!    LWORK = INT(REAL(WORK(1)))
!    DEALLOCATE(WORK)
!    ALLOCATE(WORK(LWORK))
!    CALL ZGEEV( JOBVL, JOBVR, N, H, LDA, W_SPACE, VL, LDVL, VR, LDVR,&
!      &                  WORK, LWORK, RWORK, INFO )
!    H = VR
!    IF(INFO /= 0) WRITE(*,*) "DIAG FAIL"
!    !WRITE(*,*) INFO,LWORK,abs(WORK(1))
! 
! 
!    DEALLOCATE(WORK)
!    DEALLOCATE(RWORK) 
!    
!  END SUBROUTINE LAPACK_ZGEEV

SUBROUTINE EV_QUASIENERGY(EV,QE,N,INFO)

  IMPLICIT NONE
  INTEGER,                       INTENT(IN)    :: N
  INTEGER,                       INTENT(INOUT) :: INFO
  COMPLEX*16, DIMENSION(N),      INTENT(IN)    :: EV
  DOUBLE PRECISION, DIMENSION(N),INTENT(INOUT)   :: QE
   
  DOUBLE PRECISION x,z,theta
  INTEGER i_

  DO i_ = 1,N

    z = REAL(EV(i_))
    x = AIMAG(EV(i_))
    !write(*,*) z,x
    if(x.gt.0 .and. z.gt.0) theta = atan(x/z)
    if(x.gt.0 .and. z.lt.0) theta = atan(x/z) + 4.0*ATAN(1.0)
    if(x.lt.0 .and. z.lt.0) theta = atan(x/z) + 4.0*ATAN(1.0)
    if(x.lt.0 .and. z.gt.0) theta = atan(x/z) + 8.0*ATAN(1.0)
    if(z.eq.0 .and. x.gt.0) theta = 2.0*ATAN(1.0) 
    if(z.eq.0 .and. x.lt.0) theta = 6.0*ATAN(1.0)
    if(x.eq.0 .and. z.eq.0) theta = 0.0
  
    QE(i_) = theta

  END DO
END SUBROUTINE EV_QUASIENERGY
