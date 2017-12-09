SUBROUTINE LAPACK_FULLEIGENVALUES(H,N,W_SPACE,INFO)

  IMPLICIT NONE
  INTEGER,                        INTENT(IN) :: N
  COMPLEX*16, DIMENSION(:,:),     INTENT(INOUT) :: H
!  COMPLEX*16, DIMENSION(:,:),     INTENT(OUT) :: Z
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: W_SPACE
  INTEGER,                        INTENT(OUT):: INFO
    

  !---SETTING  LAPACK VARIABLES: START ---------!
  CHARACTER         JOBZ, UPLO
  INTEGER           LWORK,LDA
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: IWORK,ISUPPZ
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: RWORK
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: WORK
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: Z,H_AUX  
  JOBZ = 'V'
  UPLO = 'L'  
 
  ALLOCATE(H_AUX(N,N))
  H_AUX = 0.0
  LDA   =  N
  LWORK = 2*N
!  ALLOCATE(W_SPACE(N))
  ALLOCATE(WORK(2*N))
  ALLOCATE(RWORK(3*N-2))
  !---- use zheev to get the optimun value of LWORK
  WORK = -1
  CALL ZHEEV(JOBZ,UPLO,N,H_AUX,LDA,W_SPACE,WORK,LWORK,RWORK,INFO)
  LWORK = INT(WORK(1))
  DEALLOCATE(WORK)
  ALLOCATE(WORK(LWORK))
  IF(INFO /= 0) WRITE(*,*) "DIAG FAIL"
  
  CALL ZHEEV(JOBZ,UPLO,N,H,LDA,W_SPACE, WORK, LWORK, RWORK,INFO)       
        
 ! DEALLOCATE(W_SPACE)
  DEALLOCATE(WORK)
  DEALLOCATE(RWORK) 
  DEALLOCATE(H_AUX)
  !DEALLOCATE(Z)

END SUBROUTINE LAPACK_FULLEIGENVALUES
      


SUBROUTINE LAPACK_ZGEEV(H,N,W_SPACE,INFO)
 
   IMPLICIT NONE
   INTEGER,                    INTENT(IN)    :: N
   COMPLEX*16, DIMENSION(:,:), INTENT(INOUT) :: H
   COMPLEX*16, DIMENSION(:),   INTENT(OUT)   :: W_SPACE
   INTEGER,                    INTENT(INOUT) :: INFO

! 
!   !---SETTING  LAPACK VARIABLES: START ---------!
!.. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            LDA, LDVL, LDVR, LWORK,i_
!*     ..
!*     .. Array Arguments ..
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
      COMPLEX*16, DIMENSION(:,:),ALLOCATABLE      :: VL,VR
      COMPLEX*16, DIMENSION(:),ALLOCATABLE        :: WORK


   JOBVL = 'N'
   JOBVR = 'N'

   LDA   =  N
   LDVL  =  1
   LDVR  =  N
   ALLOCATE(VL(LDVL,LDVL))
   ALLOCATE(VR(LDVR,LDVR))
   ALLOCATE(WORK(4*N))
   ALLOCATE(RWORK(2*N))
   !---- use zgeev to get the optimun value of LWORK
   LWORK = -1
   CALL ZGEEV( JOBVL, JOBVR, N, H, LDA, W_SPACE, VL, LDVL, VR, LDVR,&
     &                  WORK, LWORK, RWORK, INFO )
   IF(INFO /= 0) WRITE(*,*) "LWORK REQUEST FAILED"
   LWORK = INT(REAL(WORK(1)))
   DEALLOCATE(WORK)
   ALLOCATE(WORK(LWORK))
   CALL ZGEEV( JOBVL, JOBVR, N, H, LDA, W_SPACE, VL, LDVL, VR, LDVR,&
     &                  WORK, LWORK, RWORK, INFO )
   H = VR
   IF(INFO /= 0) WRITE(*,*) "DIAG FAIL"
   !WRITE(*,*) INFO,LWORK,abs(WORK(1))


   DEALLOCATE(WORK)
   DEALLOCATE(RWORK) 
   DEALLOCATE(VL)
   DEALLOCATE(VR)
   
 END SUBROUTINE LAPACK_ZGEEV

 
!!$ SUBROUTINE LAPACK_ZGEEV_GPU(H,N,W,INFO)
!!$   
!!$   USE MAGMA
!!$   IMPLICIT NONE
!!$   INTEGER,                  INTENT(IN)    :: N
!!$   COMPLEX*16, DIMENSION(:), INTENT(INOUT) :: H
!!$   COMPLEX*16, DIMENSION(:), INTENT(OUT)   :: W
!!$   INTEGER,                  INTENT(INOUT) :: INFO
!!$   
!!$   ! 
!!$   !   !---SETTING  LAPACK VARIABLES: START ---------!
!!$   !.. Scalar Arguments ..
!!$   CHARACTER          JOBVL, JOBVR
!!$   INTEGER            LDA, LDVL, LDVR, LWORK,i_
!!$   !*     ..
!!$   !*     .. Array Arguments ..
!!$   DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
!!$   COMPLEX*16, DIMENSION(:),ALLOCATABLE        :: VL,VR
!!$   COMPLEX*16, DIMENSION(:),ALLOCATABLE        :: WORK
!!$   
!!$   call cublas_init()
!!$   
!!$   JOBVL = 'N'
!!$   JOBVR = 'N'
!!$   
!!$   LDA   =  N
!!$   LDVL  =  1
!!$   LDVR  =  N
!!$   
!!$   ALLOCATE(VL(LDVL*LDVL))
!!$   ALLOCATE(VR(LDVR*LDVR))
!!$   ALLOCATE(WORK(1))
!!$   ALLOCATE(RWORK(2*N))
!!$   !---- use magma_zgeev to get the optimun value of LWORK
!!$   LWORK = -1
!!$   call magmaf_zgeev(JOBVL,JOBVR, N, H, N, W, VL,LDVL,&
!!$        & VR,LDVR, WORK, LWORK,RWORK,INFO) 
!!$   IF(INFO /= 0) WRITE(*,*) "LWORK REQUEST FAILED"
!!$   LWORK = INT(REAL(WORK(1)))
!!$   DEALLOCATE(WORK)
!!$   ALLOCATE(WORK(LWORK))
!!$   call magmaf_zgeev(JOBVL,JOBVR, N, H, N, W, VL,LDVL,&
!!$        & VR,LDVR, WORK, LWORK,RWORK,INFO) 
!!$   H = VR
!!$   IF(INFO /= 0) WRITE(*,*) "DIAG FAIL"
!!$   
!!$   
!!$   DEALLOCATE(WORK)
!!$   DEALLOCATE(RWORK) 
!!$   DEALLOCATE(VL)
!!$   DEALLOCATE(VR)
!!$   
!!$ END SUBROUTINE LAPACK_ZGEEV_GPU
