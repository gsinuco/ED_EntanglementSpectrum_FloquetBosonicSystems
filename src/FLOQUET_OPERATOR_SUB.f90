! En este programa calculamos de manera exacta el spectro de Floquet
! en un tigh-binding Hamiltonian PRL...
! Single particle!

SUBROUTINE FLOQUET_BOSE_HUBBARD_SUB(STATES,Q_T,U_FTY,U_MB,EIGENVALUES,INFO)

#ifdef CUBLAS  
  USE MAGMA
#endif
  USE physical_constants  
  USE LATTICE
  USE INTERACTION
  USE TIME_DEPENDENT
  USE INTERFACE
  USE STATE_OBJ

  IMPLICIT NONE
  TYPE(STATE),      DIMENSION(:),                INTENT(IN)    :: STATES
  INTEGER,                                       INTENT(INOUT) :: Q_T
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE, INTENT(OUT)   :: EIGENVALUES
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)   :: U_MB
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE, INTENT(OUT)   :: U_FTY
  INTEGER,                                       INTENT(INOUT) :: INFO
  

  INTEGER BAND_INDEX
  
  COMPLEX*16, DIMENSION(:),ALLOCATABLE :: EIGENVALUES_AUX
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE   :: E_MBX,E_MBY,EIGENVALUES_H,rho,norm,theta
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE   :: U_X, U_Y,U_AUX,U_FLOQUET
!  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE   :: U_F
  INTEGER,          DIMENSION(:),   ALLOCATABLE   :: EV_INDEX
  
  INTEGER BASIS_DIM,i_,k_,q,q_,D, j, i_off,j_
  DOUBLE PRECISION D_H
  COMPLEX*16 :: z_Y,z_X, alpha_,beta_

  
#ifdef CUBLAS

#ifdef CUBLAS_USE_THUNKING
!    write(*,*) "CUBLAS, USING THUNKING"
  OPEN(UNIT=1,FILE='/home/g/gs/gs249/Programs/ED-EntanglementSpectrum-FloquetV2/data/SP64x128_thunkingB.dat',ACTION='WRITE')
#else
  integer*8 devPtrUFLOQUET, devPtrUMB, devPtrUY,devPtrUX,devPtrEIGENVALUES
  OPEN(UNIT=1,FILE='/home/g/gs/gs249/Programs/ED-EntanglementSpectrum-FloquetV2/data/SP64x128_gpuB.dat',ACTION='WRITE')
#endif

  call cublas_init()

#else
  OPEN(UNIT=1,FILE='/home/g/gs/gs249/Programs/ED-EntanglementSpectrum-FloquetV2/data/SP128x128_cpuB.dat',ACTION='WRITE')
#endif
  
  BASIS_DIM =  INT(D_H(N_SITES,N_PARTICLES))   !GAMMA(1.0*N_PARTICLES+1.0*N_SITES_X*N_SITES_Y-1.0+1.0)/(GAMMA(1.0*N_PARTICLES+1)*GAMMA(1.0*N_SITES_X*N_SITES_Y-1.0+1.0))
  
  ALLOCATE(U_AUX(N_SITES_X,N_SITES_X))
  ALLOCATE(EIGENVALUES_H(N_SITES_X))
  ALLOCATE(rho(N_SITES_X))
  
  
  z_X       = DCMPLX(0.0,-1.0)*tau
  z_Y       = DCMPLX(0.0,-1.0)*(T-tau)
  
  ALPHA_ = DCMPLX(1.0,0.0)
  BETA_  = DCMPLX(0.0,0.0)
  k_     = 0
  
  
  INFO  = 0
  U_AUX     = DCMPLX(0.0,0.0)
  EIGENVALUES_H = 0.0
  !Diagonalization of H_x:
  DO i_=1,N_SITES_X-1        
     U_AUX(i_,i_+1)  = -J_X
  END DO
  U_AUX(1,N_SITES_X) = -J_X  ! add periodic boundary conditions
  U_AUX = U_AUX + transpose(conjg(U_AUX))
  CALL LAPACK_FULLEIGENVALUES(U_AUX,N_SITES_X,EIGENVALUES_H,INFO)


  !write(*,*) EIGENVALUES_h
  IF((N_PARTICLES.GT.1) .AND. (INFO.EQ.0)) THEN
     
     ALLOCATE(U_FTY(N_PARTICLES))
     q  = 1
     q_ = 2
     D  =  (N_SITES_X*N_SITES_X  + N_SITES_X)/2
     ALLOCATE(U_MB(D,D))
     U_MB = 0
     CALL MANY_BODY_SPECTRUM(STATES,U_AUX,q,q_,U_MB,D,INFO)
     U_MB = CONJG(U_MB)
     !write(*,*) abs(U_MB(1,1))
     ALLOCATE(E_MBX(D))
     ALLOCATE(E_MBY(D))
     E_MBX = 0.0
     E_MBY = 0.0
#ifdef CUBLAS
#ifdef CUBLAS_USE_THUNKING
 !   write(*,*) "CUBLAS, USING THUNKING"
#else
 !    write(*,*) "setting u_mb"
     call cublas_alloc(D*D,sizeof_complex_16,devPtrUMB)
     call cublas_set_matrix(D,D,sizeof_complex_16,U_MB,D,devPtrUMB,D)
#endif
#endif

!     DO Q_T = 1,1!N_SITES_Y!(N_SITES_Y*N_SITES_Y + N_SITES_Y)/2 !loop over blocks of total momentum_Y
     DO q  = 2,2!  2,N_SITES_Y!(N_SITES_Y*N_SITES_Y + N_SITES_Y)/2 !loop over blocks of total momentum_Y
     DO q_ = 67,67!q+1,N_SITES_Y!(N_SITES_Y*N_SITES_Y + N_SITES_Y)/2 !loop over blocks of total momentum_Y

        j = 0
        i_off = 0

        IF(q.EQ.q_) THEN
           i_off = 0
        ELSE
           i_off = (N_SITES_X*N_SITES_X  + N_SITES_X)/2
        END IF

        CALL  TWO_BODY_ENERGY_SPECTRUM(STATES,EIGENVALUES_H,D,i_off,q,q_,E_MBX,E_MBY,INFO)
        IF(INFO .NE. 0) WRITE(*,*) 'THERE WAS A PROBLEM WITH MANY_BODY_SPECTRUM: CHECK STATES OF EQUAL TOTAL MOMENTUM'

        !Definition of partial Floquet operators,
        !write(*,*) D

        !--------------------------------------------
        !------------TESTING--BEGIN------------------
        !--------------------------------------------
        ! {|eigen-basis>_i} = U_MB(:,i) {|old>}
        !
        !     ! {|new>} = U^dagger {|old>}
        ! operators in the new basis:U^dagger O U
        ! operators in the old basis:U O U^dagger
        ! U^dagger = trasnspose(U_MB)
        ! U        = conjugate(U_MB)
        !     !DO i_=1,D
        !     !   U_X(i_,i_) = E_MBX(i_) ! in its eigen-basis
        !     !END DO
        !     !     DO i_=1,D
        !!        write(*,*) REAL(U_X(i_,:))
        !!     END DO
        !!     write(*,*)
        !     U_X = MATMUL(TRANSPOSE(U_MB),CONJG(U_MB))
        !     DO i_=1,D
        !        write(*,*) abs(U_X(i_,:))
        !     END DO
        !     write(*,*)
        !     write(*,*)
        !
        !!     CALL LAPACK_ZGEEV(U_X,D,EIGENVALUES,INFO)
        !!     DO i_=1,D
        !!        WRITE(*,*)  EIGENVALUES(i_)
        !!     END DO
        !--------------------------------------------
        !------------TESTING--END--------------------
        !--------------------------------------------

        ALLOCATE(U_Y(D,D))
        U_Y = 0
        DO i_=1,D
           U_Y(i_,i_) = exp(z_Y*E_MBY(i_)) ! In its eigen-basis
        END DO

        ALLOCATE(U_FLOQUET(D,D))
        U_FLOQUET = 0

           !write(*,*)
           !write(*,*)
#ifdef CUBLAS

        ALPHA_ = DCMPLX(1.0,0.0)
        BETA_  = DCMPLX(0.0,0.0)

#ifdef CUBLAS_USE_THUNKING

  !      write(*,*) "CUBLAS, USING THUNKING"
        call cublas_zGEMM('n','n',D,D,D,alpha_,U_Y, D,U_MB,D,beta_,U_FLOQUET,D) ! U_FLOQUET IN THIS CASE IS AN AUXILIARY MATRIX
        call cublas_zGEMM('c','n',D,D,D,alpha_,U_MB,D,U_FLOQUET,  D,beta_,U_Y,      D)
!        U_Y = MATMUL(TRANSPOSE(CONJG(U_MB)),MATMUL(U_Y,U_MB)) ! transformation to the eigen-basis of U_X, |kq,k'q'>

#else
  !      write(*,*) 'no thunking',D

        call cublas_alloc(D*D,sizeof_complex_16,devPtrUY)
        call cublas_alloc(D*D,sizeof_complex_16,devPtrUFLOQUET)

!        call cublas_alloc(D*D,sizeof_complex_16,devPtrUMB)
!        call cublas_set_matrix(D,D,sizeof_complex_16,U_MB,D,devPtrUMB,D)

        call cublas_set_matrix(D,D,sizeof_complex_16,U_Y,D,devPtrUY,D)
        call cublas_set_matrix(D,D,sizeof_complex_16,U_FLOQUET,D,devPtrUFLOQUET,D)

        call cublas_zGEMM('n','n',D,D,D,alpha_,devPtrUY, D,devPtrUMB,       D,beta_,devPtrUFLOQUET,D) ! U_FLOQUET IN THIS CASE IS AN AUXILIARY MATRIX
        call cublas_zGEMM('c','n',D,D,D,alpha_,devPtrUMB,D,devPtrUFLOQUET,  D,beta_,devPtrUY,      D)
        call cublas_get_matrix(D,D,sizeof_complex_16,devPtrUY,D,U_Y,D)
!        U_Y = MATMUL(TRANSPOSE(CONJG(U_MB)),MATMUL(U_Y,U_MB)) ! transformation to the eigen-basis of U_X, |kq,k'q'>

#endif


#else
        U_Y = MATMUL(TRANSPOSE(CONJG(U_MB)),MATMUL(U_Y,U_MB)) ! transformation to the eigen-basis of U_X, |kq,k'q'>
#endif

!        ALLOCATE(U_X(D,D))
!        U_X = 0
!        DO i_=1,D
!           U_X(i_,i_) = exp(z_X*E_MBX(i_)) ! in its eigen-basis
!        END DO

        U_FLOQUET = DCMPLX(0.0,0.0)


!#ifdef CUBLAS
!        ALPHA_ = DCMPLX(1.0,0.0)
!        BETA_  = DCMPLX(0.0,0.0)
!#ifdef CUBLAS_USE_THUNKING
!        !call cublas_zGEMM('n','n',D,D,D,alpha_,U_Y,D,U_X,D,beta_,U_FLOQUET,D)
!        DO i_=1,D
!            U_FLOQUET(:,i_) = exp(z_X*E_MBX(i_))*U_Y(:,i_)
!        END DO
!
!
!#else
!   !     write(*,*) 'no thunking',D
!
!!        call cublas_alloc(D*D,sizeof_complex_16,devPtrUX)
!!        call cublas_set_matrix(D,D,sizeof_complex_16,U_X,D,devPtrUX,D)
!!        call cublas_zGEMM('n','n',D,D,D,alpha_,devPtrUY,D,devPtrUX,D,beta_,devPtrUFLOQUET,D)
!!        call cublas_free(devPtrUX)
!!        call cublas_free(devPtrUY)
!!        call cublas_alloc(D,sizeof_complex_16,devPtrEIGENVALUES)
!!        call cublas_get_matrix(D,D,sizeof_complex_16,devPtrUFLOQUET,D,U_FLOQUET,D)
!
!        DO i_=1,D
!            U_FLOQUET(:,i_) = exp(z_X*E_MBX(i_))*U_Y(:,i_)
!        END DO
!
!#endif
!
!#else
!        !U_FLOQUET = MATMUL(U_Y,U_X) ! Total Floquet operator in the eigen-basis of H_X (equiv. U_X)
!
!        DO i_=1,D
!            U_FLOQUET(:,i_) = exp(z_X*E_MBX(i_))*U_Y(:,i_)
!        END DO
!
!#endif
        DO i_=1,D
            U_FLOQUET(:,i_) = exp(z_X*E_MBX(i_))*U_Y(:,i_)
        END DO
        !write(*,*) '# Lapack starts',Q_T
        !DEALLOCATE(U_X)
        DEALLOCATE(U_Y)
        ALLOCATE(EIGENVALUES(D))
        EIGENVALUES = 0.0


#ifdef CUBLAS

#ifdef CUBLAS_USE_THUNKING
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!! USING MAGMA
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        ALLOCATE(U_F(D*D))
!        U_F = 0.0
!        DO i_=1,D
!           DO k_=1,D
!              U_F(i_ + (k_-1)*D) = U_FLOQUET(k_,i_)
!           END DO
!        END DO
!
!        CALL LAPACK_ZGEEV_GPU_THUNKING(U_F,D,EIGENVALUES,INFO)!Eigenvectors of the full Floquet operator, stored as U_MB(:,i) in the {|kq,k'q'>}
!
!        DO i_=1,D
!           DO k_=1,D
!              U_FLOQUET(k_,i_) = U_F(i_ + (k_-1)*D)
!           END DO
!        END DO
!        DEALLOCATE(U_F)

        CALL LAPACK_ZGEEV_GPU_NOTHUNKING(U_FLOQUET,D,EIGENVALUES,INFO)!Eigenvectors of the full Floquet operator, stored as U_MB(:,i) in the {|kq,k'q'>}
    ! CALL LAPACK_ZGEEV(U_FLOQUET,D,EIGENVALUES,INFO) !Eigenvectors of the full Floquet operator, stored as U_MB(:,i) in the {|kq,k'q'>}

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#else
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!! USING MAGMA

        CALL LAPACK_ZGEEV_GPU_NOTHUNKING(U_FLOQUET,D,EIGENVALUES,INFO)!Eigenvectors of the full Floquet operator, stored as U_MB(:,i) in the {|kq,k'q'>}
  !      CALL LAPACK_ZGEEV(U_FLOQUET,D,EIGENVALUES,INFO) !Eigenvectors of the full Floquet operator, stored as U_MB(:,i) in the {|kq,k'q'>}
  !      call cublas_free(devPtrUFLOQUET)
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#endif

#else
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!! USING LAPACK
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL LAPACK_ZGEEV(U_FLOQUET,D,EIGENVALUES,INFO) !Eigenvectors of the full Floquet operator, stored as U_MB(:,i) in the {|kq,k'q'>}
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#endif



        ALLOCATE(EV_INDEX(D))
        EV_INDEX = -10000
        CALL EV_QUASIENERGY(EIGENVALUES,EV_INDEX,SIZE(EIGENVALUES),INFO) ! The eigenvalues of U_MB are complex values. This function evaluates the norm of the eigenvalues
        DO i_=1,D
           WRITE(1,*)    Q_T,q,q_,REAL(EIGENVALUES(EV_INDEX(i_))),EV_INDEX(i_)
        END DO
        WRITE(1,*)
        WRITE(1,*)
        DEALLOCATE(EV_INDEX)
        DEALLOCATE(EIGENVALUES,U_FLOQUET)




!
!        !{states(i_off+1: i_off+D)} == {|nq,n'q'>}
!        !     U_MB = MATMUL(TRANSPOSE(U_FLOQUET),TRANSPOSE(U_MB))!The eigenvectors in the basis !nq,n'q'>, stored as ROWS U_MB(|nq;n'q'>, :)
!        !     U_MB = TRANSPOSE(U_MB) !The eigenvectors in the basis |nq,n'q'>, stored as COLUMNS U_X(:, |nq;n'q'>)
!        !     U_FTY(1) = q
!        !     U_FTY(2) = q_
!
!
!        !     IF (q .NE. q_) THEN
!        !         DO k_=1,N_SITES_X
!        !         DO j =1,N_SITES_X
!        !            i_ = j + (k_-1)*N_SITES_X
!        !            WRITE(*,4051) i_off+i_, STATES(i_off+i_)%n(1), STATES(i_off+i_)%m(1), STATES(i_off+i_)%n(2), &
!        !                & STATES(i_off+i_)%m(1),abs(U_MB(i_,EV_INDEX(849)))
!!         END DO
!        !            WRITE(*,*)
!        !         END DO
!        !         WRITE(*,*)
!        !         WRITE(*,*)
!        !
!        !         DO k_=1,N_SITES_X
!        !         DO j =1,N_SITES_X
!        !            i_ = j + (k_-1)*N_SITES_X
!        !            WRITE(*,4051) i_off+i_, STATES(i_off+i_)%n(1), STATES(i_off+i_)%m(1), STATES(i_off+i_)%n(2), &
!        !                & STATES(i_off+i_)%m(1),abs(U_MB(i_,EV_INDEX(850)))
!        !         END DO
!        !            WRITE(*,*)
!        !         END DO
!        !         WRITE(*,*)
!        !         WRITE(*,*)
!        !
!        !
!        !         DO k_=1,N_SITES_X
!        !         DO j =1,N_SITES_X
!        !            i_ = j + (k_-1)*N_SITES_X
!        !            WRITE(*,4051) i_off+i_, STATES(i_off+i_)%n(1), STATES(i_off+i_)%m(1), STATES(i_off+i_)%n(2), &
!        !                & STATES(i_off+i_)%m(1),abs(U_MB(i_,EV_INDEX(869)))
!        !         END DO
!        !            WRITE(*,*)
!        !         END DO
!        !         WRITE(*,*)
!        !         WRITE(*,*)
!        !
!        !         DO k_=1,N_SITES_X
!        !         DO j =1,N_SITES_X
!        !            i_ = j + (k_-1)*N_SITES_X
!        !            WRITE(*,4051) i_off+i_, STATES(i_off+i_)%n(1), STATES(i_off+i_)%m(1), STATES(i_off+i_)%n(2), &
!        !                & STATES(i_off+i_)%m(1),abs(U_MB(i_,EV_INDEX(870)))
!        !         END DO
!        !            WRITE(*,*)
!        !         END DO
!        !         WRITE(*,*)
!        !         WRITE(*,*)
!        !
!        !     ELSE
!        !
!        !         i_ = i_off
!        !         DO j  = 1,N_SITES_X
!        !         DO k_ = j,N_SITES_X
!        !            i_ = i_ + 1
!        !         !DO i_ =1,D
!        !            WRITE(*,4051) i_off+i_, STATES(i_off+i_)%n(1), STATES(i_off+i_)%m(1), STATES(i_off+i_)%n(2), &
!        !                & STATES(i_off+i_)%m(1),abs(U_MB(i_,EV_INDEX(7)))
!        !         END DO
!        !            WRITE(*,*)
!        !         END DO
!        !         WRITE(*,*)
!        !         WRITE(*,*)
!        !
!        !
!        !     END IF
!        !     CALL oneD_DENSITY(STATES,i_off,U_MB(:,EV_INDEX(849)),rho,INFO)
!        !     WRITE(*,201) rho
!        !     WRITE(*,*)
!        !     WRITE(*,*)
!        !
!        !     CALL oneD_DENSITY(STATES,i_off,U_MB(:,EV_INDEX(850)),rho,INFO)
!        !     WRITE(*,201) rho
!        !     WRITE(*,*)
!        !     WRITE(*,*)
!        !
!        !     CALL oneD_DENSITY(STATES,i_off,U_MB(:,EV_INDEX(869)),rho,INFO)
!        !     WRITE(*,201) rho
!        !     WRITE(*,*)
!        !     WRITE(*,*)
!        !
!        !     CALL oneD_DENSITY(STATES,i_off,U_MB(:,EV_INDEX(870)),rho,INFO)
!        !     WRITE(*,201) rho
!        !     WRITE(*,*)
!        !     WRITE(*,*)
!        !     DEALLOCATE(U_FLOQUET)
!        !     DEALLOCATE(EIGENVALUES)
!        !#ifdef CUBLAS
!        !     DEALLOCATE(U_F)
!        !#endif
!        !     DEALLOCATE(EV_INDEX)
!
     END DO
     END DO

     DEALLOCATE(U_FTY)
     DEALLOCATE(U_MB)
     DEALLOCATE(E_MBX)
     DEALLOCATE(E_MBY)

#ifdef CUBLAS
#ifdef CUBLAS_USE_THUNKING
#else
     call cublas_free(devPtrUMB)
#endif
     call cublas_shutdown()
#endif

  ELSE IF(INFO.EQ.0) THEN

     ALLOCATE(U_FTY(N_SITES_Y))
     ALLOCATE(U_X(N_SITES_X,N_SITES_X))
     ALLOCATE(U_Y(N_SITES_X,N_SITES_X))
     ALLOCATE(EIGENVALUES(N_SITES_X))
     ALLOCATE(EIGENVALUES_AUX(N_SITES_X))
     ALLOCATE(U_MB(N_SITES_X,N_SITES_X))
     ALLOCATE(EV_INDEX(N_SITES_X))
     ALLOCATE(norm(N_SITES_X))
     ALLOCATE(theta(N_SITES_X))

     z_X       = DCMPLX(0.0,-1.0)*tau
     z_Y       = DCMPLX(0.0,-1.0)*(T-tau)

     !Evolution operator U(0,tau) =  exp(-ihbar tau H_X)
     U_X  = DCMPLX(0.0,0.0)

     DO i_=1,N_SITES_X	
        U_X(i_,i_) = exp(z_X*EIGENVALUES_H(i_))
     END DO


     DO q=-N_SITES_Y/2 + 7,-N_SITES_Y/2 + 7!-N_SITES_Y/2,N_SITES_Y/2-1 !Q_T,Q_T!-N_SITES_Y/2,N_SITES_Y/2-1

        U_MB       = DCMPLX(0.0,0.0)
        U_Y       = DCMPLX(0.0,0.0)
        
        DO i_=1,N_SITES_X
           U_Y(i_,i_) = exp(z_Y*2.0*cos(2.0*pi*i_*alpha + 2*pi*q/N_SITES_Y))
        END DO
        U_Y = MATMUL(TRANSPOSE(U_AUX),MATMUL(U_Y,CONJG(U_AUX)))
        U_MB = MATMUL(U_Y,U_X)


        CALL LAPACK_ZGEEV(U_MB,N_SITES_X,EIGENVALUES,INFO)
        !WRITE(2,2064) U_T ! EIGENVECTORS STORED AS COLUMNS OF U_T; U_T(:,1) IN THE MOMENTUM BASIS |k,q>
        !WRITE(*,20128) REAL(EIGENVALUES)
        EIGENVALUES_AUX = EIGENVALUES
        CALL EV_QUASIENERGY(EIGENVALUES,EV_INDEX,SIZE(EIGENVALUES),INFO)
        WRITE(1,20256) REAL(EIGENVALUES(EV_INDEX))
        
        DO i_=1,N_SITES_X/3
            j_ = EV_INDEX(i_)
            IF(REAL(EIGENVALUES(j_)).GE.-PI   .AND. REAL(EIGENVALUES(j_)).LT.-PI/2) BAND_INDEX = 1
            IF(REAL(EIGENVALUES(j_)).GE.-PI/2 .AND. REAL(EIGENVALUES(j_)).LT. PI/2) BAND_INDEX = 2
            IF(REAL(EIGENVALUES(j_)).GE.PI/2  .AND. REAL(EIGENVALUES(j_)).LT. PI)   BAND_INDEX = 3
            !WRITE(*,*)N_SITES_X q,BAND_INDEX, REAL(EIGENVALUES(j_)),REAL(EIGENVALUES_AUX(j_)),AIMAG(EIGENVALUES_AUX(j_))
            DO q_=1,N_SITES_X
                CALL CARTtoEULER(REAL(U_MB(q_,j_)),AIMAG(U_MB(q_,j_)),norm(q_),theta(q_))
            END DO
            write(*,20256) theta
        END DO
        WRITE(*,*)
       WRITE(*,*)


       !U_MB   = MATMUL(U_AUX,U_MB)! EIGENVECTORS IN THE basis {|n,m>}
       !DO i_=1,N_SITES_X/3
        !j_ = EV_INDEX(i_)
        !write(*,*) abs(U_MB(:,j_))
       !END DO
!!        DO i_=1,N_SITES_Y
!!           U_FTY(i_) = exp(-DCMPLX(0.0,1.0)*2*pi*q*i_/N_SITES_Y)/SQRT(1.0*N_SITES_Y) ! phase factor for the Fourier Transform along y direction
!!        END DO
        !write(*,*) 'here', tau, T, q
     END DO
  END IF
  deallocate(u_aux)
  
  !201 format(1e15.6)
  !202 format(2e15.6)
  !203 FORMAT(3E15.6)
  !204 FORMAT(4E15.6)
  !206 FORMAT(6E15.6)
  2032 FORMAT(32E15.6)
  208 format(8e15.6)
  !2026 format(26e15.6)
  !2045 format(45e15.6)
  21552 FORMAT(1552E15.6)
  !209 format(9e15.6)
  2064 FORMAT(64E15.6)
  20128 FORMAT(128E15.6)
  20256 FORMAT(480E15.6)
  !303 format(3e15.6,3e15.6)
  !308 format(8e15.6,8e15.6)
  30128 format(128e15.6,128e15.6)
  !4021 FORMAT(2I8,1E15.6)
  !4031 FORMAT(3I8,1E15.6)
  !4051 FORMAT(5I8,1E15.6)
  
END SUBROUTINE FLOQUET_BOSE_HUBBARD_SUB

SUBROUTINE EV_QUASIENERGY(EV,EV_INDEX,N,INFO)

  USE physical_constants

  IMPLICIT NONE
  INTEGER,                       INTENT(IN)    :: N
  INTEGER,                       INTENT(INOUT) :: INFO
  INTEGER,    DIMENSION(N),      INTENT(OUT)   :: EV_INDEX
  COMPLEX*16, DIMENSION(N),      INTENT(INOUT) :: EV
  !DOUBLE PRECISION, DIMENSION(N),INTENT(INOUT)   :: QE

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

    if(theta.GE.PI) theta = theta - 2*pi

     EV(i_) = DCMPLX(theta,0.0)

  END DO

  DO i_=1,N

     EV_INDEX(i_) = i_

  END DO

  INFO  = 0

  CALL QUICK_SORT_I_T_DOUBLE(real(EV),EV_INDEX,N)
  !write(*,*) 'here',N,info
END SUBROUTINE EV_QUASIENERGY


SUBROUTINE CARTtoEULER(x,z,norm,theta)
    USE physical_constants
    implicit  none
    double precision, intent(in) :: x,z
    double precision, intent(out) :: norm, theta

     !write(*,*) z,x
     if(x.gt.0 .and. z.gt.0) theta = atan(x/z)
     if(x.gt.0 .and. z.lt.0) theta = atan(x/z) + 4.0*ATAN(1.0)
     if(x.lt.0 .and. z.lt.0) theta = atan(x/z) + 4.0*ATAN(1.0)
     if(x.lt.0 .and. z.gt.0) theta = atan(x/z) + 8.0*ATAN(1.0)
     if(z.eq.0 .and. x.gt.0) theta = 2.0*ATAN(1.0)
     if(z.eq.0 .and. x.lt.0) theta = 6.0*ATAN(1.0)
     if(x.eq.0 .and. z.eq.0) theta = 0.0

    if(theta.GE.PI) theta = theta - 2*pi

    norm = sqrt(x*x + z*z)

END SUBROUTINE CARTtoEULER

