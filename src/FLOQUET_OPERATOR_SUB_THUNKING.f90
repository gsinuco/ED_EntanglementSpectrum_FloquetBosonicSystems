! En este programa calculamos de manera exacta el spectro de Floquet
! en un tigh-binding Hamiltonian PRL...
! Single particle!

SUBROUTINE FLOQUET_BOSE_HUBBARD_SUB(STATES,Q_T,U_FTY,U_MB,EIGENVALUES,INFO)
  
  USE MAGMA
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
  

  
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE   :: E_MBX,E_MBY,EIGENVALUES_H,rho
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE   :: U_X, U_Y,U_AUX,U_FLOQUET
  INTEGER,          DIMENSION(:),   ALLOCATABLE   :: EV_INDEX
  
  INTEGER BASIS_DIM,i_,k_,q,q_,D, j, i_off
  DOUBLE PRECISION D_H
  COMPLEX*16 :: z_Y,z_X, alpha_,beta_

  
!    write(*,*) "CUBLAS, USING THUNKING"
  OPEN(UNIT=1,FILE='/home/g/gs/gs249/Programs/ED-EntanglementSpectrum-FloquetV2/data/SP_thunking.dat',ACTION='WRITE')

  call cublas_init()
  
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
  !U_AUX(1,N_SITES_X) = -J_X  ! add periodic boundary conditions
  U_AUX = U_AUX + transpose(conjg(U_AUX))
  CALL LAPACK_FULLEIGENVALUES(U_AUX,N_SITES_X,EIGENVALUES_H,INFO)

  !write(*,*) EIGENVALUES_h
  IF((N_PARTICLES.GT.1) .AND. (INFO.EQ.0)) THEN
     
     ALLOCATE(U_FTY(N_PARTICLES))
     q  = 1
     q_ = 1
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


     DO Q_T = 1,1!N_SITES_Y!(N_SITES_Y*N_SITES_Y + N_SITES_Y)/2 !loop over blocks of total momentum_Y
        
        q  = Q_T
        q_ = Q_T
        
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
        ALLOCATE(U_Y(D,D))
        U_Y = 0
        DO i_=1,D
           U_Y(i_,i_) = exp(z_Y*E_MBY(i_)) ! In its eigen-basis
        END DO

        ALLOCATE(U_FLOQUET(D,D))
        U_FLOQUET = 0
        
        ALPHA_ = DCMPLX(1.0,0.0)
        BETA_  = DCMPLX(0.0,0.0)
        call cublas_zGEMM('n','n',D,D,D,alpha_,U_Y, D,U_MB,D,beta_,U_FLOQUET,D) ! U_FLOQUET IN THIS CASE IS AN AUXILIARY MATRIX
        call cublas_zGEMM('c','n',D,D,D,alpha_,U_MB,D,U_FLOQUET,  D,beta_,U_Y,      D)

        ALLOCATE(U_X(D,D))
        U_X = 0
        DO i_=1,D
           U_X(i_,i_) = exp(z_X*E_MBX(i_)) ! in its eigen-basis
        END DO

           U_FLOQUET = DCMPLX(0.0,0.0)


        ALPHA_ = DCMPLX(1.0,0.0)
        BETA_  = DCMPLX(0.0,0.0)
        call cublas_zGEMM('n','n',D,D,D,alpha_,U_Y,D,U_X,D,beta_,U_FLOQUET,D)

        DEALLOCATE(U_X)
        DEALLOCATE(U_Y)
        ALLOCATE(EIGENVALUES(D))
        EIGENVALUES = 0.0


        CALL LAPACK_ZGEEV_GPU_NOTHUNKING(U_FLOQUET,D,EIGENVALUES,INFO)!Eigenvectors of the full Floquet operator, stored as U_MB(:,i) in the {|kq,k'q'>}

        ALLOCATE(EV_INDEX(D))
        EV_INDEX = -10000
        CALL EV_QUASIENERGY(EIGENVALUES,EV_INDEX,SIZE(EIGENVALUES),INFO) ! The eigenvalues of U_MB are complex values. This function evaluates the norm of the eigenvalues
        DO i_=1,D
           WRITE(1,*)    Q_T,q,q_,REAL(EIGENVALUES(EV_INDEX(i_))),EV_INDEX(i_)
        END DO
        WRITE(1,*)
        WRITE(1,*)
        DEALLOCATE(EV_INDEX,EIGENVALUES,U_FLOQUET)

     END DO

     DEALLOCATE(U_FTY)
     DEALLOCATE(U_MB)
     DEALLOCATE(E_MBX)
     DEALLOCATE(E_MBY)

     call cublas_shutdown()
     
  ELSE IF(INFO.EQ.0) THEN
     
     ALLOCATE(U_FTY(N_SITES_Y))
     ALLOCATE(U_X(N_SITES_X,N_SITES_X))
     ALLOCATE(U_Y(N_SITES_X,N_SITES_X))
     ALLOCATE(EIGENVALUES(N_SITES_X))
     ALLOCATE(U_MB(N_SITES_X,N_SITES_X))
     
     z_X       = DCMPLX(0.0,-1.0)*tau
     z_Y       = DCMPLX(0.0,-1.0)*(T-tau)

     !Evolution operator U(0,tau) =  exp(-ihbar tau H_X)
     DO i_=1,N_SITES_X	
        U_X(i_,i_) = exp(z_X*EIGENVALUES_H(i_))
     END DO
     
     
     DO q=-N_SITES_Y/2,N_SITES_Y/2-1 !Q_T,Q_T!-N_SITES_Y/2,N_SITES_Y/2-1
        
        !q= Q_T
        U_MB       = DCMPLX(0.0,0.0)
        U_Y       = DCMPLX(0.0,0.0)
        
        DO i_=1,N_SITES_X
           U_Y(i_,i_) = exp(z_Y*2.0*cos(2.0*pi*i_*alpha + 2*pi*q/N_SITES_Y))
        END DO
        U_Y = MATMUL(TRANSPOSE(U_AUX),MATMUL(U_Y,CONJG(U_AUX)))
        U_MB = MATMUL(U_Y,U_X)
        
        CALL LAPACK_ZGEEV(U_MB,N_SITES_X,EIGENVALUES,INFO)
        !WRITE(2,2064) U_T ! EIGENVECTORS STORED AS COLUMNS OF U_T; U_T(:,1) IN THE MOMENTUM BASIS |k,q>
        
        EIGENVALUES_H = 0.0
        CALL EV_QUASIENERGY(EIGENVALUES,EV_INDEX,SIZE(EIGENVALUES),INFO)
        !WRITE(1,30128) EIGENVALUES_H
        !WRITE(*,2064) REAL(EIGENVALUES)
        !write(*,*)
        
        U_MB   = MATMUL(U_AUX,U_MB)! EIGENVECTORS IN THE basis {|n,m>}
        DO i_=1,N_SITES_Y
           U_FTY(i_) = exp(-DCMPLX(0.0,1.0)*2*pi*q*i_/N_SITES_Y)/SQRT(1.0*N_SITES_Y) ! phase factor for the Fourier Transform along y direction
        END DO
        
     END DO
  END IF
  deallocate(u_aux)

END SUBROUTINE FLOQUET_BOSE_HUBBARD_SUB

SUBROUTINE EV_QUASIENERGY(EV,EV_INDEX,N,INFO)

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

     EV(i_) = DCMPLX(theta,0.0)

  END DO

  DO i_=1,N
     EV_INDEX(i_) = i_
  END DO

  INFO  = 0
  CALL QUICK_SORT_I_T_DOUBLE(abs(EV),EV_INDEX,N)

END SUBROUTINE EV_QUASIENERGY
