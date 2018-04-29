! Single particle case
! two particles

!BASISV2(STATES,BASIS_DIM): generates the state basis for two particles. It loops through combinations of all possible n's (first index)
!                           with fixed m's. For example: (m=1,m'=2, n=1,2,3 n'=1,2,3)
!                                                        |100 100 000>
!                                                        |100 010 000>
!                                                        |100 001 000>
!                                                        |010 100 000>
!                                                        |010 010 000>
!                                                        |010 001 000> ...
!                           STATES: a structure with information about the states
!                           BASIS_DIM: Number of states

!FLOQUET_BOSE_HUBBARD_SUB(STATES,Q_T,U_FTY,U_MB,EIGENVALUES,EIGENVALUES_H,INFO):
!   TYPE(STATES) STATES(BASIS_DIM)  : INPUT
!   INTEGER Q_T,                    : indicates the block of m,m' (named q,q')
!   COMPLEX*16, DIMENSION(:) U_FTY, : contains q,q'
!   COMPLEX*16, DIMENSION(D,D) U_MB : Eigenvectors of the Floquet matrix in the basis {|n,m ; n'm'>}
!   COMPLEX*16, DIMENSION(D) EIGENVALUES: Quasi-energies of the Floquet matrix
!   COMPLEX*16, DIMENSION(N_SITES_X) EIGENVALUES_H :
!   INTEGER INFO

PROGRAM EntanglementSpectrumFloquetSpectrum

  USE TASKS
  USE physical_constants  
  USE LATTICE
  USE INTERACTION
  USE INTERFACE
  USE SUBLATTICE
  USE MATRIX_AUX
  USE STATE_OBJ
  USE INTERFACE_MAP
  USE INTERFACE_FLOQUETSPECTRUM


  IMPLICIT NONE
  INTEGER,          DIMENSION(:,:), ALLOCATABLE :: C
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: TAG_INDEX,SITEMAP_A,SITEMAP_B,BASIS_DIM_A,BASIS_DIM_B

  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: rho_A,U_MB
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: U_FTY,EIGENVALUES
  ! DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E_S
  TYPE(STATE),      DIMENSION(:),   ALLOCATABLE :: STATES
  TYPE(C_matrix),   DIMENSION(N_PARTICLES)      :: CS



  INTEGER STATE_(N_SITES)
  LOGICAL MORE
  INTEGER BASIS_DIM,INFO,i,j,i_,k,k_,N_SITES_B,AUX,Q_T,D,j_
  INTEGER i_off,TAG_FUNCTION
  INTEGER m,m_,n,n_
  !DOUBLE PRECISION D_H
  INTEGER, DIMENSION(1) :: SEED
  DOUBLE PRECISION :: HARVEST

  SEED(1) = 35413
  CALL RANDOM_SEED(seed(1))
  CALL RANDOM_NUMBER(HARVEST)


  !WRITE(*,*) HARVEST
  MORE = .FALSE.
  INFO = 0

  STATE_ = 0


  CALL BASISV3(STATES,TAG_INDEX,BASIS_DIM)
  CALL CHECK_HASHING(STATES(:)%TAG,TAG_INDEX,BASIS_DIM,INFO)
!  DO m = 1, BASIS_DIM
!    STATES(m)%site(:)
!  END DO

  IF(N_PARTICLES.EQ.2) THEN
    CALL BASISV2(STATES,TAG_INDEX,BASIS_DIM)
    WRITE(*,*)'# BASIS_DIM (second label)', BASIS_DIM
    !CALL CHECK_HASHING(STATES(:)%TAG,TAG_INDEX,BASIS_DIM,INFO)
  END IF
  IF(FLOQUET_SPECTRUM .EQ. 1 .AND. INFO.EQ.0) THEN
     !--- FIND THE FLOQUET SPECTRUM
     !WRITE(*,*) '#',BASIS_DIM
     INFO = 0
     Q_T  = 10  !TOTAL MOMENTUM. IT CAN TAKE THE VALUES 1,(N_SITES_Y*N_SITES_Y + N_SITES_Y)/2
     CALL  FLOQUET_BOSE_HUBBARD_SUB(STATES,Q_T,U_FTY,U_MB,EIGENVALUES,INFO)
     !DO i = 3,3!1,size(U_MB,1)
     !   WRITE(*,201) (ABS(U_MB(1:size(U_MB,1),i))) ! EIGENVECTORS STORED AS COLUMNS OF U_MB; U_T(:,1) IN THE MOMENTUM BASIS |k,q>
     !END  DO
     !write(*,*)
     !write(*,*)'here'
     !write() EIGENVALUES
  ELSE
     INFO = -1
  END IF

  INFO = 0
  IF(ENTANGLEMENT_SPECTRUM.EQ.1 .AND. INFO.EQ.0) THEN

     ALLOCATE(C(5,BASIS_DIM))


     IF((N_SITES_A + N_SITES_X_B*N_SITES_Y_B) .EQ. N_SITES) THEN

        !--- DEFINE PARTITION ------------
        ! define a correspondence between sites of the subsets (A and B) and the full lattice.
        N_SITES_B = N_SITES - N_SITES_A

        ALLOCATE(SITEMAP_A(N_SITES_A))
        ALLOCATE(SITEMAP_B(N_SITES_B))

        i= 1
        DO i_=1,N_SITES_X/2
           DO k_ =1,N_SITES_Y
              SITEMAP_A(i) = (k_-1)*N_SITES_X + i_
              i = i+1
           END DO
        END DO

        i = 1
        DO i_=N_SITES_X/2+1,N_SITES_X
           DO k_ =1,N_SITES_Y
              SITEMAP_B(i) = (k_-1)*N_SITES_X + i_
              i = i+1
           END DO
        END DO

        !--- DEFINE STATES OF EACH PARTITION, CORRESPONDING TO DIFFERENT NUMBER OF PARTICLES IN EACH ONE
        !--- SET THE MAP BETWEEN THE ORIGINAL BASIS AND THE PARTITION:
        !---|\psi> = \sum_n^{BASIS_DIM} a_n \n> == \sum_i^{BASIS_DIM_A} \sum_i^{BASIS_DIM_B} C_{ij} \i>_A|j>_B

        ALLOCATE(BASIS_DIM_A(N_PARTICLES+1))
        ALLOCATE(BASIS_DIM_B(N_PARTICLES+1))

        C = 0 ! matrix containing the map a(n) to C_{ij}
        INFO = 0 ! to test: set info equal to the state to be tested.
        CALL STATEtoPARTITIONmapV2(STATES,TAG_INDEX,SITEMAP_A,SITEMAP_B,C,BASIS_DIM_A,BASIS_DIM_B,INFO)
        !write(*,*) states(2067)
        !WRITE(*,*) C(:,2067)
        DEALLOCATE(SITEMAP_A)
        DEALLOCATE(SITEMAP_B)
        !
        !---ALLOCATE MEMORY FOR AN ARRAY OF MATRICES THAT HELP TO BUILD THE DENSTIY MATRIX OF SUBSPACES CORRESPONDING TO DIFFERENT TOTAL NUMBER OF PARTICLES
        DO i=1,N_PARTICLES + 1
           ALLOCATE(CS(i)%C(BASIS_DIM_A(i),BASIS_DIM_B(N_PARTICLES - i + 2 )))
           CS(i)%C(:,:) = 0
           !WRITE(*,*) BASIS_DIM_A(i),BASIS_DIM_B(N_PARTICLES - i + 2 )
        END DO


        !--- USE THE U_MB(n)->c_{ij} MAP TO FILL THE AUXILIARY MATRICES
        !
        j = 0
        i = 0
        i_off = 0

        DO i_=1,Q_T
           j = i_off
           IF(STATES(j+1)%site(1) .EQ. STATES(j+1)%site(2)) THEN
              D = (N_SITES_X*N_SITES_X  + N_SITES_X)/2 ! EQUAL MOMENTUM q == q'
           ELSE
              D = N_SITES_X*N_SITES_X                  ! DIFFERENT MOMENTUM q != q'
           END IF
           IF(i_.LE.Q_T-1) then
              i_off = i_off + D
           END IF
        END DO


        WRITE(*,*) "#Partition basis dimension:"
        DO j=1,N_PARTICLES + 1

           WRITE(*,103) j-1,BASIS_DIM_A(j),BASIS_DIM_B(N_PARTICLES - j + 2 )!ABS(CS(j)%C(:,:))
           !!    WRITE(*,201)ABS(CS(j)%C(:,:))
        END DO
        write(*,*)
        write(*,*)
        !write(*,*) i_off,D
        i_ = 3 ! eigen state of interest
        DO m=1,N_SITES_Y
           DO m_=1,N_SITES_Y
              DO i=1,D

                 n  = STATES(i_off+i)%n(1)
                 n_ = STATES(i_off+i)%n(2)
                 IF(n.EQ.0 .OR. n_.EQ.0) WRITE(*,*) 'error, there is only one particle!'
                 j  = n  + (m  - 1)*N_SITES_X
                 j_ = n_ + (m_ - 1)*N_SITES_X
                 IF(j.GT.j_) THEN
                    AUX  = TAG_FUNCTION(m_,m,n_,n,N_SITES_X,N_SITES_Y)
                 END IF
                 IF(j.LE.j_) THEN
                    AUX  = TAG_FUNCTION(m,m_,n,n_,N_SITES_X,N_SITES_Y)
                 END IF
                 INFO = 0
                 CALL NEWTON_SEARCHING(AUX,STATES%TAG,TAG_INDEX,k,BASIS_DIM,INFO)
                 k =  tag_index(k)
                 IF(j.LE.j_) CS(C(3,k)+1)%C(C(4,k),C(5,k)) = U_MB(i,i_)*exp(-2*pi*DCMPLX(0.0,1.0)*&
                      & (U_FTY(1)*m + U_FTY(2)*m_)/N_SITES_Y)/N_SITES_Y
                 IF(j.GT.j_)CS(C(3,k)+1)%C(C(4,k),C(5,k)) = U_MB(i,i_)*exp(-2*pi*DCMPLX(0.0,1.0)*&
                      & (U_FTY(1)*m_ + U_FTY(2)*m)/N_SITES_Y)/N_SITES_Y
                 !IF((STATES(k)%m(1).EQ.m .AND. STATES(k)%m(2).EQ.m_) .OR. (STATES(k)%m(1).EQ.m_ &
                 ! & .AND. STATES(k)%m(2).EQ.m)) THEN
                 ! write(*,*) 'Problem: m,m_ do not correspond to states()%m'
                 ! WRITE(*,*)STATES(k)%m(1),STATES(k)%m(2),m,m_
                 !END IF
              END DO
           END DO
        END DO


        ! REDUCED DENSITY MATRIX:
        INFO = 0
        DO i=1,N_PARTICLES+1
           ALLOCATE(rho_A(BASIS_DIM_A(i),BASIS_DIM_A(i)))
           ALLOCATE(E_S(size(rho_A,1)))
           DO j=1,BASIS_DIM_A(i)
              DO k_=j,BASIS_DIM_A(i)
                 rho_A(j,k_) = DOT_PRODUCT(CS(i)%C(j,:),CS(i)%C(k_,:))
              END DO
           END DO
           rho_A = rho_A + transpose(conjg(rho_A))
           !DO j=1,BASIS_DIM_A(i)
           !      write(*,*) abs(rho_A(j,:))
           !END DO
           !write(*,*)
           !write(*,*)
           CALL LAPACK_FULLEIGENVALUES(rho_A,BASIS_DIM_A(i),E_S,INFO)
           DO i_=1, BASIS_DIM_A(i)
              WRITE(*,202)  1.0*i-1,E_S(i_)
           END DO
           write(*,*)
           write(*,*)
           DEALLOCATE(rho_A)
           DEALLOCATE(E_S)
        END DO

     ELSE
        WRITE(*,*) '#THE PARTITION IS NOT CORRECTLY DEFINED '
     END IF

     ! DEALLOCATE(BASIS_DIM_A)
     ! DEALLOCATE(BASIS_DIM_B)

     IF(INFO.EQ.0) WRITE(*,*) '#ENTANGLEMENT SPECTRUM DONE '
  ELSE
     WRITE(*,*) '#NO ENTANGLEMENT SPECTRUM HAS BEEN CALCULATED '
  END IF


  !  DEALLOCATE(U_MB)
  !  DEALLOCATE(EIGENVALUES)

!101 FORMAT(1I8)
!102 FORMAT(2I4)
103 FORMAT(3I4)
!105 FORMAT(5I4)
!109 FORMAT(9I6)
!201 FORMAT(1E15.6)
202 FORMAT(2E15.6)
!209 FORMAT(9E15.6)
!2010 FORMAT(10E15.6)
!2021 FORMAT(21E15.6)
!2032 FORMAT(32E15.6)
!20256 FORMAT(256E15.6)

END PROGRAM EntanglementSpectrumFloquetSpectrum



!            i_off = 0
!            DO n  = 1,4
!            DO m  = 1,3
!            DO n_ = 1,4
!            DO m_ = 1,3
!                ! n  = 1
!                ! m  = 1
!
!                ! n_ = 1
!                ! m_ = 1
!                 j  = n  + (m  - 1)*N_SITES_X
!                 j_ = n_ + (m_ - 1)*N_SITES_X
!                 IF(j.GT.j_) THEN
!                    AUX  = TAG_FUNCTION(m_,m,n_,n,N_SITES_X,N_SITES_Y)
!                 END IF
!                 IF(j.LE.j_) THEN
!                    AUX  = TAG_FUNCTION(m,m_,n,n_,N_SITES_X,N_SITES_Y)
!                 END IF
!                 i_off = i_off +1
!
!                    INFO = 0
!                    CALL NEWTON_SEARCHING(AUX,STATES%TAG,TAG_INDEX,i,BASIS_DIM,INFO)
!                    i =  tag_index(i)
!                    CALL RANDOM_NUMBER(HARVEST)
!                    CALL RANDOM_NUMBER(HARVEST)
!                    CALL RANDOM_NUMBER(HARVEST)
!                    !WRITE(*,*) HARVEST
!                    WRITE(*,*) i_off, AUX,i,C(3,i),C(4,i),C(5,i),harvest,n,m,n_,m_,j,j_
!    !                 WRITE(*,109) STATES(i)
!                     CS(C(3,i)+1)%C(C(4,i),C(5,i)) = sqrt(HARVEST)!10!U_MB(j,i_)*exp(-2*pi*DCMPLX(0.0,1.0)*&
!                      !& (U_FTY(1)*STATES(i)%m(1) + U_FTY(2)*STATEs(i)%m(2))/N_SITES_Y)/N_SITES_Y
!
!            END DO
!            END DO
!            END DO
!            END DO
