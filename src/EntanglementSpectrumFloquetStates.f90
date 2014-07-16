PROGRAM EntanglementSpectrumFloquetSpectrum

  USE TASKS
  USE physical_constants  
  USE LATTICE
  USE INTERACTION
  USE INTERFACE
  USE SUBLATTICE
  USE MATRIX_AUX
  USE INTERFACE_MAP
  USE INTERFACE_FLOQUETSPECTRUM

  IMPLICIT NONE
  INTEGER,          DIMENSION(:,:), ALLOCATABLE :: STATE,C
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: TAG,TAG_INDEX,SITEMAP_A,SITEMAP_B,BASIS_DIM_A,BASIS_DIM_B
  
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: rho_A,U_T
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: a,U_FTY
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phi
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: E_S,EIGENVALUES_H
  TYPE(C_matrix) CS(N_PARTICLES)
  
  
!!!! TEST: BEGIN
  INTEGER, DIMENSION(1) :: OLD, SEED 
  DOUBLE PRECISION, DIMENSION(N_SITES) :: HARVEST
!!!! TEST: END
 
  INTEGER STATE_(N_SITES)
  LOGICAL MORE
  INTEGER BASIS_DIM,INFO,i,j,i_,k_,N_SITES_B,AUX
  DOUBLE PRECISION D_H

!!!! TEST: BEGIN
  SEED(1) = 12345
  CALL RANDOM_SEED
  CALL RANDOM_NUMBER(HARVEST)
  !WRITE(*,*) ' Random numbers : ', HARVEST
  
!!!! TEST: END
 

  MORE = .FALSE.
  INFO = 0
  IF(FLOQUET_SPECTRUM .EQ. 1 .AND. INFO.EQ.0) THEN
     !--- FIND THE FLOQUET SPECTRUM
     !ALLOCATE(a(BASIS_DIM))
     INFO = 0
     ALLOCATE(U_T(N_SITES_X,N_SITES_X))
     ALLOCATE(EIGENVALUES_H(N_SITES_X))
     ALLOCATE(U_FTY(N_SITES_Y))   
     CALL  FLOQUET_BOSE_HUBBARD_SUB(-N_SITES_Y/2+16,U_T,U_FTY, EIGENVALUES_H,INFO)
     WRITE(*,201) EIGENVALUES_H
     write(*,*)
     write(*,*)
     ALLOCATE(phi(N_SITES_X,N_SITES_Y))
     DO j=1,N_SITES_X !ENUMERATING THE EIGEN-STATE
        !DO i_=1,N_SITES_X ! 
        DO k_ = 1,1!N_SITES_Y
           phi(:,k_) = abs(U_FTY(k_)*U_T(:,j))
           write(*,*) phi(:,k_)**2
        END DO
        !END DO
        !WRITE(*,2032) phi
     END DO

  ELSE
     INFO = -1
  END IF

  INFO = 0
  IF(ENTANGLEMENT_SPECTRUM.EQ.1 .AND. INFO.EQ.0) THEN
     BASIS_DIM =  INT(D_H(N_SITES,N_PARTICLES))
     write(*,*) '#',BASIS_DIM
     ALLOCATE(STATE(N_SITES,BASIS_DIM))
     ALLOCATE(TAG(BASIS_DIM))
     ALLOCATE(TAG_INDEX(BASIS_DIM))
     ALLOCATE(C(5,BASIS_DIM))

!!!! TEST: BEGIN
     ALLOCATE(EIGENVALUES_H(N_SITES))
     EIGENVALUES_H = HARVEST
!!!! TEST: END

     INFO  = 0
     !--- DEFINE STATES IN SITE-PARTICLE-NUMBER BASIS.
     !--- ASSIGN A TAG TO EACH STATE AND FIND THE TAG_INDEX SUCH AS TAG(TAG_INDEX) IS ORDERED
!     CALL BASIS(STATE,BASIS_DIM)  
!     CALL CHECK_ORDER(STATE,BASIS_DIM,INFO)
!     IF(INFO.NE.0) WRITE(*,*) 'WARNING: ERROR GENERATING THE BASIS. ELEMENT:', INFO

!     CALL HASHING(STATE,TAG,TAG_INDEX,BASIS_DIM)
!     CALL CHECK_HASHING(TAG,TAG_INDEX,BASIS_DIM,STATE,INFO)
!     IF(INFO.NE.0) WRITE(*,*) 'WARNING: ERROR GENERATING THE BASIS HASHING. ELEMENT:', INFO

     !WRITE(*,109) STATE
     !WRITE(*,*)

     
     STATE = 0
     i =0
     DO
        CALL COMP_NEXT(N_PARTICLES,N_SITES,STATE_,MORE)         
        !WRITE (*,109) STATE_
        i = i +1
        STATE(:,i) = STATE_
        IF ( .NOT. MORE )  THEN
           EXIT
        END IF
     END DO
     !write(*,*) i
     !WRITE(*,109) STATE
     !WRITE(*,*)
     CALL HASHING(STATE,TAG,TAG_INDEX,BASIS_DIM)
     CALL CHECK_HASHING(TAG,TAG_INDEX,BASIS_DIM,STATE,INFO)
     IF(INFO.NE.0) WRITE(*,*) 'WARNING: ERROR GENERATING THE BASIS HASHING. ELEMENT:', INFO 

    
     STATE_  = 0
     MORE = .FALSE.

     DO i=22,22!1,BASIS_DIM
        WRITE(*,109) STATE(:,i)
        DO j=1,N_SITES
           STATE_ = 0
           MORE   = .FALSE.
           IF(STATE(j,i).GE.1) THEN
              !write(*,*) j
              DO
                 CALL COMP_NEXT(STATE(j,i),N_SITES,STATE_,MORE)         
                 WRITE(*,109)  STATE_
                 CALL MULTINOMIAL_COEF2(STATE(j,i), STATE_, AUX )
                 WRITE(*,*) AUX
                 
                 IF(.NOT.MORE) THEN
                    EXIT
                 END IF
              END DO
           ELSE
              
           END IF
        END DO
     END DO
     
     
     !--- DEFINE PARTITION
     IF(N_SITES_A .GE. 1 .AND. N_SITES_A.LT.N_SITES) THEN
        
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
        
        i= 1
        DO i_=N_SITES_X/2+1,N_SITES_X
           DO k_ =1,N_SITES_Y        
              SITEMAP_B(i) = (k_-1)*N_SITES_X + i_
              i = i+1
           END DO
        END DO
        
        
!!$
        !--- DEFINE STATES OF EACH PARTITION, CORRESPONDING TO DIFFERENT NUMBER OF PARTICLES IN EACH ONE
        !--- SET THE MAP BETWEEN THE ORIGINAL BASIS AND THE PARTITION:
        !---|\psi> = \sum_n^{BASIS_DIM} a_n \n> == \sum_i^{BASIS_DIM_A} \sum_i^{BASIS_DIM_A} C_{ij} \i>_A|j>_B
        ALLOCATE(BASIS_DIM_A(N_PARTICLES+1))
        ALLOCATE(BASIS_DIM_B(N_PARTICLES+1))
        C = 0 ! matrix containing the map a(n) to C_{ij}
        
        CALL STATEtoPARTITIONmap(STATE,TAG,TAG_INDEX,SITEMAP_A,SITEMAP_B,C,BASIS_DIM_A,BASIS_DIM_B,INFO)
        !write(*,105) C
        
        DEALLOCATE(STATE)
        DEALLOCATE(TAG)
        DEALLOCATE(TAG_INDEX)
        DEALLOCATE(SITEMAP_A)
        DEALLOCATE(SITEMAP_B)
        
        !---ALLOCATE MEMORY FOR AN ARRAY OF MATRICES THAT HELP TO BUILD THE DENSTIY MATRIX OF SUBSPACES CORRESPONDING TO DIFFERENT TOTAL NUMBER OF PARTICLES
        DO i=1,N_PARTICLES + 1
           ALLOCATE(CS(i)%C(BASIS_DIM_A(i),BASIS_DIM_B(N_PARTICLES + 2 - i)))  
        END DO
        
        DEALLOCATE(BASIS_DIM_A)
        DEALLOCATE(BASIS_DIM_B)
        
        DO i=1,N_PARTICLES+1
           WRITE(*,*) "#",i-1,size(CS(i)%C,1),size(CS(i)%C,2)
        END DO
        
        
        
        !--- USE THE FLQOUET(n)->c_{ij} MAP TO FILL THE AUXILIARY MATRICES     
        
        DO j=1,N_SITES_X
           DO k_=1,N_SITES_X
              DO i_=1,N_SITES_Y
                 i= (i_-1)*N_SITES_X + k_    
                 CS(C(3,i)+1)%C(C(4,i),C(5,i))=U_FTY(i_)*U_T(k_,j) 
              END DO
           END DO
           
           !DEALLOCATE(C)
           
           INFO = 0  
           DO i=1,1!N_PARTICLES+1
              ALLOCATE(rho_A(size(CS(i)%C,1),size(CS(i)%C,1)))
              ALLOCATE(E_S(size(rho_A,1)))
              DO i_=1,size(rho_A,1)
                 DO k_=i_,size(rho_A,1)
                    !  WRITE(*,*) DOT_PRODUCT(CS(i)%C(i_,:),CS(i)%C(k_,:))
                    rho_A(i_,k_) = DOT_PRODUCT(CS(i)%C(i_,:),CS(i)%C(k_,:))
                    IF(k_.NE.i_)rho_A(k_,i_) = CONJG(rho_A(i_,k_))
                 END DO
              END DO
              CALL LAPACK_FULLEIGENVALUES(rho_A,size(rho_A,1),E_S,INFO)
              WRITE(*,202) EIGENVALUES_H(j),E_S
              !write(*,*)
              DEALLOCATE(rho_A)
              DEALLOCATE(E_S)
           END DO
        END DO
     END IF
  ELSE
     WRITE(*,*) 'PARTITION IS NOT CORRECTLY DEFINED '
  END IF

105 FORMAT(5I4)
109 FORMAT(9I4)
209 FORMAT(9E15.6)
201 FORMAT(1E15.6)
202 FORMAT(2E15.6)
2032 FORMAT(32E15.6)
20256 FORMAT(256E15.6)

END PROGRAM EntanglementSpectrumFloquetSpectrum

