SUBROUTINE STATEtoPARTITIONmap(STATE,TAG,TAG_INDEX,SITEMAP_A,SITEMAP_B,C,BASIS_DIM_A,BASIS_DIM_B,INFO)
 
  USE LATTICE  
  USE INTERFACE
  USE SUBLATTICE 
  
 
  IMPLICIT NONE
  INTEGER, DIMENSION(:,:), INTENT(IN)  :: STATE
  INTEGER, DIMENSION(:),   INTENT(IN)  :: SITEMAP_A, SITEMAP_B,TAG,TAG_INDEX
  INTEGER, DIMENSION(:,:), INTENT(OUT) :: C
  INTEGER, DIMENSION(:),   INTENT(INOUT) :: BASIS_DIM_A,BASIS_DIM_B
  INTEGER,                 INTENT(OUT)   :: INFO

  INTEGER,          DIMENSION(:,:), ALLOCATABLE :: STATES_AUX,STATES_A,STATES_B
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: TAG_A,TAG_INDEX_A
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: TAG_B,TAG_INDEX_B,TAG_AUX,TAG_INDEX_AUX
  
    
  INTEGER BASIS_DIM,i,j,i_,k_,DIM_A,DIM_B,N_SITES_B,NP_i,u,m_c,n_c
  DOUBLE PRECISION D_H
  
  BASIS_DIM = SIZE(STATE,2)  
  N_SITES_B = N_SITES - N_SITES_A
    
  ALLOCATE(TAG_AUX(1))
  ALLOCATE(TAG_INDEX_AUX(1))

  INFO  = 0
  
  BASIS_DIM_A(1) = 1
  DIM_A = 1
  DO i=1,N_PARTICLES
     BASIS_DIM_A(i+1) =  INT(D_H(N_SITES_A,i))
     DIM_A = DIM_A + BASIS_DIM_A(i+1)
  END DO
  WRITE(*,*) "#BASIS_DIM_A ",BASIS_DIM_A,DIM_A

  ALLOCATE(STATES_A(N_SITES_A,DIM_A))
  ALLOCATE(TAG_A(DIM_A))
  ALLOCATE(TAG_INDEX_A(DIM_A))
  TAG_A(1)       = 1
  TAG_INDEX_A(1) = 1
  STATES_A = 0
 
  i_ = 2
  k_ = i_+ BASIS_DIM_A(2)-1
  DO i=1,N_PARTICLES
     CALL BASIS(STATES_A(:,i_:k_),BASIS_DIM_A(i+1),i,N_SITES_A)
     CALL HASHING(STATES_A(:,i_:k_),TAG_A(i_:k_),TAG_INDEX_A(i_:k_),BASIS_DIM_A(i+1),N_SITES_A)
     CALL CHECK_HASHING(TAG_A,TAG_INDEX_A,BASIS_DIM_A(i+1),STATES_A(:,i_:k_),INFO)
     i_ = k_ + 1
     k_ = k_ + BASIS_DIM_A(i+2)
  END DO


  BASIS_DIM_B(1) = 1
  DIM_B = 1
  DO i=1,N_PARTICLES
     BASIS_DIM_B(i+1) =  INT(D_H(N_SITES_B,i))
     DIM_B = DIM_B + BASIS_DIM_B(i+1)
  END DO
  WRITE(*,*) "#BASIS_DIM_B ",BASIS_DIM_B,DIM_B

  ALLOCATE(STATES_B(N_SITES_B,DIM_B))
  ALLOCATE(TAG_B(DIM_B))
  ALLOCATE(TAG_INDEX_B(DIM_B))
  TAG_B(1)       = 1
  TAG_INDEX_B(1) = 1
  STATES_B = 0

  i_ = 2
  k_ = i_+ BASIS_DIM_B(2)-1
  DO i=1,N_PARTICLES
     CALL BASIS(STATES_B(:,i_:k_),BASIS_DIM_B(i+1),i,N_SITES_B)
     CALL HASHING(STATES_B(:,i_:k_),TAG_B(i_:k_),TAG_INDEX_B(i_:k_),BASIS_DIM_B(i+1),N_SITES_B)
     CALL CHECK_HASHING(TAG_B,TAG_INDEX_B,BASIS_DIM_B(i+1),STATES_B(:,i_:k_),INFO)     
     i_ = k_ + 1
     k_ = k_ + BASIS_DIM_B(i+2)
  END DO

  C = 0
  DO i=1,BASIS_DIM

     NP_i = 0
     ALLOCATE(STATES_AUX(N_SITES_A,1))
     DO i_=1,N_SITES_A
        NP_i = NP_i + STATE(SITEMAP_A(i_),i)     
        STATES_AUX(i_,1) = STATE(SITEMAP_A(i_),i)
     END DO
     C(3,i)= NP_i ! number of particles in region A
     i_=1
     k_=1
     IF(NP_i.GT.0) THEN
        CALL HASHING(STATES_AUX,TAG_AUX,TAG_INDEX_AUX,1,N_SITES_A)    
        !write(*,*) TAG_AUX
        !write(*,*) i_,k_
        !write(*,*) TAG_A(i_:k_)
        i_ = 2
        k_ = i_+ BASIS_DIM_A(2)-1
        DO j=1,NP_i-1
           i_ = k_ + 1
           k_ = k_ + BASIS_DIM_A(j+2)
        END DO
        IF(TAG_AUX(1).GT.0) THEN
           u = -1
           CALL NEWTON_SEARCHING(TAG_AUX(1),TAG_A(i_:k_),TAG_INDEX_A(i_:k_),u,BASIS_DIM_A(NP_i+1),INFO)
	   C(4,i) = u ! within the subset of states with number of particles NP_i in region A
           IF(INFO.EQ.0 .AND. (u.GE.1 .AND. u.LE. BASIS_DIM_A(NP_i+1)) ) THEN
              u = i_ + u - 1
           ELSE
              !INFO = -1
              !break
              !k = size(state,2)+1
              !i = N_SITES + 1
           END IF
        END IF
     ELSE
        u=1
        C(4,i) = u ! within the subset of states with number of particles (NP_i) in region A
     END IF
     m_c = u
     DEALLOCATE(STATES_AUX)
     
     
     NP_i = 0
     ALLOCATE(STATES_AUX(N_SITES_B,1))
     DO i_=1,N_SITES_B
        NP_i = NP_i + STATE(SITEMAP_B(i_),i)     
        STATES_AUX(i_,1) = STATE(SITEMAP_B(i_),i)
     END DO
     IF(NP_i.GT.0) THEN
        CALL HASHING(STATES_AUX,TAG_AUX,TAG_INDEX_AUX,1,N_SITES_B)        
        !write(*,*) TAG_AUX
        !write(*,*) i_,k_
        !write(*,*) TAG_A(i_:k_)
         i_ = 2
        k_ = i_+ BASIS_DIM_B(2)-1
        DO j=1,NP_i-1
           i_ = k_ + 1
           k_ = k_ + BASIS_DIM_B(j+2)
        END DO
        IF(TAG_AUX(1).GT.0) THEN
           u = -1
           CALL NEWTON_SEARCHING(TAG_AUX(1),TAG_B(i_:k_),TAG_INDEX_B(i_:k_),u,BASIS_DIM_B(NP_i+1),INFO)
           C(5,i)= u  ! within the subset of states with number of particles (N_PARTICLES - NP_i) in region B
           IF(INFO.EQ.0 .AND. (u.GE.1 .AND. u.LE. BASIS_DIM_B(NP_i+1)) ) THEN
              u = i_ + u - 1
           ELSE
              !INFO = -1
              !break
              !k = size(state,2)+1
              !i = N_SITES + 1
           END IF
        END IF
     ELSE
        u=1
        C(5,i) = u ! within the subset of states with number of particles (N_PARTICLES - NP_i) in region B
     END IF
     n_c = u
     DEALLOCATE(STATES_AUX)                    
     C(1,i) = m_c
     C(2,i) = n_c
  END DO
  
  !write(*,105) C
  DEALLOCATE(TAG_AUX)
  DEALLOCATE(TAG_INDEX_AUX)
  DEALLOCATE(STATES_A)
  DEALLOCATE(TAG_A)
  DEALLOCATE(TAG_INDEX_A)
  DEALLOCATE(STATES_B)
  DEALLOCATE(TAG_B)
  DEALLOCATE(TAG_INDEX_B)


105 FORMAT(5I4)
104 FORMAT(4I4)
1056 FORMAT(56I4) 
2010 FORMAT(3E15.6)
  
END SUBROUTINE STATEtoPARTITIONmap  


!  BASIS_DIM_A(i) : number of states in subregion A with total number of particles i-1
!  BASIS_DIM_B(i) : number of states in subregion B with total number of particles i-1
!  STATES_A(B), TAG_A(B), TAG_INDEX_A(B): Corresponding arrays in region A (B)). TAG_INDEX_A(B) are organized in such a way that the tags are organized in each subspace of different number of particles.
!  Then I run over all states in the system A+B and identify corresponding states m_c and n_c of subspaces A and B, respectively.
!  C(1,i) : Glogal indice in the state space of region A
!  C(2,i) : Global indice in the state space of region B
!  C(3,i) : Number of particles in A
!  C(4,i) : indice en el subspacio A, corresponding to the number of particles C(3,i) 
!  C(5,i) : indice en el subspacio B, corresponding to the number of particles N_PARTICLES - C(3,i) 
