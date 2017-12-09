!  STATES: Full set of states of the lattice
!  TAG_INDEX: A list of indices that order the tags of the states
!  SITEMAP_A
!  SITEMAP_B
!  C(1,i) : Global indice in the state space of region A (i.e. site in the full lattice)
!  C(2,i) : Global indice in the state space of region B (i.e. site in the full lattice)
!  C(3,i) : Number of particles in A
!  C(4,i) : indice en el subspacio A, corresponding to the number of particles C(3,i)
!  C(5,i) : indice en el subspacio B, corresponding to the number of particles N_PARTICLES - C(3,i)
!  BASIS_DIM_A(i) : number of states in subregion A with total number of particles i-1
!  BASIS_DIM_B(i) : number of states in subregion B with total number of particles i-1
!  STATES_A(B), TAG_A(B), TAG_INDEX_A(B): Corresponding arrays in region A (B)). TAG_INDEX_A(B) are organized in such a way that the tags are organized in each subspace of different number of particles.
!  Then I run over all states in the system A+B and identify corresponding states m_c and n_c of subspaces A and B, respectively.

SUBROUTINE STATEtoPARTITIONmapV2(STATES,TAG_INDEX,SITEMAP_A,SITEMAP_B,C,BASIS_DIM_A,BASIS_DIM_B,INFO)

  USE LATTICE
  USE INTERFACE
  USE SUBLATTICE
  USE STATE_OBJ

  IMPLICIT NONE
  TYPE(STATE),  DIMENSION(:),   INTENT(IN)    :: STATES
  INTEGER,      DIMENSION(:),   INTENT(IN)    :: SITEMAP_A, SITEMAP_B,TAG_INDEX
  INTEGER,      DIMENSION(:,:), INTENT(OUT)   :: C
  INTEGER,      DIMENSION(:),   INTENT(INOUT) :: BASIS_DIM_A,BASIS_DIM_B
  INTEGER,                      INTENT(OUT)   :: INFO

  TYPE(STATE),  DIMENSION(:),   ALLOCATABLE   :: STATES_A,STATES_B!,STATES_AUX
  INTEGER,      DIMENSION(:),   ALLOCATABLE   :: TAG_INDEX_A,TAG_INDEX_B

  INTEGER,      DIMENSION(:),   ALLOCATABLE   :: p_A
  INTEGER,      DIMENSION(:),   ALLOCATABLE   :: p_B

  INTEGER,      DIMENSION(2,2)                :: n,m
  INTEGER BASIS_DIM,i,j,i_,k_,DIM_A,DIM_B,N_SITES_B,NP_i,u,j_A,j_B

  INTEGER TAG_A, TAG_B,TAG_FUNCTION,INFO_
  DOUBLE PRECISION D_H

  BASIS_DIM = SIZE(STATES,1)

    INFO_ = 0
  N_SITES_B = N_SITES - N_SITES_A
  !INFO  = 0

  i = TAG_INDEX(1)
  i = SITEMAP_A(1)
  i = SITEMAP_B(1)
  BASIS_DIM_A(1) = 1
  DIM_A = 1
  DO i=1,N_PARTICLES
     BASIS_DIM_A(i+1) =  INT(D_H(N_SITES_A,i))
     DIM_A = DIM_A + BASIS_DIM_A(i+1)
  END DO
  !WRITE(*,*) "#BASIS_DIM_A ",BASIS_DIM_A,DIM_A

  ALLOCATE(p_A(N_SITES_A))
  ALLOCATE(p_B(N_SITES_B))


  ALLOCATE(STATES_A(DIM_A))
  ALLOCATE(TAG_INDEX_A(DIM_A))
  TAG_INDEX_A(1) = 1

  STATES_A(1)%TAG = 0
  STATES_A(1)%site(1) = 0
  STATES_A(1)%occupation(1) = 0
  IF(N_PARTICLES.EQ.2) THEN
  STATES_A(1)%site(2) = 0
  STATES_A(1)%occupation(2) = 0
  END IF
  !write(*,*) BASIS_DIM_A
  i_ = 2
  k_ = i_+ BASIS_DIM_A(2)-1
  DO i=1,N_PARTICLES
     CALL BASISV3(STATES_A(i_:k_),TAG_INDEX_A(i_:k_),BASIS_DIM_A(i+1),i,N_SITES_X_A,N_SITES_Y_A)
     CALL CHECK_HASHING(STATES_A(i_:k_)%TAG,TAG_INDEX_A(i_:k_),BASIS_DIM_A(i+1),INFO)
     !IF(i.EQ.1) write(*,109) STATES_A(i_:k_)
!     DO j= i_,k_
!        IF((STATES_A(j)%n(1) .EQ. 2) .AND. (STATES_A(j)%n(2) .EQ. 4) .AND. (STATES_A(j)%m(1) .EQ. 1) &
!            &  .AND. (STATES_A(j)%m(2) .EQ. 7)) WRITE(*,*) 'FOUND:',j,&
!            & TAG_INDEX_A(j),states_a(j)%tag
!     END DO
     i_ = k_ + 1
     k_ = k_ + BASIS_DIM_A(i+2)

  END DO


  BASIS_DIM_B(1) = 1
  DIM_B = 1
  DO i=1,N_PARTICLES
     BASIS_DIM_B(i+1) =  INT(D_H(N_SITES_B,i))
     DIM_B = DIM_B + BASIS_DIM_B(i+1)
  END DO

  ALLOCATE(STATES_B(DIM_B))
  ALLOCATE(TAG_INDEX_B(DIM_B))
  TAG_INDEX_B(1) = 1

  STATES_B(1)%TAG = 0
  STATES_B(1)%site = 0
  STATES_B(1)%occupation = 0
  STATES_B(1)%n = 0
  STATES_B(1)%m = 0

  i_ = 2
  k_ = i_+ BASIS_DIM_B(2)-1
  DO i=1,N_PARTICLES
     CALL BASISV3(STATES_B(i_:k_),TAG_INDEX_B(i_:k_),BASIS_DIM_B(i+1),i,N_SITES_X_B,N_SITES_Y_B)
     CALL CHECK_HASHING(STATES_B(i_:k_)%TAG,TAG_INDEX_B(i_:k_),BASIS_DIM_B(i+1),INFO)
 !    IF(i.EQ.2) THEN
        !WRITE(*,*) i_,k_
!        write(*,109) STATES_B(i_:k_)
 !    END IF
     i_ = k_ + 1
     k_ = k_ + BASIS_DIM_B(i+2)
  END DO


  !write(*,*) '12',STATES_B()
  !write(*,*) '474',STATES_B(474)
  C = 0
  !write(*,*) INFO
  INFO_ = INFO
  DO i=1,BASIS_DIM

     !DO j=1,N_SITES_A
     !  p_A(j) = INT(SQRT(1.0*HASH_SLOPE*SQRT(1.0*j) + 1.0*HASH_SHIFT)) ! to tag states in A
     !END DO

     !DO j=1,N_SITES_B
     !  p_B(j) = INT(SQRT(1.0*HASH_SLOPE*SQRT(1.0*j) + 1.0*HASH_SHIFT)) ! to tag states in B
     !END DO

     TAG_A = 0
     TAG_B = 0
     NP_I  = 0
     m     = 0
     n     = 0
     j_A   = 0
     j_B   = 0

     DO j=1,N_PARTICLES  ! assuming there are less particles than sites
        ! This is defined for the partition of the cilynder where A is the left half and B is the right-half
        ! different partitions requires different asigments

        IF(STATES(i)%n(j).NE.0 .AND. STATES(i)%n(j).LE.N_SITES_X/2) THEN
            j_A   = j_A + 1
            !i_=STATES(i)%n(j) + N_SITES_X_A*(STATES(i)%m(j) - 1) ! site in the sublattice A
            !TAG_A = TAG_A + p_A(i_)*STATES(i)%occupation(j)!TAG_A + INT(SQRT(1.0*HASH_SLOPE*SQRT(1.0*i_) + 1.0*HASH_SHIFT))*STATES(i)%occupation(j)
            !TAG_A = TAG_FUNCTION()
            n(1,j_A) = STATES(i)%n(j)
            m(1,j_A) = STATES(i)%m(j)
            NP_I     = NP_I + 1
            !END IF!!
            !IF(STATES(i)%n(j).NE.0 .AND. STATES(i)%n(j).GT.N_SITES_X/2) THEN
        ELSE
            j_B = j_B + 1
            n(2,j_B) = STATES(i)%n(j) - N_SITES_X/2
            m(2,j_B) = STATES(i)%m(j)
            !i_   = STATES(i)%n(j) - N_SITES_X/2 + N_SITES_X_B*(STATES(i)%m(j) - 1) ! site in the sublattice B
            !TAG_B = TAG_B + p_B(i_)*STATES(i)%occupation(j)!TAG_B + INT(SQRT(1.0*HASH_SLOPE*SQRT(1.0*i_) + 1.0*HASH_SHIFT))*STATES(i)%occupation(j)
        END IF

     END DO


     TAG_A = 0
     IF(NP_I .GT. 1) THEN

        TAG_A =     TAG_FUNCTION(m(1,1),m(1,2),n(1,1),n(1,2),N_SITES_X_A,N_SITES_Y_A)
     ELSE IF(NP_I .EQ. 1) THEN
        TAG_A = 1 + TAG_FUNCTION(m(1,1),m(1,2),n(1,1),n(1,2),N_SITES_X_A,N_SITES_Y_A)
     END IF

     TAG_B = 0
     IF(N_PARTICLES - NP_I.GT.1) THEN
        TAG_B =     TAG_FUNCTION(m(2,1),m(2,2),n(2,1),n(2,2),N_SITES_X_B,N_SITES_Y_B)
     ELSE IF (N_PARTICLES - NP_I.EQ.1) THEN
        TAG_B = 1 + TAG_FUNCTION(m(2,1),m(2,2),n(2,1),n(2,2),N_SITES_X_B,N_SITES_Y_B)
     END IF

     C(3,i) = NP_I
     IF(i.eq.info_) write(*,*) i,STATES(i),NP_I,TAG_A,TAG_B
     IF(TAG_A.NE.0) THEN
         i_ = 2
         k_ = i_+ BASIS_DIM_A(2)-1
         !write(*,*) TAG_A,i_,k_
         DO j=1,NP_I-1
            i_ = k_ + 1
            k_ = k_ + BASIS_DIM_A(j+2)
         END DO
         !write(*,*) TAG_A,i_,k_
         !write(*,101) STATES_A(i_:k_)
         CALL NEWTON_SEARCHING(TAG_A,STATES_A(i_:k_)%TAG,TAG_INDEX_A(i_:k_),u,BASIS_DIM_A(NP_i+1),INFO)
         u = TAG_INDEX_A(i_ + u - 1)! because "u" is initially found with respecto to the ordered array, and the order is given by TAG_INDEX
         IF(INFO.EQ.0 .AND. (u.GE.1 .AND. u.LE. BASIS_DIM_A(NP_i+1)) ) THEN
              C(4,i) = u ! within the subset of states with number of particles (NP_i) in region A
              C(1,i) = i_ + u - 1 ! state in subspace A
         ELSE
         write(*,*) "Error with the partition in A: ",i
              !INFO = -1
              !break
              !k = size(state,2)+1
              !i = N_SITES + 1
         END IF
     ELSE
         C(4,i) = 1! within the subset of states with number of particles = 0 in region A
         C(1,i) = 1
     END IF

     IF(TAG_B.NE.0) THEN
         i_ = 2
         k_ = i_+ BASIS_DIM_B(2)-1
         DO j=1,N_PARTICLES - NP_I - 1
            i_ = k_ + 1
            k_ = k_ + BASIS_DIM_B(j+2)
         END DO
         CALL NEWTON_SEARCHING(TAG_B,STATES_B(i_:k_)%TAG,TAG_INDEX_B(i_:k_),u,BASIS_DIM_B(N_PARTICLES - NP_i+1),INFO)
         IF(i.EQ.INFO_)write(*,101) STATES_B(i_:k_)%TAG
         !IF((i_+u) .LE. (k_-i_+1)) THEN
         !   u = TAG_INDEX_B(i_+u)
         !ELSE
         u = TAG_INDEX_B(i_+u-1)
         !END IF! because "u" is initially found with respecto to the ordered array, and the order is given by TAG_INDEX
         IF(INFO.EQ.0 .AND. (u.GE.1 .AND. u.LE. BASIS_DIM_B(N_PARTICLES - NP_i+1)) ) THEN
              C(5,i) = u! within the subset of states with number of particles (N_PARTICLES - NP_i) in region B
              C(2,i) = i_ + u - 1 ! state in subspace B
         ELSE
            write(*,*) "Error with the partition in B: ",i
              !INFO = -1
              !break
              !k = size(state,2)+1
              !i = N_SITES + 1
         END IF
    ELSE
         C(5,i) = 1! within the subset of states with number of particles = 0 in region B
         C(2,i) = 1
    END IF

    IF(i.eq.info_) write(*,*) C(:,i)

   END DO

  DEALLOCATE(STATES_A)
  DEALLOCATE(TAG_INDEX_A)
  DEALLOCATE(STATES_B)
  DEALLOCATE(TAG_INDEX_B)
  DEALLOCATE(p_A)
  DEALLOCATE(p_B)


101 FORMAT(1I8)
!104 FORMAT(4I4)
!105 FORMAT(5I10)
!109 FORMAT(9I8)
!1056 FORMAT(56I4)
!201 FORMAT(1E15.6)
!2010 FORMAT(3E15.6)

END SUBROUTINE STATEtoPARTITIONmapV2


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

  i = TAG_INDEX(1)
  i = TAG(1)
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
     !CALL CHECK_HASHING(TAG_A,TAG_INDEX_A,BASIS_DIM_A(i+1),STATES_A(:,i_:k_),INFO)
     CALL CHECK_HASHING(TAG_A,TAG_INDEX_A,BASIS_DIM_A(i+1),INFO)
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
     !CALL CHECK_HASHING(TAG_B,TAG_INDEX_B,BASIS_DIM_B(i+1),STATES_B(:,i_:k_),INFO)
     CALL CHECK_HASHING(TAG_B,TAG_INDEX_B,BASIS_DIM_B(i+1),INFO)
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


!105 FORMAT(5I4)
!104 FORMAT(4I4)
!1056 FORMAT(56I4)
!2010 FORMAT(3E15.6)
  
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
