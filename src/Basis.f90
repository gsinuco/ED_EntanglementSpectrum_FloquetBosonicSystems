  SUBROUTINE BASISV2(STATES,TAG_INDEX,BASIS_DIM,NP,NS)
  USE LATTICE
  USE STATE_OBJ
  USE A_DAGGER_

   ! ----------------- BEGINING -----------------
  ! ----- BUILD THE MANY-BODY LATTICE BASIS ----
  ! ----------------- BEGINING -----------------
  IMPLICIT NONE
  TYPE(STATE), DIMENSION(:), ALLOCATABLE, INTENT(OUT)          :: STATES
  INTEGER,     DIMENSION(:), ALLOCATABLE, INTENT(OUT)          :: TAG_INDEX
  INTEGER,                                INTENT(OUT)          :: BASIS_DIM
  INTEGER,                                INTENT(IN), OPTIONAL :: NP,NS


  !INTEGER, DIMENSION(:), ALLOCATABLE :: STATE_,TAGS,STATE_AUX
  INTEGER, DIMENSION(:), ALLOCATABLE :: p
  INTEGER i,j,j_,m,m_,n,n_
  LOGICAL MORE
  DOUBLE PRECISION D_H

  INTEGER N_PARTICLES_, N_SITES_,TAG_FUNCTION
  MORE = .FALSE.

  IF(PRESENT(NP) .AND. PRESENT(NS)) THEN
     N_PARTICLES_ = NP
     N_SITES_ = NS
  END IF

  BASIS_DIM =  INT(D_H(N_SITES,N_PARTICLES))
  !ALLOCATE(STATE_(N_SITES))
  !ALLOCATE(TAGS(BASIS_DIM))
  !ALLOCATE(STATE_AUX(N_SITES))
  BASIS_DIM = (N_SITES_X*N_SITES_X  + N_SITES_X)/2 + N_SITES_X*N_SITES_X !restricting to m=m_ and m!=m_
  ALLOCATE(STATES(BASIS_DIM))
  ALLOCATE(TAG_INDEX(BASIS_DIM))
  ALLOCATE(p(N_SITES))

  STATES(:)%TAG = 0
  STATES(:)%site(1) = 0
  STATES(:)%occupation(1) = 0
  IF(N_PARTICLES.EQ.2) THEN
    STATES(:)%site(2) = 0
    STATES(:)%occupation(2) = 0
  END IF
  MORE = .FALSE.
  i =0

  DO j=1,N_SITES
     p(j) = int(log(1.0*hash_slope*j) + hash_shift)
  END DO
  i = 0
  write(*,*) BASIS_DIM
  DO m=1,1!N_SITES_Y
     DO m_ = m,2!N_SITES_Y

!        i_off = 0
!        Q_T = 0
!        DO i_=1,m-1
!           Q_T = Q_T + N_SITES_Y - (i_-1)
!           i_off = i_off + (N_SITES_Y - (i_-1) -1)*N_SITES_X*N_SITES_X
!        END DO
!        Q_T =  Q_T + (m_ - m + 1)
!
!        i_off = i_off + (m-1)*(N_SITES_X*N_SITES_X+N_SITES_X)/2
!        if(m_.GT.m) i_off = i_off + (N_SITES_X*N_SITES_X+N_SITES_X)/2
!        DO i_=m+1,m_-1
!           i_off = i_off + N_SITES_X*N_SITES_X
!        END DO

        !write(*,*) Q_T,i_off
        DO n=1,N_SITES_X
           IF(m.EQ.m_) THEN
              DO n_= n,N_SITES_X
                 i = i + 1
                 j  = n  + (m  - 1)*N_SITES_X
                 j_ = n_ + (m_ - 1)*N_SITES_X
                 !STATE_    = 0
                 !STATE_AUX = 0
                 !CALL a_dagger(j ,STATE_,STATE_AUX)
                 !CALL a_dagger(j_,STATE_AUX,STATE_)
                 !write(*,*) i
                 !WRITE(*,1012) STATE_
!                 Q_T = 0
!                 DO i_=1,n-1
!                    Q_T = Q_T + N_SITES_X - (i_-1)
!                 END DO
!                 Q_T =  Q_T + (n_ - n + 1)
                 STATES(i)%TAG = TAG_FUNCTION(m,m_,n,n_,N_SITES_X,N_SITES_Y)!i_off+Q_T!p(j) + p(j_)! INT((1.0*HASH_SLOPE*(1.0*j) + 1.0*HASH_SHIFT)) &
                 !&+ INT((1.0*HASH_SLOPE*(1.0*j_) + 1.0*HASH_SHIFT))
                 STATES(i)%site(1)       = j
                 STATES(i)%occupation(1) = 1!STATE_(j)
                 STATES(i)%m(1)          = m
                 STATES(i)%n(1)          = n
                 IF(N_PARTICLES.EQ.2) THEN
                    STATES(i)%site(2)       = j_
                    STATES(i)%occupation(2) = 1!STATE_(j_)
                    STATES(i)%m(2)          = m_
                    STATES(i)%n(2)          = n_
                 END IF

                 !write(*,*) n,m,n_,m_, i
              END DO
           ELSE
              DO n_= 1,N_SITES_X
                 i = i +1
                 j  = n  + (m  - 1)*N_SITES_X
                 j_ = n_ + (m_ - 1)*N_SITES_X
                 !STATE_    = 0
                 !STATE_AUX = 0
                 !CALL a_dagger(j ,STATE_,STATE_AUX)
                 !CALL a_dagger(j_,STATE_AUX,STATE_)
                 !write(*,*) i
                 !WRITE(*,1012) STATE_
!                 Q_T = N_SITES_X*(n-1) + n_
                 STATES(i)%TAG =  TAG_FUNCTION(m,m_,n,n_,N_SITES_X,N_SITES_Y)!i_off + Q_T!p(j) + p(j_)!INT((1.0*HASH_SLOPE*(1.0*j) + 1.0*HASH_SHIFT)) &
                 !&+ INT((1.0*HASH_SLOPE*(1.0*j_) + 1.0*HASH_SHIFT))!p(j) + p(j_)
                 STATES(i)%site(1)       = j
                 STATES(i)%occupation(1) = 1!STATE_(j)
                 STATES(i)%m(1)          = m
                 STATES(i)%n(1)          = n
                 IF(N_PARTICLES.EQ.2) THEN
                 STATES(i)%site(2)       = j_
                 STATES(i)%occupation(2) = 1!STATE_(j_)
                 STATES(i)%m(2)          = m_
                 STATES(i)%n(2)          = n_
                 END IF
                 !write(*,*) n,m,n_,m_, i
              END DO
           END IF
        END DO
     END DO
  END DO

  WRITE(*,*) '# BULDING THE BASIS: DONE '

!!$  STATE_ = 0
!!$  i = 0
!!$  DO
!!$     CALL COMP_NEXT(N_PARTICLES,N_SITES,STATE_,MORE)
!!$     write(*,109) STATE_
!!$     i = i +1
!!$     k = 0
!!$     STATES(i)%TAG = DOT_PRODUCT(STATE_,p)
!!$     DO j=1,size(STATE_)
!!$        IF (STATE_(j).GT.0 .AND. STATE_(j).LT.N_PARTICLES) THEN
!!$           k = k + 1
!!$           STATES(i)%site(k)       = j
!!$           STATES(i)%occupation(k) = STATE_(j)
!!$           STATES(i)%m(k)          = INT((j-1)/N_SITES_X) + 1
!!$           STATES(i)%n(k)          = j - N_SITES_X*(STATES(i)%m(k) - 1)
!!$!           write(*,*) '14324123'
!!$       ELSE IF(STATE_(j).EQ.N_PARTICLES) THEN
!!$           !k = k + 1
!!$           !write(*,*)'yo'
!!$           STATES(i)%site(:)       = j
!!$           STATES(i)%occupation(:) = 1
!!$           STATES(i)%m(:)          = INT((j-1)/N_SITES_X) + 1
!!$           STATES(i)%n(:)          = j - N_SITES_X*(STATES(i)%m(1) - 1)
!!$        END IF
!!$
!!$
!!$     END DO
!!$     IF ( .NOT. MORE )  THEN
!!$        EXIT
!!$     END IF
!!$  END DO
  DO j=1,BASIS_DIM
     TAG_INDEX(j) = j ! because the basis is created in order (the tags), we save a call to quicksort. To check, checkorder.
  END DO
  !IF(BASIS_DIM.GT.1) CALL  QUICK_SORT_I_T(STATES(:)%TAG,TAG_INDEX,BASIS_DIM)
!  DO j=1,BASIS_DIM
!     write(*,*) STATES(TAG_INDEX(j))%TAG
!  END DO

!  WRITE(*,101) TAG_INDEX
!101 format(1i10)
!!$  DEALLOCATE(STATE_)
!!$  DEALLOCATE(p)


! ----------------- END -----------------
! ----- BUILD THE MANY-BODY LATTICE BASIS ----
! ----------------- END -----------------
END SUBROUTINE BASISV2

SUBROUTINE BASISV3(STATES,TAG_INDEX,BASIS_DIM,NP,NX,NY)
  USE LATTICE
  USE STATE_OBJ
  USE A_DAGGER_

   ! ----------------- BEGINING -----------------
  ! ----- BUILD THE MANY-BODY LATTICE BASIS ----
  ! ----------------- BEGINING -----------------

  IMPLICIT NONE
  TYPE(STATE), DIMENSION(:), ALLOCATABLE, INTENT(OUT)          :: STATES
  INTEGER,     DIMENSION(:), ALLOCATABLE, INTENT(OUT)          :: TAG_INDEX
  INTEGER,                                INTENT(OUT)          :: BASIS_DIM
  INTEGER,                                INTENT(IN), OPTIONAL :: NP,NX,NY


  INTEGER, DIMENSION(:), ALLOCATABLE :: STATE_
  INTEGER, DIMENSION(:), ALLOCATABLE :: p
  INTEGER i,j,k,N_PARTICLES_,N_SITES_,N_SITES_X_,N_SITES_Y_
  INTEGER TAG_FUNCTION
  LOGICAL MORE
  DOUBLE PRECISION D_H

  MORE = .FALSE.

  IF(PRESENT(NP) .AND. PRESENT(NX) .AND. PRESENT(NY)) THEN
     N_PARTICLES_ = NP
     N_SITES_     = NX*NY
     N_SITES_X_   = NX
     N_SITES_Y_   = NY
  ELSE
     N_PARTICLES_ = N_PARTICLES
     N_SITES_     = N_SITES
     N_SITES_X_   = N_SITES_X
     N_SITES_Y_   = N_SITES_Y
  END IF

  BASIS_DIM =  INT(D_H(N_SITES_,N_PARTICLES_))
  ALLOCATE(STATE_(N_SITES_))
  ALLOCATE(STATES(BASIS_DIM))
  ALLOCATE(TAG_INDEX(BASIS_DIM))
  ALLOCATE(p(N_SITES_))

  DO i=1,N_PARTICLES
    STATES(:)%TAG = 0
    STATES(:)%site(i) = 0
    STATES(:)%occupation(i) = 0
    STATES(:)%n(i) = 0
    STATES(:)%m(i) = 0
  END DO

  MORE = .FALSE.
  i =0

  DO j=1,N_SITES_
     p(j) = INT(SQRT(1.0*HASH_SLOPE*SQRT(1.0*j) + 1.0*HASH_SHIFT))
  END DO
  i = 0
  STATE_    = 0
  STATE_(1) = N_PARTICLES_
  DO
     CALL COMP_NEXT(N_PARTICLES_,N_SITES_,STATE_,MORE)
     write(*,*) i,STATE_
     i = i +1
     k = 0
    STATES(i)%TAG = DOT_PRODUCT(STATE_,p)
!
     DO j=1,size(STATE_)
        IF (STATE_(j).GT.0 .AND. STATE_(j).LT.N_PARTICLES_) THEN
           k = k + 1
           STATES(i)%site(k)       = j
           STATES(i)%occupation(k) = STATE_(j)
           STATES(i)%m(k)          = INT((j-1)/N_SITES_X_) + 1
           STATES(i)%n(k)          = j - N_SITES_X_*(STATES(i)%m(k) - 1)
       ELSE IF(STATE_(j).EQ.N_PARTICLES_) THEN

           IF(N_PARTICLES_.GT.1) THEN
                STATES(i)%site(:)       = j
                STATES(i)%occupation(:) = 1
                STATES(i)%m(:)          = INT((j-1)/N_SITES_X_) + 1
                STATES(i)%n(:)          = j - N_SITES_X_*(STATES(i)%m(1) - 1)
           ELSE
                STATES(i)%site(1)       = j
                STATES(i)%occupation(1) = 1
                STATES(i)%m(1)          = INT((j-1)/N_SITES_X_) + 1          !(along y)
                STATES(i)%n(1)          = j -  N_SITES_X_*(STATES(i)%m(1) - 1) !(along x)
           END IF
        END IF
     END DO

     IF(N_PARTICLES_.EQ.2) THEN
         STATES(i)%TAG = TAG_FUNCTION(STATES(i)%m(1),STATES(i)%m(2),STATES(i)%n(1)&
                        &,STATES(i)%n(2),N_SITES_X_,N_SITES_Y_)

     ELSE IF (N_PARTICLES_.EQ.1) THEN
          STATES(i)%TAG = 1 + TAG_FUNCTION(STATES(i)%m(1),0,STATES(i)%n(1)&
                       &,0,N_SITES_X_,N_SITES_Y_)

     END IF
!
     IF ( .NOT. MORE )  THEN
        EXIT
     END IF

  END DO

  DO j=1,BASIS_DIM
     TAG_INDEX(j) = j
  END DO

  IF(BASIS_DIM.GT.1) CALL  QUICK_SORT_I_T(STATES(:)%TAG,TAG_INDEX,BASIS_DIM)
  WRITE(*,*) BASIS_DIM
  WRITE(*,*) states(TAG_INDEX)%TAG


! ----------------- END -----------------
! ----- BUILD THE MANY-BODY LATTICE BASIS ----
! ----------------- END -----------------

!101 FORMAT(1I6)
!1024 FORMAT(24I6)

END SUBROUTINE BASISV3

SUBROUTINE BASIS(STATE,BASIS_DIM,NP,NS)

! Notese que k no esta definido en algunos casos y se usa el valor obtenido 
! en la iteracion anterior. No entiendo porque funciona asi.

  USE LATTICE
  IMPLICIT NONE
  INTEGER, DIMENSION(:,:), INTENT(OUT) :: STATE
  INTEGER, INTENT(IN) :: BASIS_DIM
  INTEGER, INTENT(IN), OPTIONAL :: NP,NS
  INTEGER i,j,N_partial,k,N_R,index_0,N_PARTICLES_,N_SITES_
  INTEGER, DIMENSION(1) :: index_min,index_aux
  
  STATE = 0
  IF(PRESENT(NP) .AND. PRESENT(NS)) THEN
     STATE(1,1) = NP
     N_PARTICLES_ = NP
     N_SITES_ = NS
     !write(*,*) "fajs",N_SITES_
  ELSE
     STATE(1,1) = N_PARTICLES
     N_PARTICLES_ = N_PARTICLES
     N_SITES_ = N_SITES
     !write(*,*) "fajs"
   END IF
  k = 0
  !WRITE(*,*) STATE(:,1)
  DO i=2,Basis_dim

     !write(*,*) Basis_dim
     index_0 = 1
     DO WHILE (STATE(index_0,i-1).EQ.0) 
        index_0 = index_0 + 1
     END DO
          
     index_min =  MINLOC(STATE(index_0:N_SITES_,i-1))
     index_min(1) =  index_min(1) + index_0 - 1 

     N_R = 0
     DO j=index_min(1)+1,N_SITES_
        N_R =  N_R + STATE(j,i-1)
     END DO
   
     DO WHILE (N_R.GT.0 .AND. index_min(1).LT.N_SITES_-1) 
        
        index_aux =  MINLOC(STATE(index_min(1)+1:N_SITES_,i-1))
        index_min(1) =  index_min(1) +  index_aux(1)
        
        N_R = 0
        DO j=index_min(1)+1,N_SITES_
           N_R =  N_R + STATE(j,i-1)
        END DO
        
     END DO
     
     index_aux(1) = index_min(1) 
    
     DO WHILE ((index_aux(1)).GE.1 .AND. STATE(index_aux(1),i-1).EQ.0)
        index_aux(1) =  index_aux(1) -1
        k = index_aux(1)
     END DO


     N_partial = 0          

     DO j=1,N_SITES_        
        
        IF(j.GE.1 .AND. j.LE.k-1) THEN
           state(j,i) = state(j,i-1)
           N_partial = N_partial +  state(j,i)
        END if
        
        IF(j.EQ.k) THEN
           state(j,i) = state(j,i-1) - 1           
           N_partial = N_partial + state(j,i)
        END IF
        
        IF(j.EQ.k+1)            state(j,i) = N_PARTICLES_ - N_partial
        IF(j.GE.k+2)            state(j,i) = 0        
     END DO

     N_R = 0
     DO j=1,N_SITES_
        N_R = N_R + STATE(j,i)
     END DO
     IF(N_R.NE.N_PARTICLES_) WRITE(*,*)  &
          & "WARNING: ERROR GENERATING SOME BASIS ELEMENTS", i,N_R,N_PARTICLES_
     !WRITE(*,*) STATE(:,i)
  END DO
!  write(*,*) "END BASIS SUBROUTINE"
 
END SUBROUTINE BASIS


SUBROUTINE CHECK_ORDER(STATE,BASIS_DIM,INFO)
  
  USE LATTICE
  IMPLICIT NONE
  INTEGER, DIMENSION(:,:), INTENT(IN) :: STATE
  INTEGER, INTENT(IN)  :: BASIS_DIM
  INTEGER, INTENT(OUT) :: INFO
  INTEGER i,j

  INFO = 0
  i=2
  DO WHILE( i.LE.BASIS_DIM .AND. INFO.EQ.0)
     
     j = 1
     DO WHILE(STATE(j,i) .EQ.STATE(j,i-1))
        j = j+1
     END DO
     IF(STATE(j,i).LT.STATE(j,i-1)) THEN
        INFO = 0
     ELSE
        INFO = j	

     END IF
     i = i+1
  END DO
 
END SUBROUTINE CHECK_ORDER

SUBROUTINE HASHING(STATE,TAG,TAG_INDEX,BASIS_DIM,NS)
  USE LATTICE
  IMPLICIT NONE
  INTEGER, DIMENSION(:,:), INTENT(IN)           :: STATE
  INTEGER, DIMENSION(:),   INTENT(OUT)          :: TAG,TAG_INDEX
  INTEGER,                 INTENT(IN), OPTIONAL :: NS
  INTEGER,                 INTENT(IN)           :: BASIS_DIM

  INTEGER, DIMENSION(:), ALLOCATABLE :: p
  INTEGER i,N_SITES_

  IF(PRESENT(NS)) THEN
     N_SITES_ = NS
     !write(*,*) NS
  ELSE
     N_SITES_ = N_SITES
  END IF

  IF(N_SITES_ .GT. 0) ALLOCATE(p(N_SITES_))

  DO i=1,N_SITES_
     p(i) = INT(SQRT(1.0*HASH_SLOPE*SQRT(1.0*i) + 1.0*HASH_SHIFT)) 
  END DO
!  write(*,*) P,N_SITES_
  TAG =  MATMUL(P,STATE)
      
  DO i=1,BASIS_DIM
     TAG_INDEX(i) = i
  END DO

  IF(BASIS_DIM.GT.1) CALL  QUICK_SORT_I_T(TAG,TAG_INDEX,BASIS_DIM)

END SUBROUTINE HASHING

!SUBROUTINE CHECK_HASHING(TAG,TAG_INDEX,BASIS_DIM,INFO)
!SUBROUTINE CHECK_HASHING(TAG,TAG_INDEX,BASIS_DIM,STATE,INFO)
SUBROUTINE CHECK_HASHING(TAG,TAG_INDEX,BASIS_DIM,INFO)
  USE LATTICE
  IMPLICIT NONE
  !INTEGER, DIMENSION(:,:),INTENT(IN) :: STATE
  INTEGER, DIMENSION(BASIS_DIM),   INTENT(IN) :: TAG
  INTEGER, DIMENSION(BASIS_DIM),   INTENT(IN) :: TAG_INDEX
!  INTEGER, DIMENSION(N_SITES,BASIS_DIM),INTENT(IN) :: STATE
  INTEGER, INTENT(IN)  :: BASIS_DIM
  INTEGER, INTENT(INOUT) :: INFO
  
  INTEGER i
  !INTEGER, DIMENSION(N_SITES) :: p

  !INFO = 0
  i    = 2
  DO WHILE(INFO.GE.0 .AND. i.LE.BASIS_DIM)
     IF(TAG(TAG_INDEX(i)).GT. TAG(TAG_INDEX(i-1))) THEN
        IF(INFO.EQ.0) INFO=0
     ELSE
        INFO = i
        IF(TAG(TAG_INDEX(i)).NE.TAG(TAG_INDEX(i-1))) THEN
           WRITE(*,*) '#HASHING ERROR:',i,TAG_INDEX(i),TAG_INDEX(i-1),TAG(TAG_INDEX(i-2)),TAG(TAG_INDEX(i-1)), &
            & TAG(TAG_INDEX(i)),TAG(TAG_INDEX(i+1)),TAG(TAG_INDEX(i+2)), "NE"
   !        WRITE(*,*) 'STATE i-1:',STATE(:,TAG_INDEX(i-1))
   !        WRITE(*,*) 'STATE i  :',STATE(:,TAG_INDEX(i))
           WRITE(*,*)
        ELSE
            WRITE(*,*) '#HASHING ERROR:',i,TAG_INDEX(i),TAG_INDEX(i-1),TAG(TAG_INDEX(i-2)),TAG(TAG_INDEX(i-1)),TAG(TAG_INDEX(i)),&
                    & TAG(TAG_INDEX(i+1)),TAG(TAG_INDEX(i+2)),(TAG(TAG_INDEX(i))+TAG(TAG_INDEX(i+1)))/2, "EQ"
                  !  TAG(TAG_INDEX(i)) = (TAG(TAG_INDEX(i))+TAG(TAG_INDEX(i+1)))/2
                  !                      info = 0
        END IF
     END IF
     i = i+1
  END DO
  
END SUBROUTINE CHECK_HASHING


SUBROUTINE NEWTON_SEARCHING(elto,array,array_index,elto_index,array_size,INFO)

  IMPLICIT NONE
  INTEGER,                        INTENT(IN)  :: elto,array_size
  INTEGER, DIMENSION(array_size), INTENT(IN)  :: array
  INTEGER, DIMENSION(array_size), INTENT(IN),OPTIONAL:: array_index
  INTEGER,                        INTENT(OUT) :: elto_index
  INTEGER,                        INTENT(INOUT) :: INFO
  INTEGER, DIMENSION(:), ALLOCATABLE :: array_
  INTEGER STOP,i_left,i_right,i_aux,i,MAXITERATIONS


  MAXITERATIONS = 50!int(LOG(1000.0*array_size)) + 5
  !INFO = 0
  i_left  = 1
  i_right = array_size
  elto_index = -1
  STOP = -1
  
  i = 0
  !write(*,*) INFO
  IF(INFO.EQ.10) THEN
  write(*,*) 'buscando a:',elto
  !write(*,*) array
  write(*,*) array_index
  write(*,*) 'en el ordered array:'
  WRITE(*,101) array(array_index)
  END IF
  ALLOCATE(array_(array_size))
  IF(PRESENT(array_index)) THEN
     array_ = array(array_index)
!     write(*,*) "ad", array_, elto
  ELSE
     array_ = array
  END IF
  IF(INFO.EQ.10) write(*,*) 'NEWTON'
  !INFO = 0
  DO WHILE(STOP.EQ.-1)
     
     i_aux   = (i_left + i_right)/2
     IF((array_(i_left) - elto)  .EQ. 0) THEN
        elto_index = i_left
        !write(*,*) array_(i_left) , elto, i_left,elto_index
        STOP = 0        
     END IF
!     write(*,*) elto_index
     IF((array_(i_right) - elto) .EQ. 0) THEN
        elto_index = i_right
        !write(*,*) "right",elto_index
        STOP = 0        
     END IF
 !    write(*,*) elto_index,STOP
     IF((array_(i_aux) - elto)   .EQ. 0) THEN
        elto_index = i_aux
        STOP = 0        
        IF(INFO.EQ.10) write(*,*) "yo",i_aux, elto,array_(i_aux), STOP
     END IF
  !   write(*,*) elto_index, stop
 
     IF((array_(i_left)-elto)*(array_(i_aux)-elto) .GT.0  .AND. STOP.NE.0) THEN
        i_left = i_aux
     END IF
   !  write(*,*) elto_index,stop

     IF((array_(i_aux)-elto)*(array_(i_right)-elto).GT.0 .AND. STOP.NE.0) THEN
        i_right = i_aux
     END IF
     IF(INFO.EQ.10) THEN
    !     write(*,*) elto_index,stop, i, MAXITERATIONS
        write(*,*) i_left,i_right,i_aux,array_(i_aux)
     END IF
     i =  i + 1
     IF(i.GT.MAXITERATIONS) STOP = -2
     if(INFO.EQ.10) WRITE(*,*) STOP
  END DO
  IF(elto_index.EQ.-1) elto_index = i_aux  
  IF(STOP.EQ.-2) THEN
     INFO = -1
     WRITE(*,*) "ERROR IN NEWTON SEARCH: MAXIMA NUMBER OF ITERATIONS REACHED"
  ELSE
     !IF(INFO.EQ.0)
     INFO = 0
!     write(*,*) 'INFO 0' 
  END IF
!  write(*,*) elto_index
  
101 FORMAT(1I8)
END SUBROUTINE NEWTON_SEARCHING


INTEGER FUNCTION  TAG_FUNCTION(m,m_,n,n_,NX,NY)

    IMPLICIT NONE
    INTEGER, INTENT(IN):: n,m,n_,m_,NX,NY

    INTEGER i_off,Q_T,i_

    i_off = 0
    DO i_=1,m-1
       i_off = i_off + (NY - (i_-1) -1)*NX*NX
    END DO
    i_off = i_off + (m-1)*(NX*NX+NX)/2

    IF(m_.GT.m) i_off = i_off + (NX*NX+NX)/2

    DO i_=m+1,m_-1
       i_off = i_off + NX*NX
    END DO

    IF(m.EQ.m_) THEN

        Q_T = 0
        DO i_=1,n-1
            Q_T = Q_T + NX - (i_-1)
        END DO
        Q_T =  Q_T + (n_ - n + 1)

    ELSE
        Q_T = NX*(n-1) + n_
    END IF

    TAG_FUNCTION = i_off + Q_T

END FUNCTION
