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
  INTEGER i,j,N_partial,k,N_R,index_0,N_R2,j_
  INTEGER, DIMENSION(1) :: index_min,index_aux

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
SUBROUTINE CHECK_HASHING(TAG,TAG_INDEX,BASIS_DIM,STATE,INFO)
  USE LATTICE
  IMPLICIT NONE
  INTEGER, DIMENSION(:,:),INTENT(IN) :: STATE
  INTEGER, DIMENSION(BASIS_DIM),   INTENT(IN) :: TAG,TAG_INDEX
!  INTEGER, DIMENSION(N_SITES,BASIS_DIM),INTENT(IN) :: STATE
  INTEGER, INTENT(IN)  :: BASIS_DIM
  INTEGER, INTENT(OUT) :: INFO
  
  INTEGER i
  !INTEGER, DIMENSION(N_SITES) :: p

  INFO = 0
  i    = 2
  DO WHILE(INFO.EQ.0 .AND. i.LE.BASIS_DIM)
     IF(TAG(TAG_INDEX(i)).GT. TAG(TAG_INDEX(i-1))) THEN
        INFO=0
     ELSE
        INFO = i
        IF(TAG(TAG_INDEX(i)).NE.TAG(TAG_INDEX(i-1))) THEN
           WRITE(*,*) 'HASHING ERROR:',i,TAG_INDEX(i),TAG_INDEX(i-1),TAG(TAG_INDEX(i)),TAG(TAG_INDEX(i-1))
           WRITE(*,*) 'STATE i-1:',STATE(:,TAG_INDEX(i-1))
           WRITE(*,*) 'STATE i  :',STATE(:,TAG_INDEX(i))
           WRITE(*,*)
        ELSE
        END IF
     END IF
     i = i+1
  END DO
  
END SUBROUTINE CHECK_HASHING


SUBROUTINE NEWTON_SEARCHING(elto,array,array_index,elto_index,array_size,INFO)

  IMPLICIT NONE
  INTEGER,                        INTENT(IN)  :: elto,array_size
  INTEGER, DIMENSION(array_size), INTENT(IN)  :: array
  INTEGER, DIMENSION(:),          INTENT(IN),OPTIONAL:: array_index
  INTEGER,                        INTENT(OUT) :: elto_index,INFO
  INTEGER, DIMENSION(:), ALLOCATABLE :: array_
  INTEGER STOP,i_left,i_right,i_aux,i,MAXITERATIONS

  MAXITERATIONS = int(LOG(10.0*array_size)) + 2
  INFO = 0
  i_left  = 1
  i_right = array_size
  elto_index = -1
  STOP = -1
  
  i = 0
  !write(*,*) elto
  !write(*,*) array
  
  ALLOCATE(array_(array_size))
  IF(PRESENT(array_index)) THEN
     array_ = array(array_index)
!     write(*,*) "ad", array_, elto
  ELSE
     array_ = array
  END IF
    
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
        !write(*,*) "yo",i_aux
        STOP = 0        
     END IF
  !   write(*,*) elto_index, stop
 
     IF((array_(i_left)-elto)*(array_(i_aux)-elto) .GT.0  .AND. STOP.NE.0) THEN
        i_left = i_aux
     END IF
   !  write(*,*) elto_index,stop

     IF((array_(i_aux)-elto)*(array_(i_right)-elto).GT.0 .AND. STOP.NE.0) THEN
        i_right = i_aux     
     END IF
    ! write(*,*) elto_index,stop, i, MAXITERATIONS
     !write(*,*) i_left,i_right,i_aux
     i =  i + 1
     IF(i.GT.MAXITERATIONS) STOP = -2
    
  END DO
  IF(elto_index.EQ.-1) elto_index = i_aux  
  IF(STOP.EQ.-2) THEN
     INFO = -1
     WRITE(*,*) "ERROR IN NEWTON SEARCH: MAXIMA NUMBER OF ITERATIONS REACHED"
  ELSE
     INFO = 0
!     write(*,*) 'INFO 0' 
  END IF
!  write(*,*) elto_index
  
 
END SUBROUTINE NEWTON_SEARCHING
