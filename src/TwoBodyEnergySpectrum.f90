! for now, Just two particle spectrum
! U_AUX:    Eigenvectors of the single particle Hamiltonian
! EIGENVALUES_H: spectrum of the single particle hamiltonian
! U_MB : eigenvectors of the many-body states built upon single particle states
! E_MB : spectrum of the many-body states


SUBROUTINE  TWO_BODY_ENERGY_SPECTRUM(STATES,EIGENVALUES_H,D,i_off,q,q_,E_MBX,E_MBY,INFO)

  USE physical_constants  
  USE LATTICE
  USE STATE_OBJ

  IMPLICIT NONE


  TYPE(STATE),       DIMENSION(:),   INTENT(IN)                 :: STATES
  DOUBLE PRECISION,  DIMENSION(:),   INTENT(IN)                 :: EIGENVALUES_H
  INTEGER,                           INTENT(IN)                 :: D, i_off,q,q_
  DOUBLE PRECISION,  DIMENSION(:),   INTENT(OUT)                :: E_MBX,E_MBY
  INTEGER,                           INTENT(INOUT)              :: INFO


  INTEGER :: k2,i_,i,j
  INTEGER :: m(2)
  
  m(1) = q
  m(2) = q_

  INFO = 0

  ! ----------------- BEGINING -----------------
  ! ----- BUILD THE TWO BODY SPECTRUM
  ! -----
  ! ----------------- BEGINING -----------------
  
  !DO j=1,BASIS_DIM ! loop over the new basis |kq k'q'>
  ! Since U_MB is block diagonal, the loop over states in the new basis can be replaced by  a loop over blocks of total momentum_Y
  ! the resulting matrix has a smaller size.
  j = 0
  i = 0

  !DO k1=Q_T,Q_T!1,(N_SITES_Y*N_SITES_Y + N_SITES_Y)/2 !loop over blocks of total momentum_Y
    ! THE SIZE OF THE BLOCK DEPENDS ON the combination of momentum_Y making up the states
    ! Also, the vector containing the states is ordered following the momentu block, therefore an index offset is needed.

    j = i_off

    !ALLOCATE(E_MBX(D))
    !ALLOCATE(E_MBY(D))

    E_MBX = 0
    E_MBY = 0

    DO k2=1,D

         j = j + 1

         DO i_=1,SIZE(STATES(j)%occupation)
            IF(STATES(j)%occupation(i_) .GT. 0) THEN
               E_MBX(k2) = E_MBX(k2) + STATES(j)%occupation(i_)*EIGENVALUES_H(STATES(j)%n(i_)) ! corresponding to the |kq,k'q'> state
               E_MBY(k2) = E_MBY(k2) + STATES(j)%occupation(i_)*2.0*cos(2.0*pi*STATES(j)%n(i_)*alpha &
                &   + 2*pi*m(i_)/N_SITES_Y) ! corresponding to the |nq,n'q'> state
                !write(*,*)k2,E_MBX(k2),E_MBY(k2), alpha,pi,N_SITES_Y,i_,m(i_)
            ELSE
            END IF
         END DO
         !WRITE(*,*) STATES(j)%n(1), STATES(j)%m(1), STATES(j)%n(2), STATES(j)%m(2),STATES(j)%occupation(1),STATES(j)%occupation(2)
         !-  The D states labelled by k2 and corresponding to j \in [i_off +1, i_off + D] are superpositions of the D states labelled by i \in [i_offf + 1 , i_off + D].-
         !    STATES(j)%n(1) : Position (x) of one particle of the state j
         !    STATES(j)%n(2) : Position (x) of another particle of the state j
         !    | 1@n(1) 1@n(2)> = sum_{l.ne.m} (U(l,1@n(1) U(m,1@n(2))+U(m,1@n(1) U(l,1@n(2))))) + sum_{m} (U(m,1@n(1) U(m,1@n(1))+U(m,1@n(1) U(l,1@n(2)))
         !     in the new basis
         !   U_AUX(STATES(i)%n(1),STATES(j)%n(1))
         !DO i=1,BASIS_DIM ! loop over the old basis |nq,n'q'>


    END DO

  !END DO

  ! ----------------- END -----------------
  ! ----- BUILD THE TWO BODY SPECTRUM
  ! ----- (new) {|k,q>} = U {|n,q>} (old) -----
  ! ----------------- END -----------------




!109 FORMAT(9I4)
!1012 FORMAT(12I4)

END SUBROUTINE TWO_BODY_ENERGY_SPECTRUM

