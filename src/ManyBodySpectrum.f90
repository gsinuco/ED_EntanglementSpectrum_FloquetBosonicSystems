! for now, Just two particle spectrum
! U_AUX:    Eigenvectors of the single particle Hamiltonian
! EIGENVALUES_H: spectrum of the single particle hamiltonian
! U_MB : eigenvectors of the many-body states built upon single particle states
! E_MB : spectrum of the many-body states

SUBROUTINE  MANY_BODY_SPECTRUM(STATES,U_AUX,q,q_,U_MB,D,INFO)

  USE TASKS
  USE physical_constants  
  USE LATTICE
  USE STATE_OBJ
  USE A_DAGGER_

  IMPLICIT NONE
  TYPE(STATE),       DIMENSION(:),   INTENT(IN)                 :: STATES
  COMPLEX*16,        DIMENSION(:,:), INTENT(IN)                 :: U_AUX
  INTEGER,                           INTENT(IN)                 :: q,q_ ! compare with previous version
  COMPLEX*16,        DIMENSION(:,:), INTENT(OUT)                :: U_MB
  INTEGER,                           INTENT(IN)                 :: D
  INTEGER,                           INTENT(INOUT)              :: INFO


  INTEGER :: k2,i_,i_off,i,j
  

  ! ----------------- BEGINING -----------------
  ! ----- BUILD THE MANY-BODY TRANSFORMATION MATRIX ---- 
  ! ----- (new) {|j==k,q>} = U_MB^t {|n,q>} (old) -----
  ! ----- U_MB(:,j)
  ! ----------------- BEGINING -----------------
  
  !DO j=1,BASIS_DIM ! loop over the new basis |kq k'q'>
  ! Since U_MB is block diagonal, the loop over states in the new basis can be replaced by  a loop over blocks of total momentum_Y
  ! the resulting matrix has a smaller size.
  j = 0
  i = 0
  i_off = 0

  !DO k1=Q_T,Q_T!1,(N_SITES_Y*N_SITES_Y + N_SITES_Y)/2 !loop over blocks of total momentum_Y
    ! THE SIZE OF THE BLOCK DEPENDS ON the combination of momentum_Y making up the states
    ! Also, the vector containing the states is ordered following the momentu block, therefore an index offset is needed.
    !
!    DO i_=1,k1
!        j = i_off
!        IF(STATES(j+1)%site(1) .EQ. STATES(j+1)%site(2)) THEN
!            D = (N_SITES_X*N_SITES_X  + N_SITES_X)/2 ! EQUAL MOMENTUM q == q'
!        ELSE
!            D = N_SITES_X*N_SITES_X                  ! DIFFERENT MOMENTUM q != q'
!        END IF
!        IF(i_.LE.k1-1) then
!             i_off = i_off + D
!        END IF
!        q  = STATES(j+1)%m(1)!?
!        q_ = STATES(j+1)%m(2)!?
!    END DO

    IF(q.eq.q_) then
        i_off = 0
        !D     = (N_SITES_X*N_SITES_X  + N_SITES_X)/2
    else
        i_off = (N_SITES_X*N_SITES_X  + N_SITES_X)/2
        !D     = N_SITES_X*N_SITES_X
    end if


    j = i_off

    !WRITE(*,*) D
    !ALLOCATE(U_MB(D,D))
    !ALLOCATE(E_MBX(1))
    !ALLOCATE(E_MBY(1))

    U_MB  = 0
    !write(*,*)
    !write(*,*) 'D, q,q_,i_off',D, q,q_,i_off
    !write(*,*)
    DO k2=1,D

         j = j + 1 ! j is shiftted
         !j=2
         !  write(*,109) STATES(j)
!         DO i_=1,SIZE(STATES(j)%occupation)
!            IF(STATES(j)%occupation(i_) .GT. 0) THEN
!               E_MBX(k2) = E_MBX(k2) + STATES(j)%occupation(i_)*EIGENVALUES_H(STATES(j)%n(i_)) ! corresponding to the |kq,k'q'> state
!               E_MBY(k2) = E_MBY(k2) + STATES(j)%occupation(i_)*2.0*cos(2.0*pi*STATES(j)%n(i_)*alpha &
!                &   + 2*pi*STATES(j)%m(i_)/N_SITES_Y) ! corresponding to the |nq,n'q'> state
!                !   write(*,*)E_MBX(k2)
!            ELSE
!            END IF
!         END DO
         !Write(*,*)
         !Write(*,*) "energies:", E_MBX(k2),E_MBY(k2)
         !Write(*,*)
         !WRITE(*,*) STATES(j)%n(1), STATES(j)%m(1), STATES(j)%n(2), STATES(j)%m(2),STATES(j)%occupation(1),STATES(j)%occupation(2)
         !-  The D states labelled by k2 and corresponding to j \in [i_off +1, i_off + D] are superpositions of the D states labelled by i \in [i_offf + 1 , i_off + D].-
         !    STATES(j)%n(1) : Position (x) of one particle of the state j
         !    STATES(j)%n(2) : Position (x) of another particle of the state j
         !    | 1@n(1) 1@n(2)> = sum_{l.ne.m} (U(l,1@n(1) U(m,1@n(2))+U(m,1@n(1) U(l,1@n(2))))) + sum_{m} (U(m,1@n(1) U(m,1@n(1))+U(m,1@n(1) U(l,1@n(2)))
         !     in the new basis
         !   U_AUX(STATES(i)%n(1),STATES(j)%n(1))
         !DO i=1,BASIS_DIM ! loop over the old basis |nq,n'q'>
         i = i_off

         DO i_=1,D

            i = i + 1

            IF(q.eq.q_) THEN

               ! write(*,109) states(i)
                !write(*,*) i,STATES(i)%site(1),STATES(i)%site(2),STATES(i)%occupation(1),STATES(i)%occupation(2)
                IF( (STATES(i)%m(1) .EQ. STATES(j)%m(1)) .AND. (STATES(i)%m(2) .EQ. STATES(j)%m(2))) THEN
                   !write(*,*) i,STATES(i)%site(1),STATES(i)%site(2),STATES(i)%occupation(1),STATES(i)%occupation(2)
                   IF(STATES(i)%n(1).NE.STATES(i)%n(2)) THEN
                        U_MB(i_,k2) = U_AUX(STATES(i)%n(1),STATES(j)%n(1))*U_AUX(STATES(i)%n(2),STATES(j)%n(2)) + &
                                    & U_AUX(STATES(i)%n(2),STATES(j)%n(1))*U_AUX(STATES(i)%n(1),STATES(j)%n(2))
                   ELSE
                        U_MB(i_,k2) = sqrt(2.0)*U_AUX(STATES(i)%n(1),STATES(j)%n(1))*U_AUX(STATES(i)%n(2),STATES(j)%n(2))
                  !      write(*,109) states(i)
                   END IF
                ELSE
                   write(*,*)'still happening: CHECK CONSERVATION OF MOMENTUM!'
                   U_MB(i_,k2) = DCMPLX(0.0,0.0)
                   INFO = 1
                END IF

            ELSE

                IF( (STATES(i)%m(1) .EQ. STATES(j)%m(1)) .AND. (STATES(i)%m(2) .EQ. STATES(j)%m(2))) THEN
                   !write(*,*) i,STATES(i)%site(1),STATES(i)%site(2),STATES(i)%occupation(1),STATES(i)%occupation(2)
                        U_MB(i_,k2) = U_AUX(STATES(i)%n(1),STATES(j)%n(1))*U_AUX(STATES(i)%n(2),STATES(j)%n(2))
                ELSE
                   write(*,*)'still happening: CHECK CONSERVATION OF MOMENTUM!'
                   U_MB(i_,k2) = DCMPLX(0.0,0.0)
                   INFO = 1
                END IF

            END IF
            !   write(*,*) i_,k2,U_MB(i_,k2)
         END DO
         IF((STATES(j)%m(1) .eq. STATES(j)%m(2)) .AND. (STATES(j)%n(1) .eq. STATES(j)%n(2)) ) then
             U_MB(:,k2) = U_MB(:,k2)/SQRT(2.0)
         END IF
     END DO

  !END DO

  !DEALLOCATE(STATES)
  !DEALLOCATE(TAGS)

  ! ----------------- END -----------------
  ! ----- BUILD THE MANY-BODY TRANSFORMATION MATRIX ---- 
  ! ----- (new) {|k,q>} = U {|n,q>} (old) -----
  ! ----------------- END -----------------




!109 FORMAT(9I4)
!1012 FORMAT(12I4)
END SUBROUTINE MANY_BODY_SPECTRUM
