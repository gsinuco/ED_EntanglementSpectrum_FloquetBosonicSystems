SUBROUTINE a_dagger(i,state_in,state_out)

  USE LATTICE

  IMPLICIT NONE
  INTEGER, DIMENSION(:), INTENT(IN)  :: state_in
  INTEGER, DIMENSION(:), INTENT(OUT) :: state_out
  INTEGER,               INTENT(IN)  :: i


  IF ((SIZE(STATE_in) .EQ. N_SITES)) THEN
     state_out = 0
     IF(state_in(i).LT.N_PARTICLES) THEN
        IF(i.GE.2) THEN
           state_out(1:i-1) = state_in(1:i-1)
        END IF
        state_out(i) = state_in(i) + 1
        IF(i.LT.N_SITES) THEN
           state_out(i+1:N_SITES) = state_in(i+1:N_SITES)
        END IF
     ELSE
        state_out = 0
     END IF
  ELSE
     state_out = 0
  END IF

END SUBROUTINE a_dagger

!!$function [state_] = a_dagger(i,state)
!!$global N L;
!!$
!!$if size(state,2) == L
!!$    state_ = zeros(1,L);
!!$    if(state(i)<N) % state(i) occupation number of site i
!!$        if(i>=2)
!!$            state_(1:i-1) = state(1:i-1);
!!$        end
!!$        state_(i)     = state(i)+1;
!!$        if(i<L)
!!$            state_(i+1:L) = state(i+1:L);
!!$        end
!!$    else
!!$        state_ = 0;
!!$    end
!!$else
!!$    state_ = 0;
!!$end
!!$    
!!$
!!$
!!$end
