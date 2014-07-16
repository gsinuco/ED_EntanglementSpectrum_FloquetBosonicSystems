DOUBLE PRECISION FUNCTION  D_H(N_SITES,N_PARTICLES) 
IMPLICIT NONE

INTEGER, INTENT(IN):: N_SITES,N_PARTICLES
!DOUBLE PRECISION D_H

INTEGER          :: N_,K_,i
DOUBLE PRECISION :: Num,Den



N_ = N_SITES+N_PARTICLES-1;
K_ = N_PARTICLES;

Num = 1.0;
do i = 1,k_
    Num = Num*n_;
    n_ = n_-1;
end do
Den = 1.0;
do i = 1,k_-1
    Den = Den*k_;
    k_ = k_-1;
end do

D_H = Num/Den;

END FUNCTION
