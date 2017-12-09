SUBROUTINE fromKQtoNM(U_MB,U_AUX,U_NM)
    IMPLICIT NONE
    COMPLEX*16, DIMENSION(:,:), INTENT(IN)  :: U_MB
    COMPLEX*16, DIMENSION(:,:), INTENT(IN)  :: U_AUX
    COMPLEX*16, DIMENSION(:,:), INTENT(OUT) :: U_NM

    INTEGER k,k_,n,n_

    DO n=1,N_SITES_X
        DO n_=1,N_SITES_X
            U_NM(n,n_) = DCMPLX(0.0,0.0)
            DO k=1,N_SITES_X
                DO k_=1,N_SITES_X

                U_NM(,n,n_) = U_NM(n,n_) + U_MB()*U_AUX(n,k)*U_AUX(n_,k_)

                END DO
            END DO


        END DO
    END DO


    j = i_off

    U_NM = DCMPLX(0.0,0.0)


    END DO

END SUBROUTINE fromKQtoNM
