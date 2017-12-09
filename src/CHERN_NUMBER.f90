SUBROUTINE ChernNumber()

  USE MAGMA
  USE physical_constants
  USE LATTICE
  USE INTERACTION
  USE TIME_DEPENDENT
  USE INTERFACE
  USE STATE_OBJ

  IMPLICIT NONE
  TYPE(S    TATE),      DIMENSION(:),                INTENT(IN)    :: STATES
  INTEGER,                                           INTENT(INOUT) :: Q_T
  COMPLEX*16,           DIMENSION(:),   ALLOCATABLE, INTENT(OUT)   :: EIGENVALUES
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE, INTENT(OUT)   :: U_MB
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE, INTENT(OUT)   :: U_FTY
  INTEGER,                                       INTENT(INOUT) :: INFO





  END SUBROUTINE
