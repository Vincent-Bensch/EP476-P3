MODULE fea_init !Module to set up variables
IMPLICIT NONE
    !Variable type definitions
    INTEGER, PARAMETER :: rknd=SELECTED_REAL_KIND(14,14)
    INTEGER, PARAMETER :: iknd=SELECTED_INT_KIND(14)
    
    !Initializing input variables
    REAL(KIND=rknd) :: rho_o, rho_L, E, L !rho_0 and rho_L are the densities at the fixed and free end of the beam respectivly. E is the modulus of elasticity, and L is the length of the beam.
    INTEGER(KIND=iknd) :: num_elem !Intiger to store the number of elements the beam is to be devided up into.
    
    !Initializing process variables
    REAL(KIND=rknd) :: rho_m !The density slope heading from fixed to free
    REAL(KIND=rknd) :: dx !element x step allong L
    REAL(KIND=rknd) :: rho_b !The density y offset accounting for the first element not being #0
    REAL(KIND=rknd) :: tmp !dummy computation variable
    REAL(KIND=rknd), DIMENSION(:, :), ALLOCATABLE :: M, K !Initialize M and K matricies
    REAL(KIND=rknd), DIMENSION(2, 2) :: element, K_elem, M_elem !Sub matrix for current computation

    INTEGER(KIND=iknd) :: i,j !integer for loop iteration

    !Variabled to catch DSBGV output
    REAL(KIND=rknd), DIMENSION(:), ALLOCATABLE :: eig_out, work_out
    REAL(KIND=rknd), DIMENSION(1, 1) :: z_out
    INTEGER(KIND=iknd) :: info_out, one
    
    !INTERFACE
    !    SUBROUTINE DSBGV(JOBZ, UPLO, N, KA, KB, AB, LDAB, BB, LDBB, W, Z, LDZ, WORK, INFO)
    !        INTEGER, PARAMETER :: rknd=SELECTED_REAL_KIND(14,14)
    !        INTEGER, PARAMETER :: iknd=SELECTED_INT_KIND(14)

    !        INTEGER(KIND=iknd), INTENT(IN) :: N, KA, KB, LDAB, LDBB, LDZ
    !        REAL(KIND=rknd), DIMENSION(LDAB, *), INTENT(INOUT) :: AB
    !        REAL(KIND=rknd), DIMENSION(LDBB, *), INTENT(INOUT) :: BB
    !        CHARACTER, INTENT(IN) :: JOBZ, UPLO
    !        REAL(KIND=rknd), DIMENSION(*), INTENT(OUT) :: W, WORK
    !        REAL(KIND=rknd), DIMENSION(LDZ, *), INTENT(OUT) :: Z
    !        INTEGER(KIND=iknd), INTENT(OUT) :: INFO
    !    END SUBROUTINE DSBGV  
    !END INTERFACE

END MODULE fea_init


PROGRAM fea_driver !Main program
USE fea_init
IMPLICIT NONE

    !Get parameters from user
    PRINT *,"Enter the density at x=0 (rho_0):"
    READ *,rho_o
    PRINT *,"Enter the density at x=L (rho_L):"
    READ *,rho_L
    PRINT *,"Enter the beam length (L):"
    READ *,L
    PRINT *,"Enter the modulus of elasticity (E):"
    READ *,E
    PRINT *,"Enter the desired number of elements:"
    READ *,num_elem
    
    rho_m = (rho_L - rho_o)/(num_elem - 1) !Computing density slope
    rho_b = rho_L - rho_m !Computing density y offset
    dx = L/num_elem !Compute dx
    
    ALLOCATE ( M(num_elem,num_elem) ) !Alocate size for M and K matricies
    ALLOCATE ( K(num_elem,num_elem) )
    
    M=0 !Set matricies to all 0s
    K=0
    
    K_elem = 1/6 !Set K element contributions
    K_elem(1,1) = 1/3 
    K_elem(2,1) = 1/3
    
    M_elem = -1  !Set M element contributions
    M_elem(1,1) = 1
    M_elem(2,2) = 1
    
    M(1,1) = dx * rho_o / 3 !Set top left matrix element from first beam element
    DO i = 2,num_elem !Loop to create M matrix
        j = i-1
        tmp = dx * ( ( rho_m * i ) + rho_b)
        M(j:i , j:i) = M_elem * tmp
    END DO
    
    tmp = E/dx !Set contribution for each element in the K matrix
    K(1,1) = tmp !Set top left matrix element from first beam element
    DO i = 2,num_elem !Loop to create K matrix
        j = i-1
        M(j:i , j:i) = K_elem * tmp
    END DO
    
    ALLOCATE ( eig_out( num_elem ) )
    ALLOCATE ( work_out( 3*num_elem ) )
    one = 1
    CALL dsbgv('N', 'L', num_elem, one, one, K, num_elem, M, num_elem, eig_out, z_out, one, work_out, info_out)
    
    DO i =1, SIZE(eig_out) 
        PRINT *, eig_out(1)
    END DO    
END PROGRAM fea_driver
