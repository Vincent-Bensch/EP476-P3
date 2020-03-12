MODULE fea_init !Module to set up variables
IMPLICIT NONE
    !Variable type definitions
    INTEGER, PARAMETER :: rknd=SELECTED_REAL_KIND(14,14)
    INTEGER, PARAMETER :: iknd=SELECTED_INT_KIND(14)
    
    !Initializing input variables
    REAL(KIND=rknd) :: rho_0, rho_L, E, L !rho_0 and rho_L are the densities at the fixed and free end of the beam respectivly. E is the modulus of elasticity, and L is the length of the beam.
    INTEGER(KIND=iknd) :: num_elem !Intiger to store the number of elements the beam is to be devided up into.
    
    !Initializing process variables
    REAL(KIND=rknd) :: rho_m !The density slope heading from fixed to free
    REAL(KIND=rknd), DIMENSION(:, :), ALLOCATABLE :: M, K !Initialize M and K matricies
    REAL(KIND=rknd), DIMENSION(2, 2), :: element !Pointer to current matrix
    
END MODULE fea_init


PROGRAM fea_driver !Main program
USE fea_init
IMPLICIT NONE
    !Get parameters from user
    PRINT *,"Enter the density at x=0 (rho_0):"
    READ *,rho_0
    PRINT *,"Enter the density at x=L (rho_L):"
    READ *,rho_L
    PRINT *,"Enter the beam length (L):"
    READ *,L
    PRINT *,"Enter the modulus of elasticity (E):"
    READ *,E
    PRINT *,"Enter the desired number of elements:"
    READ *,num_elem
    
    rho_m = (rho_L - rho_0)/L !Computing density slope
    
    ALLOCATE ( M(num_elem,num_elem) ) !Alocate size for M and K matricies
    ALLOCATE ( K(num_elem,num_elem) )
    
    M=0 !Set matricies to all 0s
    K=0
    
    
    
    print '(*(A))', NEW_LINE('a'), NEW_LINE('a')
    print "(es10.3)",rho_0
    print "(es10.3)",rho_L
    print "(es10.3)",L
    print "(es10.3)",E
    print "(i7)",num_elem


CONTAINS
    
    SUBROUTINE place_elem(elem_in, elem_num_in)
        REAL(KIND=rknd), INTENT(IN) :: elem_in(2,2) !The incoming sub-matrix will always be shape (2x2)
        INTEGER(KIND=iknd), INTENT(IN) :: elem_num_in !Element number associated with sub-matrix
        INTEGER :: i, j !Integers for loop iteration
        INTEGER :: pos_i, pos_j !Integers for position tracking in input matrix
        REAL(KIND=rknd) :: tmp !Real for holding intermediate product
        
        DO i = 1, 2
            DO j = 1,2
                pos_i = elem_num_in - 2 + i
                pos_j = elem_num_in - 2 + i
                
                tmp = ptr(pos_i,pos_j) + elem_in(i,j)
                ptr(pos_i,pos_j) = tmp
            END DO
        END DO
    END SUBROUTINE place_elem
    
END PROGRAM fea_driver
