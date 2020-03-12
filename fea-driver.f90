MODULE fea_init
IMPLICIT NONE
    !Type parameters
    INTEGER, PARAMETER :: rknd=SELECTED_REAL_KIND(14)
    INTEGER, PARAMETER :: iknd=SELECTED_INT_KIND(7)
    
    !Initializing variables
    REAL(KIND=rknd) :: rho_0, rho_L, E, L
    INTEGER(KIND=iknd) :: num_elem
END MODULE fea_init


PROGRAM fea_driver
USE fea_init
IMPLICIT NONE

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
    
    print '(*(A))', NEW_LINE('a'), NEW_LINE('a')
    print "(es10.3)",rho_0
    print "(es10.3)",rho_L
    print "(es10.3)",L
    print "(es10.3)",E
    print "(i7)",num_elem
          
END PROGRAM fea_driver
    
