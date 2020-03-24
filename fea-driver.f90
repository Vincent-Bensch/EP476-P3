
MODULE fea_init !Module to set up variables
IMPLICIT NONE
    !Variable type definitions
    INTEGER, PARAMETER :: rknd=SELECTED_REAL_KIND(14,14)
    INTEGER, PARAMETER :: iknd=SELECTED_INT_KIND(14)
    
    !Initializing input variables
    REAL(KIND=rknd) :: rho_0, rho_1, rho_i, E, L !rho_0 and rho_1 are the density at the fixed end of the beam, and the slope of the density line respectively. E is the  modulus of elasticity, and L is the length of the beam. rho_i is the density at the point of iteration
    INTEGER(KIND=iknd) :: num_elem !Intiger to store the number of elements the beam is to be devided up into.
    
    !Initializing process variables
    REAL(KIND=rknd) :: dx !element x step allong L
    REAL(KIND=rknd) :: tmp !dummy computation variable
    REAL(KIND=rknd) :: K_mul
    REAL(KIND=rknd), DIMENSION(:, :), ALLOCATABLE :: M, K !Initialize M and K matricies
    REAL(KIND=rknd), DIMENSION(:, :), ALLOCATABLE :: MB, KB !Banded M and K
    REAL(KIND=rknd), DIMENSION(2, 2) :: K_elem, K_ins, M_elem !Sub matrix for current computation

    INTEGER(KIND=iknd) :: i,j !integer for loop iteration

    !Variabled to catch DSBGV output
    REAL(KIND=rknd), DIMENSION(:), ALLOCATABLE :: eig_out, work_out
    REAL(KIND=rknd), DIMENSION(1, 1) :: z_out
    INTEGER(KIND=iknd) :: info_out, one
    
END MODULE fea_init


PROGRAM fea_driver !Main program
USE fea_init
IMPLICIT NONE

    !Get parameters from user
    PRINT *,"Enter the density at x=0 (rho_0):"
    READ *,rho_0
    PRINT *,"Enter the density slope (rho_1):"
    READ *,rho_1
    PRINT *,"Enter the beam length (L):"
    READ *,L
    PRINT *,"Enter the modulus of elasticity (E):"
    READ *,E
    PRINT *,"Enter the desired number of elements:"
    READ *,num_elem
    
    dx = REAL(L)/REAL(num_elem) !Compute dx
    
    ALLOCATE ( M(num_elem,num_elem) ) !Alocate size for M and K matricies
    ALLOCATE ( K(num_elem,num_elem) )
    ALLOCATE ( MB(num_elem,2) ) !Alocate size for MB and KB matricies
    ALLOCATE ( KB(num_elem,2) )
    
    M=0 !Set matricies to all 0s
    K=0
    
    M_elem = 1/6. !Set M element contributions
    M_elem(1,1) = 1/3.
    M_elem(2,2) = 1/3.
    
    K_elem = -1  !Set K element contributions
    K_elem(1,1) = 1
    K_elem(2,2) = 1
    
    M(1,1) = dx * rho_0 / 3 !Set top left matrix element from first beam element
    j = 1 ! J counter should always be one less than i
    DO i = 2,num_elem !Loop to create M matrix
        rho_i = rho_0 + (REAL(i) / REAL(num_elem)) * rho_1 !Compute density at i
        M(j:i , j:i) = M(j:i, j:i) + M_elem * ( dx * rho_i ) !Compute M element and add to M matrix
        j = i !Set j counter for next loop
    END DO
    
    K_mul = E/dx !Set contribution for each element in the K matrix
    K_ins = K_mul * K_elem !Set insertion element for K matrix
    K(1,1) = k_mul !Set top left matrix element comtribution from first beam element
    j=1
    DO i = 2,num_elem !Loop to create K matrix
        K(j:i , j:i) = K(j:i, j:i) +  K_ins
        j = i 
    END DO
	
    MB = bandify(num_elem, M)
    KB = bandify(num_elem, K)
    
    ALLOCATE ( eig_out( num_elem ) )
    ALLOCATE ( work_out( 3*num_elem ) )
    one = 1
    
    !dsbgv(   JOBZ, UPLO,       N,  KA,  KB, AB,    LDAB, BB,    LDBB,       W,     Z, LDZ,     WORK, INFO )
    CALL dsbgv('N', 'L', num_elem, one, one, KB, num_elem, MB, num_elem, eig_out, z_out, one, work_out, info_out)
    
    DO i =1, SIZE(eig_out) 
        PRINT *, eig_out(i)
    END DO    

    PRINT *, "INFO:"
    PRINT *, info_out

CONTAINS
    FUNCTION bandify(mat_size, mat_in)
        INTEGER(KIND=iknd) :: mat_size !the n by n size of mat_in
    	REAL(KIND=rknd), INTENT(IN), DIMENSION(mat_size, mat_size) :: mat_in !Ordinary matrix in
    	REAL(KIND=rknd), DIMENSION(mat_size, 2) :: bandify !Banded matrix out
    	INTEGER(KIND=iknd) :: ib !Iteration counter
	
	    DO ib = 1, (mat_size - 1)
	        bandify(ib,1) = mat_in(ib,ib)
	        bandify(ib,2) = mat_in(ib,ib+1)
    	END DO
	    bandify(mat_size,1) = mat_in(mat_size,mat_size)	
    END FUNCTION bandify
END PROGRAM fea_driver
