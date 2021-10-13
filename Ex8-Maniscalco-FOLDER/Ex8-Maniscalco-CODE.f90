MODULE USEFUL
  USE MATRICES
  IMPLICIT NONE
  CONTAINS

SUBROUTINE DENSITY_MATRIX(vector, DENS_MAT)
  INTEGER :: nrow, ii, jj
  TYPE(DMATRIX):: DENS_MAT
  TYPE(DVECTOR) :: vector

  nrow = vector%N
  DENS_MAT%N(1) = vector%n
  DENS_MAT%N(2) = vector%n
  ALLOCATE (DENS_MAT%elem(nrow,nrow))

  DO ii=1, nrow
    DO jj=1, nrow
      DENS_MAT%elem(ii,jj) = conjg(vector%elem(ii))*(vector%elem(jj))
    END DO
  END DO

  DENS_MAT%tr = TRACE(DENS_MAT)

END SUBROUTINE

SUBROUTINE REDUCED_A(DENS_MAT, RED_A, DD)
  TYPE(DMATRIX) :: DENS_MAT
  TYPE(DMATRIX) :: RED_A 
  INTEGER :: ii, jj, kk ,index1, index2, NN
  INTEGER:: DD

  NN = 2
 
  RED_A%N(1) = DD
  RED_A%N(2) = DD
  ALLOCATE (RED_A%elem(DD,DD))
  RED_A%elem = complex(0d0,0d0)

  DO ii=0, DD-1
    DO jj=0, DD-1
      DO kk=0, DD-1
        index1 = 1 + ii*DD + kk
        index2 = 1 + jj*DD + kk
         RED_A%elem(ii+1,jj+1) = RED_A%elem(ii+1,jj+1) + DENS_MAT%elem(index1,index2) 
      END DO
    END DO
  END DO

  RED_A%tr = TRACE(RED_A)

  END SUBROUTINE

SUBROUTINE REDUCED_B(DENS_MAT, RED_B, DD)
  TYPE(DMATRIX) :: DENS_MAT
  TYPE(DMATRIX) :: RED_B 
  INTEGER :: ii, jj, kk ,index1, index2, NN
  INTEGER :: DD

  NN = 2

  RED_B%N(1) = DD
  RED_B%N(2) = DD
  ALLOCATE (RED_B%elem(DD,DD))
  RED_B%elem = complex(0d0,0d0)

  DO ii=0, DD-1
    DO jj=0, DD-1
      DO kk=0, DD-1
        index1 = 1 + kk*DD + ii
        index2 = 1 + kk*DD + jj
         RED_B%elem(ii+1,jj+1) = RED_B%elem(ii+1,jj+1) + DENS_MAT%elem(index1,index2) 
      END DO
    END DO
  END DO

  RED_B%tr = TRACE(RED_B)
END SUBROUTINE

SUBROUTINE VEC_DEBUGGING(psi, debug,name)
  LOGICAL :: debug
  DOUBLE COMPLEX :: norm
  INTEGER :: dim, ii 
  TYPE(DVECTOR) :: psi
  CHARACTER(len=*) :: name

  norm = complex(0d0,0d0)

  IF (debug .eqv. .true.) THEN
    DO ii=1, psi%N
       norm = norm + (psi%elem(ii))*conjg(psi%elem(ii))
    END DO
    print*, 'Dimension of ',name,' is: ',dim
    print*, 'Norm of ',name,' is: ',norm
END IF
  
  dim = psi%n

 END SUBROUTINE

SUBROUTINE MATRIX_DEBUGGING(rho, debug,name)
  LOGICAL :: debug
  DOUBLE COMPLEX :: trace
  INTEGER :: dim1, dim2
  TYPE(DMATRIX) :: rho
  CHARACTER(len=*) :: name

  trace = rho%tr
  dim1 = rho%N(1)
  dim2 = rho%N(2)

print*, 'Trace of ',name,' is: ', trace
print*, 'Dimensions of ',name, ' are: ', dim1,'',dim2
END SUBROUTINE

END MODULE

!********************************************************************************************************
! main program
!********************************************************************************************************
PROGRAM MAIN
  USE MATRICES        
  USE USEFUL
  IMPLICIT NONE
  
  TYPE(DVECTOR) :: psi, psi_sep, psi2, checkvec1
  TYPE(DMATRIX) :: rho, rhoA, rhoB, checkmatrix, checkA, checkB
  INTEGER*4 :: NN, dd, ii, ss, reason
  LOGICAL :: debug
  REAL*8 :: T1, T2, T3, T4

  debug = .false.

!*******************************************************************************************************
! First part, general states
!*******************************************************************************************************
print*,  '1) GENERAL WAVEFUNCTIONS, calculating efficiency for inputdim.dat data.'!,char(10),' Select the number of bodies N: '
!READ(*,*)  NN
!print*, 'Select the dimension D: '
!READ(*,*) DD

OPEN(UNIT=77,FILE='inputdim.dat',status='old')
ss = 1                           
DO
READ(77,*,iostat=reason) DD, NN
  IF(reason<0) THEN               
  EXIT
  ELSE

  CALL CPU_TIME(T1)
  psi_sep = RANDOM_INITV(NN*dd)
  CALL CPU_TIME(T2)
  CALL CPU_TIME(T3)
  psi = RANDOM_INITV(dd**NN)
  CALL CPU_TIME(T4)

CALL  VEC_DEBUGGING(psi, debug,'general state')
CALL  VEC_DEBUGGING(psi_sep,debug,'separable state')

IF (ss==1) THEN
  OPEN(unit=10,file='separable_times.txt',status='replace')
  OPEN(unit=20,file='general_times.txt',status='replace')
ELSE
  OPEN(unit=10,file='separable_times.txt',status='old')
  OPEN(unit=20,file='general_times.txt',status='old')
END IF
write(10,*) DD, NN, NN*DD, T2-T1
write(20,*) DD, NN, DD**NN, T4-T3

  ss = 2
END IF
END DO

CLOSE(10)
CLOSE(20)

DEALLOCATE(psi%elem)
DEALLOCATE(psi_sep%elem)

!*************************************************************************
  !Second part, denisty matrix and reduced, NN = 2 fixed
!*************************************************************************

  print*, char(10),'2)DENSITY AND REDUCED MATRICES, PURE STATE, N=2. &
            &Select the dimension D:'
NN = 2
READ(*,*) DD

psi2 = RANDOM_INITV(dd**NN)

CALL DENSITY_MATRIX(psi2, rho)
CALL MATRIX_DEBUGGING(rho, debug, 'density matrix')
CALL REDUCED_A(rho,rhoA, dd)
CALL REDUCED_B(rho,rhoB, dd)

print*, 'Printing files: reducedA_random.txt, reducedB_random.txt,random.txt'
CALL PRINTFILE(rhoA, 'reducedA_random.txt')
CALL PRINTFILE(rhoB, 'reducedB_random.txt')
CALL PRINTFILE(rho,'random.txt')

!*****************************************************************************************
! Third part, check two qubits
!*****************************************************************************************

ALLOCATE(checkvec1%elem(4))
checkvec1%N = 4

checkvec1%elem(1) = complex(0.5,0.0)
checkvec1%elem(2) = complex(0.5,0.0)
checkvec1%elem(3) = complex(0.0,0.0)
checkvec1%elem(4) = complex(0.0,0.0)

checkvec1 = VEC_NORMALIZE(checkvec1, 4)

CALL DENSITY_MATRIX(checkvec1, checkmatrix)
CALL REDUCED_A(checkmatrix,checkA,2)
CALL REDUCED_B(checkmatrix,checkB, 2)

print*, char(10),'3)TEST ON TWO SPINS. ',char(10), 'Printing files: reducedA_check.txt, reducedB_check.txt, check.txt'
CALL PRINTFILE(checkA, 'reducedA_check.txt')
CALL PRINTFILE(checkB, 'reducedB_check.txt')
CALL PRINTFILE(checkmatrix,'check.txt')

DEALLOCATE(rho%elem, rhoA%elem, rhoB%elem, checkmatrix%elem, checkA%elem, checkB%elem)
DEALLOCATE(psi2%elem, checkvec1%elem)

STOP
END PROGRAM MAIN

