MODULE USEFUL
  USE MATRICES
  IMPLICIT NONE
  CONTAINS

  !***********************************************************************
  ! diagonaliztion
  !***********************************************************************
  SUBROUTINE DIAGONALIZE(ham, eigenvalues)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK, WORK, eigenvalues
    TYPE(RMATRIX) :: ham
    INTEGER :: LWORK, INFO, ii, NN
    LOGICAL debug

    NN = ham%N(1)
  debug = .false.
  LWORK = NN*3+NN/2
  
  ALLOCATE(WORK(LWORK), RWORK(3*NN-2)) 
  ALLOCATE(eigenvalues(NN))

  CALL DSYEV('N','U',NN,ham%elem,NN,eigenvalues,WORK,LWORK,INFO)

  !OPEN(unit=10,file='eigenvalues.dat',status="replace")
   ! DO ii=1,NN
    !    write(10,*) ii, ' ',eigenvalues(ii)
    !END DO
  !CLOSE(10)

  IF (debug .eqv..true.) THEN
  print*, 'First and last three eigenvalues: '
  print*, eigenvalues(1),eigenvalues(2),eigenvalues(3)
  print*, eigenvalues(NN-2),eigenvalues(NN-1),eigenvalues(NN)
  END IF

  DEALLOCATE(WORK, RWORK)

END SUBROUTINE

  !***********************************************************************
  ! ising hamiltonian
  !***********************************************************************

  FUNCTION HAM_ISING(NN, lambda)
    TYPE(RMATRIX) :: prev_id, next_id, pauli_zz, pauli_xx, ham_z, ham_x, temp, temp2, temp3,ham_ising
    INTEGER :: NN, pp
    REAL*4 :: lambda

    ALLOCATE(pauli_xx%elem(2,2))
    ALLOCATE(pauli_zz%elem(2,2))
    pauli_zz%N(1) = 2
    pauli_zz%N(2) = 2
    pauli_xx%N(1) = 2
    pauli_zz%N(2) = 2

    pauli_xx%elem(1,1) = 0
    pauli_xx%elem(1,2) = 1
    pauli_xx%elem(2,1) = 1
    pauli_xx%elem(2,2) = 0

    pauli_zz%elem(1,1) = 1
    pauli_zz%elem(1,2) = 0
    pauli_zz%elem(2,1) = 0
    pauli_zz%elem(2,2) = -1

    ALLOCATE(ham_z%elem(2**NN,2**NN))
    ALLOCATE(ham_x%elem(2**NN,2**NN))
    ALLOCATE(ham_ising%elem(2**NN,2**NN))
    ham_x%N(1) = 2**NN
    ham_x%N(2) = 2**NN
    ham_z%N(1) = 2**NN
    ham_z%N(2) = 2**NN
    ham_ising%N(1) = 2**NN
    ham_ising%N(2) = 2**NN   

    ham_z%elem = 0.0
    ham_x%elem = 0.0 
    ham_ising%elem = 0.0

    DO pp=1, NN
      prev_id = IDENTITY(2**(pp-1))
      next_id = IDENTITY(2**(NN-pp))

      temp = TENSOR_PRODUCT(prev_id, pauli_zz)
      temp2 = TENSOR_PRODUCT(temp, next_id)

      ham_z%elem = ham_z%elem + temp2%elem

      DEALLOCATE(prev_id%elem, next_id%elem, temp%elem, temp2%elem)
    END DO

    DO pp=1, NN-1
      prev_id = IDENTITY(2**(pp-1))
      next_id = IDENTITY(2**(NN-pp-1))

      temp = TENSOR_PRODUCT(prev_id,pauli_xx)
      temp2 = TENSOR_PRODUCT(temp,pauli_xx)
      temp3 = TENSOR_PRODUCT(temp2,next_id)

      ham_x%elem = ham_x%elem + temp3%elem

      DEALLOCATE(prev_id%elem, next_id%elem, temp%elem, temp2%elem, temp3%elem)
    END DO
    
    ham_z%elem = lambda*ham_z%elem
    ham_ising%elem = ham_x%elem + ham_z%elem
    ham_ising%tr = .TR.(ham_ising)

    RETURN
    DEALLOCATE(ham_z%elem, ham_x%elem, pauli_xx%elem, pauli_zz%elem)
  END FUNCTION
  
  !*********************************************************************
  ! TENSOR PRODUCT
  !*********************************************************************
FUNCTION TENSOR_PRODUCT(AA, BB)
  INTEGER :: dima, dimb
  INTEGER :: ib, jb
  TYPE(RMATRIX):: AA, BB
  TYPE(RMATRIX) :: TENSOR_PRODUCT

  dima = AA%N(1)
  dimb = BB%N(1)

  TENSOR_PRODUCT%N(1) = dima*dimb
  TENSOR_PRODUCT%N(2) = dima*dimb
  
  ALLOCATE(TENSOR_PRODUCT%elem(TENSOR_PRODUCT%N(1),TENSOR_PRODUCT%N(2)))

  DO ib=1, dima     !dima = #of blocks; dimb = dimension of blocks
    DO jb=1, dima
      TENSOR_PRODUCT%elem((ib-1)*dimb+1:(ib+1)*dimb,(jb-1)*dimb+1:(jb+1)*dimb) = AA%elem(ib,jb)*BB%elem
    END DO
  END DO

  TENSOR_PRODUCT%tr = .TR.(TENSOR_PRODUCT)
  RETURN
  DEALLOCATE(TENSOR_PRODUCT%elem)
END FUNCTION

!**********************************************************************
! IDENTITY
!**********************************************************************

FUNCTION IDENTITY(dim)
  TYPE(RMATRIX) :: IDENTITY
  INTEGER :: dim, ii

  IDENTITY%N(1) = dim
  IDENTITY%N(2) = dim

  ALLOCATE(IDENTITY%elem(dim,dim))
  IDENTITY%elem = 0.0

  DO ii=1, dim
    IDENTITY%elem(ii,ii) = 1.0
  END DO

  IDENTITY%tr = float(dim)
  RETURN
  DEALLOCATE(IDENTITY%elem)
  END FUNCTION

!*******************************************************************
! debugging
!*******************************************************************
SUBROUTINE TENSOR_DEBUGGING()
  LOGICAL :: debug
  TYPE(RMATRIX) :: trial, supertrial, tensor
  INTEGER :: dim_t, dim_st

  dim_t = 2
  dim_st = 3

  ALLOCATE(trial%elem(dim_t,dim_t))
ALLOCATE(supertrial%elem(dim_st,dim_st))

trial%N(1) = dim_t
trial%N(2) = dim_t
supertrial%N(1) = dim_st
supertrial%N(2) = dim_st


trial%elem(1,1) = 1.0
trial%elem(1,2) = 0.5
trial%elem(2,1) = 0.0
trial%elem(2,2) = 0.0

supertrial%elem(1,1) = 2.0
supertrial%elem(1,2) = 4.0
supertrial%elem(1,3) = 8.0
supertrial%elem(2,1) = 2.0
supertrial%elem(2,2) = 4.0
supertrial%elem(2,3) = 8.0
supertrial%elem(3,1) = 2.0
supertrial%elem(3,2) = 4.0
supertrial%elem(3,3) = 8.0

tensor = TENSOR_PRODUCT(trial, supertrial)

CALL PRINTFILE(tensor, 'tensor_product_debugging.txt')

DEALLOCATE(trial%elem, supertrial%elem, tensor%elem)

END SUBROUTINE

!************************************************************************
! several lambdas
!************************************************************************
SUBROUTINE SEVERAL_LAMBDAS()
INTEGER :: NN, reason, kk
REAL*4 :: lambda
TYPE(RMATRIX) :: hamilt
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eigens

print*, 'Choose the dimension of the lattice N: '
READ(*,*) NN
print*, 'Choose the number k of eigenvalues to be plotted: '
READ(*,*) kk

IF(kk > 2**NN) THEN
  print*, 'ERROR: received ',kk,' eigenvalues for ',NN,'lattice'
  stop
end if


OPEN(UNIT=77,FILE='lambdas.dat',status='old')
OPEN(UNIT=20,FILE='for_graphics_eigens.dat',status='replace')

OPEN(UNIT=30,FILE='N.dat',status='replace')
write(30,*) NN, kk            ! just needed by the python script for graphics
CLOSE(30)


DO
READ(77,*,iostat=reason) lambda
  IF(reason<0) THEN      !reason<0 means that it reached the end of the file
  EXIT
  ELSE

    hamilt = HAM_ISING(NN,lambda)
    CALL DIAGONALIZE(hamilt, eigens)
    write(20,*) lambda, eigens(1:kk)
    DEALLOCATE(hamilt%elem, eigens)
  END IF
END DO

CLOSE(77)
CLOSE(20)

END SUBROUTINE

END MODULE

!********************************************************************************************************
! main program
!********************************************************************************************************
PROGRAM MAIN
  USE MATRICES        
  USE USEFUL
  IMPLICIT NONE
  
  TYPE(RMATRIX) :: hamilt
  LOGICAL :: debug
  REAL*4 :: lambda
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eigens

  debug = .false.

  IF (debug.eqv..true.) THEN
    CALL TENSOR_DEBUGGING()
  END if

  IF (debug.eqv..true.) THEN
      lambda = 0.0
      hamilt = HAM_ISING(2,lambda)
      CALL DIAGONALIZE(hamilt, eigens)
      CALL PRINTFILE(hamilt, 'debugging_hamiltonian.txt')
      DEALLOCATE(hamilt%elem, eigens)
  END IF

  CALL SEVERAL_LAMBDAS()

 STOP
END PROGRAM MAIN

