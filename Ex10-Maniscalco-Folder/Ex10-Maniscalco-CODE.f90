MODULE RG
    USE MATRICES
    IMPLICIT NONE
    CONTAINS

  !***********************************************************************
  ! diagonaliztion
  !***********************************************************************
  SUBROUTINE DIAGONALIZE(ham, eigenvalues)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
    COMPLEX(8), DIMENSION(:), ALLOCATABLE :: WORK
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eigenvalues
    TYPE(DMATRIX) :: ham
    INTEGER :: LWORK, INFO, NN

    NN = ham%N(1)
    LWORK = 2*NN - 1

  
  ALLOCATE(WORK(LWORK), RWORK(3*NN-2)) 
  ALLOCATE(eigenvalues(NN))

  CALL ZHEEV('V','U',NN,ham%elem,NN,eigenvalues,WORK,LWORK,RWORK,INFO)

  DEALLOCATE(WORK, RWORK)

END SUBROUTINE

SUBROUTINE DIAGONALIZE_ONLYVEC(ham)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
  COMPLEX(8), DIMENSION(:), ALLOCATABLE :: WORK
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eigenvalues
  TYPE(DMATRIX) :: ham
  TYPE(DMATRIX) :: debug_mat
  INTEGER :: LWORK, INFO, NN
  LOGICAL debug
 
  debug = .false.

  IF(debug.eqv..true.) THEN

    NN = 2

    debug_mat%N(1) = NN
    debug_mat%N(2) = NN
    ALLOCATE(debug_mat%elem(NN,NN))
    debug_mat%elem = complex(0.0,0.0)
    debug_mat%elem(1,1) = complex(-3.0,0.0)
    debug_mat%elem(1,2) = complex(1.0,0.0)
    debug_mat%elem(2,2)= complex(1.0,0.0)
  
    LWORK = 2*NN
    ALLOCATE(WORK(LWORK), RWORK(LWORK)) 
    ALLOCATE(eigenvalues(NN))
    CALL ZHEEV('V','U',NN,debug_mat%elem,NN,eigenvalues,WORK,LWORK,RWORK,INFO)
    CALL PRINTFILE(debug_mat,'zheevdebug.dat')
    DEALLOCATE(WORK, RWORK, eigenvalues)
    
  END IF

  NN = ham%N(1)
  LWORK = 2*NN-1
  ALLOCATE(WORK(LWORK), RWORK(3**NN-2)) 
  ALLOCATE(eigenvalues(NN))

  CALL ZHEEV('V','U',NN,ham%elem,NN,eigenvalues,WORK,LWORK,RWORK,INFO)
  DEALLOCATE(WORK, RWORK, eigenvalues)

END SUBROUTINE

    !***************************************************************************
    ! pauli matrices
    !***************************************************************************
FUNCTION pauli_xx()
   TYPE(DMATRIX) :: pauli_xx
   ALLOCATE(pauli_xx%elem(2,2))

   pauli_xx%N(1) = 2
   pauli_xx%N(2) = 2

   pauli_xx%elem(1,1) = complex(0.0,0.0)
   pauli_xx%elem(1,2) = complex(1.0,0.0)
   pauli_xx%elem(2,1) = complex(1.0,0.0)
   pauli_xx%elem(2,2) = complex(0.0,0.0)

   RETURN
END FUNCTION

FUNCTION pauli_zz()
    TYPE(DMATRIX) :: pauli_zz

    ALLOCATE(pauli_zz%elem(2,2))
    
    pauli_zz%N(1) = 2
    pauli_zz%N(2) = 2

    pauli_zz%elem(1,1) = complex(1.0,0.0)
    pauli_zz%elem(1,2) = complex(0.0,0.0)
    pauli_zz%elem(2,1) = complex(0.0,0.0)
    pauli_zz%elem(2,2) = complex(-1.0,0.0)

    RETURN
 END FUNCTION
  
!***********************************************************************
  ! ising hamiltonian
  !***********************************************************************

  FUNCTION HAM_ISING(NN, lambda)
    TYPE(DMATRIX) :: prev_id, next_id, sigma_z, sigma_x, ham_z, ham_x, temp, temp2, temp3,ham_ising
    INTEGER :: NN, pp
    REAL*4 :: lambda

    sigma_x = pauli_xx()
    sigma_z = pauli_zz()

    ALLOCATE(ham_z%elem(2**NN,2**NN))
    ALLOCATE(ham_x%elem(2**NN,2**NN))
    ALLOCATE(ham_ising%elem(2**NN,2**NN))
    ham_x%N(1) = 2**NN
    ham_x%N(2) = 2**NN
    ham_z%N(1) = 2**NN
    ham_z%N(2) = 2**NN
    ham_ising%N(1) = 2**NN
    ham_ising%N(2) = 2**NN   

    ham_z%elem = complex(0.0,0.0)
    ham_x%elem = complex(0.0,0.0) 
    ham_ising%elem = complex(0.0,0.0)

    DO pp=1, NN
      prev_id = IDENTITY(2**(pp-1))
      next_id = IDENTITY(2**(NN-pp))

      temp = TENSOR_PRODUCT(prev_id, sigma_z)
      temp2 = TENSOR_PRODUCT(temp, next_id)

      ham_z%elem = ham_z%elem + temp2%elem

      DEALLOCATE(prev_id%elem, next_id%elem, temp%elem, temp2%elem)
    END DO

    DO pp=1, NN-1
      prev_id = IDENTITY(2**(pp-1))
      next_id = IDENTITY(2**(NN-pp-1))

      temp = TENSOR_PRODUCT(prev_id,sigma_x)
      temp2 = TENSOR_PRODUCT(temp,sigma_x)
      temp3 = TENSOR_PRODUCT(temp2,next_id)

      ham_x%elem = ham_x%elem + temp3%elem

      DEALLOCATE(prev_id%elem, next_id%elem, temp%elem, temp2%elem, temp3%elem)
    END DO
    
    ham_z%elem = lambda*ham_z%elem
    ham_ising%elem = ham_x%elem + ham_z%elem
    ham_ising%tr = .TR.(ham_ising)

    RETURN
    DEALLOCATE(ham_z%elem, ham_x%elem, sigma_x%elem, sigma_z%elem)
  END FUNCTION
  
  !*********************************************************************
  ! TENSOR PRODUCT
  !*********************************************************************
FUNCTION TENSOR_PRODUCT(AA, BB)
  INTEGER :: dima, dimb
  INTEGER :: ib, jb
  TYPE(DMATRIX):: AA, BB
  TYPE(DMATRIX) :: TENSOR_PRODUCT

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
    TYPE(DMATRIX) :: IDENTITY
    INTEGER :: dim, ii
  
    IDENTITY%N(1) = dim
    IDENTITY%N(2) = dim
  
    ALLOCATE(IDENTITY%elem(dim,dim))
    IDENTITY%elem = complex(0.0,0.0)
  
    DO ii=1, dim
      IDENTITY%elem(ii,ii) = complex(1.0,0.0)
    END DO
  
    IDENTITY%tr = complex(dim,0.0)
    RETURN
    DEALLOCATE(IDENTITY%elem)
    END FUNCTION

!********************************************************************************************
    ! RG SUBROUTINE
!********************************************************************************************

SUBROUTINE RG_ALGORITHM(debug, lambda, n_iter, groundstate)
    REAL*4 :: lambda
    TYPE(DMATRIX) :: ham_is, ham_int, adjpp, temp1, temp2, temp3, temp4, left, right, sigma_x
    TYPE(DMATRIX) :: ham_doubled, copy, pp, ham_loop ! da allocare
    TYPE(DMATRIX) :: id_2_nm1, id_2_n
    INTEGER :: NN, n_iter, ii
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eigenvalues
    LOGICAL :: debug
    REAL*8 :: groundstate

    NN = 2

    !here below i must allocate properly all matrices that are not results of operations
    ! that allocate the matrices inside (eg tensor product, matadjoint)

    ham_doubled%N(1) = 2**(2*NN)
    ham_doubled%N(2) = 2**(2*NN)
    ALLOCATE(ham_doubled%elem(ham_doubled%N(1),ham_doubled%N(2)))

    copy%N(1) = 2**(2*NN)
    copy%N(2) = 2**(2*NN)
    ALLOCATE(copy%elem(copy%N(1), copy%N(2)))

    ham_loop%N(1) = 2**NN
    ham_loop%N(2) = 2**NN
    ALLOCATE(ham_loop%elem(ham_loop%N(1), ham_loop%N(2)))

    pp%N(1) = 2**(2*NN)
    pp%N(2) = 2**(NN)
    ALLOCATE(pp%elem(pp%N(1),pp%N(2)))

    id_2_n = IDENTITY(2**NN)
    id_2_nm1 = IDENTITY(2**(NN-1))
    sigma_x = pauli_xx()

    ham_is = HAM_ISING(NN,lambda)   !function calculates starting hamiltonian
    IF(debug.eqv..true.) THEN
      !print*, 'ham ising: ',char(10),ham_is%elem
    END IF

  !*************************************************** first iteration **********************
    temp1 = TENSOR_PRODUCT(ham_is, id_2_n)
    temp2 = TENSOR_PRODUCT(id_2_n, ham_is)
    left = TENSOR_PRODUCT(id_2_nm1, sigma_x)
    right = TENSOR_PRODUCT(sigma_x, id_2_nm1)
    ham_int = TENSOR_PRODUCT(left, right)

    IF(debug.eqv..true.) THEN   !debugging
    CALL PRINTFILE(left, 'leftdebugging.dat')
    END IF

    ham_doubled%elem = temp1%elem + temp2%elem + ham_int%elem
    copy = ham_doubled

    IF(debug.eqv..true.) THEN   !debugging
      CALL PRINTFILE(ham_doubled, 'hamdoubled_debugging.dat')
    END IF

    CALL DIAGONALIZE_ONLYVEC(copy)

    IF(debug.eqv..true.) THEN   !debugging
      CALL PRINTFILE(copy, 'copy_debugging.dat')
    END IF

    pp%elem(1:,1:) = copy%elem(1:,1:2**NN)
    adjpp = MATADJOINT(pp)

    IF(debug.eqv..true.) THEN    !debugging
      CALL PRINTFILE(pp,'pp_debugging.dat')
      CALL PRINTFILE(adjpp, 'ppstar_debugging.dat')
    END IF
    
    ALLOCATE(temp3%elem(1,1)) !stratagemma, o non riesce a deallocare nel ciclo alla 1st iter
    ALLOCATE(temp4%elem(1,1))

    !****************************************** loop ****************************************
   DO ii=1, n_iter
      
      ham_loop%elem = matmul(adjpp%elem,matmul(ham_doubled%elem,pp%elem))
     
      DEALLOCATE(ham_int%elem,temp1%elem,temp2%elem)
      DEALLOCATE(temp3%elem,temp4%elem)

      IF(ii == n_iter) THEN
        exit
      END IF
      
      temp1 = TENSOR_PRODUCT(ham_loop, id_2_n)
      temp2 = TENSOR_PRODUCT(id_2_n, ham_loop)
      temp3 = TENSOR_PRODUCT(id_2_n, left)  !per aggiornare left e right servono i old p e p+
      temp4 = TENSOR_PRODUCT(right, id_2_n)

      IF(debug.eqv..true.) THEN
        OPEN(unit=10,file='dimension_debug.dat',status='replace')
        WRITE(10,*) 'adjoint pp: ', adjpp%N(1), adjpp%N(2), char(10)
        WRITE(10,*) 'temp3: ', temp3%N(1), temp3%N(2)
        WRITE(10,*) 'temp4: ', temp4%N(1), temp4%N(2)
        write(10,*) 'pp: ',pp%N(1), pp%N(2), char(10)
        write(10,*) 'left: ', left%N(1), left%N(2)
        write(10,*) 'right: ', right%N(1), right%N(2)
        CLOSE(10)
      END IF

      left%elem = matmul(adjpp%elem,matmul(temp3%elem,pp%elem))
      right%elem = matmul(adjpp%elem, matmul(temp4%elem,pp%elem))
    
     !print*, ii    DEBUGGING
    
      ham_int = TENSOR_PRODUCT(left, right)

      ham_doubled%elem = temp1%elem + temp2%elem + ham_int%elem
      copy%elem = ham_doubled%elem

      CALL DIAGONALIZE_ONLYVEC(copy)

      pp%elem(1:,1:) = copy%elem(1:,1:2**NN)
      DEALLOCATE(adjpp%elem)
      adjpp = MATADJOINT(pp)

 END DO
DEALLOCATE(copy%elem, ham_doubled%elem, pp%elem, adjpp%elem)

    CALL DIAGONALIZE(ham_loop, eigenvalues)

    IF(debug.eqv..true.) THEN
      print*, 'first three eigenvalues debugging: ', eigenvalues
    END IF

    groundstate = eigenvalues(1)!/(NN*(2**n_iter))

    !print*, 'Ground state: ', groundstate
    DEALLOCATE(ham_loop%elem, eigenvalues)
END SUBROUTINE

!************************************************************************
! several lambdas
!************************************************************************
SUBROUTINE SEVERAL_LAMBDAS(debug)
  INTEGER :: reason, niter
  REAL*4 :: lambda
  REAL*8 groundst
  LOGICAL :: debug
  
  print*, 'Choose the number of iterations n of the RG algorithm. Lattice dimension will be 2^(n+1). '
  READ(*,*) niter
  
  OPEN(UNIT=77,FILE='lambdas.dat',status='old')
  OPEN(UNIT=20,FILE='for_graphics_eigens.dat',status='replace')
  
  OPEN(UNIT=30,FILE='niter.dat',status='replace')
  write(30,*) niter            ! just needed by the python script for graphics
  CLOSE(30)
  
DO
  READ(77,*,iostat=reason) lambda
    IF(reason<0) THEN      !reason<0 means that it reached the end of the file
    EXIT
    ELSE
  
      CALL RG_ALGORITHM(debug, lambda, niter, groundst)
      write(20,*) lambda, groundst
      
    END IF
  END DO
  
  CLOSE(77)
  CLOSE(20)
  
  END SUBROUTINE



END MODULE

PROGRAM MAIN
    USE MATRICES        
    USE RG
    IMPLICIT NONE
    LOGICAL :: debug

    debug = .false.

    CALL SEVERAL_LAMBDAS(debug)

END PROGRAM

