MODULE MATRICES
IMPLICIT NONE
  
!********************************** TYPES ****************
  
  TYPE DMATRIX
     INTEGER, DIMENSION(2)::N
     DOUBLE COMPLEX det
     DOUBLE COMPLEX tr
     DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE:: elem !this is the whole matrix
  END TYPE DMATRIX
  
!*******************************
  
  TYPE DVECTOR
     INTEGER :: N
     DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE:: elem !this is the whole vector
  END TYPE DVECTOR

!********************************
  TYPE RMATRIX
	INTEGER, DIMENSION(2) :: N
	REAL*8 det
	REAL*8 tr
	REAL*8, DIMENSION(:,:), ALLOCATABLE:: elem
  END TYPE RMATRIX

!*********************************
TYPE DDMATRIX
     INTEGER, DIMENSION(2)::N
     COMPLEX*16 det
     COMPLEX*16 tr
     COMPLEX*16, DIMENSION(:,:), ALLOCATABLE:: elem
  END TYPE DDMATRIX
  
! ******************************INTERFACES***************************
  
  INTERFACE OPERATOR (.ADJ.)
   MODULE PROCEDURE MATADJOINT,CONJUGATE,VECADJOINT
END INTERFACE OPERATOR (.ADJ.)

INTERFACE OPERATOR (.TR.)
   MODULE PROCEDURE TRACE, RTRACE
END INTERFACE OPERATOR (.TR.)

!******************************
  
CONTAINS

! ******************************  FUNCTION INITV, INITIALIZES THE VECTOR ****************
  FUNCTION INITV(length)
    INTEGER :: ii, length
    TYPE(DVECTOR) :: V, INITV
    V%N = length
    ALLOCATE(V%elem(length))
    DO ii=1,length
       V%elem(ii) = complex(1d0*ii,1d0*ii)
    END DO
    INITV = V
    RETURN
  END FUNCTION INITV

! ******************************  FUNCTION RANDOM_INITV, INITIALIZES THE VECTOR ****************
  FUNCTION RANDOM_INITV(length)
    INTEGER, INTENT(IN) :: length
    INTEGER :: ii
    REAL*8 xx,yy
    DOUBLE COMPLEX :: temp

    TYPE(DVECTOR) :: RANDOM_INITV
    RANDOM_INITV%N = length
    ALLOCATE(RANDOM_INITV%elem(length))
    temp = 0.0

    DO ii=1,length
    CALL RANDOM_NUMBER(xx)
    CALL RANDOM_NUMBER(yy)
       RANDOM_INITV%elem(ii) = complex(xx,yy)
       temp = temp + (RANDOM_INITV%elem(ii))*conjg(RANDOM_INITV%elem(ii))
    END DO

    RANDOM_INITV%elem = RANDOM_INITV%elem/cdsqrt(temp)
 
    RETURN
  END FUNCTION RANDOM_INITV

!**********************************  FUNCTION VEC_NORMALIZE  ******************************
 FUNCTION VEC_NORMALIZE(vector, length)
  INTEGER, INTENT(IN) :: length
  TYPE(DVECTOR), INTENT(IN) :: vector
  TYPE(DVECTOR) :: vec_normalize
  INTEGER :: ii
  DOUBLE COMPLEX :: temp

  temp = complex(0d0,0d0)
  DO ii=1, length
     temp = temp + (vector%elem(ii))*conjg(vector%elem(ii))
  END DO

  vec_normalize%elem = vector%elem/cdsqrt(temp)
  vec_normalize%N = length
  RETURN
  END FUNCTION VEC_NORMALIZE

!*****************************   FUNCTION INIT, INITIALIZES THE DMATRIX ***************
  
  FUNCTION INIT(nrow,ncol)
    INTEGER, INTENT(IN) :: nrow,ncol   !stuff MUST be declared AS FIRST. Otherwise, ERROR
    INTEGER :: ii,jj,x
    TYPE(DMATRIX) :: INIT
    INIT%N(1)=nrow
    INIT%N(2)=ncol
    ALLOCATE (INIT%elem(nrow,ncol))
     DO ii=1,nrow
         DO jj=1,ncol
            INIT%elem(ii,jj) = complex(1d0*ii,1d0*jj)
         ENDDO
      ENDDO
      INIT%tr = .TR.(INIT)    ! this calls the trace function to calculate it
      x=0
      INIT%det = complex(0.0/x,0.0/x) ! To be implemented. Meanwhile returns NAN
    RETURN
  END FUNCTION INIT

!********************** FUNCTION REAL TRACE

FUNCTION RTRACE(A)
  TYPE(RMATRIX), INTENT(IN) :: A
  INTEGER :: ii,x
  REAL :: RTRACE
  RTRACE=(0d0)
  
  IF (A%N(1) /= A%N(2)) THEN
     x=0
     RTRACE = 0.0/x
     print*, 'Warning: not squared matrix given in input for TRACE calculus'
     
  ELSE
  DO ii=1,A%N(1)
     RTRACE = RTRACE + A%elem(ii,ii)
     
  END DO
  RETURN
END IF
END FUNCTION RTRACE


!***********************************  FUNCTION TRACE **************
  
FUNCTION TRACE(A)
  TYPE(DMATRIX), INTENT(IN) :: A
  INTEGER :: ii,x
  DOUBLE COMPLEX :: TRACE
  TRACE=(0d0,0d0)
  
  IF (A%N(1) /= A%N(2)) THEN
     x=0
     TRACE = complex(0.0/x,0.0/x)
     !print*, 'Warning: not squared matrix given in input for TRACE calculus'
     
  ELSE
  DO ii=1,A%N(1)
     TRACE = TRACE + A%elem(ii,ii)
     
  END DO
  RETURN
END IF
END FUNCTION TRACE

! ********************************** FUNCTION MATADJOINT *************************

FUNCTION MATADJOINT(A)
  INTEGER :: ii,jj,x
  TYPE(DMATRIX), INTENT(IN):: A   !I MUST declare what is the INPUT, with INTENT(IN)
  TYPE(DMATRIX) ::  MATADJOINT
  ALLOCATE(MATADJOINT%elem(A%N(2),A%N(1))) !output matrix has nrow & ncol swapped

  MATADJOINT%N(1) = A%N(2)      ! these must be RE-INITIALIZED!!
  MATADJOINT%N(2) = A%N(1)
  
  DO ii=1,A%N(2)
     DO jj=1,A%N(1)
        MATADJOINT%elem(ii,jj) = conjg(A%elem(jj,ii))
     END DO
  END DO
  
  MATADJOINT%tr = .TR.(MATADJOINT)    ! TRACE MUST BE RE-CALCULATED!!
  x=0
  MATADJOINT%det = complex(0.0/x,0.0/x) ! To be implemented. Meanwhile returns NAN

  RETURN
END FUNCTION MATADJOINT

! ******************************** FUNCTION CONJUGATE ************************

FUNCTION CONJUGATE(A)
  DOUBLE COMPLEX, INTENT(IN) :: A
  DOUBLE COMPLEX :: CONJUGATE
  CONJUGATE = conjg(A)
  RETURN
END FUNCTION CONJUGATE

! ***************************** FUNCTION VECADJOINT **************************

FUNCTION VECADJOINT(A)
  INTEGER ::  ii
  TYPE(DVECTOR), INTENT(IN) :: A
  TYPE(DVECTOR) ::  VECADJOINT
  VECADJOINT%N = A%N              !it must be set again
  ALLOCATE(VECADJOINT%elem(A%N))
  DO ii=1,A%N
     VECADJOINT%elem(ii) = conjg(A%elem(ii))
  END DO
  RETURN
END FUNCTION VECADJOINT

!*********************************  PRINTING SUBROUTINE ********************

SUBROUTINE PRINTFILE(OBJECT,filename)
  
  CHARACTER(LEN=*), PARAMETER :: FMT2 = "(A,I0,'X',I0)"
  CHARACTER(LEN=*) :: filename
  TYPE(RMATRIX):: OBJECT
  INTEGER ii,jj
  OPEN(UNIT=50,FILE=filename,status='unknown')
  WRITE(50,*) 'Trace: ', OBJECT%tr
  WRITE(50,*) 'Determinant: ',OBJECT%det
  WRITE(50,FMT2) 'Dimensions: ',OBJECT%N(1),OBJECT%N(2)
  WRITE(50,*) 'Matrix: '
  DO ii=1,OBJECT%N(1)
        WRITE(50,*) OBJECT%elem(ii,:) !this prints an entire row
  END DO
  CLOSE(50)
  
  RETURN
 
END SUBROUTINE PRINTFILE
!***************************************
! *******     "ISING MODULE"    **********************************
 !***********************************************************************
  ! diagonalization
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

  CALL DSYEV('V','U',NN,ham%elem,NN,eigenvalues,WORK,LWORK,INFO)

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

!**************************************************************
SUBROUTINE DIAGONALIZE_ONLYVEC(ham)
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK, WORK, eigenvalues
    TYPE(RMATRIX) :: ham
    INTEGER :: LWORK, INFO, ii, NN
    LOGICAL debug

    NN = ham%N(1)
  debug = .false.
  LWORK = NN*3+NN/2
  
  ALLOCATE(WORK(LWORK), RWORK(3*NN-2)) 
  ALLOCATE(eigenvalues(NN))

  CALL DSYEV('V','U',NN,ham%elem,NN,eigenvalues,WORK,LWORK,INFO)

DEALLOCATE(WORK, RWORK, eigenvalues)

END SUBROUTINE

    !***************************************************************************
    ! pauli matrices
    !***************************************************************************
FUNCTION pauli_xx()
   TYPE(RMATRIX) :: pauli_xx
   ALLOCATE(pauli_xx%elem(2,2))

   pauli_xx%N(1) = 2
   pauli_xx%N(2) = 2

   pauli_xx%elem(1,1) = 0.0
   pauli_xx%elem(1,2) = 1.0
   pauli_xx%elem(2,1) = 1.0
   pauli_xx%elem(2,2) = 0.0

   RETURN
END FUNCTION

!****************************************************************
FUNCTION pauli_zz()
    TYPE(RMATRIX) :: pauli_zz

    ALLOCATE(pauli_zz%elem(2,2))
    
    pauli_zz%N(1) = 2
    pauli_zz%N(2) = 2

    pauli_zz%elem(1,1) = 1.0
    pauli_zz%elem(1,2) = 0.0
    pauli_zz%elem(2,1) = 0.0
    pauli_zz%elem(2,2) = -1.0

    RETURN
 END FUNCTION
  
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

END MODULE MATRICES

!*********************************     MAIN     ***********************************

  !PROGRAM MAIN
  !USE MATRICES        !it gives an error if you swap USE.. and IMPLICIT NONE
  !IMPLICIT NONE
  
  !TYPE(DMATRIX) :: mat,admat
  
  !mat = INIT(2,2)        !init function takes as an input the dimensions
  !admat = .ADJ.(mat)
  
  !CALL PRINTFILE(mat,'matrix') !takes as input a DMATRIX object and the FILENAME
  !CALL PRINTFILE(admat,'adjointmatrix')
  
 !****************************** TEST AND OTHERS *****************

!  TYPE(DVECTOR) :: vec 
!  test = complex(1.0,2.0)     !Number used just to test .ADJ. operator
!  vec = INITV(4)              !Vector used just to test .ADJ. operator
!  test = .ADJ.(test)
!  vec = .ADJ.(vec)
  
  ! REMEMBER: fortran prints matrices BY COLUMN
  ! REMEMBER: to put the '%elem' when printing a VECTOR or a MATRIX
  
! print*,mat%elem
! print*, admat%elem
  ! print*,test
 ! print*, vec%elem

!****************************  END PROGRAM  *************************
  
!STOP
!END PROGRAM MAIN




!*****************************************************************************
! Failed output format attempts. It seems that you have to know a priori the length of the output in order to have a nice print. It compiles, but the graphical result isn't better from the one implemented in the code.
!CHARACTER(LEN=*), PARAMETER :: FMT1 = "(A,'(',F18.9,',',F18.9,')')"
!CHARACTER(LEN=*), PARAMETER :: FMT3 = "('(',F18.9,',',F18.9,')',$)"
! CHARACTER(LEN=*), PARAMETER :: FMT4 = "('(',*,',',*,')',$)"  
