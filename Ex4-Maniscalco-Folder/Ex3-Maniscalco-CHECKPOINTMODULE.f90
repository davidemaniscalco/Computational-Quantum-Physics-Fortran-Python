!  Checkpoint subroutines for scalars, arrays, and 2d matrices, for all
!  integers, real, complex double. All the subroutines print a string and
! the variable to be checked.
! ############################################################
MODULE CHECKPOINTS


  INTERFACE CHECKPOINT
     MODULE PROCEDURE LINE, INTEGER_CHECK, REAL_CHECK, COMPLEX_CHECK ,integer_ARRAY_CHECK
     MODULE PROCEDURE real_ARRAY_CHECK,doublecomplex_ARRAY_CHECK,integer_MATRIX_CHECK,real_MATRIX_CHECK,doublecomplex_MATRIX_CHECK
  END INTERFACE CHECKPOINT

CONTAINS     !this must be put AFTER the interface definition

!----------------------------------------------------------
SUBROUTINE LINE(debug,string)
    LOGICAL  debug
    CHARACTER(LEN=*) :: string
    IF(debug.eqv..false.) THEN
       STOP
    ELSE
       print*,string
    END IF
  END SUBROUTINE LINE

!---------------------------------------------------------- 
SUBROUTINE INTEGER_CHECK(debug,string,intvar)
    LOGICAL debug
    CHARACTER(LEN=*) :: string
    INTEGER :: intvar
  IF(debug.eqv..false.) THEN
     RETURN
  ELSE
     print*, string, intvar
  END IF
END SUBROUTINE INTEGER_CHECK

!-------------------------------------------------------------
SUBROUTINE REAL_CHECK(debug,string,realvar)
    LOGICAL debug
    CHARACTER(LEN=*) :: string
    REAL :: realvar
  IF(debug.eqv..false.) THEN
     RETURN
  ELSE
     print*, string, realvar
  END IF
END SUBROUTINE REAL_CHECK

!--------------------------------------------------------------
SUBROUTINE COMPLEX_CHECK(debug,string,complexvar)
    LOGICAL debug
    CHARACTER(LEN=*) :: string
    DOUBLE COMPLEX :: complexvar
  IF(debug.eqv..false.) THEN
     RETURN
  ELSE
     print*, string, complexvar
  END IF
END SUBROUTINE COMPLEX_CHECK

!-------------------------------------------------------------
SUBROUTINE integer_ARRAY_CHECK(debug,string,arrayvar,dim)
  INTEGER dim
  LOGICAL debug
  CHARACTER(LEN=*) :: string
  INTEGER, DIMENSION(dim) :: arrayvar

   IF(debug.eqv..false.) THEN
     RETURN
  ELSE
    print*, string, arrayvar
  END IF
END SUBROUTINE integer_ARRAY_CHECK

!-------------------------------------------------------------
SUBROUTINE real_ARRAY_CHECK(debug,string,arrayvar,dim)
  INTEGER dim
  LOGICAL debug
  CHARACTER(LEN=*) :: string
  REAL, DIMENSION(dim) :: arrayvar

   IF(debug.eqv..false.) THEN
     RETURN
  ELSE
    print*, string, arrayvar
  END IF
END SUBROUTINE real_ARRAY_CHECK

!-------------------------------------------------------------
SUBROUTINE doublecomplex_ARRAY_CHECK(debug,string,arrayvar,dim)
  INTEGER dim
  LOGICAL debug
  CHARACTER(LEN=*) :: string
  DOUBLE COMPLEX, DIMENSION(dim) :: arrayvar

   IF(debug.eqv..false.) THEN
     RETURN
  ELSE
    print*, string, arrayvar
  END IF
END SUBROUTINE doublecomplex_ARRAY_CHECK

!------------------------------------------------------------------
SUBROUTINE integer_MATRIX_CHECK(debug,string,matrixvar,dim1,dim2)
  INTEGER dim1, dim2
  LOGICAL debug
  CHARACTER(LEN=*) :: string
  INTEGER, DIMENSION(dim1,dim2) :: matrixvar

   IF(debug.eqv..false.) THEN
     RETURN
  ELSE
    print*, string
    DO ii=1,dim1
	print*, matrixvar(ii,:)
    END DO
    print*,
  END IF
END SUBROUTINE integer_MATRIX_CHECK

!-------------------------------------------------------------
SUBROUTINE real_MATRIX_CHECK(debug,string,matrixvar,dim1,dim2)
  INTEGER dim1,dim2
  LOGICAL debug
  CHARACTER(LEN=*) :: string
  REAL, DIMENSION(dim1,dim2) :: matrixvar

   IF(debug.eqv..false.) THEN
     RETURN
  ELSE
    print*, string
    DO ii=1,dim1
	print*, matrixvar(ii,:)
    END DO
    print*, 
END IF
END SUBROUTINE real_MATRIX_CHECK

!------------------------------------------------------------------
SUBROUTINE doublecomplex_MATRIX_CHECK(debug,string,matrixvar,dim1,dim2)
  INTEGER dim1,dim2
  LOGICAL debug
  CHARACTER(LEN=*) :: string
  DOUBLE COMPLEX, DIMENSION(dim1,dim2) :: matrixvar

   IF(debug.eqv..false.) THEN
     RETURN
  ELSE
    print*, string
    DO ii=1,dim1
	print*, matrixvar(ii,:)
    END DO
    print*,
  END IF
END SUBROUTINE doublecomplex_MATRIX_CHECK


END MODULE CHECKPOINTS

!############################################################
! small testing program
! ###########################################################

!PROGRAM TEST
!USE CHECKPOINTS

 !IMPLICIT NONE
  !LOGICAL  debug
 ! INTEGER :: a
 ! REAL, DIMENSION(5) :: vector

  ! a = 2
  ! vector(1) = 1
  ! vector(5) = 5
  !debug=.true.

  !CALL CHECKPOINT(debug,'line test')
  !CALL CHECKPOINT(debug,'variable a value is: ',a)
  !CALL CHECKPOINT(debug,'array vector is: ',vector,5)

  !STOP
!END PROGRAM TEST

