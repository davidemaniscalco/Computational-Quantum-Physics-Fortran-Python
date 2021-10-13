MODULE MATRICES
IMPLICIT NONE
  
!********************************** TYPES ****************

  TYPE RMATRIX
  INTEGER, DIMENSION(2)::N
  REAL :: det
  REAL ::  tr
  REAL, DIMENSION(:,:), ALLOCATABLE:: elem !this is the whole matrix
  END TYPE RMATRIX
  

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
  
! ******************************INTERFACES***************************
  
  INTERFACE OPERATOR (.ADJ.)
   MODULE PROCEDURE MATADJOINT,CONJUGATE,VECADJOINT
END INTERFACE OPERATOR (.ADJ.)

INTERFACE OPERATOR (.TR.)
   MODULE PROCEDURE TRACE
END INTERFACE OPERATOR (.TR.)

!******************************
  
CONTAINS

! ******************************  FUNCTION REAL_INITV, INITIALIZES THE VECTOR with real randoms ****************
  FUNCTION REAL_INITV(length)
    INTEGER :: ii, length
    TYPE(DVECTOR) :: REAL_INITV
    REAL*8 xx

    REAL_INITV%N = length
    ALLOCATE(REAL_INITV%elem(length))

    DO ii=1,length
       CALL RANDOM_NUMBER(xx)
       REAL_INITV%elem(ii) = xx
    END DO
    
    RETURN
  END FUNCTION REAL_INITV


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

!******************************* FUNCTION DIAG_INIT *******************************
 FUNCTION DIAG_INIT(nrow,ncol)
    REAL*8 :: xx,yy
    INTEGER, INTENT(IN) :: nrow,ncol   !stuff MUST be declared AS FIRST. Otherwise, ERROR
    
    INTEGER :: ii,jj,zero
    TYPE(DMATRIX) :: DIAG_INIT
    DIAG_INIT%N(1)=nrow
    DIAG_INIT%N(2)=ncol
    ALLOCATE (DIAG_INIT%elem(nrow,ncol))
    yy = 1d0
     DO ii=1,nrow
         DO jj=ii,ncol
             CALL RANDOM_NUMBER(xx)
             	IF (ii==jj) THEN
            	DIAG_INIT%elem(ii,jj) = complex(xx,yy)
                ELSE
		DIAG_INIT%elem(ii,jj) = complex(yy,yy)
		END IF
         ENDDO
      ENDDO
      DIAG_INIT%tr = .TR.(DIAG_INIT)    ! this calls the trace function to calculate it
      zero=0
      DIAG_INIT%det = complex(0.0/zero,0.0/zero) ! To be implemented. Meanwhile returns NAN
    RETURN
  END FUNCTION DIAG_INIT

!*****************************   FUNCTION INIT, INITIALIZES THE DMATRIX ***************
  
  FUNCTION INIT(nrow,ncol)
    REAL*8 :: xx,yy
    INTEGER, INTENT(IN) :: nrow,ncol   !stuff MUST be declared AS FIRST. Otherwise, ERROR
    
    INTEGER :: ii,jj,zero
    TYPE(DMATRIX) :: INIT
    INIT%N(1)=nrow
    INIT%N(2)=ncol
    ALLOCATE (INIT%elem(nrow,ncol))
     DO ii=1,nrow
         DO jj=ii,ncol
             CALL RANDOM_NUMBER(xx)
             CALL RANDOM_NUMBER(yy)
             	IF (ii==jj) THEN
             	yy = 0
            	INIT%elem(ii,jj) = complex(xx,yy)
                ELSE
		INIT%elem(ii,jj) = complex(xx,yy)
                INIT%elem(jj,ii) = conjg(INIT%elem(ii,jj))
		END IF
         ENDDO
      ENDDO
      INIT%tr = .TR.(INIT)    ! this calls the trace function to calculate it
      zero=0
      INIT%det = complex(0.0/zero,0.0/zero) ! To be implemented. Meanwhile returns NAN
    RETURN
  END FUNCTION INIT

!***********************************  FUNCTION TRACE **************************
  
FUNCTION TRACE(A)
  TYPE(DMATRIX), INTENT(IN) :: A
  INTEGER :: ii,x
  DOUBLE COMPLEX :: TRACE
  TRACE=(0d0,0d0)
  
  IF (A%N(1) /= A%N(2)) THEN
     x=0
     TRACE = complex(0.0/x,0.0/x)
     print*, 'Warning: not squared matrix given in input for TRACE calculus'
     
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
  TYPE(DMATRIX):: OBJECT
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
END MODULE MATRICES
