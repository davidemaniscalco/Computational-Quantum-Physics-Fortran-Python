MODULE EDWIN
    IMPLICIT NONE
    CONTAINS
FUNCTION SHROD_INITIALIZE(dim,dd,xmin,xmax)
    REAL*8, DIMENSION(:,:), ALLOCATABLE :: SHROD_INITIALIZE
    INTEGER :: ii, jj, dim, xmin, xmax
    REAL*8 :: dd,omega,m,hbar

ALLOCATE(SHROD_INITIALIZE(dim+1,dim+1))

SHROD_INITIALIZE = 0
omega = 1d0
hbar = 1d0
 m = 0.5


    DO ii= 1, dim+1
        DO jj=ii, ii+1
            IF(jj==(dim+2)) then
              exit
            END IF

            IF(jj == ii) THEN
           SHROD_INITIALIZE(ii,jj) = (2.0*hbar/2/m + omega**2*dd**2*(xmin+(ii-1)*dd)**2)
           ELSE
               SHROD_INITIALIZE(ii,jj) = -1.0*hbar/2/m
            END IF
     END DO
    END DO
    
END FUNCTION

SUBROUTINE DIAGONALIZE(NN,ham,dd,xmin)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK, WORK,lambdas,norm_lambdas
    REAL*8, DIMENSION(NN+1,NN+1) :: ham
    INTEGER :: LWORK, INFO, NN,ii,jj,xmin
    REAL*8 dd
    LOGICAL debug

  debug = .true.
  LWORK = (NN+1)*(3+(NN+1)/2)
  ALLOCATE(lambdas(NN+1),norm_lambdas(NN+1))
  ALLOCATE(WORK(LWORK), RWORK(3*(NN+1)-2)) 

  CALL DSYEV('V','U',NN+1,ham,NN+1,lambdas,WORK,LWORK,INFO)
  norm_lambdas = lambdas/dd**2

  OPEN(unit=10,file='eigenvalues.dat',status="replace")
    DO ii=1,NN+1
        write(10,*) ii, ' ',norm_lambdas(ii)
    END DO
  CLOSE(10)
  
  OPEN(unit=20,file='eigenvectors.dat',status="replace")
  DO ii=1,NN+1
      write(20,*) xmin+(ii-1)*dd, ham(ii,:)/sqrt(dd)
  END DO
CLOSE(20)


  IF (debug .eqv..true.) THEN
  print*, 'First and last three eigenvalues: '
  print*, norm_lambdas(1),norm_lambdas(2),norm_lambdas(3)
  print*, norm_lambdas(NN-1),norm_lambdas(NN),norm_lambdas(NN+1)
  END IF

END SUBROUTINE

END MODULE

PROGRAM MAIN
    USE EDWIN
    IMPLICIT NONE
    REAL*8 ::  dd
    INTEGER NN, xmin, xmax,reason
    REAL*8, DIMENSION(:,:),ALLOCATABLE :: ham
   ! LOGICAL debug

   ! debug = .true.

OPEN(UNIT=77,file='input.dat',status='old')

READ(77,*,iostat=reason) xmax, NN
 ! IF(reason<0) THEN               !reason<0 means that it reached the end of the file
  !EXIT
  !ELSE

    print*, 'Input values N = ',NN,', xmax = ',xmax
    xmin = -xmax
    dd = 2.0*xmax/NN

    ham = SHROD_INITIALIZE(NN,dd, xmin,xmax)
    CALL DIAGONALIZE(NN,ham,dd,xmin)

  !END IF
STOP
END
