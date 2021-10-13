MODULE USEFUL_FUNCTIONS
CONTAINS

FUNCTION SPACED_EIGENS(eigens,ssize)

INTEGER :: ii,ssize,number_bins
REAL*8  mean
REAL*8, DIMENSION(ssize) :: eigens
REAL*8, DIMENSION(ssize-1) :: spaced_eigens

DO ii=1, ssize-1
	spaced_eigens(ii) = (eigens(ii+1)-eigens(ii))
END DO

mean = sum(spaced_eigens)/(ssize-1)

DO ii=1, ssize-1
	spaced_eigens(ii) = spaced_eigens(ii)/mean
END DO

CLOSE(10)
RETURN
END FUNCTION

END MODULE 
!**********************************************************************
! MAIN
!**********************************************************************

PROGRAM MAIN
  USE MATRICES        
  USE USEFUL_FUNCTIONS
  IMPLICIT NONE
!*******************************************************************
! Declarations 
!*******************************************************************  
  TYPE(DMATRIX) :: mat,diag_mat
  REAL*8, DIMENSION(:), ALLOCATABLE :: RWORK, spaced, diag_spaced, diag_lambdas, lambdas
  INTEGER :: NN, LWORK, INFO, ii, nndim, reason,ss,threshold_index
  DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: WORK
  REAL*8 xx
  REAL*8, DIMENSION(:), ALLOCATABLE :: sps,hist,half_points
  INTEGER :: number_bins,jj
  REAL*8, DIMENSION(:,:),ALLOCATABLE :: to_print,dto_print
  CHARACTER(len=20) filename,file_id
  REAL*8  threshold,step

!*******************************************************************
! Read from file
!*******************************************************************

OPEN(UNIT=77,file='inputdim.dat',status='old')
DO

READ(77,*,iostat=reason) nndim
  IF(reason<0) THEN               !reason<0 means that it reached the end of the file
  EXIT
  ELSE
!************************************************************************  
! Allocation and initialization
!*************************************************************************
  NN = nndim
  LWORK = 2*NN-1
  ALLOCATE(WORK(LWORK), RWORK(3*NN-2))   !for zheev
  ALLOCATE(lambdas(NN),diag_lambdas(NN),spaced(NN),diag_spaced(NN))     

  mat = INIT(NN,NN) 
  diag_mat = DIAG_INIT(NN,NN)

!**************************************************************************
! Calculations
!**************************************************************************
  CALL ZHEEV('N','U',NN,mat%elem,NN,lambdas,WORK,LWORK,RWORK,INFO)
  CALL ZHEEV('N','U',NN,diag_mat%elem,NN,diag_lambdas,WORK,LWORK,RWORK,INFO)

  spaced = SPACED_EIGENS(lambdas,NN)
  diag_spaced = SPACED_EIGENS(diag_lambdas,NN)

!**************************************************************************
! HISTOGRAM
! *************************************************************************
threshold_index = int(size(spaced)*98/100+1)
threshold = spaced(threshold_index)
step = 0.01
sps = pack(spaced,spaced < threshold)
number_bins = int(threshold/step)+1

print*, threshold

ALLOCATE(hist(number_bins))
ALLOCATE(half_points(number_bins))

hist = 0d0

do ii=1,number_bins
  do jj=1,size(sps)
     if(step*(ii-1) <= sps(jj) .and. sps(jj)<step*ii)  then 
     hist(ii) = hist(ii) + 1
	end if
 end do
half_points(ii) = step*(ii-1) + step/2
end do

ALLOCATE(to_print(number_bins,2))

hist = hist/sum(hist)
hist = hist/step

to_print(:,1) = half_points
to_print(:,2) = hist

DEALLOCATE(hist,half_points)

!**************************************************************************
!DIAG HISTOGRAM
!**************************************************************************
threshold_index = int(size(diag_spaced)*96/100+1)
threshold = diag_spaced(threshold_index)
step = 0.001
sps = pack(diag_spaced,diag_spaced < threshold)
number_bins = int(threshold/step)+1

ALLOCATE(hist(number_bins))
ALLOCATE(half_points(number_bins))

hist = 0d0

do ii=1,number_bins
  do jj=1,size(sps)
     if(step*(ii-1) <= sps(jj) .and. sps(jj)<step*ii)  then 
     hist(ii) = hist(ii) + 1
	end if
 end do
half_points(ii) = step*(ii-1) + step/2
end do

ALLOCATE(dto_print(number_bins,2))

hist = hist/sum(hist)
hist = hist/step

dto_print(:,1) = half_points
dto_print(:,2) = hist

!*************************************************************************
! Printings
!************************************************************************
write(file_id, '(i0)') nndim
filename = './FITS/'//trim(adjustl(file_id)) // '.dat'

open(unit=30,file=filename,status="replace") 
            do ii=1,size(to_print,1)
                  write(30,*) (to_print(ii,jj), jj=1,size(to_print,2))
            end do
            close(30)

write(file_id, '(i0)') nndim
filename = './FITS/diag_' // trim(adjustl(file_id)) // '.dat'

open(unit=40,file=filename,status="replace") 
            do ii=1,size(dto_print,1)
                  write(40,*) (dto_print(ii,jj), jj=1,size(dto_print,2))
            end do
            close(40)
       

!*************************************************************************
! Loops closing, deallocations
!*************************************************************************
END IF

DEALLOCATE(WORK, RWORK)
DEALLOCATE(lambdas,diag_lambdas,spaced, diag_spaced)
DEALLOCATE(to_print,hist,half_points,dto_print)

END DO

STOP
END PROGRAM MAIN

