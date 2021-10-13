!########################################################################
! POINT 2 EXERCISE 3 INFORMATION THEORY AND COMPUTATION course,
! master degree in physics of data, university of Padova, October 2019
! Davide Maniscalco
!--------------------------------------------------------------------
! MATRIX MULTIPLICATION with CHECKPOINTS, COMMENTS, PRE/POST CONDS., DOCUMENTATION
! ------------------------------------------------------------------
! Matrix-matrix multiplication in 3 diff ways with CPU time monitoring
!--------------------------------------------------------------------
!
! ##################################################################
! IMPORTANT NOTE ABOUT COMPILATION:
! This code uses, for the checkpoints, an external module, named 'checkpoints'.
! Therefore the external module file must be compiled in the same folder of this code
! with the -c flag. A .mod file and a .o file must deliver after that compilation.
! Moreover, this code must be compiled writing the name of the .o file at the end, e.g:
! gfortran Ex3-Maniscalco-CODE.f90 Ex3-Maniscalco-CHECKPOINTMODULE.o
!---------------------------------------------------------------
! VARIABLES DECLARATION:
!
! real allocatable arrays for the matrices
! integers for the do loops
! reals for the cpu_times
! integer arrays for matrices dimensions
! debugging logical variable
! -------------------------------------------------------------
! DIMENSIONS INITIALIZATION:

! Dimensions must me handled by hand in the code: nn,mm,ll, meant for the matrix-
!matrix product (nn x mm)(mm x ll) = (nn x ll)
! shape function returns an integer array containing the dimension of an input matrix
! -------------------------------------------------------------
! MATRICES INITIALIZATION:

! made in the code arbitarily
! first_matrix(ii,jj) = ii (its row number)
! second_matrix(ii,jj) = jj (its column number)
! output matrices are initialized to 0
!----------------------------------------------------------------
! PRODUCTS:

! prodmat1: output matrix product from by-column loop
! prodmat2: output matrix product from by-row loop
! prodmat3: output matrix product from matmul
!----------------------------------------------------------------
! CHECKPOINTS:
! Checkpoints are meant to print matrices and matrices dimension. Set debug.eqv..false.
! to deactivate checks. 
! --------------------------------------------------------------
! CALL_CPU_TIME:
!
! Will be used to store istant ina real variable. Istants differences will be
! made in the end
! ----------------------------------------------------------------
! TIMES PRINTING:
!
! The print of the final matrix dimensions is based on the matmul one
!##################################################################
!
! ---------------------------------------------------------------
! Variables declaration
! -----------------------------------------------------------------
PROGRAM MATRIX_MULTIPLICATION
      USE CHECKPOINTS
      IMPLICIT NONE
      
      REAL,DIMENSION(:,:),ALLOCATABLE ::mat1,mat2
      REAL,DIMENSION(:,:),ALLOCATABLE ::prodmat1,prodmat2,prodmat3
      INTEGER ii,jj,kk,nn,mm,ll
      REAL T1,T2,T3,T4,T5,T6
      
      INTEGER, DIMENSION(2) :: dim1,dim2,dimp1,dimp2,dimp3
      LOGICAL debug
!      CHARACTER ( LEN =*) , PARAMETER :: FMT2 = " (A , I0 , ’X ’ , I0 ) "
! --------------------------------------------------------------------------
      !set to false to ignore all debugs
! --------------------------------------------------------------------------
      debug = .true.
! ---------------------------------------------------------------------------
      ! Dimension initialization      
      ! Product between matrices of dimensions: (nn x mm)(mm x ll)=(nn x ll)
! ----------------------------------------------------------------------------      
      nn = 5
      mm = 3
      ll = 2

      ALLOCATE(mat1(nn,mm))
      ALLOCATE(mat2(mm,ll))
      ALLOCATE(prodmat1(nn,ll))
      ALLOCATE(prodmat2(nn,ll))
      ALLOCATE(prodmat3(nn,ll))

      dim1 = SHAPE(mat1)
      dim2 = SHAPE(mat2)
      dimp1 = SHAPE(prodmat1)
      dimp2 = SHAPE(prodmat2)
      dimp3 = SHAPE(prodmat3)
!--------------------------------------------------------------
      !checkpoint, matrices dimensions
!--------------------------------------------------------------
 
       IF(debug.eqv..true.) THEN
         print 1, '1st matrix dimensions: ',dim1
         print 1, '2nd matrix dimensions: ',dim2
         print 1, 'by-column matrix dimensions: ',dimp1
         print 1, 'by-row matrix dimensions: ',dimp2
         print 1, 'matmul matrix dimensions: ',dimp3
      END IF


     ! matrices dimensions precondition
     IF (dim1(2) /= dim2(1)) THEN
      print*, 'Error: first input matrix has ',dim1(2), 'columns, while second input matrix has ',dim2(1), 'rows'
      STOP
   END IF 
!----------------------------------------------------------
   ! matrices initialization
! ---------------------------------------------------------

   ! first matrix
       DO ii=1,nn
         DO jj=1,mm
            mat1(ii,jj)=ii
         ENDDO
      ENDDO

      !second matrix
      DO ii=1,mm
         DO jj=1,ll
            mat2(ii,jj)=jj
         ENDDO
      ENDDO
      
   ! initialization to 0 of prodmatrices
   DO ii=1,nn
         DO jj=1,ll
            prodmat1(ii,jj)=0
            prodmat2(ii,jj)=0
            prodmat3(ii,jj)=0
         ENDDO
      ENDDO   
!--------------------------------------------------------
      ! Input matrices checkpoint
!--------------------------------------------------------
 CALL real_MATRIX_CHECK(debug,'1st matrix is: ',mat1,dim1(1),dim1(2))
 CALL real_MATRIX_CHECK(debug,'2nd matrix is: ',mat2,dim2(1),dim2(2))
!--------------------------------------------------------
 ! By column loop, with post condition
!--------------------------------------------------------
      CALL CPU_TIME(T3)
      DO kk=1,mm
         DO jj=1,ll
            DO ii=1,nn
              prodmat1(ii,jj)=prodmat1(ii,jj)+mat1(ii,kk)*mat2(kk,jj)
            ENDDO
         ENDDO
      ENDDO
      CALL CPU_TIME(T4)
      
!post condition      
IF ((dim1(1) /= dimp1(1)) .or. (dim2(2) /= dimp1(2))) THEN
 print*, 'Error: by-column product matrix has not compatible dimensions with input ones'
 STOP
 END IF
!---------------------------------------------------------
 ! By row loop, with post condition
!---------------------------------------------------------
      CALL CPU_TIME(T1)
      DO ii=1,nn
         DO jj=1,ll
            DO kk=1,mm
               prodmat2(ii,jj)=prodmat2(ii,jj)+mat1(ii,kk)*mat2(kk,jj)
            ENDDO
         ENDDO
      ENDDO
      CALL CPU_TIME(T2)
      
!post condition      
IF ((dim1(1) /= dimp2(1)) .or. (dim2(2) /= dimp2(2))) THEN
 print*, 'Error: by-row product matrix has not compatible dimensions with input ones'
 STOP
END IF
!---------------------------------------------------------
 ! Fortran function, with post condition
!---------------------------------------------------------
      CALL CPU_TIME(T5)
      prodmat3 = MATMUL(mat1,mat2)
      CALL CPU_TIME(T6)

!post condition
IF ((dim1(1) /= dimp3(1)) .or. (dim2(2) /= dimp3(2))) THEN
 print*, 'Error: matmul product matrix has not compatible dimensions with input ones'
 STOP
 END IF
!----------------------------------------------------------
 ! Output matrices checkpoint
 ! --------------------------------------------------------
 CALL real_MATRIX_CHECK(debug,'by-column matrix is: ',prodmat1,dimp1(1),dimp1(2))
CALL real_MATRIX_CHECK(debug,'by-row matrix is: ',prodmat2,dimp2(1),dimp2(2))
CALL real_MATRIX_CHECK(debug,'matmul matrix is: ',prodmat3,dimp3(1),dimp3(2))
!----------------------------------------------------------
! Times printing
! ---------------------------------------------------------
      print 1, 'Final matrix dimensions: ',dimp3(1),dimp3(2) 
      print*, 'By row time: ',T2-T1
      print*, 'By col time: ',T4-T3
      print*, 'Matmul time: ',T6-T5
!-----------------------------------------------------------
! Deallocations
!-----------------------------------------------------------
DEALLOCATE(mat1)
DEALLOCATE(mat2)
DEALLOCATE(prodmat1)
DEALLOCATE(prodmat2)
DEALLOCATE(prodmat3)
!-----------------------------------------------------------
! Printing formats
!-----------------------------------------------------------   
1	FORMAT(A,I0,' ',I0)
!-----------------------------------------------------------
  
      STOP
    END PROGRAM MATRIX_MULTIPLICATION
