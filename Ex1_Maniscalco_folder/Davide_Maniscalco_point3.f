      PROGRAM SHEET1PART2
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(500,500) ::mat1
      DOUBLE PRECISION,DIMENSION (500,500) :: mat2
      DOUBLE PRECISION,DIMENSION(500,500) ::prodmat1,prodmat2,prodmat3

      INTEGER ii,jj,kk,nn,mm,ll
      REAL T1,T2,T3,T4,T5,T6
      
! Product between matrices of dimensions: (nn x mm)(mm x ll)=(nn x ll)
      nn = 500
      mm = 500
      ll = 500
      
! definition of the first matrix
       DO ii=1,nn
         DO jj=1,mm
            mat1(ii,jj)=ii
         ENDDO
      ENDDO

! definition of the second matrix
      DO ii=1,mm
         DO jj=1,ll
            mat2(ii,jj)=jj
         ENDDO
      ENDDO
      
! By column loop
      CALL CPU_TIME(T3)
      DO kk=1,mm
         DO jj=1,ll
            DO ii=1,nn
              prodmat1(ii,jj)=prodmat1(ii,jj)+mat1(ii,kk)*mat2(kk,jj)
            ENDDO
         ENDDO
      ENDDO
      CALL CPU_TIME(T4)

! By row loop
      CALL CPU_TIME(T1)
      DO ii=1,nn
         DO jj=1,ll
            DO kk=1,mm
               prodmat2(ii,jj)=prodmat2(ii,jj)+mat1(ii,kk)*mat2(kk,jj)
            ENDDO
         ENDDO
      ENDDO
      CALL CPU_TIME(T2)

! Fortran function
      CALL CPU_TIME(T5)
      prodmat3 = MATMUL(mat1,mat2)
      CALL CPU_TIME(T6)

! print of output matrices, uncomment to check correctness of the result
!      print*, prodmat1
!      print*, prodmat2
!      print*, prodmat3

      print*, 'Input Matrices dimensions: ',nn,'',mm,';',mm,'',ll
      print*, 'By row time: ',T2-T1
      print*, 'By col time: ',T4-T3
      print*, 'Matmul time: ',T6-T5
      
      STOP
      END
