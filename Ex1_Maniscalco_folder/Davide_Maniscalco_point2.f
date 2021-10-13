      SUBROUTINE POINT2A
      INTEGER*4 x,y
      INTEGER*2 q,w
      x=1
      y=2000000
      q=1
      w=2000000
      print*,'Sum with integer*4: ' , x+y
      print*,'Sum with integer*2: ', q+w
      END SUBROUTINE POINT2A

      SUBROUTINE POINT2B
      REAL*8 x,y
      REAL*4 q,w
      x=4*ATAN(1.0)*1E32
      y=sqrt(2.0)*1E21
      q=4*ATAN(1.0)*1E32
      w=sqrt(2.0)*1E21
      print*,'Sum with double precision: ',x+y
      print*,'Sum with single precision: ',q+w
      print*,'(pi)E32 in single precision: ',q
      END SUBROUTINE POINT2B

      PROGRAM SHEET1PART1
      IMPLICIT NONE
      call POINT2A()
      call POINT2B()
      STOP
      END
