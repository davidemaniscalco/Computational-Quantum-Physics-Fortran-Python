MODULE MATRICES
    IMPLICIT NONE

    TYPE DMATRIX
       INTEGER, DIMENSION(2)::N   !2 dimensional array
       DOUBLE COMPLEX, DIMENSION(:,:),ALLOCATABLE :: Elem !dimension(:,:) because a complex number is a pair of constants
       DOUBLE COMPLEX TRACE,DET
    END TYPE DMATRIX

    TYPE REALMATRIX
        INTEGER, DIMENSION(2)::N   !2 dimensional array
        DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE :: Elem !dimension(:,:) because a complex number is a pair of constants
        DOUBLE PRECISION TRACE,DET
    END TYPE REALMATRIX
  
  CONTAINS
  
  
  !initializing an harmonic oscilation Hamiltonian
    FUNCTION INIT(N,xmax,hbar,w,mass)
      INTEGER ii,jj,N,xmax,xmin
      REAL*8 d,hbar,w,mass
      TYPE(REALMATRIX)::A,INIT
      A%N(1)=N+1  !creating matrix with given number of row and columns
      A%N(2)=N+1
      xmin=-xmax
      d=2.0*xmax/N
      ALLOCATE (A%elem(N+1,N+1))
      A%elem=0
      DO ii=1,N+1
         DO jj=ii,N+1
            IF (ii==jj) THEN
                A%elem(ii,jj)=((hbar**2)/(mass))+(d**2)*((xmin+(ii-1)*d)**2)*w**2
            ELSEIF (jj==ii+1) THEN
                A%elem(ii,jj)=-1.0*(hbar**2)/(2.0*mass)
            END IF
         END DO
      END DO
      INIT=A
      DEALLOCATE(A%elem)
      RETURN
    END FUNCTION INIT


    SUBROUTINE WRITEFILE(OBJECT)
      TYPE(REALMATRIX)::OBJECT
      INTEGER ii
      OPEN(unit=3, file='matrix_out.txt', ACTION="write")
      DO ii=1,OBJECT%N(1)
        WRITE (3,*) OBJECT%elem(ii,:)
      END DO
      CLOSE(3)
    END SUBROUTINE WRITEFILE


    FUNCTION POTENTIAL(ii,jj,xmax,T,tt,d)
      INTEGER ii,jj,xmax,xmin,tmin,tt
      REAL*8 POTENTIAL,d, deltat, T
      
      deltat = T/float(tt)
      tmin=0.0
      xmin=-xmax
      POTENTIAL = ((xmin)+(ii-1)*d-((tmin+(jj-1)*deltat)/T))**2
      
      RETURN
    END FUNCTION

    FUNCTION EXP_POTENTIAL(jj,prevpsi,N,xmax,T,tt,d)
      INTEGER *4 N,xmax,tt
      DOUBLE COMPLEX, DIMENSION(N+1):: psi,prevpsi
      
      DOUBLE PRECISION ::  V,d
      DOUBLE COMPLEX, DIMENSION(N+1)::EXP_POTENTIAL
      INTEGER ii,jj
      REAL*8 deltat, T

      deltat = T/float(tt)

     !                         to print potentials
      !OPEN(unit=22,file='potentials.txt', ACTION='write',status='unknown')
      !WRITE (22,*) -xmax+(ii-1)*d
      !CLOSE(unit=22)
      !END IF

      DO ii=1,N+1
        V=POTENTIAL(ii,jj,xmax,T,tt,d)
        psi(ii)=cdexp(complex(0.0,-V*deltat/float(2)))*prevpsi(ii)
        !WRITE (22,*) -xmax+(ii-1)*d, V !prevpsi(ii), psi(ii)
      END DO
      EXP_POTENTIAL=psi
      !CLOSE(unit=22)
      RETURN 
    END FUNCTION


    FUNCTION EXP_MOMENTA(prevpsi,N,T,tt,d,xmax)
      INTEGER *4 N,tt
      DOUBLE COMPLEX, DIMENSION(N+1):: psi,prevpsi
      DOUBLE PRECISION d
      
      DOUBLE COMPLEX, DIMENSION(N+1)::EXP_MOMENTA,x
      REAL*8 p,deltat, deltap, T, range, width_k, range_k
      REAL*8, DIMENSION(:), ALLOCATABLE :: ks
      INTEGER ii,jj,xmax

      ALLOCATE(ks(N+1))
      deltat = T/float(tt)
      
      !*******************************************************************************
      !                                             vecchio
      !DO ii=1, N/2
      !p=2.0*acos(-1.0)*ii/float(xmax)
      !  psi(ii)=cdexp(complex(0.0,-p*p*deltat*0.5))*prevpsi(ii)
        
      !END DO

      !DO ii=N/2+1, N+1
      !p=2.0*acos(-1.0)*ii/float(xmax) - 2.0*acos(-1.0)/d
      ! psi(ii)=cdexp(complex(0.0,-p*p*deltat*0.5))*prevpsi(ii)
        
      !END DO
      !******************************************************************************
      !             fs

      range = float(xmax)
      width_k = acos(-1.0)/range
      range_k = acos(-1.0)*(N)/(2.0*range)
      ks = [(width_k*ii, ii=0,(N/2 + MOD(N,2))),(-range_k + width_k*ii, ii=0, (N/2-1))]
      
                               
      DO ii=1, N+1
        psi(ii)=cdexp(complex(0.0,-ks(ii)*ks(ii)*deltat*0.5))*prevpsi(ii)
      END DO
    !********************************************************************************  
    EXP_MOMENTA=psi
      RETURN 
    END FUNCTION EXP_MOMENTA
!***********************************************************************************

    FUNCTION FOURIER_TRANSFORM(psi,N)
      use, intrinsic :: iso_c_binding
      include "fftw3.f03"
      INTEGER*4 N
      DOUBLE COMPLEX, DIMENSION(N+1):: psi,FOURIER_TRANSFORM
      INTEGER*8 plan!,FFTW_FORWARD,FFTW_ESTIMATE

      call dfftw_plan_dft_1d(plan,N+1,psi,FOURIER_TRANSFORM,FFTW_FORWARD,FFTW_ESTIMATE)
      CALL DFFTW_EXECUTE_DFT(PLAN,psi,FOURIER_TRANSFORM)
      CALL DFFTW_DESTROY_PLAN(PLAN)
      FOURIER_TRANSFORM=FOURIER_TRANSFORM/sqrt(float(N+1))
      RETURN

    END FUNCTION FOURIER_TRANSFORM


    FUNCTION ANTIFOURIER_TRANSFORM(psi,N)
      use, intrinsic :: iso_c_binding
      include "fftw3.f03"
      INTEGER*4 N
      DOUBLE COMPLEX, DIMENSION(N+1):: psi,ANTIFOURIER_TRANSFORM
      INTEGER*8 plan_backward!,FFTW_FORWARD,FFTW_ESTIMATE

      call dfftw_plan_dft_1d(plan_backward,N+1,psi,ANTIFOURIER_TRANSFORM,FFTW_BACKWARD,FFTW_ESTIMATE)
      CALL DFFTW_EXECUTE_DFT(plan_backward,psi,ANTIFOURIER_TRANSFORM)
      CALL DFFTW_DESTROY_PLAN(plan_backward)
      ANTIFOURIER_TRANSFORM=ANTIFOURIER_TRANSFORM/sqrt(float(N+1))

      RETURN

    END FUNCTION ANTIFOURIER_TRANSFORM

    !initializes an harmonic oscillator matrix
    !diagonalizes it
    
    SUBROUTINE DIAGONALIZE(N,kk,xmax,w,hbar,mass,DEBUG,T,tt)
      TYPE (REALMATRIX)::H
      INTEGER*4 N,INFO,LWORK,xmax,dim,kk,ii,jj,tt
      REAL*8 d,w,hbar,mass,T
      REAL*8,DIMENSION(:),ALLOCATABLE::eigen_hermit,work
      DOUBLE COMPLEX, DIMENSION(N+1)::psi0
      DOUBLE COMPLEX, DIMENSION(:,:),ALLOCATABLE::psi
      LOGICAL DEBUG
      
      ALLOCATE(psi(N+1,tt+1))
      !bin size

      d=(2.0*xmax)/N
      !effective matrix dimension
      dim=N+1
      LWORK=dim*(3+dim/2)
      ALLOCATE(eigen_hermit(dim),WORK(LWORK))
      !initialize matrix
      H=INIT(N,xmax,hbar,w,mass)
      !check matrix is correct 
      !remember only the upper triangular part is assigned
      !lower triangular part (exluding diagonal) is zero
      if (DEBUG .eqv. .true.) THEN
        CALL WRITEFILE(H)
      END IF
      !diagonalize matrix
      CALL dsyev( 'V', 'U', dim, H%Elem,dim, eigen_hermit, WORK, LWORK,  INFO)

      !multipliyng by required constants
      eigen_hermit=hbar*w*eigen_hermit/d**2
      !OPEN(unit=2, file='eigenvalues.txt', ACTION="write")
      OPEN(unit=11,file='eigenstates.txt', ACTION='write',status='replace')

      !ground state
      psi0=H%Elem(:,1)
      psi(:,1)= psi0/sqrt(d)

      !j is the temporal index
      !i is the spatial index
      DO jj = 2, tt+1
          psi(:,jj) = EXP_POTENTIAL(jj,psi(:,jj-1),N,xmax,T,tt,d)
          psi(:,jj) = FOURIER_TRANSFORM(psi(:,jj),N)
          psi(:,jj) = EXP_MOMENTA(psi(:,jj),N,T,tt,d,xmax)
          psi(:,jj) = ANTIFOURIER_TRANSFORM(psi(:,jj),N)
          psi(:,jj) = EXP_POTENTIAL(jj,psi(:,jj),N,xmax,T,tt,d)    
          !psi(:,jj) = psi(:,jj)/sqrt(d)             
      END DO

      !write eigenstates
      DO ii=1, H%N(1)
        WRITE (11,*) -xmax+(ii-1)*d,(abs(psi(ii, jj)), jj = 1, tt+1)
      END DO

      !CLOSE(2)
      CLOSE(11)

      DEALLOCATE(eigen_hermit,WORK,psi)
    END SUBROUTINE

  END MODULE MATRICES
  
PROGRAM HARMONIC_OSCILLATOR
    use, intrinsic :: iso_c_binding
    USE MATRICES
    IMPLICIT NONE
    include "fftw3.f03"

    INTEGER*4 N,xmax,i,kk, tt
    TYPE(REALMATRIX)::M
    REAL*8 hbar,w,mass, T
    LOGICAL DEBUG

    print*, 'Write the value for xmax (natural): '
    READ(*,*)  xmax
    print*, 'Write the value for the number N of spacings (natural): '
    READ(*,*)  N
    print*, 'Write the value for the total time T: '
    READ(*,*)  T
    print*, 'Write the value for the number of temporal spacings (natural): '
    READ(*,*)  tt

   ! print*, 'Running for deltax = ',2.0*xmax/N,', deltat = ',T/tt
    kk = 1
    !reduced planck constant, angular frequency and mass
    hbar=1.0
    w=1.0
    mass=0.5
    DEBUG=.false.
    CALL DIAGONALIZE(N,kk,xmax,w,hbar,mass,DEBUG,T,tt)
    STOP
END PROGRAM HARMONIC_OSCILLATOR
