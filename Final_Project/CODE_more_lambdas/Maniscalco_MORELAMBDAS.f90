!############################################################################
! FINAL EXERCISE INFORMATION THEORY AND COMPUTATION COURSE
! Master degree in Physics of Data, University of Padova, January 2020
! Davide Maniscalco
!--------------------------------------------------------------------------
! ADIABATIC QUANTUM COMPUTING for RANDOM FIELD ISING MODEL
! -------------------------------------------------------------------------
!  SUGGESTED COMPILATION (NEEDED MATRIX MODULE)
! gfortran MANISCALCO-RMATRIX.f90 Maniscalco_MORELAMBDAS.f90 -O -llapack -o Exfinale.out
!---------------------------------------------------------------------------
! DIFFERENCES RESPECT TO MAIN CODE
! Here the values of lambda are taken as input form file 'input_list.dat'
!
! ############################################################################
! LINEAR SYSTEM SOLVING
! Input = A,X,B such that  Ax = b linear system. ZGESV solves the system
!----------------------------------------------------------------------------
! HAM_0
! Input: dimension N, Output: RMATRIX(ham_0). Ising 1d chain, field along z axis
!---------------------------------------------------------------------------
!MF_HAM_P
! Input: dimension N, strength of interaction. Output: Ising 1d chain hamilt , nn interaction
!---------------------------------------------------------------------------
! HAM_P
! Input: dimension N. Output: Ising 1d chain hamilt with NN interaction (Type dmatrix)
!---------------------------------------------------------------------------
! TOTAL_MF(ham0, hamp, tt, bigt, total) 
! Input: start hamiltonian, end hamiltonian (RMATRIX), instant, final instant, dmatrix where store result
! Output: total = (1-tt/bigt)*ham0 + tt/bigt*hamp
!---------------------------------------------------------------------------------
! SPECTRUM
! Input: hamiltonian (RMATRIX) !Output: file containing spectrum
!--------------------------------------------------------------------------------
! GROUNDSTATE
! Input: hamiltonian (RMATRIX) !Output: real vector containing the ground state
! --------------------------------------------------------------------------------
! TIME EVOLUTION
! Obsolete: use CN_TIME_EVOLUTION
!---------------------------------------------------------------------------------
! CN_TIME_EVOLUTION
! Performs the time evolution of the given state with Crank Nicholson, for given input
! start and end hamiltonian. Total ham done with TOTAL_MF
!
! Input: start and end hamilt (RMATRIX); start state (DOUBLE COMPLEX), total time, #of steps,
  ! vector (DOUBLE COMPLEX) where to save evolution
! Output: file(30): evolved eigenvector at each step; file(40): respective eigenvalue at each 
  ! time step.
! --------------------------------------------------------------------------------
! TRUE_RESULT_DURING_TIME(ham0, hamp, bigt, nsteps, NN)
! Input: start and end hamilt (RMATRIX); total time, #of steps, dimension of the system
! At each step matrix is diagonalized to compute eigenvalues and the epsilon factor
! Output: file(80): complete spectrum at each iteration; file(30): groundstate vector at each iter,
         ! file(10): in each row, parity expectation value for each eigenstate
         ! at a given time.
! Debugging output: file(70): eigenvectors print, to be handled in the code; file(20): 
                   !groundstate at each step; file(60) diagonalized ham at each steo
!-----------------------------------------------------------------------------------
! PARITY HAMILTONIAN
! Input: Dimension of the system
! Output: parity matrix P=\prod\sigma_z^i with i=1,N
!-----------------------------------------------------------------------------------
! GS_PARITY
! Input: dimension of the system
! Output: the last element of the parity matrix
!-----------------------------------------------------------------------------------
! PARITY_CHECK
! Input: start and end hamiltonian, dimension of the system
! Output: parity matrix, start hamilt, end hamilt; commutators
!##############################################################################

MODULE AQC
    USE MATRICES
    IMPLICIT NONE
    CONTAINS

    !********************************************************************
    ! LINEAR SYSTEM SOLVING
    !********************************************************************
    SUBROUTINE SYSTEM_SOLVE(AA, xx, bb)
        TYPE(DMATRIX) :: AA
        DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: xx, bb
        INTEGER :: dim, INFO
        INTEGER, DIMENSION(:), ALLOCATABLE :: IPIV

        dim = AA%N(1)
        ALLOCATE(IPIV(dim))

        CALL ZGESV(dim, 1, AA%elem, dim, IPIV, bb, dim, INFO )

        xx = bb  !zgesv after excut puts the sol, i.e. what I want is xx, in bb

        DEALLOCATE(IPIV, bb)

    END SUBROUTINE SYSTEM_SOLVE

    !*********************************************************************
    ! HAM 0 
    !*********************************************************************
    FUNCTION HAM_0(NN)
        TYPE(RMATRIX) :: sigma_z, prev_id, next_id, temp, temp2, ham_0
        INTEGER :: NN, pp

        sigma_z = pauli_zz()

        ham_0%N(1) = 2**NN
        ham_0%N(2) = 2**NN  
        ALLOCATE(ham_0%elem(2**NN,2**NN))
        ham_0%elem = 0.0

        DO pp=1, NN
            prev_id = IDENTITY(2**(pp-1))
            next_id = IDENTITY(2**(NN-pp))
      
            temp = TENSOR_PRODUCT(prev_id, sigma_z)
            temp2 = TENSOR_PRODUCT(temp, next_id)
      
            ham_0%elem = ham_0%elem + temp2%elem
                 DEALLOCATE(prev_id%elem, next_id%elem, temp%elem, temp2%elem)
          END DO
          RETURN
          DEALLOCATE(sigma_z%elem)
    END FUNCTION HAM_0

!**********************************************************************************
    ! MF HAM P
!*********************************************************************************
    FUNCTION MF_HAM_P(NN)
            TYPE(RMATRIX) :: sigma_x, prev_id, next_id, temp, temp2, temp3, mf_ham_p
            INTEGER :: NN, pp
            REAL*4 :: lambda

            print*, 'Insert lambda for mean field: '
            READ(*,*) lambda
    
            sigma_x = pauli_xx()
    
            mf_ham_p%N(1) = 2**NN
            mf_ham_p%N(2) = 2**NN  
            ALLOCATE(mf_ham_p%elem(2**NN,2**NN))
            mf_ham_p%elem = 0.0

            DO pp=1, NN-1
                prev_id = IDENTITY(2**(pp-1))
                next_id = IDENTITY(2**(NN-pp-1))
          
                temp = TENSOR_PRODUCT(prev_id,sigma_x)
                temp2 = TENSOR_PRODUCT(temp,sigma_x)
                temp3 = TENSOR_PRODUCT(temp2,next_id)
          
                mf_ham_p%elem = mf_ham_p%elem + temp3%elem
          
                DEALLOCATE(prev_id%elem, next_id%elem, temp%elem, temp2%elem, temp3%elem)
            END DO

              mf_ham_p%elem = lambda*mf_ham_p%elem
    
              RETURN
              DEALLOCATE(sigma_x%elem)
    END FUNCTION MF_HAM_P

!**********************************************************************************
    ! HAM P WITHOUT MEAN FIELD
!*********************************************************************************
 FUNCTION HAM_P(NN)
        TYPE(RMATRIX) :: sigma_x, prev_id, next_id, temp, temp2, temp3, ham_p
        INTEGER :: NN, pp
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lambdas

        ALLOCATE(lambdas(NN-1))             !  one lambda for each interaction
	OPEN(unit=77, file='input_list.dat')
        READ(77,*) lambdas
        print*, 'lambdas debugging: ',lambdas
	CLOSE(77)

        ham_p%N(1) = 2**NN
        ham_p%N(2) = 2**NN  
        ALLOCATE(ham_p%elem(2**NN,2**NN))
        ham_p%elem = 0.0

        DO pp=1, NN-1
            sigma_x = pauli_xx()
            sigma_x%elem = sigma_x%elem 
            prev_id = IDENTITY(2**(pp-1))
            next_id = IDENTITY(2**(NN-pp-1))
      
            temp = TENSOR_PRODUCT(prev_id,sigma_x)
            temp2 = TENSOR_PRODUCT(temp,sigma_x)
            temp3 = TENSOR_PRODUCT(temp2,next_id)
      
            ham_p%elem = ham_p%elem + lambdas(pp)*temp3%elem
      
            DEALLOCATE(prev_id%elem, next_id%elem, temp%elem, temp2%elem, temp3%elem)
        END DO

          RETURN
          DEALLOCATE(sigma_x%elem)
END FUNCTION HAM_P

!*******************************************************************************
    ! TOTAL HAMILTONIAN MF
!******************************************************************************
    SUBROUTINE TOTAL_MF(ham0, hamp, tt, bigt, total) 
        REAL*8 :: tt, bigt
        TYPE(RMATRIX) :: ham0, hamp, total
        total%elem = (1-tt/bigt)*ham0%elem + tt/bigt*hamp%elem
        RETURN
    END SUBROUTINE

!*******************************************************************************
    ! SPECTRUM
!*******************************************************************************
    SUBROUTINE SPECTRUM(hamiltonian)
        TYPE(RMATRIX) :: hamiltonian, copy
        DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: eigenvalues

        copy%N(1) = hamiltonian%N(1)
        copy%N(2) = hamiltonian%N(2)
        copy%elem = hamiltonian%elem

        CALL DIAGONALIZE(copy, eigenvalues)

        OPEN(UNIT=20,FILE='spectrum.dat',status='replace')
        WRITE(20,*) eigenvalues(1:)
        CLOSE(20)

        DEALLOCATE(eigenvalues, copy%elem)
    END SUBROUTINE SPECTRUM

!***************************************************************************************
! GROUND STATE
!***************************************************************************************
FUNCTION GROUNDSTATE(hamiltonian)
    TYPE(RMATRIX) :: hamiltonian, copy
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: groundstate

    copy%N(1) = hamiltonian%N(1)
    copy%N(2) = hamiltonian%N(2)
    copy%elem = hamiltonian%elem

    CALL DIAGONALIZE_ONLYVEC(copy)

    groundstate = copy%elem(:,1)

    DEALLOCATE(copy%elem)
    RETURN
END FUNCTION GROUNDSTATE

!***************************************************************************************
! TIME EVOLUTION
!**************************************************************************************
SUBROUTINE TIME_EVOLUTION(ham0, hamp, gs_start, bigt, nsteps, gs)
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: gs, application
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gs_start
    REAL*8 :: deltat, bigt, tt, norm
    INTEGER :: nsteps, step, dim, ii
    TYPE(RMATRIX) :: ham0, hamp, hamiltonian
    LOGICAL :: debug

    debug = .true.

    deltat = bigt/float(nsteps)
    dim = ham0%N(1)

    DO ii=1, dim
        gs(ii) = complex(gs_start(ii),0.0)
    END DO

    print*, 'deltat= ',deltat

    hamiltonian%N(1) = ham0%N(1)
    hamiltonian%N(2) = ham0%N(2)
    ALLOCATE(hamiltonian%elem(hamiltonian%N(1),hamiltonian%N(2)))

    ALLOCATE(application(dim))

    IF(debug.eqv..true.) THEN
        OPEN(unit=20, file='timeevoldebug.dat',status='replace')
        WRITE(20,*) 'gs_start: ', gs_start, char(10)
        WRITE(20,*) 'gs: ', gs, char(10)
        WRITE(20,*) dim, hamp%N(1), hamiltonian%N(1), hamiltonian%N(2)
        write(20,*) bigt, nsteps, deltat
        CLOSE(20)
    END IF

    OPEN(unit=30, file='evol.dat',status='replace')
    OPEN(unit=40, file='eigenvalues_evol.dat',status='replace')
    WRITE(30,*) gs, char(10)

    DO step=1, nsteps
        tt = deltat*float(step-1)           !need the ham at previous t
        CALL TOTAL_MF(ham0, hamp, tt, bigt, hamiltonian)
        gs = gs - complex(0.0,deltat)*matmul(hamiltonian%elem,gs)
        norm = dot_product(gs,gs)   !automathically conjugates 1st argument
        gs = gs/norm
        application = matmul(hamiltonian%elem,gs)  !gs should be eigv, then i should find gs=appl
        
        WRITE(30,*) gs, char(10)
        WRITE(40,*) tt/bigt, real(dot_product(gs, application))
    END DO

    CLOSE(30)
    CLOSE(40)

END SUBROUTINE TIME_EVOLUTION
!**********************************************************************************************
! TRUE RESULT DURING TIME
!*********************************************************************************************
SUBROUTINE TRUE_DURING_TIME(ham0, hamp, bigt, nsteps, NN)
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gs, temp_gs, epsilon_gs, eigenvalues, parity_vec
    REAL*8 :: deltat, bigt, tt, gsp
    INTEGER :: nsteps, step, dim, ii, num_epsilon, NN
    TYPE(RMATRIX) :: ham0, hamp, hamiltonian, ham_parity
    LOGICAL :: debug, verbose

    debug = .false.
    verbose = .false.

    deltat = bigt/float(nsteps)
    dim = ham0%N(1)
    num_epsilon = dim             ! how many gaps to calculate

    ! parity
    gsp = GS_PARITY(NN)
    ham_parity = parity_hamiltonian(NN)

    !print*, 'deltat= ',deltat
    !print*, 'ground state parity: ', gsp

    hamiltonian%N(1) = ham0%N(1)
    hamiltonian%N(2) = ham0%N(2)
    ALLOCATE(hamiltonian%elem(hamiltonian%N(1),hamiltonian%N(2)))

    ALLOCATE(gs(ham0%N(1)), temp_gs(ham0%N(1)))
    ALLOCATE(epsilon_gs(num_epsilon))
    ALLOCATE(parity_vec(num_epsilon))

   ! debugging files
    IF(debug.eqv..true.) THEN
        OPEN(unit=20, file='true_evol_gs.dat',status='replace')
        OPEN(unit=60, file='true_vectors_evol.dat',status='replace')
        OPEN(unit=70, file='eigenvectors_debugging.dat', status='replace') 
    END IF

    OPEN(unit=80, file='true_spectrum.dat',status='replace')
    OPEN(unit=30, file='epsilon_gs.dat',status='replace')

    DO step=0, nsteps
        IF(verbose.eqv..true.) THEN
            print*, step
        END IF
        tt = deltat*float(step)
        CALL TOTAL_MF(ham0, hamp, tt, bigt, hamiltonian)
        CALL DIAGONALIZE(hamiltonian, eigenvalues)

        gs =  hamiltonian%elem(:,1)
        temp_gs = matmul(hamp%elem-ham0%elem,gs)

        DO ii=1, num_epsilon
            epsilon_gs(ii) = abs(dot_product(hamiltonian%elem(:,ii),temp_gs))
        END DO
    
        WRITE(30,*) tt/bigt, epsilon_gs(1:)
        WRITE(80,*) tt/bigt, eigenvalues(1:)
      
        !debugging
        IF(debug.eqv..true.) THEN
            IF(step == 150) THEN
                WRITE(70,*) hamiltonian%elem(:,4), char(10), hamiltonian%elem(:,5)
                WRITE(70,*) char(10), 'temp gs: ', char(10), temp_gs
                WRITE(70,*) char(10), 'gs: ' ,char(10), gs
            END IF

            WRITE(20,*) gs, char(10)
            WRITE(60,*) hamiltonian%elem, char(10)
        END IF

        DEALLOCATE(eigenvalues)
    END DO

    ! debugging
    IF(debug.eqv..true.) THEN
        CLOSE(20)
        CLOSE(70)
        CLOSE(60)
    END IF

    CLOSE(30)
    CLOSE(80)

    DEALLOCATE(gs, epsilon_gs, temp_gs)

END SUBROUTINE

!************************************************************************************
! PARITY
!*************************************************************************************
SUBROUTINE PARITY_CHECK(ham0, hamp, NN)
    TYPE(RMATRIX) :: sigma_z,ham0, hamp, ham_parity, hamiltonian, commutator
    INTEGER :: ii, NN

    hamiltonian%N(1) = ham0%N(1)
    hamiltonian%N(2) = ham0%N(2)
    ALLOCATE(hamiltonian%elem(hamiltonian%N(1),hamiltonian%N(2)))
    commutator%N(1) = ham0%N(1)
    commutator%N(2) = ham0%N(2)
    ALLOCATE(commutator%elem(hamiltonian%N(1),hamiltonian%N(2)))

    sigma_z = pauli_zz()

    ham_parity%elem = sigma_z%elem
    ham_parity%N(1) = 2
    ham_parity%N(2) = 2

    DO ii=1, NN-1
        ham_parity = TENSOR_PRODUCT(ham_parity,sigma_z)
    END DO

    CALL PRINTFILE(ham_parity, 'parity_matrix.dat')
    CALL PRINTFILE(ham0, 'ham0.dat')
    CALL PRINTFILE(hamp, 'hamp.dat')
    commutator%elem = matmul(hamp%elem,ham0%elem) - matmul(ham0%elem,hamp%elem)
    CALL PRINTFILE(commutator,'comm_ham_0_p')
    commutator%elem = matmul(hamp%elem, ham_parity%elem) - matmul(ham_parity%elem,hamp%elem)
    CALL PRINTFILE(commutator,'comm_hamp_parity')
    commutator%elem = matmul(ham0%elem, ham_parity%elem) - matmul(ham_parity%elem,ham0%elem)
    CALL PRINTFILE(commutator,'comm_ham0_parity')

    hamiltonian%elem = (1-0.5)*hamp%elem - 0.5*ham0%elem
    CALL PRINTFILE(hamiltonian, 'hamiltonian_check.dat')

    commutator%elem = matmul(hamiltonian%elem, ham_parity%elem) - matmul(ham_parity%elem, hamiltonian%elem)
    CALL PRINTFILE(commutator, 'comm_parity_ham.dat')

    DEALLOCATE(ham_parity%elem, hamiltonian%elem, commutator%elem)

END SUBROUTINE PARITY_CHECK

!*************************************************************************************
! GS PARITY CHECK
!*************************************************************************************
FUNCTION GS_PARITY(NN)
    TYPE(RMATRIX) :: ham_parity
    INTEGER :: NN
    REAL*8 :: gs_parity

    ham_parity = PARITY_HAMILTONIAN(NN)
    gs_parity = ham_parity%elem(2**NN,2**NN)

    OPEN(unit=10,file='gs_parity_vs_dim.dat',status='replace')
    WRITE(10,*) NN, ' ', gs_parity, char(10)

    RETURN

END FUNCTION GS_PARITY

!**************************************************************************************
! PARITY HAMILOTNIAN
!**************************************************************************************
FUNCTION PARITY_HAMILTONIAN(NN)
    TYPE(RMATRIX) :: parity_hamiltonian, sigma_z
    INTEGER :: ii, NN

    sigma_z = pauli_zz()

    parity_hamiltonian%elem = sigma_z%elem
    parity_hamiltonian%N(1) = 2
    parity_hamiltonian%N(2) = 2

    DO ii=1, NN-1
        parity_hamiltonian = TENSOR_PRODUCT(parity_hamiltonian,sigma_z)
    END DO
    RETURN
END FUNCTION PARITY_HAMILTONIAN

!**************************************************************************************
! CRANCK NICHOLSON TIME EVOLUTION
!*************************************************************************************
SUBROUTINE CN_TIME_EVOLUTION(ham0, hamp, gs_start, bigt, nsteps, gs)
    
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: gs, application, phi_vec
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gs_start
    REAL*8 :: deltat, bigt, tt
    INTEGER :: nsteps, step, dim, ii
    TYPE(RMATRIX) :: ham0, hamp, hamiltonian, id
    TYPE(DMATRIX) :: temp1, temp2
    LOGICAL :: debug, verbose
    
    debug = .false.
    verbose = .true.
    
    deltat = bigt/float(nsteps)
    dim = ham0%N(1)
    
    DO ii=1, dim
        gs(ii) = complex(gs_start(ii),0.0)
    END DO
    
    print*, 'deltat= ',deltat
    
    hamiltonian%N(1) = ham0%N(1)
    hamiltonian%N(2) = ham0%N(2)
    ALLOCATE(hamiltonian%elem(hamiltonian%N(1),hamiltonian%N(2)))
    temp1%N(1) = ham0%N(1)

    !THE TWO HAMILTONIANS AND THE GS
    temp1%N(2) = ham0%N(2)
    ALLOCATE(temp1%elem(hamiltonian%N(1),temp1%N(2)))
    temp2%N(1) = ham0%N(1)
    temp2%N(2) = ham0%N(2)
    ALLOCATE(temp2%elem(hamiltonian%N(1),temp1%N(2)))

    id = IDENTITY(dim)
    
    ALLOCATE(phi_vec(dim))
    ALLOCATE(application(dim))
    
    OPEN(unit=30, file='cn_evol.dat',status='replace')
    OPEN(unit=40, file='cn_eigenvalues_evol.dat',status='replace')
    tt = 0
    CALL TOTAL_MF(ham0, hamp, tt, bigt, hamiltonian)
    WRITE(30,*) gs, char(10)
    WRITE(40,*) 0.0, real(dot_product(gs,matmul(hamiltonian%elem,gs)))
    
    DO step=1, nsteps
        IF(verbose.eqv..true.) THEN
            print*, step
        END IF
            tt = deltat*float(step-1)           
        CALL TOTAL_MF(ham0, hamp, tt, bigt, hamiltonian)
        temp1%elem = id%elem + complex(0.0,0.5*deltat)*hamiltonian%elem
        temp2%elem = id%elem - complex(0.0,0.5*deltat)*hamiltonian%elem
        phi_vec = matmul(temp2%elem, gs)
        CALL SYSTEM_SOLVE(temp1, gs, phi_vec)

        application = matmul(hamiltonian%elem,gs)  
        WRITE(30,*) gs, char(10)
        WRITE(40,*) (tt+deltat)/bigt, real(dot_product(gs, application))
    END DO
    
    CLOSE(30)
    CLOSE(40)

    DEALLOCATE(temp1%elem, temp2%elem, application, hamiltonian%elem)
END SUBROUTINE CN_TIME_EVOLUTION

END MODULE

!*************************************************************************************
!       MAIN
!*************************************************************************************
PROGRAM MAIN

    USE MATRICES
    USE AQC
    IMPLICIT NONE

    TYPE(RMATRIX) :: ham1, ham2
    INTEGER :: NN, nsteps, decision
    REAL*8 :: bigt, par
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE :: gs_evol, gs_exp
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gs_start

    NN = 3
    bigt = 1000
    nsteps = 1000

    !THE TWO HAMILTONIANS AND THE GS
    ham1 = ham_0(NN)
    !ham2 = mf_ham_p(NN)     ! mean field ising
    ham2 = ham_p(NN)

    ALLOCATE(gs_start(2**NN))
    ALLOCATE(gs_exp(2**NN))
    ALLOCATE(gs_evol(2**NN))
    gs_start = GROUNDSTATE(ham1)

    CALL TRUE_DURING_TIME(ham1, ham2, bigt, nsteps, NN)
    
    DEALLOCATE(gs_evol,gs_exp,gs_start, ham1%elem, ham2%elem)

END PROGRAM

