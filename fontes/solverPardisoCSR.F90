  module mSolverPardiso

        use mEstruturasDadosSistEq
        
        logical :: primeiravezVel

      contains

!=======================================================================
!    
      subroutine solverPardisoEsparso(estrutSistEq_, parte, label)

      use mEstruturasDadosSistEq 

      implicit none 
!
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(LEN=4), intent(in) :: parte
      character(LEN=*), intent(in) :: label

! 
        call solverPardisoPPD_Nodal(estrutSistEq_%Ap, estrutSistEq_%Ai, estrutSistEq_%alhs, estrutSistEq_%brhs, &
                     estrutSistEq_%neq, estrutSistEq_%nalhs, estrutSistEq_%pt, estrutSistEq_%iparm,           &
                     estrutSistEq_%dparm, estrutSistEq_%simetria, parte)

      end subroutine solverPardisoEsparso
!
!=======================================================================
!  

      subroutine solverPardisoPPD_Nodal(ia, ja, a, b, neq, nonzeros, pt, iparm, dparm, simetria, parte)
      implicit none

        INTEGER*4, intent(in)  :: ia(neq+1), ja(nonzeros)
        REAL*8, INTENT(IN)     :: a(nonzeros)
        REAL*8, INTENT(INOUT)        :: b(neq)
        INTEGER*4, INTENT(IN)  :: NEQ, NONZEROS

        INTEGER*4, INTENT(INOUT) :: pt(64), iparm(64)        
        real*8, INTENT(INOUT)    :: dparm(64)
        LOGICAL, INTENT(IN)          :: simetria
        character(LEN=4), INTENT(IN) :: parte

        REAL*8,  save :: ddum
        INTEGER, save :: maxfct, mnum, mtype, phase, nrhs, error,  msglvl
        INTEGER, save :: idum, solver
        INTEGER :: i
        REAL*8  :: t1, t2 , tt1, tt2
        integer :: omp_get_num_threads
        
        real*8, allocatable :: x(:)
        
! 
#ifdef withPardiso
!
!  .. Setup Pardiso control parameters und initialize the solvers     
!     internal adress pointers. This is only necessary for the FIRST   
!     call of the PARDISO solver.                                     
!     
! !$OMP PARALLEL

      if(parte=='reor'.or. parte=='full') then

      iparm=0
      !iparm(0)=0
      mtype     = -2   ! real and symmetric matrix, inefinite
      iparm(11) = 0    !Do not use (symmetric matrices).
      if(simetria.eqv..false.)  mtype     = 11   ! unsymmetric matrix, indefinite

      solver    = 0    ! use 0 for sparse direct method or 1 multi-recursive iterative solver
      msglvl    = 0    ! with statistical information
      iparm(33) = 0    ! compute determinant 
      iparm(52) = 1    !For OpenMP-threaded solver
      iparm(2)  = 0 !ou 0?      !Fill-In reduction reordering.
      iparm(27) = 1
      
      nrhs      = 1
      mnum      = 1
      pt       = 0
      idum      = 0
      ddum      = 0.0
      dparm    = 0.0
      maxfct    = 1
      error     = 0

!
!  .. Numbers of Processors ( value of OMP_NUM_THREADS )
!
       iparm(3) = 1
! 
#ifdef withOMP
!$OMP PARALLEL
       iparm(3) =   omp_get_num_threads()
!$OMP END PARALLEL
#endif
       print*, "em pardiso G com, numThreads", iparm(3) 

!        
       call timing(tt1)
!
!  .. PARDISO license check and initialize solver
!       call pardisoinit(pt, mtype, solver, iparm, dparm, error)
!
      IF (error .NE. 0) THEN
        IF (error.EQ.-10 ) WRITE(*,*) 'No license file found'
        IF (error.EQ.-11 ) WRITE(*,*) 'License is expired'
        IF (error.EQ.-12 ) WRITE(*,*) 'Wrong username or hostname'
        STOP ' (error .NE. 0) '
      ELSE
  !      WRITE(*,*) '[PARDISO]: License check was successful ... '
      END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!..   Reordering and Symbolic Factorization, This step also allocates
!     all memory that is necessary for the factorization
!
         call timing(t1)
         phase     = 11     ! only reordering and symbolic factorization

         CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                    idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm)
         call timing(t2)

#ifdef mostrarTempos
         write(*,*) "reordering, tempo de parede = ", t2 - t1
#endif

      IF (error .NE. 0) THEN
        WRITE(*,*) 'The following ERROR was detected: ', error
        STOP
      END IF


      endif       !if(parte=='reor'.or. parte=='full') then

      if(parte=='fact'.or. parte=='full') then
!
!.. Factorization.
!
      WRITE(*,*) 'Begining factorization  ... '
       call timing(t1)
      phase     = 22  ! only factorization
      CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja,  &
                    idum, nrhs, iparm, msglvl, ddum, ddum, error, dparm) 
      call timing(t2)
#ifdef mostrarTempos
        write(*,*) "factorization, tempo de parede = ", t2 - t1
#endif

    IF (error .NE. 0) THEN
       WRITE(*,*) 'The following ERROR was detected: ', error
      STOP
    ENDIF 

      endif  !     if(parte=='fact'.or. parte=='full') then

! 
      if(parte=='back'.or. parte=='full') then
      
      allocate(x(neq)); x=0.d0

!.. Back substitution and iterative refinement
      WRITE(*,*) 'Begining backsubstitution  ... '
       call timing(t1)
      iparm(8)  = 1   ! max numbers of iterative refinement steps
      phase     = 33  ! only solve
      CALL pardiso (pt, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
                   idum, nrhs, iparm, msglvl, b, x, error, dparm) 
       call timing(t2)
#ifdef mostrarTempos
      write(*,*) "backsubstitution, tempo de parede = ", t2 - t1
#endif

       call timing(tt2)
! #ifdef mostrarTempos
!       WRITE(*,*) 'Solve completed ...  ',label, ", tempo =", tt2-tt1
! #endif

      b(1:neq) = x(1:neq)
      deallocate(x)

      endif !if(parte=='back'.or. parte=='full') then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!       if(parte=='full') then
! 
! !       desalocando memoria
!       WRITE(*,*) 'Begining memory desalocation  ... '
!       call timing(t1)
!       phase = -1
! 
!       CALL pardiso (ptP, maxfct, mnum, mtype, phase, neq, a, ia, ja, &
!                      idum, nrhs, iparmP, msglvl, b, xP, error, dparmP)
!          
!       call timing(t2)
! #ifdef mostrarTempos
! !       write(*,*) "memory desalocation: ", label, ", tempo = ", t2 - t1
! #endif
! 
!       endif !if(parte=='full') then

#endif
      end subroutine 

!

! **** new *********************************************************************
!
       subroutine addlhsCSR(alhs,eleffm,Ap,Ai,lm,nee)
! 
! .... program to add element left-hand-side matrix to
!         global left-hand-side matrix
! 
! 
       implicit none
! 
! .... deactivate above card(s) for single-precision operation
! 
       real*8  :: alhs(*),eleffm(nee,*)
       integer :: lm(*), Ap(*), Ai(*)
       integer :: nee
!
       integer :: i, j, linha, coluna
       integer :: inicio, fim, jj
       character(40) :: formatoALHS
       formatoALHS = "(10f8.4)"
! 
     !write(*,*) "em addlhsCSR, lm =  ", lm(1:nee)
     !write(*,formatoALHS) eleffm(1,1:nee)
     !write(*,formatoALHS) eleffm(2,1:nee)
     !write(*,formatoALHS) eleffm(3,1:nee)
     !if(nee==4) write(*,formatoALHS) eleffm(4,1:nee)
          do   j=1,nee
          coluna = abs(lm(j))
          if (coluna.gt.0) then
             do i=1,nee
             linha = abs(lm(i))
             if (linha.gt.0) then
                inicio=Ap(coluna)
                fim=Ap(coluna+1)-1
                 do jj=inicio, fim
                      if(Ai(jj)==linha) then
                          alhs(jj)= alhs(jj)+eleffm(i,j)
                          !write(*,'(3(1x,i0), 2e10.3)') i, j, jj, eleffm(i,j), alhs(jj)
                      endif
                 enddo
             endif
           enddo
! 
          endif
      enddo
!     write(*,formatoAlhs) "em addlhsCSR, alhs =  ",alhs(1:6)
! 
       return
       end subroutine addlhsCSR
!     
     end module mSolverPardiso
