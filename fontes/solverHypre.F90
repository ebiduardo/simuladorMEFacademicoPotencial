 module mSolverHypre
 
      use mEstruturasDadosSistEq
      implicit none
#ifdef withHYPRE
      include 'mpif.h'
      !use  mpi
#endif
         
      integer*4 :: myid, num_procs
      integer*4 :: mpi_comm

      integer*4 :: posPonteiro, contPonteiro, posColunas, posCoef
      integer*4 :: nonzeros, nonzerosEst
      logical   :: primeiravezVel, primeiravezGeo
      integer*8, parameter :: HYPRE_PARCSR=5555

INTERFACE criarVetorHypre
      MODULE PROCEDURE criarVetorHypre_A, criarVetorHypre_B
END INTERFACE criarVetorHypre 

INTERFACE criarMatrizHYPRE
      MODULE PROCEDURE criarMatrizHYPRE_A, criarMatrizHYPRE_B, criarMatrizHYPRE_C
END INTERFACE criarMatrizHYPRE 
INTERFACE resolverSistemaAlgHYPRE
      MODULE PROCEDURE resolverSistemaAlgHYPRE_A, resolverSistemaAlgHYPRE_B, & 
       resolverSistemaAlgHYPRE_C!, resolverSistemaAlg_HYPRE
END INTERFACE resolverSistemaAlgHYPRE 
INTERFACE fecharMatrizHypre
      MODULE PROCEDURE fecharMatrizHypre_A, fecharMatrizHypre_B
END INTERFACE fecharMatrizHypre 
INTERFACE fecharVetorHYPRE
      MODULE PROCEDURE fecharVetorHYPRE_A, fecharVetorHYPRE_B, fecharVetorHYPRE_C
END INTERFACE fecharVetorHYPRE 

INTERFACE addnslHYPRE
      MODULE PROCEDURE addnslHYPRE_A, addnslHYPRE_B
END INTERFACE addnslHYPRE

      contains

      subroutine lixoB()
      end subroutine lixoB


      subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)
      integer*4, intent(out) :: myid_, num_procs_
      integer*4 :: mpi_comm_
      integer*4 ::  ierr  
!      write(*,'(a)') " ++++  em subroutine inicializarMPI(myid_, num_procs_, mpi_comm_)"
       myid_=-1; num_procs_=-1; mpi_comm_=-1
#ifdef withHYPRE
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid_,      ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs_, ierr)
!   Convert the Fortran MPI communicator to a C version that can be used in hypre.
!   Uncomment the line below if you are using MPICH
    !myid_      = myid
    !num_procs_ = num_procs
    mpi_comm_  = MPI_COMM_WORLD
!   Uncomment the line below if you are using LAM-MPI
      !call HYPRE_MPI_Comm_f2c(mpi_comm, MPI_COMM_WORLD, ierr)
#endif
      write(*,'(3(a,i0))')  " em subroutine inicializarMPI, myid= ", myid_, &
                   " numprocs= ", num_procs_,", mpi_comm= ", mpi_comm_
!=============================================================================
      end subroutine inicializarMPI
      subroutine finalizarMPI                 ()
      integer*4 :: ierr=-1
#ifdef withHYPRE
         call MPI_Finalize(ierr)
#endif
      end subroutine finalizarMPI   
!=============================================================================
      subroutine resolverSistemaAlgHYPRE_C  (estrutSistEq_) 
      implicit none
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      !write(*,*)" em  resolverSistemaAlgHYPRE_C  (estrutSistEq_),", estrutSistEq_%A_HYPRE, estrutSistEq_%b_HYPRE, estrutSistEq_%u_HYPRE , estrutSistEq_%solver
      call resolverSistemaAlgHYPRE_A( estrutSistEq_%A_HYPRE, estrutSistEq_%parcsr_A, &
                       estrutSistEq_%b_HYPRE, estrutSistEq_%par_b, &
                       estrutSistEq_%u_HYPRE, estrutSistEq_%par_u, &
                       estrutSistEq_%solver, & 
                      estrutSistEq_%solver_id, estrutSistEq_%precond_id, & 
                      estrutSistEq_%tol, estrutSistEq_%num_iterations, estrutSistEq_%final_res_norm, &
                      estrutSistEq_%initialGuess, estrutSistEq_%brhs, estrutSistEq_%rows, &
                      estrutSistEq_%neq, myid, mpi_comm)
!      call resolverSistemaAlgHYPRE_B( estrutSistEq_%A_HYPRE, estrutSistEq_%b_HYPRE, estrutSistEq_%u_HYPRE, &
!                      estrutSistEq_%solver, estrutSistEq_%solver_id, estrutSistEq_%precond_id, & 
!                      estrutSistEq_%tol, estrutSistEq_%num_iterations, estrutSistEq_%final_res_norm, &
!                      estrutSistEq_%initialGuess, estrutSistEq_%brhs, estrutSistEq_%rows, &
!                      estrutSistEq_%neq, myid, mpi_comm)
      end subroutine resolverSistemaAlgHYPRE_C

!=============================================================================
      subroutine resolverSistemaAlgHYPRE_B  (A_, b_, u_, solver_,  solver_id_, precond_id_, &
                                             tol_, num_iterations_, final_res_norm_, &
                                             initialGuess_, brhs_, rows_, neq_, myid_, mpi_comm_)  

      implicit none
      integer*8, intent(in)         :: A_, b_, u_
      integer*8, intent(out)        :: solver_
      integer*4, intent(in)         :: solver_id_, precond_id_
      double precision,  intent(in) :: tol_ ! = 1.0e-08 
      integer*4, intent(out)        :: num_iterations_
      double precision, intent(out) :: final_res_norm_
      double precision, intent(in)  :: initialGuess_(neq_), brhs_(neq_)
      integer*4, intent(in)         :: rows_(neq_), neq_, myid_, mpi_comm_ 

      !integer*8         :: parcsr_A, par_b, par_u

      !write(*,*)" em  resolverSistemaAlgHYPRE_B  (estrutSistEq_),", A_, b_, u_, solver_ 
      solver_=0;num_iterations_=0; final_res_norm_=0;
      !call resolverSistemaAlg_HYPRE (A_, parcsr_A, b_, par_b, u_, par_u, solver_, solver_id_, precond_id_, &
      !           tol_, num_iterations_, final_res_norm_, initialGuess_, brhs_, rows_, neq_, myid_, mpi_comm_)  
      !call resolverSistemaAlgHYPRE_A (A_, parcsr_A, b_, par_b, u_, par_u, solver_, solver_id_, precond_id_, &
      !           tol_, num_iterations_, final_res_norm_, initialGuess_, brhs_, rows_, neq_, myid_, mpi_comm_, 0)  

      end subroutine resolverSistemaAlgHYPRE_B 


!=============================================================================
      subroutine resolverSistemaAlgHYPRE_A  (A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_, solver_id_, precond_id_, &
                                             tol_, num_iterations_, final_res_norm_, initialGuess_, brhs_, rows_, &
                                             neq_, myid_, mpi_comm_) !, lixo)  

      implicit none
      integer*8, intent(in)         :: A_, b_, u_
      integer*8, intent(inout)      :: parcsr_A_, par_b_, par_u_, solver_
      integer*4, intent(in)         :: solver_id_, precond_id_
      double precision,  intent(in) :: tol_ ! = 1.0e-08 
      integer*4, intent(out)        :: num_iterations_
      double precision, intent(out) :: final_res_norm_
      double precision, intent(in)  :: initialGuess_(neq_), brhs_(neq_)
      integer*4, intent(in)         :: rows_(neq_), neq_, myid_, mpi_comm_ 
      integer*4 :: lixo

      integer*4         :: ierr, printLevel=0, numMaxIterations = 10000
      real*8            :: t1, t2, t3, t4
      real*8,    save   :: elapsedT 
      integer*8, save   :: solverAnt=-1
      character(LEN=80) :: mensagem

      !block
      character(len=12) :: nomeA= "Alhs00.out."
      character(len=12) :: nomeB= "brhs00.out."
      character(len=15) :: nomeU= "solucao00.out."
      !end block
      write(*,*)" em  resolverSistemaAlgHYPRE_A  (estrutSistEq_),", A_, b_, u_ , solver_
      write(*,*) "solver_id_ = ", solver_id_, ", precond_id_ = ", precond_id_, ", tol_ = ", tol_
      final_res_norm_=0.0; num_iterations_=0
#ifdef withHYPRE
!       call HYPRE_IJMatrixPrint( A_, nomeA, ierr)
!       call HYPRE_IJVectorPrint( b_, nomeB, ierr)
      !write(*,*) " set initialGuess_, " !, initialGuess_(1:5)
      call HYPRE_IJVectorSetValues (u_, neq_, rows_, initialGuess_, ierr)
!      Get parcsr matrix object
      call HYPRE_IJMatrixGetObject( A_, parcsr_A_, ierr)
!      Get parcsr brhs vector 
      call HYPRE_IJVectorGetObject (b_, par_b_, ierr)
!      Get parcsr solution vector 
      call HYPRE_IJVectorGetObject (u_, par_u_, ierr)

        mensagem=".                   :"
      write(*,'(2a,i5,a,e13.5)') trim(mensagem), "  num MAX it = ", numMaxIterations, " , tol = ", tol_
      call timing(t1)
      if     (solver_id_==0) then
        mensagem="BoomerAMG solver_id 0"
        call solverOptionA ()
      elseif (solver_id_==1) then
        mensagem="         PCG with AMG"
        call solverOptionB ()
      elseif (solver_id_==2) then
        mensagem="   PCG with ParaSails"
        call solverOptionC ()
      else
         if (myid_ .eq. 0) then 
           print *,'Invalid solver id specified'
           stop
         endif  
      endif
      write(*,'(a)', ADVANCE="NO") trim(mensagem)
      write(*,'(a,i5,a,e13.5)') ", num it = ", num_iterations_+1, " , residual norm = ", final_res_norm_
      call timing(t2)
      elapsedT = elapsedT+t2 - t1
      solverAnt = solver_
      if( num_iterations_ == numMaxIterations) then 
           write(*,*) "num_iterations_=",num_iterations_, ", ==  numMaxIterations =",  numMaxIterations
           write(*,*) "verifique suas escolhas!!! "
      end if
      call flush(6)
!
      contains

      subroutine solverOptionA
!        Create solver
         call HYPRE_BoomerAMGCreate          (solver_, ierr)
!        Set some parameters (See Reference Manual for more parameters)
!        print solve info + parameters 
         call HYPRE_BoomerAMGSetPrintLevel   (solver_, printLevel,   ierr)  
!        Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (solver_, 6,    ierr) 
!        G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (solver_, 6,    ierr)     
!        Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (solver_, 1,    ierr)  
!        maximum number of levels 
         call HYPRE_BoomerAMGSetMaxLevels    (solver_, 200,  ierr) 
!        conv. tolerance
         call HYPRE_BoomerAMGSetTol          (solver_, tol_, ierr)    
         call HYPRE_BoomerAMGSetTol          (solver_, 0.0d-6, ierr)    
!        Now setup and solve!
         call HYPRE_BoomerAMGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr )
         call HYPRE_BoomerAMGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr )
!        Run info - needed logging turned on 
         call HYPRE_BoomerAMGGetNumIterations(solver_, num_iterations_, ierr)
         call HYPRE_BoomerAMGGetFinalReltvRes(solver_, final_res_norm_, ierr)
!        Destroy solver_
         call HYPRE_BoomerAMGDestroy         (solver_, ierr )
      end subroutine solverOptionA

      subroutine solverOptionB ()
         integer*8:: precond
         call HYPRE_ParCSRPCGCreate        (MPI_COMM_WORLD, solver_, ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (solver_, printLevel,   ierr)  
!        Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter    (solver_, numMaxIterations,   ierr)
         call HYPRE_ParCSRPCGSetTol        (solver_, tol_, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (solver_, 1,      ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (solver_, printLevel,      ierr)
!        Now set up the AMG preconditioner and specify any parameters
         call HYPRE_BoomerAMGCreate        (precond, ierr)
         call HYPRE_BoomerAMGSetPrintLevel (precond, printLevel,   ierr)  
!        Set some parameters (See Reference Manual for more parameters)
!  Relaxation Parameters:
!   Visiting Grid:                     down   up  coarse
!            Number of sweeps:            1    1     1 
!   Type 0=Jac, 3=hGS, 6=hSGS, 9=GE:      6    6     6 
!        print less solver_ info since a preconditioner
         call HYPRE_BoomerAMGSetPrintLevel   (precond, printLevel,     ierr); 
!        Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (precond, 6,     ierr) 
!        SYMMETRIC G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (precond, 6,     ierr)     
!        Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (precond, 1,     ierr)  
!        conv. tolerance
         call HYPRE_BoomerAMGSetTol          (precond, tol_, ierr)     
!        do only one iteration! 
         call HYPRE_BoomerAMGSetMaxIter      (precond, 1,     ierr)
!        set amg as the pcg preconditioner
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id_, precond, ierr)
!        Now setup and solve!
         call HYPRE_ParCSRPCGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr)
!        Run info - needed logging turned on 
         call HYPRE_ParCSRPCGGetNumIterations (solver_, num_iterations_, ierr)
         call HYPRE_ParCSRPCGGetFinalRelative (solver_, final_res_norm_, ierr)
!       Destroy precond and solver
         call HYPRE_BoomerAMGDestroy          (precond, ierr )
         call HYPRE_ParCSRPCGDestroy          (solver_, ierr)
      end subroutine solverOptionB

      subroutine solverOptionC()

!        Create solver
         integer*8:: precond
         call HYPRE_ParCSRPCGCreate          (MPI_COMM_WORLD, solver_, ierr)

!        Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter      (solver_, 1000, ierr)
         call HYPRE_ParCSRPCGSetTol          (solver_, tol_, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm      (solver_, 1, ierr)
         call HYPRE_ParCSRPCGSetPrintLevel   (solver_, printLevel, ierr)
         call HYPRE_ParCSRPCGSetLogging      (solver_, 1, ierr)
!        Now set up the Parasails preconditioner and specify any parameters
         call HYPRE_ParaSailsCreate          (MPI_COMM_WORLD, precond,ierr)
         call HYPRE_ParaSailsSetParams       (precond, 0.1d0, 1, ierr)
         call HYPRE_ParaSailsSetFilter       (precond, 0.05d0, ierr)
         call HYPRE_ParaSailsSetSym          (precond, 1)
         call HYPRE_ParaSailsSetLogging      (precond, 3, ierr)
!        set parsails as the pcg preconditioner
         !precond_id_ = 4
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id_, precond, ierr)
!        Now setup and solve!
         call HYPRE_ParCSRPCGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr)
!        Run info - needed logging turned on 
        call HYPRE_ParCSRPCGGetNumIterations (solver_, num_iterations_, ierr)
        call HYPRE_ParCSRPCGGetFinalRelative (solver_, final_res_norm_, ierr)
!       Destroy precond and solver
        call HYPRE_ParaSailsDestroy          (precond, ierr )
        call HYPRE_ParCSRPCGDestroy          (solver_, ierr)
      end subroutine solverOptionC

#endif

      end subroutine resolverSistemaAlgHYPRE_A


      subroutine resolverSistemaAlg_HYPRE  (A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_, solver_id_, precond_id_, &
                                             tol_, num_iterations_, final_res_norm_, &
                                             initialGuess_, brhs_, rows_, &
                                             neq_, myid_, mpi_comm_)  

      implicit none
      integer*8, intent(in)         :: A_, b_, u_
      integer*8, intent(inout)      :: parcsr_A_, par_b_, par_u_, solver_
      integer*4, intent(in)         :: solver_id_, precond_id_
      double precision,  intent(in) :: tol_ ! = 1.0e-08 
      integer*4, intent(out)        :: num_iterations_
      double precision, intent(out) :: final_res_norm_
      integer*4, intent(in)         :: neq_
      integer*4, intent(in)         :: rows_(neq_)
      double precision, intent(in)  :: initialGuess_(neq_), brhs_(neq_)
      integer*4, intent(in)         :: myid_
      integer*4, intent(in)         :: mpi_comm_ 

      integer*8        :: precond
      integer*4        :: precond_id;
      integer*4        :: ierr
      integer*4, save  :: contador = 1
      integer*4        :: local_size , numMaxIterations = 10000
      integer*4        :: printLevel, solver_id
      integer*4        :: i, Flower, Fupper
      real*8           :: t1, t2, t3, t4, elapsedT 

      !block
      character(len=15) :: nomeU 
      character(len=12) :: nomeB
      character(len=12) :: nomeA
          
      nomeU = "solucao00.out."
      nomeA = "Alhs00.out."
      nomeB = "brhs00.out."
    !  call HYPRE_IJVectorPrint( b_, nomeB, ierr)
    !  call HYPRE_IJMatrixPrint( A_, nomeA, ierr)
      !end block
      !write(*,*) "em subroutine resolverSistemaAlg_HYPRE  (A_, parcsr_A_, b_, par_b_, "
!      write(*,'(3(a,i0))') "s+++, em resolverSistemaAlgHYPRE, myid= ", myid,", &
!                    numprocs= ", num_procs,", mpi_comm= ", mpi_comm
      local_size = neq_
      printLevel = 0
#ifdef withHYPRE

      call HYPRE_IJVectorSetValues (u_, neq_, rows_, initialGuess_(1), ierr)
 !     call HYPRE_IJVectorSetValues (b_, neq_, rows_, brhs_, ierr)
 !     call HYPRE_IJVectorPrint( b_, nomeB, ierr)
 !     call HYPRE_IJVectorPrint( A_, nomeA, ierr)

      !write(*,*) rows_(1:neq_)
      !write(*,*) brhs_(1:neq_)
!     write(*,*) "solver_id_=", solver_id_
!     write(*,*) "precond =",  precond
!     write(*,*) "precond_id_ =",  precond_id_
!     write(*,*) "neq_ =", neq_
!     write(*,*) "++ precond_id_ =", precond_id_,", tol =", tol_

      call timing(t1)
      if     (solver_id_==0) then
        write(*,*) " ... BoomerAMGCreate , BD"
        call solverOptionA ()

      elseif (solver_id_==1) then
        write(*,*) " ... PCG with AMG preconditioner, BD_dez";
        call solverOptionB ()

      elseif (solver_id_==2) then
        write(*,*) " ... PCG with ParaSails"
        call solverOptionC ()
      else
         if (myid_ .eq. 0) then 
           print *,'Invalid solver id specified'
           stop
         endif  
      endif
      call timing(t2)
      elapsedT = t2 - t1

      contains

      subroutine solverOptionA
!        Create solver
         call HYPRE_BoomerAMGCreate          (solver_, ierr)

!        Set some parameters (See Reference Manual for more parameters)

!        print solve info + parameters 
         call HYPRE_BoomerAMGSetPrintLevel   (solver_, printLevel,      ierr)  
!        Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (solver_, 6,      ierr) 
!        G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (solver_, 6,      ierr)     
!        Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (solver_, 1,      ierr)  
!        maximum number of levels 
         call HYPRE_BoomerAMGSetMaxLevels    (solver_, 200,     ierr) 
!        conv. tolerance
         call HYPRE_BoomerAMGSetTol          (solver_, tol_, ierr)    

!        Now setup and solve!
         call HYPRE_BoomerAMGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr )
         call HYPRE_BoomerAMGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr )

!        Run info - needed logging turned on 
         call HYPRE_BoomerAMGGetNumIterations(solver_, num_iterations_, ierr)
         call HYPRE_BoomerAMGGetFinalReltvRes(solver_, final_res_norm_, ierr)

!        Destroy solver_
         call HYPRE_BoomerAMGDestroy         (solver_, ierr )

!     PCG with AMG preconditioner
      end subroutine solverOptionA

      subroutine solverOptionB ()
!        Create solver_
         call HYPRE_ParCSRPCGCreate        (MPI_COMM_WORLD, solver_, ierr)
!        Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter    (solver_, numMaxIterations,   ierr)
         call HYPRE_ParCSRPCGSetTol        (solver_, tol_, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm    (solver_, 1,      ierr)
!         call HYPRE_ParCSRPCGSetPrintLevel (solver_, 2,      ierr)
         call HYPRE_ParCSRPCGSetPrintLevel (solver_, printLevel,      ierr)

!        Now set up the AMG preconditioner and specify any parameters

         call HYPRE_BoomerAMGCreate          (precond, ierr)

!        Set some parameters (See Reference Manual for more parameters)

!  Relaxation Parameters:
!   Visiting Grid:                     down   up  coarse
!            Number of sweeps:            1    1     1 
!   Type 0=Jac, 3=hGS, 6=hSGS, 9=GE:      6    6     6 

!        print less solver_ info since a preconditioner
         !call HYPRE_BoomerAMGSetPrintLevel   (precond, printLevel,     ierr); 
!        Falgout coarsening
         call HYPRE_BoomerAMGSetCoarsenType  (precond, 6,     ierr) 
!        SYMMETRIC G-S/Jacobi hybrid relaxation 
         call HYPRE_BoomerAMGSetRelaxType    (precond, 6,     ierr)     
!        Sweeeps on each level
         call HYPRE_BoomerAMGSetNumSweeps    (precond, 1,     ierr)  
!        conv. tolerance
         call HYPRE_BoomerAMGSetTol          (precond, tol_, ierr)     
!        do only one iteration! 
         call HYPRE_BoomerAMGSetMaxIter      (precond, 1,     ierr)

!        set amg as the pcg preconditioner
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id_, precond, ierr)

!        Now setup and solve!
         call HYPRE_ParCSRPCGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr)
!        Run info - needed logging turned on 
        call HYPRE_ParCSRPCGGetNumIterations (solver_, num_iterations_, ierr)
        call HYPRE_ParCSRPCGGetFinalRelative (solver_, final_res_norm_, ierr)
        write(*,*) "num_iterations= ", num_iterations_; 
        write(*,*) "final_res_norm_= ", final_res_norm_; 

!       Destroy precond and solver
        call HYPRE_BoomerAMGDestroy          (precond, ierr )
        call HYPRE_ParCSRPCGDestroy          (solver_, ierr)
        stop
      end subroutine solverOptionB

      subroutine solverOptionC()

!        Create solver
         call HYPRE_ParCSRPCGCreate          (MPI_COMM_WORLD, solver_, ierr)

!        Set some parameters (See Reference Manual for more parameters) 
         call HYPRE_ParCSRPCGSetMaxIter      (solver_, 1000, ierr)
         call HYPRE_ParCSRPCGSetTol          (solver_, tol_, ierr)
         call HYPRE_ParCSRPCGSetTwoNorm      (solver_, 1, ierr)
         call HYPRE_ParCSRPCGSetPrintLevel   (solver_, printLevel, ierr)
!         call HYPRE_ParCSRPCGSetLogging      (solver_, 1, ierr)

!        Now set up the Parasails preconditioner and specify any parameters
         call HYPRE_ParaSailsCreate          (MPI_COMM_WORLD, precond,ierr)
         call HYPRE_ParaSailsSetParams       (precond, 0.1d0, 1, ierr)
         call HYPRE_ParaSailsSetFilter       (precond, 1.05d0, ierr)
         call HYPRE_ParaSailsSetSym          (precond, 1)
!         call HYPRE_ParaSailsSetLogging      (precond, 3, ierr)

!        set parsails as the pcg preconditioner
         !precond_id_ = 4
         call HYPRE_ParCSRPCGSetPrecond      (solver_, precond_id_, precond, ierr)


!        Now setup and solve!
         call HYPRE_ParCSRPCGSetup           (solver_, parcsr_A_, par_b_, par_u_, ierr)
         call HYPRE_ParCSRPCGSolve           (solver_, parcsr_A_, par_b_, par_u_, ierr)

!        Run info - needed logging turned on 

        call HYPRE_ParCSRPCGGetNumIterations (solver_, num_iterations_, ierr)
        call HYPRE_ParCSRPCGGetFinalRelative (solver_, final_res_norm_, ierr)

!       Destroy precond and solver
        call HYPRE_ParaSailsDestroy          (precond, ierr )
        call HYPRE_ParCSRPCGDestroy          (solver_, ierr)
      end subroutine solverOptionC

#endif

      end subroutine resolverSistemaAlg_HYPRE
!=============================================================================
!
      subroutine addnslHYPRE_geo(A_, eleffm,lmT,nee,nel)
!
!         program to add element left-hand-side matrix to          
!                global left-hand-side matrix                      
!                                                                  
!        ldiag = .true.,  add diagonal element matrix              
!                                                                  
!        ldiag = .false, then                                     
!        add full nonsymmetric element matrix                   
!                                                                  
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*8, intent(inout) :: A_
      integer*4:: nee, nel
      real*8  :: eleffm(nee,*)
      integer, pointer :: lmT(:,:,:)
      integer*4:: lm(nee)
!
      integer*4, parameter:: numCoefPorLinha=100 
      integer*4:: i,j,k,l, eq, c, nnz, ierr 
      integer*4::  cols(numCoefPorLinha)
      real*8  :: values(numCoefPorLinha) 
      integer, save :: contador = 0
#ifdef withHYPRE

      lm(:)=reshape(lmT(:,:,nel),(/nee/));
   
      contador = contador + 1
      !write(*,*) contador, "lm = ", lm(:)
         do 400 j=1,nee
            nnz = 0
            eq = abs(lm(j))
            if(eq > 0) then 
               do 200 i=1,nee
                   c = abs(lm(i))
                   if(c > 0) then 
                       nnz = nnz + 1
                       cols(nnz) = c
                     values(nnz) = eleffm(j,i)
                   end if
  200          continue
              ! cols = cols - 1
              ! eq = eq - 1
              !write(*,'(a,i5,a,10i5)') " eq = ", eq,  ", ", cols(1:nnz)
               call HYPRE_IJMatrixAddToValues(A_, 1, nnz, eq, cols, values, ierr)
            end if
  400    continue
#endif
!
      return
      end subroutine addnslHYPRE_geo

      subroutine addnsl_HYPRE(A_, eleffm,lmT,nee,nel)
!
!         program to add element left-hand-side matrix to          
!                global left-hand-side matrix                      
!                                                                  
!        ldiag = .true.,  add diagonal element matrix              
!                                                                  
!        ldiag = .false, then                                     
!        add full nonsymmetric element matrix                   
!                                                                  
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*8, intent(inout) :: A_
      integer*4:: nee, nel
      real*8  :: eleffm(nee,*)
      integer, pointer :: lmT(:,:,:)
      integer*4:: lm(nee)
!
      integer*4, parameter:: numCoefPorLinha=100 
      integer*4:: i,j,k,l, eq, c, nnz, ierr 
      integer*4::  cols(numCoefPorLinha)
      real*8  :: values(numCoefPorLinha) 
      integer, save :: contador = 0

      if( nel <=5) write(*,*) "em subroutine addnsl_HYPRE(A_, eleffm,lmT,nee,nel)"
#ifdef withHYPRE

      lm(:)=reshape(lmT(:,:,nel),(/nee/));
   
      contador = contador + 1
      !write(*,*) contador, "lm = ", lm(:)
         do 400 j=1,nee
            nnz = 0
            eq = abs(lm(j))
            if(eq > 0) then 
               do 200 i=1,nee
                   c = abs(lm(i))
                   if(c > 0) then 
                       nnz = nnz + 1
                       cols(nnz) = c
                     values(nnz) = eleffm(j,i)
                   end if
  200          continue
               cols = cols - 1  + 1 !Bjan23
               eq   =   eq - 1  + 1 !Bjan23
               !if(eq <= 5) write(*,'(a,i5,a,10i5)') " eq = ", eq,  ", ", cols(1:nnz)
               call HYPRE_IJMatrixAddToValues(A_, 1, nnz, eq, cols, values, ierr)
            end if
  400    continue
#endif
!
      return
      end subroutine addnsl_HYPRE
!*************************************************************************************************
      subroutine addnslHYPRE_A( A_HYPRE_,eleffm_, idiag_, lm_,  nee_, diag_, lsym_)
!     call addnslHYPRE(estrutSistEqP_%A_HYPRE, eleffm, lmLocal, estrutSistEqP_%idiag, nee, diag, lsym)
!
!         program to add element left-hand-side matrix to          
!                global left-hand-side matrix                      
!        ldiag = .true.,  add diagonal element matrix              
!        ldiag = .false, then                                     
!        add full nonsymmetric element matrix                   
!
      implicit none
      
       integer*8:: A_HYPRE_ 
       real*8  ::  eleffm_(:,:)
       integer*4:: lm_(:) , idiag_(:), nee_
       logical :: diag_, lsym_
       
!
      integer*4, parameter:: numCoefPorLinha=100 
      integer*4 :: i,j,k,l, eq, c, nnz, ierr
      integer*4 ::  cols(numCoefPorLinha)
      real*8  :: values(numCoefPorLinha) 
      integer*4, save :: contador = 0

      contador = contador + 1
         do 400 j=1,nee_
            nnz = 0
            eq = abs(lm_(j))
            if(eq > 0) then 
               do 200 i=1,nee_
                   c = abs(lm_(i))
                   if(c > 0) then 
                       nnz         = nnz + 1
                       cols(nnz)   = c
                       values(nnz) = eleffm_(j,i)
                   end if
  200          continue
               !if(eq <= 5) write(*,'(a,i5,a,10i5)') " eq = ", eq,  ", ", cols(1:nnz)
#ifdef withHYPRE
               call HYPRE_IJMatrixAddToValues(A_HYPRE_, 1, nnz, eq, cols, values, ierr)
#endif
            end if
  400    continue
!
      return
      end subroutine addnslHYPRE_A
!*************************************************************************************************
      subroutine addnslHYPRE_B(estrutSistEq_, eleffm, lm, nee)
!
!         program to add element left-hand-side matrix to          
!                global left-hand-side matrix                      
!                                                                  
!        ldiag = .true.,  add diagonal element matrix              
!                                                                  
!        ldiag = .false, then                                     
!        add full nonsymmetric element matrix                   
!                                                                  
!
      implicit none
!
!.... remove above card for single-precision operation
!
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_  
      
       real*8  :: eleffm(:,:)
       integer*4:: lm(:), nee
!
      integer*4, parameter:: numCoefPorLinha=100 
      integer*4 :: i,j,k,l, eq, c, nnz, ierr
      integer*4 ::  cols(numCoefPorLinha)
      real*8  :: values(numCoefPorLinha) 
      integer*4, save :: contador = 0

      write(*,*) "em addnslHYPRE_B"
 !     integer*4, allocatable:: lmLocal(:)
 !      nee = size(estrutSistEq_%lm,1)*size(estrutSistEq_%lm,2)
 !      allocate (lmLocal(nee))
 !      lmLocal(:)=reshape(estrutSistEq_%lm(:,:,nel),(/nee/))
   
 !     contador = contador + 1
     ! write(*,*) nel, "lm = ", lm(:)
         do 400 j=1,nee
            nnz = 0
            eq = abs(lm(j))
            if(eq > 0) then 
               do 200 i=1,nee
                   c = abs(lm(i))
                   if(c > 0) then 
                       nnz         = nnz + 1
                       cols(nnz)   = c
                       values(nnz) = eleffm(j,i)
                   end if
  200          continue
               cols = cols - 1 + 1
               eq   = eq   - 1 + 1
               !if(eq <= 5) write(*,'(a,i5,a,10i5)') " eq = ", eq,  ", ", cols(1:nnz)
#ifdef withHYPRE
              !write(*,'(i0, a,i3,a,10i3)') estrutSistEq_%A_HYPRE, ",  eq = ", eq,  ", ", cols(1:nnz)
!              write(*,'(i0, a,i3,a,10e12.3)') estrutSistEq_%A_HYPRE, ",  eq = ", eq,  ", ", values(1:nnz)
             !  write(*,*) "em ne addnslHYPRE, estrutSistEq_%A_HYPRE=", estrutSistEq_%A_HYPRE, nnz, eq, cols(1:nnz), values(1:nnz) 
               call HYPRE_IJMatrixAddToValues(estrutSistEq_%A_HYPRE, 1, nnz, eq, cols, values, ierr)
#endif
            end if
  400    continue
!
      return
      end subroutine addnslHYPRE_B
!
!**** new **********************************************************************
      subroutine addrhsHYPRE (b_, brhs_, elresf, lm, nee_)
!
!.... program to add element residual-force vector to
!        global right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*8 , intent(inout) ::  b_
      real*8, intent(in)  :: brhs_(*)
      integer*4:: nee_
      real*8  :: elresf(nee_)
      integer*4:: lm(:)
!
      integer*4:: k, j, nnz, ierr, lowerInd=1
      integer*4:: rows(nee_+20), i
      real*8   :: values(nee_+20)
!     integer, save :: nel = 0

      write(*,*) " em subroutine addrhsHYPRE (b_, brhs_, elresf,lm,nee_)"
      values = 0.0
!     nel = nel + 1
!       write(*,'(i5,a,16i4)') nel, "lm = ", lm(:)
      nnz = 0     
      do 100 j=1,nee_
      k = lm(j)
      if (k.gt.0) then
            nnz = nnz + 1
            values(nnz) = elresf(j) + brhs_(k)  
            rows  (nnz) = k-1 +1
      end if
  100 continue
!      rows = rows - 1
!     write(*,'(a,i3,a,16i6)') " nel = ", nel,  ", ", rows(1:nnz)
!     write(*,'(a,i3,a,i9,16e15.4)') " nel = ", nel,  ", ", b_, values(1:nnz)
      !call HYPRE_IJVectorSetValues (b_, nnz, rows, elresf, ierr)
#ifdef withHYPRE
     ! write(*,*) "rows =", rows(1:nnz)
     ! write(*,*) "values =", values(1:nnz)
      call HYPRE_IJVectorAddToValues (b_, nnz, rows, values(lowerInd), ierr)
#endif
!

!     write(*,'(a,4i5)') "+++++  nee_ = ", nee_
!    write(*,'(a,4i5)') "+++++",   lm(1:nee_)
!    write(*,'(a,4i5)') "+++++",   rows(1:nnz)
!    write(*,'(a,4f8.4)') "+++++", elresf(1:nee_)
!    write(*,'(a,4f8.4)') "+++++", values(1:nnz)
      return
      end subroutine addrhsHYPRE

      subroutine lerValoresSistemaAlgHYPRE (alhs_ , Ap_, Ai_, rhs_, x_, rows_, &
                                   neq_, nonzerosT_, lowerInd_, upperInd_,  &
                                A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_,  &
                                 myid_, mpi_comm_) 
      implicit none
      integer*4, intent(in)  :: nonzerosT_
      integer*4, intent(in)  :: neq_
      integer*4, intent(in)  :: lowerInd_, upperInd_
      double precision       :: alhs_(nonZerosT_)  
      integer*4              :: Ap_(neq_), Ai_(nonZerosT_)  
      double precision       :: rhs_(neq_), x_(neq_)
      integer*4              :: rows_(neq_)  
      integer*8              :: A_, parcsr_A_, b_, par_b_, u_, par_u_, solver_
      integer*4, intent(in)  :: myid_
      integer*4, intent(in)  :: mpi_comm_

      integer*4     :: nnz, i, j, k, eq, ierr, printLevel
      integer*4     :: local_size, Flower, Fupper
      integer*4     :: cols(90)
      integer*4     :: colsB(90), eqB
      double precision ::  values(90)
     
#ifdef withHYPRE
 !   write(*,*) " em subroutine lerValoresSistemaAlgCSR_HYPRE ( A_, parcsr_A_, b_, par_b_, u_, par_u_, .."
    print*, " atribuindo valores da matrix para o HYPRE_IJMatrix "
    Flower    = lowerInd_+1; Fupper    = upperInd_+1;
   !   HYPRE_IJMatrixRead( <filename>, MPI_COMM_WORLD, HYPRE_PARCSR, &A );
   !   HYPRE_IJMatrixWrite( "matrixHypre.dat", MPI_COMM_WORLD, HYPRE_PARCSR, &A );
    local_size = upperInd_ - lowerInd_ + 1

       !HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD, HYPRE_PARCSR, &b ); 
    print*, " ++++ B atribuindo valores de RHS para o HYPRE_IJVector ", Flower, Fupper
      !call HYPRE_IJVectorSetValues (b_, local_size, rows_, rhs_, ierr )
   !   call HYPRE_IJVectorRead( <filename>, MPI_COMM_WORLD, par_b_, b_ ); 

      print*, " +++++ atribuindo valores de guess para o HYPRE_IJVector "
      !call HYPRE_IJVectorSetValues (u_, local_size, rows_, x_, ierr)
   !   call HYPRE_IJVectorRead(  <filename>,  , MPI_COMM_WORLD, par_u_, u_ ); 

#endif
 end subroutine lerValoresSistemaAlgHYPRE

!=========================================
      subroutine criarVetorHYPRE_A(v_, lowerInd_, upperInd_, mpi_comm_)
      integer*8, intent(out) :: v_
      integer*4, intent(in)  :: lowerInd_, upperInd_
      integer*4, intent(in)  :: mpi_comm_
      integer*4 :: ierr 
      v_ = -9
      !write(*,'(2(a,i0))') "em criarVetorHYPRE, i, v_ = ", v_, lowerInd_, upperInd_
#ifdef withHYPRE
      call HYPRE_IJVectorCreate        (mpi_comm_, lowerInd_, upperInd_, v_, ierr )
      !write(*,*) "call HYPRE_IJVectorSetObjectType,  HYPRE_PARCSR = ", HYPRE_PARCSR
      call HYPRE_IJVectorSetObjectType (v_, HYPRE_PARCSR, ierr)
      !write(*,*) "call HYPRE_IJVectorInitialize,     HYPRE_PARCSR = ", HYPRE_PARCSR
      call HYPRE_IJVectorInitialize    (v_, ierr)
#endif
      !write(*,'(2(a,i0))') "em criarVetorHYPRE, f, v_ = ", v_, lowerInd_, upperInd_
      end subroutine criarVetorHYPRE_A
!=========================================
      subroutine criarVetorHYPRE_B  (estrutSistEq_) !, mpi_comm_) 
       implicit none
       type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      !write(*,*) " em subroutine criarMatrizHYPRE_C, M_ = ", estrutSistEq_%A_HYPRE
       call criarVetorHYPRE_A (estrutSistEq_%B_HYPRE, estrutSistEq_%flower, estrutSistEq_%fupper, mpi_comm)
      end subroutine criarVetorHYPRE_B  
!=========================================
      subroutine criarMatrizHYPRE_C  (estrutSistEq_) !, mpi_comm_) 
       implicit none
       type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      !write(*,*) " em subroutine criarMatrizHYPRE_C, M_ = ", estrutSistEq_%A_HYPRE
       call criarMatrizHYPRE_A (estrutSistEq_, mpi_comm)
      end subroutine criarMatrizHYPRE_C  
!=============================================================================
      subroutine criarMatrizHYPRE_A  (estrutSistEq_, mpi_comm_) 
       implicit none
       integer*4, intent(in)  :: mpi_comm_
       type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
     !write(*,*) " em subroutine criarMatrizHYPRE_A, M_ = ", estrutSistEq_%A_HYPRE, estrutSistEq_%flower, estrutSistEq_%fupper 
       call criarMatrizHYPRE_B (estrutSistEq_%A_HYPRE, estrutSistEq_%flower, estrutSistEq_%fupper, mpi_comm_)
      end subroutine criarMatrizHYPRE_A  
!====================================
      subroutine criarMatrizHYPRE_B(M_, lowerInd_, upperInd_,  mpi_comm_)
      integer * 8, intent(inout)  :: M_
      integer * 4, intent(in)  :: lowerInd_, upperInd_
      integer * 4, intent(in)  :: mpi_comm_
      integer :: ierr 
      !write(*,'(/,3(a,i0))') "em criarMatrizHYPRE_B, i, M_ ", M_, ",  lowerInd_ = ",  lowerInd_, ',upperInd_=',  upperInd_
#ifdef withHYPRE
!     Create the matrix.
!     Note that this is a square matrix, so we indicate the row partition
!     size twice (since number of rows = number of cols)
      call HYPRE_IJMatrixCreate        (mpi_comm, lowerInd_, upperInd_, lowerInd_, upperInd_, M_, ierr )
!     Choose a parallel csr format storage (see the User's Manual)
      !write(*,*) "call HYPRE_IJMatrixSetObjectType,  parcsr_ = ", parcsr_
      !call HYPRE_IJMatrixSetObjectType (M_, parcsr_, ierr)
      !call HYPRE_IJMatrixSetObjectType (M_, parcsr_, ierr)
      call HYPRE_IJMatrixSetObjectType (M_, HYPRE_PARCSR, ierr)
!     call HYPRE_StructMatrixSetSymmetric (M_, 1)
!     Initialize before setting coefficients
      call HYPRE_IJMatrixInitialize    (M_, ierr)
#endif
!      write(*,'(7(a,i0))') "cM+++, ", M_ , ", lowerInd_= ", lowerInd_, ",  upperInd_=" , upperInd_
      if(M_<=0) then
        write(*,*) " em subroutine criarMatrizHYPRE_B, erro, M_ = ", M_
        stop
      end if
       !write(*,'(/,2(a,i0))') "em criarMatrizHYPRE_B, f, M_ ", M_
      end subroutine criarMatrizHYPRE_B
!=============================================================================
      subroutine inicializarMatrizHYPRE(M_)
      integer * 8              :: M_
      integer :: ierr 
!      write(*,'(7(a,i0))') "iM+++, ", M_ 

#ifdef withHYPRE
!     Initialize before setting coefficients
      call HYPRE_IJMatrixInitialize    (M_, ierr)
#endif
!      write(*,'(7(a,i0))') "iM+++, ", M_ 
      if(M_==0) then
        write(*,*) " em subroutine inicializarMatrizHYPRE, erro, M_ = ", M_
        stop
      end if
      end subroutine inicializarMatrizHYPRE
      subroutine inicializarVetorHYPRE(v_)
      integer * 8              :: v_
      integer :: ierr 
!      write(*,'(7(a,i0))') "iv+++, ", v_ 

#ifdef withHYPRE
!     Initialize before setting coefficients
      call HYPRE_IJVectorInitialize( v_, ierr)
#endif
!      write(*,'(7(a,i0))') "iv+++, ", v_ 
      if(v_==0) then
        write(*,*) " em subroutine inicializarVetorHYPRE, erro, v_ = ", v_
        stop
      end if
      end subroutine inicializarVetorHYPRE
!=============================================================================
      subroutine fecharMatriz_HYPRE  (estrutSistEq_) 
       implicit none
       type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
       call fecharMatrizHYPRE_A(estrutSistEq_)
      end subroutine fecharMatriz_HYPRE
!=============================================================================
      subroutine fecharMatrizHYPRE_A  (estrutSistEq_) 
       implicit none
       type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
       call fecharMatrizHYPRE_B(estrutSistEq_%A_HYPRE, estrutSistEq_%parcsr_A )
      end subroutine fecharMatrizHYPRE_A
!=============================================================================
      subroutine fecharMatrizHYPRE_B (A_, parcsr_A_) !
       implicit none
       integer*8, intent(in)  :: A_
       integer*8, intent(in) :: parcsr_A_
       integer*4 :: ierr
#ifdef withHYPRE
!     Assemble after setting the coefficients
       call HYPRE_IJMatrixAssemble( A_, ierr)
!      Get parcsr matrix object
       call HYPRE_IJMatrixGetObject( A_, parcsr_A_, ierr)
#endif
      end subroutine fecharMatrizHYPRE_B
!=============================================================================
!      subroutine fecharVetor_HYPRE (v_, par_v_)
!       implicit none
!       integer*8, intent(in)  :: v_
!       integer*8, intent(inout) :: par_v_
!       integer*4 :: ierr
!       !write(*,'(2(a,i0))') "fecharVetor_HYPRE, i, v_ = ", v_, ", par_v_=",par_v_
!      call fecharVetorHYPRE_A (v_,par_v_)
!       !write(*,'(2(a,i0))') "fecharVetor_HYPRE, f, v_ = ", v_, ", par_v_=",par_v_
!end subroutine fecharVetor_HYPRE
!=============================================================================
!=============================================================================
      subroutine fecharVetorHYPRE_C  (estrutSistEq_) 
       implicit none
       type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
       call fecharVetorHYPRE_B(estrutSistEq_%B_HYPRE)
      end subroutine fecharVetorHYPRE_C
!=============================================================================
      subroutine fecharVetorHYPRE_A (v_, par_v_)
       implicit none
       integer*8, intent(in)  :: v_
       integer*8, intent(inout) :: par_v_
       integer*4 :: ierr
       !write(*,'(2(a,i0))') "fecharVetor_HYPRE_A, i, v_ = ", v_, ", par_v_=",par_v_
#ifdef withHYPRE
!      Assemble after setting the coefficients
       call HYPRE_IJVectorAssemble  (v_, ierr)
       call HYPRE_IJMatrixGetObject( v_, par_v_, ierr)
#endif
       !write(*,'(2(a,i0))') "fecharVetor_HYPRE_A, f, v_ = ", v_, ", par_v_=",par_v_
!      write(*,'(2(a,i0))') "fV+++, v_ ", v_, ", par_v =", par_v_
end subroutine fecharVetorHYPRE_A
!=============================================================================
      subroutine fecharVetorHYPRE_B (v_)
       implicit none
       integer*8, intent(in)  :: v_
       integer*4 :: ierr
       !write(*,'(2(a,i0))') "fecharVetor_HYPRE_B, i, v_ = ", v_!, ", par_v_=",par_v_
#ifdef withHYPRE
!      Assemble after setting the coefficients
       call HYPRE_IJVectorAssemble  (v_, ierr)
       !call HYPRE_IJMatrixGetObject( v_, par_v_, ierr)
#endif
       !write(*,'(2(a,i0))') "fecharVetor_HYPRE_B, f, v_ = ", v_!, ", par_v_=",par_v_
!      write(*,'(2(a,i0))') "fV+++, v_ ", v_, ", par_v =", par_v_
end subroutine fecharVetorHYPRE_B
!=============================================================================
      subroutine destruirVetor_HYPRE   ( v_)
      integer*8 , intent(inout) ::  v_
      integer*4 ::  ierr  
!      write(*,'(a,i0)') "dV+++, v=",  v_ 
#ifdef withHYPRE
         call destruirVetorHYPRE(v_)
#endif
      end subroutine destruirVetor_HYPRE       
!=============================================================================
      subroutine destruirVetorHYPRE   ( v_)
      integer*8 , intent(inout) ::  v_
      integer*4 ::  ierr  
!      write(*,'(a,i0)') "dV+++, v=",  v_ 
#ifdef withHYPRE
         call HYPRE_IJVectorDestroy(v_, ierr)
#endif
      end subroutine destruirVetorHYPRE       
!=============================================================================
      subroutine destruirMatriz_HYPRE   (M_)
      integer*8 , intent(inout) :: M_
      integer*4 ::  ierr  
!       write(*,'(a,i0)') "dM+++, M_= ", M_ 
#ifdef withHYPRE
         call destruirMatrizHYPRE(M_)
#endif
      end subroutine destruirMatriz_HYPRE       
!=============================================================================
      subroutine destruirMatrizHYPRE   (M_)
      integer*8 , intent(inout) :: M_
      integer*4 ::  ierr  
!       write(*,'(a,i0)') "dM+++, M_= ", M_ 
#ifdef withHYPRE
         call HYPRE_IJMatrixDestroy(M_, ierr)
#endif
      end subroutine destruirMatrizHYPRE       
!=============================================================================
      subroutine extrairValoresVetor_HYPRE(v_HYPRE_, lowerInd_, upperInd_, rows_, destino_) 
      integer*8, intent(in) :: v_HYPRE_
      integer*4, intent(in)  :: lowerInd_, upperInd_
      integer*4, intent(in)  :: rows_(:)
      real*8, intent(out)    :: destino_(:)
      integer*4 :: ierr 
!     write(*,'(3(a,i0))') "extV+++, ", v_HYPRE_, ", lowerInd_= ", lowerInd_, ", upperInd_= ", upperInd_
       !destino_=0; 
#ifdef withHYPRE
      call extrairValoresVetorHYPRE (v_HYPRE_, lowerInd_, upperInd_, rows_, destino_)
#endif
      end subroutine extrairValoresVetor_HYPRE

      subroutine extrairValoresVetorHYPRE(v_HYPRE_, lowerInd_, upperInd_, rows_, destino_) 
      integer*8, intent(in) :: v_HYPRE_
      integer*4, intent(in)  :: lowerInd_, upperInd_
      integer*4, intent(in)  :: rows_(:)
      real*8, intent(out)    :: destino_(:)
      integer*4 :: ierr 
!     write(*,'(3(a,i0))') "extV+++, ", v_HYPRE_, ", lowerInd_= ", lowerInd_, ", upperInd_= ", upperInd_
       !destino_=0; 
#ifdef withHYPRE
      call HYPRE_IJVectorGetValues (v_HYPRE_, upperInd_- lowerInd_ + 1, rows_, destino_(lowerInd_), ierr)
#endif
      end subroutine extrairValoresVetorHYPRE

      subroutine atribuirValoresVetor_HYPRE(v_HYPRE_, lowerInd_, upperInd_, rows_, values_) 
      integer*8, intent(in) :: v_HYPRE_
      integer*4, intent(in)  :: lowerInd_, upperInd_
      integer*4, intent(in)  :: rows_(:)
      real*8, intent(in)    :: values_(:)
      integer*4 :: ierr 
!      write(*,'(3(a,i0))') "atrV+++, ", v_HYPRE_, ", lowerInd_= ", lowerInd_, ", upperInd_= ", upperInd_
#ifdef withHYPRE
      call atribuirValoresVetorHYPRE(v_HYPRE_, lowerInd_, upperInd_, rows_, values_)
#endif
      end subroutine atribuirValoresVetor_HYPRE

      subroutine atribuirValoresVetorHYPRE(v_HYPRE_, lowerInd_, upperInd_, rows_, values_) 
      integer*8, intent(in) :: v_HYPRE_
      integer*4, intent(in)  :: lowerInd_, upperInd_
      integer*4, intent(in)  :: rows_(:)
      real*8, intent(in)    :: values_(:)
      integer*4 :: ierr 
!      write(*,'(3(a,i0))') "atrV+++, ", v_HYPRE_, ", lowerInd_= ", lowerInd_, ", upperInd_= ", upperInd_
#ifdef withHYPRE
      call HYPRE_IJVectorSetValues (v_HYPRE_, upperInd_- lowerInd_ + 1, rows_, values_(lowerInd_), ierr)
#endif
      end subroutine atribuirValoresVetorHYPRE

      subroutine atribuirValorVetorHYPRE(v_HYPRE_, i_, valor_) 
      integer*8, intent(in) :: v_HYPRE_
      integer*4, intent(in)  :: i_
      real*8    :: valor_(1)
      integer*4 :: ierr 
      integer*4, parameter :: ntermos=1  
      integer*4 :: rows(nTermos)  
      rows(1)=i_
!      write(*,'(3(a,i0))') "atrV+++, ", v_HYPRE_, ", lo  rInd_= ", lowerInd_, ", upperInd_= ", upperInd_
#ifdef withHYPRE
      call HYPRE_IJVectorSetValues (v_HYPRE_, nTermos, rows, valor_(1) , ierr)
#endif
      end subroutine atribuirValorVetorHYPRE

      subroutine atribuirValorMatrizIJHYPRE(A_HYPRE_, i_, j_, valor_) 
      integer*8, intent(in) :: A_HYPRE_
      integer*4, intent(in)  :: i_, j_
      real*8, intent(in)    :: valor_(1)
      integer*4 :: ierr 
!      write(*,'(3(a,i0))') "atrV+++, ", v_HYPRE_, ", lowerInd_= ", lowerInd_, ", upperInd_= ", upperInd_M
#ifdef withHYPRE
       call HYPRE_IJMatrixSetValues(A_HYPRE_, 1, 1,  i_, j_, valor_, ierr)
#endif
      end subroutine atribuirValorMatrizIJHYPRE

      subroutine adicionarValoresVetorHYPRE(v_HYPRE_, lowerInd_, upperInd_, rows_, destino_) 
      integer*8, intent(in) :: v_HYPRE_
      integer*4, intent(in)  :: lowerInd_, upperInd_
      integer*4, intent(in)  :: rows_(:)
      real*8, intent(in)    :: destino_(:)
      integer*4 :: ierr 
!     write(*,'(3(a,i0))') "addV+++, ", v_HYPRE_, ", lowerInd_= ", lowerInd_, ", upperInd_= ", upperInd_
!     write(*,*) rows_(lowerInd_:upperInd_) 
#ifdef withHYPRE
      call HYPRE_IJVectorAddToValues (v_HYPRE_, upperInd_- lowerInd_ + 1, rows_, destino_(lowerInd_), ierr)
#endif
     !write(*,'(3(a,i0))') "addV+++, ", v_HYPRE_, ", lowerInd_= ", lowerInd_, ", upperInd_= ", upperInd_
      end subroutine adicionarValoresVetorHYPRE

      subroutine adicionarValorVetorHYPRE(v_HYPRE_, i_, valor_) 
      integer*8, intent(in) :: v_HYPRE_
      integer*4, intent(in)  :: i_
      real*8    :: valor_
     ! real*8, intent(inout)   :: valor_
      integer*4 :: ierr 
      integer*4 :: lowerInd=1 , rows(1)
      rows(1)=i_
!      write(*,'(3(a,i0))') "atrV+++, ", v_HYPRE_, ", lowerInd_= ", lowerInd_, ", upperInd_= ", upperInd_
#ifdef withHYPRE
      call HYPRE_IJVectorAddToValues (v_HYPRE_, lowerInd, rows, valor_(lowerInd), ierr)
#endif
      end subroutine adicionarValorVetorHYPRE

subroutine  escreverMatrizHYPRE(A_HYPRE_, matrixFile_)
  implicit none
  integer*8, intent(in) :: A_HYPRE_
  character(len=*) :: matrixFile_
#ifdef withHYPRE
           call HYPRE_IJMatrixPrint(A_HYPRE_, matrixFile_)
#endif
end subroutine  escreverMatrizHYPRE

      subroutine escreverResultadosHYPRE (x_, num_iterations_, final_res_norm_, &
                                  elapsedT_, myId_, print_solution_)
      implicit none
      integer*8, intent(in) :: x_
      integer*4, intent(in) :: num_iterations_
      double precision, intent(in) :: final_res_norm_, elapsedT_
      integer*4, intent(in) :: myId_
      integer*4, intent(in) :: print_solution_ 
!
      integer*8 ::  ierr, ps
      character (LEN=16)  :: SolverStatus="incomplete"

    write(*,*) 'em subroutine escreverResultadosHYPRE'
    SolverStatus="complete"
    if(SolverStatus=="complete" .and. myid_ .eq. 0) then
        print *
        print *, "Final Relative Residual Norm = ", final_res_norm_
        print *, "Iterations                   = ", num_iterations_
        print *, 'Elapsed real time            = ', elapsedT_
        print *
    endif

!     Print the solution
      ps=print_solution_
      ps=0
      if ( ps .ne. 0 ) then
#ifdef withHYPRE
         call HYPRE_IJVectorPrint( x_, "ij.out.x", ierr)
#endif
      endif
      end subroutine escreverResultadosHYPRE       
end module mSolverHypre

