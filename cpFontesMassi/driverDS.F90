!=================================================================================
!
!         programa de elementos finitos em fortran 90 ++
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!
!         novo projeto: 
!         Eduardo Garcia,        bidu@lncc.br [4]
!         Tuane Lopes,           tuane@lncc.br [5]
!
!         LNCC/MCT
!         Petropolis, 09.2015
!
!         propriedades e funcionalidades:
!         0. calculo do potencial nos pontos nodais pelo metodo de galerkin
!           e do fluxo nodal pela tecnica de posprocessamento global (Loula, 1995) 
!         1. resultados identicos ao programa original
!         2. mantem parte das rotinas originais do livro do hughes: 
!             shlt, shlq, shgt, shgq, local, 
!             addrhs, addlhs, kdbc, load, btod, ftod, fact, back
!         3  eliminacao do vetor A estatico e da funcao mpoint
!         4. alocacao dinamica de memoria: stack (automaticas) + heap (allocatable)
!         5. novos nomes, + significativoes e maiores, de variaveis e procedimentos  
!         6. utilizacao "module":
!             em substituicao aos commons: exigencia fortran 90 em conhecer as interfaces dos procedimentos
!             agrupamento de variaveis e procedimentos
!             
!         7. obrigatoriedade de declaracao de variaveis: implicit none
!         8. composto de varios arquivos
!               variaveisGlobais.F90 malha.F90 utilitarios.F90
!               leituraEscrita.F90 mInputReader.F90
!               funcoesDeForma.F90 estruturasDadosSistEq.F90
!               potencial.F90 fluxo.F90
!               solverGaussSkyline.F90 solverHypre.F90 solverPardisoCSR.F90
!               utilSistemaEquacoes.F90 driverDS.F90
!           e varios diretorios
!               fontes, include, bin 
!         9. extensao "F90" dos arquivos com maiuscula sendo importante para uso de opcao de preprocessamento
!        10. 3 opcoes de solvers: 2 diretos (gauss original e pardiso) e 1 iterativo (HYPRE)
!        11. leitura de dados usando palavras chaves (desenvolvido pela DeepSoft)
!        10. utilizacao do comando make no linux
!              porem com possibilidade de unir os arquivos em um unico
!              colocando-os na ordem que aparecem no item 8 
!              (exigencia por existencia de dependencia entre "module"
!=================================================================================
!
  program poisson
!
      use mGlobaisEscalares, only: exec
      use mLeituraEscrita,   only: abrirArquivos, fecharArquivos

      implicit none

      logical  :: comDS = .true.
      character(LEN=20)  :: optSolverP='PardisoEsparso'
      character(LEN=20)  :: optSolverF='HYPREEsparso'
!
!
!-----------------------------------------------------------------------
!
       optSolverP='GaussSkyline'
       optSolverF='GaussSkyline'

#ifdef withPardiso
       optSolverP='PardisoEsparso'
       optSolverF='PardisoEsparso'
#endif

#ifdef withHYPRE
      optSolverP='HYPREEsparso'
      optSolverF='HYPREEsparso'
#endif

#ifdef withTwoPH
      optSolverP='PardisoEsparso'
      optSolverF='HYPREEsparso'
#endif

#ifdef withTwoHP
      optSolverP='HYPREEsparso'
      optSolverF='PardisoEsparso'
#endif

      write(*,*) "optSolverP= ", optSolverP
      write(*,*) "optSolverF= ", optSolverF

      print*, "      call abrirArquivos  ()"
!
     comDS = .false.
     comDS = .true.

      call abrirArquivos  (comDS)

     if(comDS) then
      print*, " call preprocessadorDS ()"; call preprocessamentoDS()  
     else
      print*, " call preprocessador ()";   call preprocessamento()  
     endif
!
      if(exec==1)  then
          print*, "      call processamento  ()"
          call processamento  (optSolverP, optSolverF)
      end if
!
      call fecharArquivos ()

!
end program poisson

!**** new **********************************************************************
      subroutine preprocessamento ()
      use mLeituraEscrita,   only: leituraGeracaoCoordenadas,  leituraCodigosCondContorno, prntel
      use mLeituraEscrita,   only: leituraValoresCondContorno, leituraGeracaoConectividades
      
      use mGlobaisArranjos,  only: etime, title, mat
      use mGlobaisEscalares, only: exec, iprtin, npint
      use mMalha,            only: numnp, numel, nen, nsd
      use mMalha,            only: x, conecNodaisElem
      use mLeituraEscrita,   only: iin, iecho, icoords, echo
      use mLeituraEscrita,   only: prntel
        
      use mPotencial,        only: estrutSistEqP
      use mFluxo,            only: estrutSistEqF
!
      implicit none

!
      logical :: simetria
!
!.... input phase
!
      call echo

      etime = 0.0

      read(iin,1000) title
      if (title(1).eq.'*end') return

      read(iin,'(3i10)') exec,iprtin,nsd
      read(iin,'(4i10)') numnp, numel, nen, npint
      read(iin,'(2i10)') estrutSistEqP%nlvect, estrutSistEqF%nlvect
!
      write(iecho,1000) title 
      write(iecho,3000) exec, iprtin, nsd
      write(iecho,4000) numnp, numel, estrutSistEqP%nlvect, estrutSistEqF%nlvect
!
!.... initialization phase
!
!
!....    set memory pointers for static data arrays,&
!        and call associated input routines 
!
      estrutSistEqP%ndof=1    
      estrutSistEqF%ndof=nsd
!
!.... input coordinate data
!
      allocate(          x(nsd,numnp));              x=0.0
      write(*,*) "call leituraGeracaoCoordenadas"
      call leituraGeracaoCoordenadas(x,nsd,numnp,iin, icoords,iprtin)

      allocate(mat(numel))
      allocate(conecNodaisElem(nen,numel))
      write(*,*) "call leituraGeracaoConectividades"
      call leituraGeracaoConectividades(conecNodaisElem,mat,nen,iin) ! genel
      if (iprtin.eq.0) call prntel(mat,conecNodaisElem,nen,numel,1_4)
!
!.... input boundary condition data and establish equation numbers
!
      allocate(estrutSistEqP%u (estrutSistEqP%ndof,numnp));  estrutSistEqP%u=0.0
      allocate(estrutSistEqP%id(estrutSistEqP%ndof,numnp)); estrutSistEqP%id=0
      write(*,*) "call leituraCodigosCondContorno(idPotencial"
      
      call leituraCodigosCondContorno(estrutSistEqP%id,estrutSistEqP%ndof,numnp,estrutSistEqP%neq, iin, iecho, iprtin)
      if (estrutSistEqP%nlvect.gt.0) then
        allocate(estrutSistEqP%f(estrutSistEqP%nlvect,estrutSistEqP%ndof,numnp))
        estrutSistEqP%f = 0.0
           write(*,*) ' call leituraValoresCondContorno(fPotencial,ndofP,numnp,0,nlvectP,iprtin)'
           call leituraValoresCondContorno(estrutSistEqP%f,estrutSistEqP%ndof,numnp,1_4,estrutSistEqP%nlvect,iprtin)
      end if
!
!.... input nodal force and prescribed kinematic boundary-value data
!
      allocate(estrutSistEqF%u(estrutSistEqF%ndof,numnp));   estrutSistEqF%u=0.0 
      allocate(estrutSistEqF%id(estrutSistEqF%ndof,numnp)); estrutSistEqF%id=0
      write(*,*) "call leituraCodigosCondContorno(idFluxo," 
      call leituraCodigosCondContorno(estrutSistEqF%id, estrutSistEqF%ndof,numnp,estrutSistEqF%neq,iin, iecho, iprtin)

      if (estrutSistEqF%nlvect.gt.0)  then
        allocate(estrutSistEqF%f(estrutSistEqF%nlvect,estrutSistEqF%ndof,numnp))
        estrutSistEqF%f = 0.0
          write(*,*) 'call leituraValoresCondContorno(fFluxo,ndofF,numnp,0,nlvectF,iprtin)'
          call leituraValoresCondContorno(estrutSistEqF%f,estrutSistEqF%ndof,numnp,1_4,estrutSistEqF%nlvect,iprtin)
      end if
!
!.... input element data
!
      write (*,*) "call leituraParamNumericosPropFisica()"
      call leituraParamNumericosPropFisica()
!     

 1000 format(20a4)
 2000 format(4i10,i10,11i10)
 3000 format(5x,&
     ' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,&
     ' execution code  . . . . . . . . . . . . . . (exec  ) = ',i10//5x,&
     '    eq. 0, data check                                   ',   //5x,&
     '    eq. 1, execution                                    ',  //5x,&
     ' input data print code . . . . . . . . . . . (iprtin ) = ',i10//5x,&
     '    eq. 0, print nodal and element input data           ',   //5x,&
     '    eq. 1, do not print nodal and element input data    ',   //5x, &
     ' number of space dimensions  . . . . . . . . (nsd    ) = ',i10)
 4000 format(5x,&
     ' number of nodal points  . . . . . . . . . . (numnp  ) = ',i10//5x,&
     ' number of elements      . . . . . . . . . . (numel  ) = ',i10//5x,&
     ' number of potencial load vectors   . . . . . (nlvectP) = ',i10//5x,&
     ' number of fluxos   load vectors   . . . . . (nlvectF) = ',i10//5x)
!
      end subroutine preprocessamento ! (optSolver_)
      
!
!**** new **********************************************************************
!
      subroutine leituraParamNumericosPropFisica ()
!
      use mGlobaisEscalares, only: iprtin, optSolver
      use mGlobaisEscalares, only: numParElem, nrowsh, numat, npint, nicode
      use mGlobaisArranjos,  only: npar,c, mat, grav
      use mMalha,            only: numel,nsd,numnp,nen,conecNodaisElem
      use mLeituraEscrita,   only: iin, iecho

!
!.... program to read, generate and write element data
!
      implicit none
!
!.... remove above card for single precision operation
!
!
      integer*4 :: m, n, i
      character(len=80) :: formatoLeitura

      nrowsh = 3
      if (nsd==3) nrowsh=nrowsh+1

      allocate(npar(numParElem))
      allocate(grav(3))
!
      formatoLeitura='(16I10)'
      read(iin, formatoLeitura) (npar(i),i=1,numParElem)

      nicode = npint
      if (nicode.eq.0) nicode=nen
!
      numat  = npar( 1)
      allocate(c(6,numat))
!
      write(iecho,1000) numel,numat,nen,npint
!
!      read material properties
!
      do 400 n=1,numat
      if (mod(n,50).eq.1) write(iecho,4000) numat
      read(iin,5000) m,(c(i,m),i=1,3)
      write(iecho,6000) m,(c(i,m),i=1,3)
  400 continue
!
!     constant body forces
!
      read (iin,7000) (grav(i),i=1,3)
      write (iecho,8000) (grav(i),i=1,3)
!
      return
!
 1000 format(//,&
     ' two/three-n o d e    e l e m e n t s ',//,5x,&
     ' number of elements  . . . . . . . . . . . (numel ) = ',i8,//5x,&
     ' number of element material sets . . . . . (numat ) = ',i5,//5x,&
     ' number of element nodes . . . . . . . . . (nen   ) = ',i5,//5x,&
     ' number of integration points. . . . . . . (npint  ) = ',i5)
 2000  format(16i5)
 4000  format(///,&
     ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
     ' number of material sets . . . . . . . . . . (numat ) = ', i5///,&
       2x,'set',4x,'Kx ', 10x,'Ky',10x,'Kz')
 5000  format(i10,5x,5f10.0)
 6000  format(2x,i3,1x,5(1x,1pe14.4))
 7000 format(8f10.0)
 8000 format(///,&
     ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
     ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
     ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
     ' exemplo 3............................ = ',      1pe15.8,//)
 9000  format(i5)
!
      end subroutine leituraParamNumericosPropFisica
 
!**** new **********************************************************************
      subroutine processamento(optSolverP_, optSolverF_)
     
      use mMalha,            only: nen, nsd, numnp, numConexoesPorElem
      use mMalha,            only: x, conecNodaisElem
      use mLeituraEscrita,   only: ignuplotPotencial, ignuplotFluxo, iparaview
      use mLeituraEscrita,   only: prtgnup, escreverArqParaview, escreverArqParaviewIntermed 
!
      use mPotencial,        only: estrutSistEqP
      use mFluxo,            only: estrutSistEqF
!      
      use mSolverHypre,      only: inicializarMPI, finalizarMPI
      use mSolverHypre,      only: myid, num_procs, mpi_comm
!
      use mUtilSistemaEquacoes, only : montarEstrutDadosSistEqAlq
      use mPotencial,           only: montarSistEqAlgPotencial
      use mFluxo,               only: montarSistEqAlgFluxo
  
      implicit none
      character(LEN=*), intent(inout)  :: optSolverP_, optSolverF_
      logical :: firstP=.true., firstF=.true.
!
!.... solution driver program 
!
      real*8 :: t1, t2, t3, t4
      logical :: simetria
      character(LEN=12) :: label, etapa
      integer * 4 :: r, numRepeticoes = 1 

      logical :: escreverSistP=.false., escreverSolP=.false.
      logical :: escreverSistF=.false., escreverSolF=.false.

!      escreverSistP=.true.; escreverSolP=.true.
!      escreverSistF=.true.; escreverSolF=.true.


      if(optSolverP_=='HYPREEsparso'.or.optSolverF_=='HYPREEsparso') then
         call inicializarMPI                 (myid, num_procs, mpi_comm)
      end if
      
        numConexoesPorElem=nen !  usado para criar lista de visinhos dos nOs
      if (optSolverP_=='PardisoEsparso') then
        if(nsd==2)estrutSistEqP%numCoefPorLinha=9    ! usado para criar LMStencilEq
        if(nsd==3)estrutSistEqP%numCoefPorLinha=27
      end if
      write(*,*) " call montarEstrutDadosSistEqAlq(optSolverP_ =", optSolverP_ 
      call montarEstrutDadosSistEqAlq(optSolverP_, estrutSistEqP) ; 

      if (optSolverF_=='PardisoEsparso') then
        if(nsd==2)estrutSistEqF%numCoefPorLinha=18
        if(nsd==3)estrutSistEqF%numCoefPorLinha=81
      end if
      write(*,*) " call montarEstrutDadosSistEqAlq(optSolverF_ =", optSolverF_ 
      call montarEstrutDadosSistEqAlq(optSolverF_, estrutSistEqF) ; 


      do r = 1, numRepeticoes
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULO DA PRESSAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
        call timing(t1)
        if(firstP) then 
          call montarSistEqAlgPotencial(optSolverP_, estrutSistEqP)
      !firstP=.false.
        endif
        call timing(t2)
        if(escreverSistP) call escreverSistema_MTX(optSolverP_, estrutSistEqP, "coefSistEqAlgEsparsoP.mtx"); 
        call timing(t3)
!stop
        label='potencial'
        call solver(optSolverP_, estrutSistEqP,label)        
        call timing(t4)
        if(escreverSolP)call escreverSolSistema_MTX(optSolverP_, estrutSistEqP, "coefSistEqAlgEsparsoP.sol"); 
         estrutSistEqP%brhs = 0.0 ;! stop! ladoB
!
!       call prtgnup(label, x,potencial,nsd,ndofP,numnp,ignuplotPotencial)
!                             2, label, len(label)) !1=por elemento, 2=nodal
!

                              ! nao lineares ou transientes
        write(*,'(a,f12.4,a)') ' tempo de parede (potencial) =',(t4-t3)+(t2-t1), ", segundos"
        write(*,'(a,f12.4,a)') ' .............. montagem matriz   =', t2 - t1, ", segundos"
        write(*,'(a,f12.4,a)') ' .............. solver            =', t4 - t3, ", segundos"
      !stop "em driver"
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULO DO FLUXO  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

        call timing(t1)
    !  if(firstF) then 
         call montarSistEqAlgFluxo (optSolverF_,estrutSistEqF,estrutSistEqP)
      !  firstF=.false.
 !     endif
        call timing(t2)
        if(escreverSistF) call escreverSistema_MTX(optSolverF_, estrutSistEqF, "coefSistEqAlgEsparsoF.mtx")
        call timing(t3)

        label='fluxo';
        call solver(optSolverF_, estrutSistEqF,label)        
        call timing(t4) ; 
        if(escreverSolF) call escreverSolSistema_MTX(optSolverF_, estrutSistEqF, "coefSistEqAlgEsparsoF.sol"); 
        estrutSistEqF%brhs = 0.0 ;! stop! ladoB
        estrutSistEqP%u = 0.0 ! ladoB
!
!     call prtgnup(label, x,fluxo,nsd,ndofF,numnp,ignuplotFluxo)
!
      write(*,'(a,f12.4,a)') ' tempo de parede (fluxo)     =',(t4-t3)+(t2-t1), ", segundos"
      write(*,'(a,f12.4,a)') ' .............. montagem matriz   =', t2 - t1, ", segundos"
      write(*,'(a,f12.4,a)') ' .............. solver            =', t4 - t3, ", segundos"
!
      end do

      if(optSolverP_=='HYPREEsparso'.or.optSolverF_=='HYPREEsparso') then
        call finalizarMPI                 ()
      end if

      return
      contains
!
!**** new **********************************************************************
!
      subroutine solver(optSolver_ , estrutSistEq_, label_) 
!
      use mSolverGaussSkyline, only: solverGaussSkyline
      use mSolverPardiso,      only: solverPardisoEsparso
      use mSolverHypre,        only: solverHYPRE, escreverResultadosHYPRE  
      use mSolverHypre,        only: resolverSistemaAlgHYPRE
      use mSolverHypre,        only: resolverSistemaAlgHYPRE01
      use mSolverHypre,        only: extrairValoresVetor_HYPRE 
  !    use mSolverHypre,        only: destruirSistemaAlgHYPRE, destruirU_HYPRE
      use mUtilSistemaEquacoes,only: btod
!
      use mEstruturasDadosSistEq  !, only : estruturasArmazenamentoSistemaEq

      implicit none
      character(len=*), intent(in) :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*), intent(in) :: label_
      
      integer * 4 ::  num_iterations, printSol = 2
      real*8 ::  final_res_norm, elapsedT, tol, solutionNorm
      integer*4   :: i, ierr
      character(LEN=12) :: etapa
      logical :: simetria=.true.
      etapa='back'; 
      etapa='full'; 

      write(*,*) 'solver ', optSolver_ ,', ',  label_ 
      if (optSolver_=='GaussSkyline') then
          call solverGaussSkyline(estrutSistEq_, simetria, etapa, label_) 
          deallocate(estrutSistEq_%alhs) ! pode ser desnecessario para problemas 
      endif
      if (optSolver_=='PardisoEsparso') then
          simetria=.true.; etapa='full'; 
          write(*,*) "simetria=.true.; etapa='full';"
          call solverPardisoEsparso(estrutSistEq_, simetria, etapa, label_)
          deallocate(estrutSistEq_%alhs) ! pode ser desnecessario para problemas 
      endif

      if(optSolver_=='HYPREEsparso') then
            simetria=.false.;
            estrutSistEq_%initialGuess=>null()
            if(.not.associated(estrutSistEq_%initialGuess)) then
               write(*,'(a)') ', allocate(initialGuess_(neqP)); initialGuessP=0.0 '
               allocate(estrutSistEq_%initialGuess(estrutSistEq_%neq));  estrutSistEq_%initialGuess=0.0
            endif

        estrutSistEq_%solver_id=1
!    if(estrutSistEq_%ndof == 1) then   ! potencial
     if(label_ == "potencial") then   ! potencial
        estrutSistEq_%precond_id=2 ! cg + AMG preconditioner
                            tol = 1.0e-08
     endif
!    if(estrutSistEq_%ndof == nsd) then ! fluxo 
     if(label_ == "fluxo") then   ! fluxo
        estrutSistEq_%precond_id=1 ! cg + jacobi preconditioner
                            tol = 1.0e-06
     endif

            !call resolverSistemaAlgHYPRE01 (estrutSistEq_, tol, num_iterations, final_res_norm)
     !call solverHYPRE    (estrutSistEq_%A_HYPRE, estrutSistEq_%parcsr_A, &
     call resolverSistemaAlgHYPRE  (estrutSistEq_%A_HYPRE, estrutSistEq_%parcsr_A, &
                          estrutSistEq_%b_HYPRE, estrutSistEq_%par_b, &
                          estrutSistEq_%u_HYPRE, estrutSistEq_%par_u, &
                          estrutSistEq_%solver, estrutSistEq_%solver_id, estrutSistEq_%precond_id,&
                          tol, num_iterations, final_res_norm,  &
                          estrutSistEq_%initialGuess, estrutSistEq_%brhs, estrutSistEq_%rows, estrutSistEq_%neq, &
                          myId, mpi_comm)
      printSol = 0
        call escreverResultadosHYPRE (estrutSistEq_%u_HYPRE, num_iterations, final_res_norm, elapsedT, myId, printSol) 

        call extrairValoresVetor_HYPRE(estrutSistEq_%u_HYPRE, 1, estrutSistEq_%neq, estrutSistEq_%rows, estrutSistEq_%brhs)
        do i = 1, estrutSistEq_%neq
             estrutSistEq_%initialGuess(i)=estrutSistEq_%brhs(i)
        end do
      endif

      write(*,*) " valores nos extremos do vetor solucao,  ", label, r
      write(*,'(6e16.8)') estrutSistEq_%brhs(1    :6)
      write(*,'(6e16.8)') estrutSistEq_%brhs(estrutSistEq_%neq-5: estrutSistEq_%neq)
            solutionNorm = 0.0 
            do i = 1, estrutSistEq_%neq
            solutionNorm = solutionNorm + estrutSistEq_%brhs(i)**2
            end do
            solutionNorm = sqrt(solutionNorm)
       write(*,'(a,1e16.8)') "Euclid norms of computed solution: ", solutionNorm


      call btod(estrutSistEq_%id,estrutSistEq_%u, estrutSistEq_%brhs,estrutSistEq_%ndof,numnp)
 

      end subroutine solver
!
!
!**** new **********************************************************************
!
      subroutine escreverSistema_MTX(optSolver_, estrutSistEq_, nomeArq_)
      use mMalha,            only: numnp, numel, nen
      use mSolverPardiso,      only: escreverSistemaAlgCSRemMTX
      use mSolverGaussSkyline, only: escreverSistemaSkylineEmMTX
      use mEstruturasDadosSistEq  !, only : estruturasArmazenamentoSistemaEq

      character(LEN=*), intent(in)  :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*) :: nomeArq_
      integer*4 :: nee

      nee = estrutSistEq_%ndof * nen
            
      if (optSolver_=='GaussSkyline') then
          call escreverSistemaSkylineEmMTX (estrutSistEq_%alhs, estrutSistEq_%brhs, estrutSistEq_%idiag, &
                estrutSistEq_%lm,estrutSistEq_%id, &
                conecNodaisElem,nee,nen, numel,numnp, estrutSistEq_%neq, nomeArq_)
       end if


      if (optSolver_=='PardisoEsparso') then
          call escreverSistemaAlgCSRemMTX    (estrutSistEq_%Alhs, estrutSistEq_%brhs, estrutSistEq_%Ap, estrutSistEq_%Ai, &
                                   estrutSistEq_%nAlhs, estrutSistEq_%neq, nomeArq_)
       end if
      end subroutine escreverSistema_MTX
      

!**** new **********************************************************************
!
      subroutine escreverSolSistema_MTX(optSolver_, estrutSistEq_, nomeArq_)
      use mMalha,            only: numnp, numel, nen
      use mSolverPardiso,      only: escreverSistemaAlgCSRemMTX
      use mSolverGaussSkyline, only: escreverSistemaSkylineEmMTX
      use mEstruturasDadosSistEq  !, only : estruturasArmazenamentoSistemaEq

      character(LEN=*), intent(in)  :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*) :: nomeArq_
      integer*4 :: nee

          call escreverBRHS (estrutSistEq_%brhs, estrutSistEq_%neq, nomeArq_)

      end subroutine escreverSolSistema_MTX
   
      subroutine escreverBRHS (brhs_, neq_, nomeArq_)
      real*8, pointer :: brhs_(:)
      integer :: neq_
      character (LEN=*) :: nomeArq_

      integer :: i

      open (unit=1111, file=trim(nomeArq_))
      write(1111, *) 1, neq_
      do i = 1, neq_
      write(1111, *) i, brhs_(i)
      end do
      close(1111)
      end subroutine 
      
      

  end subroutine processamento
!**** new **********************************************************************
!
   subroutine preprocessamentoDS ()
      use mInputReader,      only: readInputFileDS
      use mInputReader,      only: leituraGeracaoCoordenadasDS
      use mInputReader,      only: leituraCodigosCondContornoDS
      use mInputReader,      only: leituraValoresCondContornoDS
      use mInputReader,      only: leituraGeracaoConectividadesDS
      
      use mGlobaisArranjos,  only: etime, title, mat
      use mGlobaisEscalares, only: exec, iprtin, npint
      use mMalha,            only: numnp, numel, nen, nsd
      use mMalha,            only: x, conecNodaisElem
      use mLeituraEscrita,   only: iin, iecho, icoords, echo
      use mLeituraEscrita,   only: prntel
     
      use mPotencial,        only: estrutSistEqP
      use mFluxo,            only: estrutSistEqF
      
      implicit none

        character(len=50) keyword_name
!
        logical :: simetria

!.... input phase

!        call echo
        call readInputFileDS()
        call readSetupPhaseDS()
        
        etime = 0.0

!.... initialization phase
!
!
!....    set memory pointers for static data arrays,&
!        and call associated input routines
!
      estrutSistEqP%ndof=1
      estrutSistEqF%ndof=nsd
!
!.... input coordinate data
      allocate(          x(nsd,numnp));              x=0.0
      write(*,*) "call leituraGeracaoCoordenadasDS"
      call leituraGeracaoCoordenadasDS(x,nsd,numnp, icoords, iprtin)

!    generation of conectivities
      allocate(mat(numel))
      allocate(conecNodaisElem(nen,numel))
      write(*,*) "call leituraGeracaoConectividadesDS"
      keyword_name = "conectividades_nodais"
      call leituraGeracaoConectividadesDS(keyword_name, conecNodaisElem,mat,nen)
      
      if (iprtin.eq.0) call prntel(mat,conecNodaisElem,nen,numel,1_4)
!
!.... input boundary condition data and establish equation numbers

      allocate(estrutSistEqP%u (estrutSistEqP%ndof,numnp));  estrutSistEqP%u=0.0
      allocate(estrutSistEqP%id(estrutSistEqP%ndof,numnp)); estrutSistEqP%id=0
        write(*,*) "call leituraCodigosCondContornoDS(idPotencial"
        keyword_name = "codigos_cond_contorno_potencial"
        call leituraCodigosCondContornoDS(keyword_name, estrutSistEqP%id, estrutSistEqP%ndof, numnp, &
                                          estrutSistEqP%neq, iecho, iprtin)
                                          

      if (estrutSistEqP%nlvect.gt.0) then
        allocate(estrutSistEqP%f(estrutSistEqP%nlvect,estrutSistEqP%ndof,numnp))
        estrutSistEqP%f = 0.0
        write(*,*) "call leituraValoresCondContornoDS(fPotencial,ndofP,numnp,0,nlvectP,iprtin)"
        keyword_name = "valores_cond_contorno_potencial"
        call leituraValoresCondContornoDS(keyword_name, estrutSistEqP%f, estrutSistEqP%ndof, numnp, 1_4,&
                                         estrutSistEqP%nlvect, iprtin)
      end if

!.... input nodal force and prescribed kinematic boundary-value data
      allocate(estrutSistEqF%u (estrutSistEqF%ndof,numnp));  estrutSistEqF%u=0.0
      allocate(estrutSistEqF%id(estrutSistEqF%ndof,numnp)); estrutSistEqF%id=0
       write(*,*) "call leituraCodigosCondContornoDS(idFluxo,"
      keyword_name = "codigos_cond_contorno_fluxo"
       call leituraCodigosCondContornoDS(keyword_name, estrutSistEqF%id, estrutSistEqF%ndof, numnp, &
                                          estrutSistEqF%neq, iecho, iprtin)

      if (estrutSistEqF%nlvect.gt.0)  then
        allocate(estrutSistEqF%f(estrutSistEqF%nlvect,estrutSistEqF%ndof,numnp))
        estrutSistEqF%f = 0.0
        write(*,*) 'call leituraValoresCondContornoDS(fFluxo,ndofF,numnp,0,nlvectF,iprtin)'
        keyword_name = "valores_cond_contorno_fluxo"
        call leituraValoresCondContornoDS(keyword_name, estrutSistEqF%f, estrutSistEqF%ndof, numnp, 1_4, &
                                           estrutSistEqF%nlvect, iprtin)
      end if
!
!.... input element data
!
     write (*,*) "call leituraParamNumericosPropFisicaDS()"
     call leituraParamNumericosPropFisicaDS()

 1000 format(20a4)
 2000 format(4i10,i10,11i10)
 3000 format(5x,&
     ' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,&
     ' execution code  . . . . . . . . . . . . . . (exec  ) = ',i10//5x,&
     '    eq. 0, data check                                   ',   //5x,&
     '    eq. 1, execution                                    ',  //5x,&
     ' input data print code . . . . . . . . . . . (iprtin ) = ',i10//5x,&
     '    eq. 0, print nodal and element input data           ',   //5x,&
     '    eq. 1, do not print nodal and element input data    ',   //5x, &
     ' number of space dimensions  . . . . . . . . (nsd    ) = ',i10)
 4000 format(5x,&
     ' number of nodal points  . . . . . . . . . . (numnp  ) = ',i10//5x,&
     ' number of elements      . . . . . . . . . . (numel  ) = ',i10//5x,&
     ' number of potencial load vectors  . . . . . (nlvectP) = ',i10//5x,&
     ' number of fluxos   load vectors   . . . . . (nlvectF) = ',i10//5x)

    end subroutine preprocessamentoDS
!
!**** new **********************************************************************
!
    !> Efetua a leitura completa dos dados da etapa de setup.
    subroutine readSetupPhaseDS
        use mInputReader,      only:readStringKeywordValue, readIntegerKeywordValue, readIntegerKeywordValue
        use mGlobaisArranjos,  only: title
        use mGlobaisEscalares, only: exec, iprtin, npint
        use mMalha,            only: nsd, numnp, numel, nen
        use mPotencial,        only: estrutSistEqP
        use mFluxo,            only: estrutSistEqF

        implicit none
        character(len=50) keyword_name
        character*4  default_title_value(20)

        !Reads the title
        keyword_name = "title"
        default_title_value = "unknown title"
        call readStringKeywordValue(keyword_name, title, default_title_value)
        !print*, "title: ", title
        !Reads exec
        keyword_name = "exec"
        call readIntegerKeywordValue(keyword_name, exec, 0_4)
        !print*, "exec: ", exec
        !Reads iprtin
        keyword_name = "iprtin"
        call readIntegerKeywordValue(keyword_name, iprtin, 0_4)
        !print*, "iprtin: ", iprtin
        !Reads nsd
        keyword_name = "nsd"
        call readIntegerKeywordValue(keyword_name, nsd, 0_4)
        !print*, "nsd: ", nsd
        !Reads numnp
        keyword_name = "numnp"
        call readIntegerKeywordValue(keyword_name, numnp, 0_4)
        !print*, "numnp: ", numnp
        !Reads numel
        keyword_name = "numel"
        call readIntegerKeywordValue(keyword_name, numel, 0_4)
        !print*, "numel: ", numel
        !Reads nen
        keyword_name = "nen"
        call readIntegerKeywordValue(keyword_name, nen, 0_4)
        !print*, "nen: ", nen
        !Reads npint
        keyword_name = "npint"
        call readIntegerKeywordValue(keyword_name, npint, 0_4)
        !print*, "npint: ", npint
        !Reads nlvectP
        keyword_name = "nlvectP"
        call readIntegerKeywordValue(keyword_name, estrutSistEqP%nlvect, 0_4)
        !print*, "nlvectP: ", nlvectP
        !Reads nlvectF
        keyword_name = "nlvectF"
        call readIntegerKeywordValue(keyword_name, estrutSistEqF%nlvect, 0_4)
        !print*, "nlvectF: ", nlvectF
        return
    end subroutine readSetupPhaseDS 
!
!**** new **********************************************************************
!
    !> Faz a leitura de constant body forces.
    !> @param keyword_name  Keyword especifica para constant body forces.
    subroutine leituraParamNumericosPropFisicaDS()
        use mGlobaisEscalares, only: iprtin, optSolver
        use mGlobaisEscalares, only: numParElem, nrowsh, numat, npint, nicode
        use mGlobaisArranjos,  only: npar,c, mat, grav
        use mMalha,            only: numel,nsd,numnp,nen,conecNodaisElem
        use mLeituraEscrita,   only: iecho
        use mInputReader,      only: findKeyword, file_lines


        !
        !.... program to read, generate and write element data
        !
        implicit none
        !
        !.... remove above card for single precision operation
        !
        !
        integer*4 :: m, n, i, keyword_line, nLinhaArqInput
        character(len=80) :: formatoLeitura
        character(len=50) keyword_name

        keyword_name = "nummat"
        keyword_line = findKeyword(keyword_name)
        nLinhaArqInput = keyword_line
        if (keyword_line.eq.-1) return

        nrowsh = 3
        if (nsd==3) nrowsh=nrowsh+1

        allocate(npar(numParElem))
        allocate(grav(3))
        !
        formatoLeitura='(16I10)'
!        read(iin, formatoLeitura) (npar(i),i=1,numParElem)
        !print *, " BD em driver, nLinhaArqInput=", nLinhaArqInput; stop
        read(file_lines(nLinhaArqInput:), formatoLeitura) (npar(i),i=1,numParElem)
        nLinhaArqInput = nLinhaArqInput + 1

        nicode = npint
        if (nicode.eq.0) nicode=nen
        !
        numat  = npar( 1)
        allocate(c(6,numat))
        !
        write(iecho,1000) numel,numat,nen,npint
        !
        !      read material properties
        !
        keyword_name = "prop_fisica_meio"
        keyword_line = findKeyword(keyword_name)
        nLinhaArqInput = keyword_line
        do 400 n=1,numat
        if (mod(n,50).eq.1) write(iecho,4000) numat
!        read(iin,5000) m,(c(i,m),i=1,3)
        read(file_lines(nLinhaArqInput:),5000) m,(c(i,m),i=1,3)
        nLinhaArqInput = nLinhaArqInput + 1
        write(iecho,6000) m,(c(i,m),i=1,3)
        400 continue
        !
        !     constant body forces
        !
        keyword_name = "grav"
        keyword_line = findKeyword(keyword_name)
        nLinhaArqInput = keyword_line
        read (file_lines(nLinhaArqInput:),7000) (grav(i),i=1,3)
        nLinhaArqInput = nLinhaArqInput + 1
        write (iecho,8000) (grav(i),i=1,3)

        return

        1000 format(//,&
        ' two/three-n o d e    e l e m e n t s ',//,5x,&
        ' number of elements  . . . . . . . . . . . (numel ) = ',i8,//5x,&
        ' number of element material sets . . . . . (numat ) = ',i5,//5x,&
        ' number of element nodes . . . . . . . . . (nen   ) = ',i5,//5x,&
        ' number of integration points. . . . . . . (npint  ) = ',i5)
        2000  format(16i5)
        4000  format(///,&
        ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
        ' number of material sets . . . . . . . . . . (numat ) = ', i5///,&
        2x,'set',4x,'Kx ', 10x,'Ky',10x,'Kz')
        5000  format(i10,5x,5f10.0)
        6000  format(2x,i3,1x,5(1x,1pe14.4))
        7000 format(8f10.0)
        8000 format(///,&
        ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
        ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
        ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
        ' exemplo 3............................ = ',      1pe15.8,//)
        9000  format(i5)
        !
    end subroutine leituraParamNumericosPropFisicaDS

