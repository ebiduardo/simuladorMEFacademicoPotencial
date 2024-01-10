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
   
      write(*,*) "exec =", exec
      if(exec==1)  then
          print*, "      call processamento  ()"
          call processamento  ()
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
!.... input phase
!
      call echo

      etime = 0.0

      read(iin,1000) title
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
      use mGlobaisEscalares, only: iprtin
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
      subroutine processamento()
     
      use mMalha,            only: nen, nsd, numnp, numConexoesPorElem
      use mMalha,            only: x, conecNodaisElem, numel
      use mLeituraEscrita,   only: ignuplotPotencial, ignuplotFluxo, iparaviewPotencial, iparaviewFluxo
      !use mLeituraEscrita,   only: prtgnup, escreverArquivosSaida
      use mLeituraEscrita,   only: prtgnup, escreverArquivosSaida_Potencial

!
      use mPotencial,        only: estrutSistEqP
      use mFluxo,            only: estrutSistEqF
!      
      use mSolverHypre,      only: inicializarMPI, finalizarMPI
      use mSolverHypre,      only: myid, num_procs, mpi_comm
      use mSolverGaussSkyline, only: escreverSistemaSkylineEmMTX
!
      use mUtilSistemaEquacoes, only: montarEstrutDadosSistEqAlq
      use mPotencial,           only: montarSistEqAlgPotencial
      use mFluxo,               only: montarSistEqAlgFluxo
  
      implicit none

      logical :: firstP=.true., firstF=.true.
!
!.... solution driver program 
!
      real*8 :: t1, t2, t3, t4
      character(LEN=12) :: label, etapa
      integer * 4 :: r, numRepeticoes = 1 

      logical :: escreverSistP=.false., escreverSolP=.false.
      logical :: escreverSistF=.false., escreverSolF=.false.
      logical :: escreverSaidaVTK=.false., escreverSaidaVTK_F=.false.
      character(LEN=40) ::  nomeArqSist
      integer :: passoTempo
      real*8 :: tempo
      escreverSistF=.false.; escreverSolF=.false.
      escreverSistF=.true. ; escreverSolF=.true.
      escreverSistP=.false.; escreverSolP=.false.
      escreverSistP=.true. ; escreverSolP=.true.
      escreverSaidaVTK=.false.
      escreverSaidaVTK=.true.
      if(estrutSistEqP%optSolver(1:12)=='HYPREEsparso'.or.estrutSistEqF%optSolver(1:12)=='HYPREEsparso') then
         call inicializarMPI                 (myid, num_procs, mpi_comm)
      end if
      numConexoesPorElem=nen !  usado para criar lista de visinhos dos nOs
      !Jan22if (estrutSistEqP%optSolver(1:14)=='PardisoEsparso') then
        if(nsd==2)estrutSistEqP%numCoefPorLinha=9    ! usado para criar LMStencilEq
        if(nsd==3)estrutSistEqP%numCoefPorLinha=27
      !Jan22end if

      write(*,*) "em processamento, numCoefPorLinhaP = ", estrutSistEqP%numCoefPorLinha
      write(*,*) " call montarEstrutDadosSistEqAlq(optSolverP_ =", estrutSistEqP%optSolver 
      call montarEstrutDadosSistEqAlq(estrutSistEqP) ; 

   !   if (estrutSistEqF%optSolver(1:14)=='PardisoEsparso') then
   !     if(nsd==2)estrutSistEqF%numCoefPorLinha=18
   !     if(nsd==3)estrutSistEqF%numCoefPorLinha=81
   !   end if
   !   write(*,*) "em processamento, numCoefPorLinhaF = ", estrutSistEqF%numCoefPorLinha
      write(*,*) " call montarEstrutDadosSistEqAlq(optSolverF_ =", estrutSistEqF%optSolver 
      call montarEstrutDadosSistEqAlq(estrutSistEqF) ; 


      do r = 1, numRepeticoes
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULO DA PRESSAO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     
        call timing(t1)
        if(firstP) then 
          call montarSistEqAlgPotencial(estrutSistEqP)
      !firstP=.false.
        endif
        call timing(t2)
                                   nomeArqSist="coefSistEqAlgEsparsoP_NSYM"
        if(estrutSistEqP%simetria) nomeArqSist="coefSistEqAlgEsparsoP_SYM"
        if(escreverSistP) call escreverSistema_MTX(estrutSistEqP, nomeArqSist); 
        call timing(t3)


!       write(*,'(a,20f10.5)') 'alhs =', estrutSistEqP%alhs
!       write(*,'(a,20f10.5)') 'brhs =', estrutSistEqP%brhs
        label='potencial'
        !write(*,'(a,21E14.5)')  "estrutSistEqP%brhs =",  estrutSistEqP%brhs; 
        !write(*,'(a,21E14.5)', advance="NO")  "estrutSistEqP%alhs =" 
        !write(*,'(21E14.5)')   estrutSistEqP%alhs; 
        call solver(estrutSistEqP,label)        
        !write(*,'(a,21E14.5)')  "estrutSistEqP%u    =",  estrutSistEqP%u   ; 
        !stop 1001 
        call timing(t4)
!                                  nomeArqSist="coefSistEqAlgEsparsoP_NSYM"
!       if(estrutSistEqP%simetria) nomeArqSist="coefSistEqAlgEsparsoP"
        if(escreverSolP)call escreverSolSistema_MTX(estrutSistEqP, nomeArqSist); 
         estrutSistEqP%brhs = 0.0 ;! stop! ladoB
!
!       call prtgnup(label, x,potencial,nsd,ndofP,numnp,ignuplotPotencial)
!                             2, label, len(label)) !1=por elemento, 2=nodal
        !if(escreverSaidaVTK) call escreverArquivosSaida(estrutSistEqP, estrutSistEqF)
!

                              ! nao lineares ou transientes
        write(*,'(a,f12.4,a)') ' tempo de parede (potencial) =',(t4-t3)+(t2-t1), ", segundos"
        write(*,'(a,f12.4,a)') ' .............. montagem matriz   =', t2 - t1, ", segundos"
        write(*,'(a,f12.4,a)') ' .............. solver            =', t4 - t3, ", segundos"

        passoTempo=0; tempo=0;
        if(escreverSaidaVTK) call escreverArquivosSaida_Potencial(estrutSistEqP, passoTempo, tempo)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CALCULO DO FLUXO  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call timing(t1)
    !  if(firstF) then 
         call montarSistEqAlgFluxo (estrutSistEqF,estrutSistEqP)
      !  firstF=.false.
 !     endif
        call timing(t2)
                                   nomeArqSist="coefSistEqAlgEsparsoF_NSYM"
        if(estrutSistEqF%simetria) nomeArqSist="coefSistEqAlgEsparsoF_SYM"
        if(escreverSistF) call escreverSistema_MTX(estrutSistEqF, nomeArqSist); 
        call timing(t3)
        label='fluxo';
        call solver(estrutSistEqF,label)        
        call timing(t4) ; 
        if(escreverSolF) call escreverSolSistema_MTX(estrutSistEqF, nomeArqSist); 

        !if(escreverSaidaVTK) call escreverArquivosSaida(estrutSistEqP, estrutSistEqF)

        estrutSistEqP%brhs = 0.0 ! 
        estrutSistEqP%u    = 0.0 ! 

!     call prtgnup(label, x,fluxo,nsd,ndofF,numnp,ignuplotFluxo)
!
      write(*,'(a,f12.4,a)') ' tempo de parede (fluxo)     =',(t4-t3)+(t2-t1), ", segundos"
      write(*,'(a,f12.4,a)') ' .............. montagem matriz   =', t2 - t1, ", segundos"
      write(*,'(a,f12.4,a)') ' .............. solver            =', t4 - t3, ", segundos"
!
      end do

      if(estrutSistEqP%optSolver(1:12)=='HYPREEsparso'.or.estrutSistEqF%optSolver(1:12)=='HYPREEsparso') then
        call finalizarMPI                 ()
      end if

      return
      contains
!
!**** new **********************************************************************
!
      subroutine solver(estrutSistEq_, label_) 
!
      use mSolverGaussSkyline, only: solverGaussSkyline
      use mSolverPardiso,      only: solverPardisoEsparso
      use mSolverHypre,        only: escreverResultadosHYPRE  
      use mSolverHypre,        only: resolverSistemaAlgHYPRE
      use mSolverHypre,        only: extrairValoresVetorHYPRE 
      use mUtilSistemaEquacoes,only: btod
!
      use mEstruturasDadosSistEq  !, only : estruturasArmazenamentoSistemaEq

      implicit none
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*), intent(in) :: label_
      
      integer * 4 ::  num_iterations, printSol = 2
      real*8 ::  final_res_norm, tol
      real*8 ::  elapsedT
      real*8 :: tA, tB
      integer*4   :: ierr
      integer*4 :: i
      character(LEN=12) :: etapa
       elapsedT =0.0

     !write(*,'(4a)') 'solver ', trim(estrutSistEq_%optSolver) ,', ',  trim(label_)
     ! stop

      if (estrutSistEq_%optSolver(1:12)=='GaussSkyline') then
          etapa='full'; 
          call solverGaussSkyline(estrutSistEq_, etapa, label_) 
      !stop 2222 
          deallocate(estrutSistEq_%alhs) ! pode ser desnecessario para problemas 
      endif
      if (estrutSistEq_%optSolver(1:14)=='PardisoEsparso') then
          etapa='full'; 
          call solverPardisoEsparso(estrutSistEq_, etapa, label_)         
          deallocate(estrutSistEq_%alhs) ! pode ser desnecessario para problemas 
      endif

      if(estrutSistEq_%optSolver(1:12)=='HYPREEsparso') then
            estrutSistEq_%initialGuess=>null()
            if(.not.associated(estrutSistEq_%initialGuess)) then
               !write(*,'(a)') ', allocate(initialGuess_(neqP)); initialGuessP=0.0 '
               allocate(estrutSistEq_%initialGuess(estrutSistEq_%neq));  estrutSistEq_%initialGuess=0.0
            endif

        estrutSistEq_%solver_id=1
!    if(estrutSistEq_%ndof == 1) then   ! potencial
     if(label_ == "potencial") then   ! potencial
        estrutSistEq_%precond_id=2
 ! cg + AMG preconditioner
        estrutSistEq_%tol = 1.0e-08
     endif
!    if(estrutSistEq_%ndof == nsd) then ! fluxo 
     if(label_ == "fluxo") then   ! fluxo
        estrutSistEq_%precond_id=1
 ! cg + jacobi preconditioner
        estrutSistEq_%tol = 1.0e-08
     endif

            !call resolverSistemaAlgHYPRE01 (estrutSistEq_, tol, num_iterations, final_res_norm)
     !call solverHYPRE    (estrutSistEq_%A_HYPRE, estrutSistEq_%parcsr_A, &
     call timing(tA)
     call resolverSistemaAlgHYPRE  (estrutSistEq_)
!     call resolverSistemaAlgHYPRE  (estrutSistEq_%A_HYPRE, estrutSistEq_%parcsr_A, &
!                          estrutSistEq_%b_HYPRE, estrutSistEq_%par_b, &
!                          estrutSistEq_%u_HYPRE, estrutSistEq_%par_u, &
!                          estrutSistEq_%solver, estrutSistEq_%solver_id, estrutSistEq_%precond_id,&
!                          tol, num_iterations, final_res_norm,  &
!                          estrutSistEq_%initialGuess, estrutSistEq_%brhs, estrutSistEq_%rows, estrutSistEq_%neq, &
!                          myId, mpi_comm)
     call timing(tB)
      printSol = 0
      elapsedT=elapsedT+tB-tA; 
        !call escreverResultadosHYPRE (estrutSistEq_%u_HYPRE, estrutSistEq_%num_iterations, estrutSistEq_%final_res_norm, elapsedT, myId, printSol) 
        call extrairValoresVetorHYPRE(estrutSistEq_%u_HYPRE, 1, estrutSistEq_%neq, estrutSistEq_%rows, estrutSistEq_%brhs)
        do i = 1, estrutSistEq_%neq
             estrutSistEq_%initialGuess(i)=estrutSistEq_%brhs(i)
        end do
      endif

      call btod(estrutSistEq_%id,estrutSistEq_%u, estrutSistEq_%brhs,estrutSistEq_%ndof,numnp)
      ! call escreverSolucaoB(estrutSistEq_) 
      if (label_ == "potencial") then
          call escreverSolucaoP(estrutSistEq_) 
      else
          call escreverSolucaoF(estrutSistEq_) 
      end if
      end subroutine solver

  subroutine escreverSolucaoP(estrutSistEq_)
      use mMalha,            only: numnp
      use mEstruturasDadosSistEq  , only : estruturasArmazenamentoSistemaEq
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      integer*4 :: i, n=6+1
      real*8    :: solutionNorm

      write(*,'(2a,i3)') " valores nos extremos do vetor solucao completo,  ", trim(label), r
      write(*,'(9e12.4)') estrutSistEq_%u(1, 1    :n)
      write(*,'(9e12.4)') estrutSistEq_%u(1, numnp-(n-1): numnp)
      solutionNorm = 0.0 
        do i = 1, numnp
           solutionNorm = solutionNorm + estrutSistEq_%u(1, i)**2
        end do
      solutionNorm = sqrt(solutionNorm)
      write(*,'(a,1e15.8)') "Euclid norms of computed solution: ", solutionNorm
  end subroutine escreverSolucaoP

  subroutine escreverSolucaoF(estrutSistEq_)
      use mMalha,            only: numnp, nsd
      use mEstruturasDadosSistEq  , only : estruturasArmazenamentoSistemaEq
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      integer*4 :: i, j, ndof
      real*8    :: solutionNorm
      ndof = estrutSistEq_%ndof
      write(*,'(2a,i3)') " valores nos extremos do vetor solucao completo,  ", trim(label), r
      write(*,'(9e12.4)') estrutSistEq_%u(1:ndof, 1    :3)
      write(*,'(9e12.4)') estrutSistEq_%u(1:ndof, numnp-2: numnp)
      solutionNorm = 0.0 
        do j = 1, nsd
         do i = 1, numnp
           solutionNorm = solutionNorm + estrutSistEq_%u(j, i)**2
         end do
        end do
      solutionNorm = sqrt(solutionNorm)
      write(*,'(a,1e15.8)') "Euclid norms of computed solution: ", solutionNorm
  end subroutine escreverSolucaoF

  subroutine escreverSolucaoB(estrutSistEq_)
      use mEstruturasDadosSistEq  , only : estruturasArmazenamentoSistemaEq
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      integer*4 :: i
      real*8    :: solutionNorm
      write(*,*) " valores nos extremos do vetor solucao BRHS,  ", label, r
      write(*,'(6e12.4)') estrutSistEq_%brhs(1    :6)
      write(*,'(6e12.4)') estrutSistEq_%brhs(estrutSistEq_%neq-5: estrutSistEq_%neq)
      solutionNorm = 0.0 
        do i = 1, estrutSistEq_%neq
           solutionNorm = solutionNorm + estrutSistEq_%brhs(i)**2
        end do
       solutionNorm = sqrt(solutionNorm)
       write(*,'(a,1e15.8)') "Euclid norms of computed solution: ", solutionNorm
  end subroutine escreverSolucaoB
!
!
!**** new **********************************************************************
!
      subroutine escreverSistema_MTX(estrutSistEq_, nomeArq_)
      use mMalha,              only: numnp, numel, nen
      use mSolverPardiso,      only: escreverSistemaAlgCSRemMTX
      use mSolverGaussSkyline, only: escreverSistemaSkylineEmMTX
      use mEstruturasDadosSistEq  !, only : estruturasArmazenamentoSistemaEq

      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*) :: nomeArq_
      character(len=80) :: nomeArqCompleto
      integer*4 :: nee

!      write(*,'(a)') "subroutine escreverSistema_MTX(estrutSistEq_, nomeArq_)"
!      write(*,'(2a)') "optSolver_,", estrutSistEq_%optSolver

      nee = estrutSistEq_%ndof * nen
            
      if (estrutSistEq_%optSolver(1:12)=='GaussSkyline') then
     !     nomeArq_="coefSistEqAlgEsparsoP_SKYLINE_SYM.mtx"
     !     nomeArq_="coefSistEqAlgEsparsoP_SKYLINE_SYM"
          nomeArqCompleto=trim(nomeArq_)//"_SKYLINE.mtx"
     !     if(.not. estrutSistEq_%simetria) nomeArq_="coefSistEqAlgEsparsoP_SKYLINE_NSYM.mtx"
          if(.not. estrutSistEq_%simetria) nomeArqCompleto=trim(nomeArq_)//"_SKYLINE.mtx"
     !    nomeArqCompleto=trim(nomeArq_)//".mtx"
          call escreverSistemaSkylineEmMTX (estrutSistEq_%alhs, estrutSistEq_%brhs, &
                estrutSistEq_%idiag, estrutSistEq_%lm,estrutSistEq_%id, &
                conecNodaisElem,nee,nen, numel,numnp, estrutSistEq_%neq, nomeArqCompleto)
       end if

      if (estrutSistEq_%optSolver(1:14)=='PardisoEsparso') then
          !nomeArq_="coefSistEqAlgEsparsoP_CSR_SYM.mtx"
          nomeArqCompleto=trim(nomeArq_)//"_CSR.mtx"
          !if(.not. estrutSistEq_%simetria) nomeArq_=trim(nomeArq_)//"_CSR_NSYM.mtx"
          if(.not. estrutSistEq_%simetria) nomeArqCompleto=trim(nomeArq_)//"_CSR.mtx"
       !   nomeArqCompleto=trim(nomeArq_)//".mtx"
          call escreverSistemaAlgCSRemMTX    (estrutSistEq_%Alhs, estrutSistEq_%brhs, &
              estrutSistEq_%Ap, estrutSistEq_%Ai, estrutSistEq_%nAlhs, estrutSistEq_%neq, nomeArqCompleto)
       end if
      end subroutine escreverSistema_MTX
      

!**** new **********************************************************************
!
      subroutine escreverSolSistema_MTX(estrutSistEq_, nomeArq_)
      use mMalha,            only: numnp, numel, nen
      use mSolverPardiso,      only: escreverSistemaAlgCSRemMTX
      use mSolverGaussSkyline, only: escreverSistemaSkylineEmMTX
      use mEstruturasDadosSistEq  !, only : estruturasArmazenamentoSistemaEq

      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
      character(len=*) :: nomeArq_
      character(len=80) :: nomeArqCompleto
      integer*4 :: nee

        nomeArqCompleto=trim(nomeArq_)//".sol"
        call escreverBRHS (estrutSistEq_%brhs, estrutSistEq_%neq, nomeArqCompleto)

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
      end subroutine escreverBRHS


  end subroutine processamento
!**** new **********************************************************************
!
   subroutine preprocessamentoDS ()
      use mInputReader,      only: readInputFileDS
      use mInputReader,      only: leituraGeracaoCoordenadasDS
      use mInputReader,      only: leituraCodigosCondContornoDS
      use mInputReader,      only: leituraCodigosCondContornoDSF
      use mInputReader,      only: leituraValoresCondContornoDS
      use mInputReader,      only: leituraGeracaoConectividadesDS
      
      use mGlobaisArranjos,  only: etime, title, mat, listaSolverDisponivel
      use mGlobaisEscalares, only: exec, iprtin, npint
      use mMalha,            only: numnp, numel, nen, nsd
      use mMalha,            only: x, conecNodaisElem
      use mLeituraEscrita,   only: iin, iecho, icoords, echo
      use mLeituraEscrita,   only: prntel
     
      use mPotencial,        only: estrutSistEqP
      use mFluxo,            only: estrutSistEqF
      
      implicit none

        character(len=50) keyword_name
        integer :: ierr
!
!.... input phase
!
!     call echo
      call readInputFileDS(iin)
      call readSetupPhaseDS()
        
      call identificaSolversDisponiveis(listaSolverDisponivel)
      call verificarSolver(estrutSistEqP%optSolver, listaSolverDisponivel)
      call verificarSolver(estrutSistEqF%optSolver, listaSolverDisponivel) 

      write(*,*) "optSolverP= ", estrutSistEqP%optSolver
      write(*,*) "optSolverF= ", estrutSistEqF%optSolver
      !stop
      estrutSistEqF%simetria=.true.
      if(estrutSistEqP%optSolver=="GaussSkyline") estrutSistEqP%simetria=.true.
      if(estrutSistEqF%optSolver=="GaussSkyline") estrutSistEqF%simetria=.true.
      if(estrutSistEqP%optSolver=="PardisoEsparso") estrutSistEqP%simetria=.true.
      if(estrutSistEqF%optSolver=="PardisoEsparso") estrutSistEqF%simetria=.true.
      estrutSistEqP%simetria=.false.
      estrutSistEqP%simetria=.true.
      if(estrutSistEqP%optSolver=="HYPREEsparso")   estrutSistEqP%simetria=.false.
      if(estrutSistEqF%optSolver=="HYPREEsparso")   estrutSistEqF%simetria=.false.
!     estrutSistEqP%simetria=.false.
!        
      etime = 0.0
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
      allocate(          x(nsd,numnp));              x=0.0
      write(*,*) "call leituraGeracaoCoordenadasDS"
      !biduOut2018call leituraGeracaoCoordenadasDS(x,nsd,numnp, icoords, iprtin)
      call leituraGeracaoCoordenadasDS(x,nsd,numnp, iprtin, icoords)

!    generation of conectivities
      allocate(mat(numel))
      allocate(conecNodaisElem(nen,numel))
      write(*,*) "call leituraGeracaoConectividadesDS"
      keyword_name = "conectividades_nodais"
      call leituraGeracaoConectividadesDS(keyword_name, conecNodaisElem,mat,nen, ierr)
      
      write(*,*) "nen = ", nen; !stop

    !  if (iprtin.eq.0)
      call prntel(mat,conecNodaisElem,nen,numel,1_4)
!
!.... input boundary condition data and establish equation numbers

      allocate(estrutSistEqP%u (estrutSistEqP%ndof,numnp));  estrutSistEqP%u=0.0
      allocate(estrutSistEqP%id(estrutSistEqP%ndof,numnp)); estrutSistEqP%id=0
        write(*,*) "call leituraCodigosCondContornoDS(idPotencial"
        keyword_name = "codigos_cond_contorno_potencial"
        if(estrutSistEqP%eliminate(1:3)=='YES') then
           call leituraCodigosCondContornoDS(keyword_name, estrutSistEqP%id, estrutSistEqP%ndof, numnp, &
                                          estrutSistEqP%neq, iecho, iprtin)
        else
           call leituraCodigosCondContornoDSF(keyword_name, estrutSistEqP%id, estrutSistEqP%ndof, numnp, &
                                          estrutSistEqP%neq, iecho, iprtin)
        end if
                                          

      if (estrutSistEqP%nlvect.gt.0) then
        allocate(estrutSistEqP%f(estrutSistEqP%nlvect,estrutSistEqP%ndof,numnp))
        estrutSistEqP%f = 0.0
        write(*,*) "call leituraValoresCondContornoDS(fPotencial,ndofP,numnp,0,nlvectP,iprtin)"
        keyword_name = "valores_cond_contorno_potencial"
        call leituraValoresCondContornoDS(keyword_name, estrutSistEqP%f, estrutSistEqP%ndof, numnp, &
                 1_4, estrutSistEqP%nlvect, iprtin,  iecho )
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
        call leituraValoresCondContornoDS(keyword_name, estrutSistEqF%f, estrutSistEqF%ndof, numnp, &
                1_4, estrutSistEqF%nlvect, iprtin,iecho)
      end if
!
!.... input element data
!
     write (*,*) "call leituraParamNumericosPropFisicaDS()"
     call leituraParamNumericosPropFisicaDS()
     write (*,*) "apos call leituraParamNumericosPropFisicaDS()"

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
        character(len=80):: stringDefault
        character(len=80):: keyword_name
        character(len=80):: default_title_value
        integer :: ierr

        !Reads the title
        keyword_name        = "title"
        default_title_value = "unknown title"
        call readStringKeywordValue(keyword_name, title, default_title_value, ierr)
        !print*, "title: ", title
        !Reads exec
        keyword_name = "exec"
        call readIntegerKeywordValue(keyword_name, exec, 0_4, ierr)
        !print*, "exec: ", exec
        !Reads iprtin
        keyword_name = "iprtin"
        call readIntegerKeywordValue(keyword_name, iprtin, 0_4, ierr)
        !print*, "iprtin: ", iprtin
        !Reads nsd
        keyword_name = "nsd"
        call readIntegerKeywordValue(keyword_name, nsd, 0_4, ierr)
        !print*, "nsd: ", nsd
        !Reads numnp
        keyword_name = "numnp"
        call readIntegerKeywordValue(keyword_name, numnp, 0_4, ierr)
        !print*, "numnp: ", numnp
        !Reads numel
        keyword_name = "numel"
        call readIntegerKeywordValue(keyword_name, numel, 0_4, ierr)
        !print*, "numel: ", numel
        !Reads nen
        keyword_name = "nen"
        call readIntegerKeywordValue(keyword_name, nen, 0_4, ierr)
        !print*, "nen: ", nen
        !Reads npint
        keyword_name = "npint"
        call readIntegerKeywordValue(keyword_name, npint, 0_4, ierr)
        !print*, "npint: ", npint
        !Reads nlvectP
        keyword_name = "nlvectP"
        call readIntegerKeywordValue(keyword_name, estrutSistEqP%nlvect, 0_4, ierr)
        !print*, "nlvectP: ", nlvectP
        !Reads nlvectF
        keyword_name = "nlvectF"
        call readIntegerKeywordValue(keyword_name, estrutSistEqF%nlvect, 0_4, ierr)
        !print*, "nlvectF: ", nlvectF
        
#ifdef withPardiso
        stringDefault='PardisoEsparso'
#endif
#ifdef withHYPRE
        stringDefault='HYPREEsparso'
#endif
        stringDefault='GaussSkyline'
        
        keyword_name = "solver_pressao"
        call readStringKeywordValue(keyword_name, estrutSistEqP%optSolver, stringDefault, ierr)
        !write(*,*) "++", trim(keyword_name(:)), estrutSistEqP%optSolver; stop
        
        keyword_name = "solver_fluxo"
        call readStringKeywordValue(keyword_name, estrutSistEqF%optSolver, stringDefault, ierr)

        keyword_name ="eliminateEq"
        stringDefault='YES'
        call readStringKeywordValue(keyword_name, estrutSistEqP%eliminate, stringDefault, ierr)
        
        return
    end subroutine readSetupPhaseDS 
!
!**** new **********************************************************************
!
    !> Faz a leitura de constant body forces.
    !> @param keyword_name  Keyword especifica para constant body forces.
    subroutine leituraParamNumericosPropFisicaDS()
        use mGlobaisEscalares, only: iprtin
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
        write(*,6000) m,(c(i,m),i=1,3)
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
!
!**** new *******************************************************************
!     
      subroutine verificarSolver(optSolver_, listaSolverDisponivel_) 
!
         character(len=*), intent(in) :: optSolver_
         logical, intent(in) :: listaSolverDisponivel_(*)
!
         if(optSolver_(1:14)=='PardisoEsparso')then
            if (listaSolverDisponivel_(2).eqv..false.) then
               print*, "O Solver escolhido ...,  ", optSolver_,",  não está disponível"
               !stop
            endif
         endif

         if(optSolver_(1:12)=='HYPREEsparso') then
            if(listaSolverDisponivel_(3).eqv..false.) then
               print*, "O Solver escolhido ...,  ", optSolver_, ",  não está disponível"
               stop
            endif
         endif
    
      end subroutine   
!
!**** new *******************************************************************
!     
      subroutine identificaSolversDisponiveis(listaSolverDisponivel_)
      
      logical, intent(inout) :: listaSolverDisponivel_(*)
         

      listaSolverDisponivel_(1)=.true. !skyline
#ifdef withPardiso
      listaSolverDisponivel_(2)=.true. !pardiso
#else
      listaSolverDisponivel_(2)=.false. !pardiso
#endif
#ifdef withHYPRE
      listaSolverDisponivel_(3)=.true. !hypre
#else
      listaSolverDisponivel_(3)=.false. !hypre
#endif

      print*, "Solvers disponiveis:"
      print*, "                      SKYLINE"
      if(listaSolverDisponivel_(2).eqv..true.) print*, "                      PARDISO"
      if(listaSolverDisponivel_(3).eqv..true.) print*, "                      HYPRE"

      end subroutine
