!=================================================================================
!
!         programa de elementos finitos em fortran 90, 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!         + implementacoes de Abimael Loula
!
!         novo projeto: 
!         Eduardo Garcia,        bidu@lncc.br [4]
!         Tuane Lopes,           tuane@lncc.br [5]
!
!         LNCC/MCT
!         Petropolis, 02.2014
!=================================================================================
!     
      module mPotencial
        use mEstruturasDadosSistEq
        type (estruturasArmazenamentoSistemaEq) :: estrutSistEqP
        contains 
!**** new **********************************************************************
      subroutine montarSistEqAlgPotencial(estrutSistEqP_ )
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqP_
      call  montarSistEqAlgPotencial0(estrutSistEqP_%optsolver, estrutSistEqP_ )
      end subroutine montarSistEqAlgPotencial
!**** new **********************************************************************
      subroutine montarSistEqAlgPotencial0(optsolver_, estrutSistEqP_ )
      use mMalha,               only: x, conecNodaisElem, numel, nen, nsd, numnp
      use mUtilSistemaEquacoes, only: load, dirichletConditions
      implicit none
      character(len=*), intent(in)                           :: optSolver_
      type (estruturasArmazenamentoSistemaEq), intent(inout) :: estrutSistEqP_
      real*8 :: t1, t2, t3

      write(*,*) "em subroutine montarSistEqAlgPotencial0"
      print*, " ..... montando sistema de equacoes de potencial"
      if (estrutSistEqP_%nlvect.gt.0) call load                (estrutSistEqP_%id,estrutSistEqP_%f,estrutSistEqP_%brhs,&
                                                estrutSistEqP_%ndof,numnp,estrutSistEqP_%nlvect)
    
      if (estrutSistEqP_%nlvect.gt.0) call dirichletConditions (estrutSistEqP_%id,estrutSistEqP_%u,estrutSistEqP_%f, &
                                                 estrutSistEqP_%ndof,numnp,estrutSistEqP_%nlvect)
      estrutSistEqP_%uTempoAnt=estrutSistEqP_%u
      call timing(t1)
      call calcCoefSistAlgPotencial (optsolver_, estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd )
      call timing(t2) 
      write(*,*) " calculo dos coeficientes, tempo = ", t2 - t1
      end subroutine montarSistEqAlgPotencial0
!
!**** new **********************************************************************
!
      subroutine calcCoefSistAlgPotencial( optSolver_, estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd )
!
      use mGlobaisEscalares, only: nrowsh, npint, pTempo
      use mGlobaisArranjos,  only: mat, c, grav
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq, shlt3D
      use mMalha,            only: local
      use mUtilSistemaEquacoes,  only: kdbc, kdbcF
      use mSolverGaussSkyline,   only: addrhs, addlhs, addlhsN
      use mSolverPardiso,        only: addlhsCSR!, addlhsCSR01
      use mSolverHypre,          only: addnslHYPRE, atribuirValoresVetorHYPRE, adicionarValoresVetorHYPRE
      use mSolverHypre,       only: fecharMatrizHYPRE, fecharVetorHYPRE, HYPRE_PARCSR, mpi_comm
      use mSolverHypre,       only: escreverMatrizHYPRE
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!     
      character(len=*), intent(in) :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqP_  
      real*8,    intent(in)    :: x(nsd, numnp)      
      integer*4, intent(in)    :: numnp, numel, nen, nsd
      integer*4, intent(in)    :: conecNodaisElem(nen,numel)
!
      real*8 :: xl(nsd,nen), pl(1,nen), plTempoAnt(1,nen), pTA_int
      real*8 :: shg(nrowsh,nen,npint), shl(nrowsh,nen,npint)
      real*8 :: det(npint), w(npint)
!
      integer*4 :: nee
      real*8   :: elresf(estrutSistEqP_%ndof*nen), eleffm(estrutSistEqP_%ndof*nen,estrutSistEqP_%ndof*nen)
!
      integer*4:: nel, m, l, i, j, k, ni, nj
      integer*4, parameter :: um = 1
      real*8  :: pi, Kx, Ky, Kz
      real*8  :: temp1, gf1, gf2, gf3
      real*8  :: pss
      real*8  :: djx, djy, djz, djn, din,  dix, diy, diz
      logical :: diag,zerodl,quad,lsym
      character(len=30) :: nomeA, nomeB
      integer*4, allocatable:: lmLocal(:)
      character(len=200) ::  matrixFile ="matrixHypre.dat"
!
      write(*,*) "Om subroutine calcCoefSistAlgPotencial"
      nee = nen*estrutSistEqP_%ndof
       allocate (lmLocal(nee))
      diag = .false.
      pi=4.d00*datan(1.d00)
!
      w=0.0
      shl=0.0


      if(nen==2) call oneshl(shl,w,npint,nen) 
      if(nen==3) call   shlt(shl,w,npint,nen)
      if(nen==4.and.nsd==2) call   shlq(shl,w,npint,nen)
      if(nen==4.and.nsd==3) call shlt3D(shl,w,npint,nen)
      if(nen==8) call shlq3d(shl,w,npint,nen)

      estrutSistEqP_%brhs= pTempo*estrutSistEqP_%brhs
      !print* , "materiais:...", c(:,1);stop

      do 500 nel=1,numel
!
!      clear stiffness matrix and force array
!
      eleffm=0.0
      elresf=0.0
!
!      LOCALIZE COORDINATes and Dirichlet b.c.
!
      call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
      call local(conecNodaisElem(1,nel),estrutSistEqP_%u,pl,nen,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
      call local(conecNodaisElem(1,nel),estrutSistEqP_%uTempoAnt,pltempoAnt,nen,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
!
      m = mat(nel)
      quad = .true.
      if (nen.eq.4.and.conecNodaisElem(3,nel).eq.conecNodaisElem(4,nel)) quad = .false.
!
      if(nen==2) call oneshg(xl,det,shl,shg,nen,npint,nsd,um,nel,um)
      if(nen==3) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)
      if(nen==4.and.nsd==2) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)
      if(nen==4.and.nsd==3) call shg3d (xl,det,shl,shg,npint,nel,nen)
      if(nen==8)            call shg3d (xl,det,shl,shg,npint,nel,nen)
!
!....... form stiffness matrix
!
!... length of the element
!
      Kx=c(1,m)
      Ky=c(2,m)
      if(nsd==3)Kz=c(3,m)
!
!.... loop on integration points
!
      do 400 l=1,npint

      temp1 = w(l)*det(l)
! 
!.... vetor de carga - RHS - f  =  gf0  
!
                 gf1 = grav(1)
                 gf2 = grav(2)
      if(nsd==3) gf3 = grav(3)
      pss = 0.0 
!
      pTA_int=0.0 
      DO J=1,NEN
             pTA_int  = pTA_int + SHG(nrowsh,J,L)*plTempoAnt(1,J)
      ENDDO

      do j=1,nen
           nj=estrutSistEqP_%ndof*j
                 djx=shg(1,j,l)*temp1
                 djy=shg(2,j,l)*temp1
      if(nsd==3) djz=shg(3,j,l)*temp1
                 djn=shg(nrowsh,j,l)*temp1
!.... source terms      
!
      elresf(nj)=elresf(nj)+pTA_int*djn
!
      do  i=1,nen
      ni = estrutSistEqP_%ndof*i
                 dix=shg(1,i,l)
                 diy=shg(2,i,l) 
      if(nsd==3) diz=shg(3,i,l) 
                 din=shg(nrowsh,i,l)
!
                 eleffm(ni,nj)=eleffm(ni,nj)+ pTempo*Kx*dix*djx 
                 eleffm(ni,nj)=eleffm(ni,nj)+ pTempo*Ky*diy*djy
      if(nsd==3) eleffm(ni,nj)=eleffm(ni,nj)+ pTempo*Kz*diz*djz
      eleffm(ni,nj)=eleffm(ni,nj)+ din*djn

      end do
      end do
  400 continue  
!
!      computation of Dirichlet b.c. contribution
!      write(*,*) "call addnslHYPRE(estrutSistEqP_%A_HYPRE, eleffm, estrutSistEqP_%idiag, lmLocal, nee, diag, lsym)"
!
      lsym=.true.
      lmLocal(:)=reshape(estrutSistEqP_%lm(:,:,nel),(/nee/))

      call ztest(pl,nee,zerodl)

      if(.not.zerodl) then
        if(estrutSistEqP_%eliminate(1:3)=='YES') then
          call kdbc(eleffm,elresf,pl,lmLocal, nee)
        else
          call kdbcF(eleffm,elresf,pl,lmLocal, nee)
        endif
      endif
!
!.... assemble element stifness matrix and force array into global
!        left-hand-side matrix and right-hand side vector

      !write(*,*) " lmLocal =", lmLocal(:)
      if (optSolver_=='GaussSkyline')   then
         call addlhsN   (estrutSistEqP_%alhs, eleffm, estrutSistEqP_%idiag, estrutSistEqP_%lm,  nee, nel, diag, lsym) 
         !call addlhsN   (estrutSistEqP_%alhs, eleffm, estrutSistEqP_%idiag, estrutSistEqP_%lm(:,:,nel),  nee, nel, diag, lsym) 
         !       subroutine addlhsN(alhs,eleffm, idiag, lmT, nee, nel, ldiag,lsym)
         !call addlhs   (estrutSistEqP_%alhs, eleffm, estrutSistEqP_%idiag, lmLocal,  nee, diag, lsym) 
      endif
      if (optSolver_=='PardisoEsparso') then
        ! call addlhsCSR01  (estrutSistEqP_, eleffm, lm, nee)
      !   write(*,*) "call addlhsCSR, nee= ", nee 
         call addlhsCSR  (estrutSistEqP_%alhs, eleffm, estrutSistEqP_%Ap, estrutSistEqP_%Ai, lmLocal,  nee)
      endif
      if (optSolver_=='HYPREEsparso')   then
       !write(*,*) nel, "  (optSolver_==HYPREEsparso) then"
!      call addnslHYPRE(estrutSistEqP_, eleffm, nel)
       call addnslHYPRE(estrutSistEqP_%A_HYPRE, eleffm, estrutSistEqP_%idiag, lmLocal, nee, diag, lsym)
      endif
      call addrhs     (estrutSistEqP_%brhs, elresf, lmLocal, nee)
  500 continue

    print*, estrutSistEqP_%id;! stop

      if (optSolver_=='HYPREEsparso')   then
       do i = 1, estrutSistEqP_%neq
           estrutSistEqP_%rows(i) = i
       end do
      call adicionarValoresVetorHYPRE(estrutSistEqP_%b_HYPRE, 1, estrutSistEqP_%neq, estrutSistEqP_%rows, estrutSistEqP_%brhs)
      write(*,*) "call adicionarValoresVetorHYPRE(estrutSistEqP_%b_HYPRE," 

      call fecharMatrizHYPRE            (estrutSistEqP_)
      call fecharVetorHYPRE             (estrutSistEqP_)
      call fecharVetorHYPRE             (estrutSistEqP_)!

      endif
      return
      end subroutine calcCoefSistAlgPotencial

      end module
!
