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
!
!**** new **********************************************************************
!
      subroutine montarSistEqAlgPotencial(optsolver_, estrutSistEqP_ )
!
      use mMalha,              only: x, conecNodaisElem, numel, nen, nsd, numnp
      use mUtilSistemaEquacoes, only: load, dirichletConditions
!
      implicit none

      character(len=*), intent(in) :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqP_
        
      real*8 :: t1, t2, t3
!
      print*, " ..... montando sistema de equacoes da potencial"
      write(*,*) estrutSistEqP_%ndof, estrutSistEqP_%nalhs, estrutSistEqP_%neq 
     
      write(*,*) " +++ apos call montarEstrutDadosSistEqAlq("
      if (estrutSistEqP_%nlvect.gt.0) call load                (estrutSistEqP_%id,estrutSistEqP_%f,estrutSistEqP_%brhs,&
                                                 estrutSistEqP_%ndof,numnp,estrutSistEqP_%nlvect)
    
      if (estrutSistEqP_%nlvect.gt.0) call dirichletConditions (estrutSistEqP_%id,estrutSistEqP_%u,estrutSistEqP_%f, &
                                                 estrutSistEqP_%ndof,numnp,estrutSistEqP_%nlvect)
      write(*,*) " +++ apos call dirichletConditions("
! 
      call timing(t1)
      call calcCoefSistAlgPotencial (optsolver_, estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd )
      write(*,*) " +++ apos call calcCoefSistAlgPotencial("
      call timing(t2) 
      write(*,*) " calculo dos coeficientes, tempo = ", t2 - t1

      end subroutine montarSistEqAlgPotencial
!
!**** new **********************************************************************
!
      subroutine calcCoefSistAlgPotencial( optSolver_, estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd )
!
      use mGlobaisEscalares, only: nrowsh, npint
      use mGlobaisArranjos,  only: mat, c, grav
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq
      use mMalha,            only: local

      use mUtilSistemaEquacoes,  only: kdbc

      use mSolverGaussSkyline,   only: addrhs, addlhs
      use mSolverPardiso,        only: addlhsCSR, addlhsCSR01
      use mSolverHypre,          only: addnslHYPRE, atribuirValoresVetor_HYPRE, adicionarValoresVetor_HYPRE
      use mSolverHypre,       only: fecharMatriz_HYPRE, fecharVetor_HYPRE
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
      real*8 :: xl(nsd,nen), pl(1,nen)
      real*8 :: shg(nrowsh,nen,npint), shl(nrowsh,nen,npint)
      real*8 :: det(npint), w(npint)
!
      integer*4 :: nee
      real*8   :: elresf(estrutSistEqP_%ndof*nen), eleffm(estrutSistEqP_%ndof*nen,estrutSistEqP_%ndof*nen)
!
      integer*4:: nel, m, l, i, j, ni, nj
      integer*4, parameter :: um = 1
      real*8  :: pi, Kx, Ky, Kz
      real*8  :: temp1, gf1, gf2, gf3
      real*8  :: pss
      real*8  :: djx, djy, djz, djn, dix, diy, diz
      logical :: diag,zerodl,quad,lsym
      character(len=30) :: nomeA, nomeB
      integer*4, allocatable:: lmLocal(:)
!
      nee = nen*estrutSistEqP_%ndof
       allocate (lmLocal(nee))
      diag = .false.
      pi=4.d00*datan(1.d00)
!
      w=0.0
      shl=0.0

      if(nen==2) call oneshl(shl,w,npint,nen) 
      if(nen==3) call   shlt(shl,w,npint,nen)
      if(nen==4) call   shlq(shl,w,npint,nen)
      if(nen==8) call shlq3d(shl,w,npint,nen)

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
!
      m = mat(nel)
      quad = .true.
      if (nen.eq.4.and.conecNodaisElem(3,nel).eq.conecNodaisElem(4,nel)) quad = .false.
!
      if(nen==2) call oneshg(xl,det,shl,shg,nen,npint,nsd,um,nel,um)
      if(nen==3) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)
      if(nen==4) call shgq  (xl,det,shl,shg,npint,nel,quad,nen)
      if(nen==8) call shg3d (xl,det,shl,shg,npint,nel,nen)
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
      do 300 j=1,nen
           nj=estrutSistEqP_%ndof*j
                 djx=shg(1,j,l)*temp1
                 djy=shg(2,j,l)*temp1
      if(nsd==3) djz=shg(3,j,l)*temp1
      djn=shg(nrowsh,j,l)*temp1
!     
!.... source terms      
!
      elresf(nj)=elresf(nj)+pss*djn
!
      do 300 i=1,nen
      ni = estrutSistEqP_%ndof*i
                 dix=shg(1,i,l)
                 diy=shg(2,i,l) 
      if(nsd==3) diz=shg(3,i,l) 
!
                 eleffm(ni,nj)=eleffm(ni,nj)+ Kx*dix*djx 
                 eleffm(ni,nj)=eleffm(ni,nj)+ Ky*diy*djy
      if(nsd==3) eleffm(ni,nj)=eleffm(ni,nj)+ Kz*diz*djz

!
  300 continue
  400 continue  
!
!      computation of Dirichlet b.c. contribution
!
       call ztest(pl,nee,zerodl)
!
      if(.not.zerodl) then
          call kdbc(eleffm,elresf,pl,nee)
      endif
!
!.... assemble element stifness matrix and force array into global
!        left-hand-side matrix and right-hand side vector
      lsym=.true.

      lmLocal(:)=reshape(estrutSistEqP_%lm(:,:,nel),(/nee/))

      if (optSolver_=='GaussSkyline')   then
         call addlhs   (estrutSistEqP_%alhs, eleffm, lmLocal, estrutSistEqP_%idiag, nee, diag, lsym) 
      endif
      if (optSolver_=='PardisoEsparso') then
        ! call addlhsCSR01  (estrutSistEqP_, eleffm, lm, nee)
         call addlhsCSR  (estrutSistEqP_%alhs, eleffm, lmLocal, estrutSistEqP_%Ap, estrutSistEqP_%Ai,  nee)
      endif
      if (optSolver_=='HYPREEsparso')   then
!      call addnslHYPRE(estrutSistEqP_, eleffm, nel)
       call addnslHYPRE(estrutSistEqP_%A_HYPRE, eleffm, lmLocal, estrutSistEqP_%idiag, nee, diag, lsym)
      endif
      call addrhs     (estrutSistEqP_%brhs, elresf, lmLocal, nee)
  500 continue

      if (optSolver_=='HYPREEsparso')   then
       do i = 1, estrutSistEqP_%neq
           estrutSistEqP_%rows(i) = i-1 
       end do
      !call atribuirValoresVetor_HYPRE(estrutSistEqP_%b_HYPRE, 1, estrutSistEqP_%neq, estrutSistEqP_%rows, estrutSistEqP_%brhs)
      call adicionarValoresVetor_HYPRE(estrutSistEqP_%b_HYPRE, 1, estrutSistEqP_%neq, estrutSistEqP_%rows, estrutSistEqP_%brhs)
      call fecharMatriz_HYPRE            (estrutSistEqP_%A_HYPRE, estrutSistEqP_%parcsr_A)
      call fecharVetor_HYPRE             (estrutSistEqP_%b_HYPRE, estrutSistEqP_%par_b   )
      call fecharVetor_HYPRE             (estrutSistEqP_%u_HYPRE, estrutSistEqP_%par_u   )

      endif
!         write(*,*) estrutSistEqP_%alhs ! estrutSistEqP_%Ap, estrutSistEqP_%Ai,  nel)

      

      return
      end subroutine calcCoefSistAlgPotencial

      end module
!
