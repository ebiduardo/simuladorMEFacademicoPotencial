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
      module mFluxo

         use mEstruturasDadosSistEq
              
        type (estruturasArmazenamentoSistemaEq) :: estrutSistEqF
              
      contains
!
!**** new **********************************************************************
!
      subroutine montarSistEqAlgFluxo(estrutSistEqF_, estrutSistEqP_)
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqF_, estrutSistEqP_
      call montarSistEqAlgFluxo0(estrutSistEqF_%optsolver, estrutSistEqF_, estrutSistEqP_)
      end subroutine montarSistEqAlgFluxo

      subroutine montarSistEqAlgFluxo0(optsolver_, estrutSistEqF_, estrutSistEqP_)
!
      use mMalha,              only: x, conecNodaisElem, numel, nen, nsd, numnp
      use mUtilSistemaEquacoes,only: load, dirichletConditions
!
      implicit none
!
      character(len=*), intent(in                          ) :: optSolver_
      type (estruturasArmazenamentoSistemaEq), intent(inout) :: estrutSistEqF_, estrutSistEqP_

      real*8 :: t1, t2, t3
      print*, " ..... montando sistema de equacoes de fluxo"
      write(*,*) estrutSistEqF_%ndof, estrutSistEqF_%nalhs, estrutSistEqF_%neq 
      if (estrutSistEqF_%nlvect.gt.0) call load                (estrutSistEqF_%id,estrutSistEqF_%f,estrutSistEqF_%brhs,&
                                estrutSistEqF_%ndof,numnp,estrutSistEqF_%nlvect)
      if (estrutSistEqF_%nlvect.gt.0) call dirichletConditions (estrutSistEqF_%id,estrutSistEqF_%u,estrutSistEqF_%f,&
                                estrutSistEqF_%ndof,numnp,estrutSistEqF_%nlvect)
! 
      call timing(t1)
      call calcCoefSistAlgFluxo(optSolver_, estrutSistEqF_, estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd )
      call timing(t2)
      write(*,*) " calculo dos coeficientes, tempo = ", t2 - t1
      end subroutine montarSistEqAlgFluxo0
!
!**** new **********************************************************************
!
      subroutine calcCoefSistAlgFluxo(optSolver_, estrutSistEqF_, estrutSistEqP_, x, conecNodaisElem, numnp, numel, nen, nsd)

      use mGlobaisEscalares,  only: zero, one, two, four, npint, nrowsh
      use mGlobaisArranjos,   only: mat, c, grav
      use mfuncoesDeForma,   only: oneshl, oneshg, shlt, shlq3d, shg3d, shgq, shlq, shlt3D
      use mmalha,             only: local

      use mUtilSistemaEquacoes,      only: kdbcF

      use mSolverGaussSkyline, only: addrhs, addlhs
      use mSolverPardiso,    only: addlhsCSR!, addlhsCSR01
      use mSolverHypre,      only: addnslHYPRE, atribuirValoresVetorHYPRE, adicionarValoresVetorHYPRE
      use mSolverHypre,       only: fecharMatrizHYPRE, fecharVetorHYPRE


      implicit none
!
      character(len=*), intent(in) :: optSolver_
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEqF_ , estrutSistEqP_
      real*8,  intent(in)    :: x(nsd,numnp)
      integer*4, intent(in)  :: numnp, numel, nen, nsd
      integer*4, intent(in)  :: conecNodaisElem(nen,numel)
           
!
      real*8 :: xl(nsd,nen), pl(estrutSistEqP_%ndof,nen)
      real*8 :: shl(nrowsh,nen,npint), shg(nrowsh,nen,npint)
      real*8 :: w(npint), det(npint)
      real*8 :: elref(estrutSistEqF_%ndof*nen), eleff(estrutSistEqF%ndof*nen,estrutSistEqF%ndof*nen)
! 
      logical :: diag,quad,lsym, zerodl
      integer*4:: nel, m, l, i, j, nee, ni, nj
      real*8  :: pi, dpi, rmu, rmu0, rmm, rkx, rky, rkz, coefx, coefy, coefz
      real*8  :: delta, delta1, delta2, deltax, deltay, deltaz, h, h2, grtx, grty, grtz
      real*8  :: temp1, dix, diy, diz, di, dj, djx, djy, djz
      integer*4, parameter :: um = 1
      integer*4, allocatable:: lmLocal(:)


!
      print*, " em calcCoefSistAlgFluxo, nen =", nen
      !write(*,*) (i,  brhs(i), i=1,neq); stop
      pi  =four*datan(one)
      dpi =two*pi
      nee = nen*estrutSistEqF%ndof

      allocate (lmLocal(nee))
!
      w   = 0.0
      det = 0.0
      shl = 0.0
!
      if(nen==2) call oneshl(shl,w,npint,nen) 
      if(nen==3) call shlt  (shl,w,npint,nen)
      if(nen==4) call shlq  (shl,w,npint,nen)
      if(nen==4.and.nsd==3) call shlt3D(shl,w,npint,nen)
      if(nen==8) call shlq3d(shl,w,npint,nen)
!
!      consistent matrix
!
      diag = .false.
!
      do 500 nel=1,numel
      
!
!      set material properties
!
      m     = mat(nel)
      rmu0  = 1.00
      rmm   = 1.00
      rkx = c(1,m)
      rky = c(2,m)
      if(nsd==3) rkz = c(3,m)
!
!     clear stiffness matrix and force array
!
      eleff=0.0
      elref=0.0
!
!      localize coordinates and dirichlet b.c.
!
      !write(*,*) "nel=", nel
      call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd) 
      call local(conecNodaisElem(1,nel),estrutSistEqP_%u,pl,nen,estrutSistEqP_%ndof,estrutSistEqP_%ndof)
!      call local(conecNodaisElem(1,nel),v,vl,nen,ndof,ndof)
!
      quad   = .true.
!
      if(nen==2) call oneshg(xl,det,shl,shg,nen,npint,nsd,um,nel,um)
      if(nen==3) call shgq (xl,det,shl,shg,npint,nel,quad,nen)
      if(nen==4.and.nsd==2) call shgq (xl,det,shl,shg,npint,nel,quad,nen)
      if(nen==4.and.nsd==3) call shg3d (xl,det,shl,shg,npint,nel,nen)
      if(nen==8) call shg3d(xl,det,shl,shg,npint,nel,nen)
!
      h2=zero
      if(nen==2) h2=h2+(xl(1,2)-xl(1,1))**2+(xl(2,2)-xl(2,1))**2

      if(nsd==2.and.nen==4) then
         h2=xl(1,2)*xl(2,3)+xl(1,1)*xl(2,2)+xl(1,3)*xl(2,1) &
           -xl(1,1)*xl(2,3)-xl(1,2)*xl(2,1)-xl(1,3)*xl(2,2)
      endif

      if(nsd==3.and.nen==4) then
         h2=h2+(xl(1,1)-xl(1,2))**2+(xl(2,1)-xl(2,4))**2
         h2=h2+(xl(3,1)-xl(3,4))**2 
      endif

      if(nen==8) then
         h2=h2+(xl(1,1)-xl(1,2))**2+(xl(2,1)-xl(2,4))**2
         h2=h2+(xl(3,1)-xl(3,5))**2 
      endif

      h=dsqrt(h2)
!
      delta = delta1*(h)**delta2 ! parametro do posprocessamento
!
      delta1 = 1.0
      delta2 = 1.0
      delta = delta1*(h)**delta2
      delta = 1.0
      deltax = delta
      deltay = delta
      if(nsd==3) deltaz = delta

!.... loop on integration points
!
      do 400 l=1,npint
!
      temp1 = w(l)*det(l)
! 
!.... source terms 
!
      grtx=zero
      grty=zero
      if(nsd==3) grtz=zero
!
      do 301 j=1,nen 
      grtx=grtx+shg(1,j,l)*pl(1,j)
      grty=grty+shg(2,j,l)*pl(1,j)
      if(nsd==3) grtz=grtz+shg(3,j,l)*pl(1,j)
301   continue
!
      rmu = 1.0
!      
      do 300 j=1,nen
      nj=estrutSistEqF%ndof*j
     
!
      djx =shg(1,j,l)*temp1
      djy =shg(2,j,l)*temp1
      if(nsd==3) djz =shg(3,j,l)*temp1
      dj  =shg(nrowsh,j,l)*temp1
 
!
!.....body forces due to the regular part of the pressure
!
      if(nsd==2) then
         elref(nj-1)= elref(nj-1) - rkx * grtx*dj
         elref(nj)  = elref(nj)   - rky * grty*dj
      else
         elref(nj-2)= elref(nj-2) - rkx * grtx*dj
         elref(nj-1) =elref(nj-1) - rky * grty*dj
         elref(nj  ) =elref(nj  ) - rkz * grtz*dj
        
      endif

      do 311 i=1,nen
      ni = estrutSistEqF%ndof*i
!
      dix =shg(1,i,l) 
      diy =shg(2,i,l)
      if(nsd==3)diz =shg(3,i,l)
      di  =shg(nrowsh,i,l)
!
      coefx = one      
      coefy = one 
      if(nsd==3) coefz = one
!
      if(nsd==2) then
         eleff(ni-1,nj-1)=eleff(ni-1,nj-1)+coefx*di*dj+deltax*dix*djx
         eleff(ni-1,nj)=eleff(ni-1,nj)      +delta*dix*djy

         eleff(ni,nj-1)=eleff(ni,nj-1)      +delta*diy*djx
         eleff(ni,nj)=eleff(ni,nj)+coefy*di*dj+deltay*diy*djy

      else
         eleff(ni-2,nj-2)=eleff(ni-2,nj-2)+coefx*di*dj+deltax*dix*djx
         eleff(ni-2,nj-1)=eleff(ni-2,nj-1)      +delta*dix*djy
         eleff(ni-2,nj  )=eleff(ni-2,nj  )      +delta*dix*djz

         eleff(ni-1,nj-2)=eleff(ni-1,nj-2)      +delta*diy*djx
         eleff(ni-1,nj-1)=eleff(ni-1,nj-1)+coefy*di*dj+deltay*diy*djy
         eleff(ni-1,nj  )=eleff(ni-1, nj )      +delta*diy*djz

         eleff(ni,nj-2)=eleff(ni,nj-2)          +delta*diz*djx
         eleff(ni,nj-1)=eleff(ni,nj-1)          +delta*diz*djy
         eleff(ni,nj  )=eleff(ni,nj  )    +coefz*di*dj+deltaz*diz*djz
      endif
!
  311 continue
!
  300 continue
!
  400 continue

!      computation of dirichlet b.c. contribution
!
!      if(.not.zerodl) call kdbc(eleff,elref,vl,nee)!,lmF(1,1,nel))
!
!.... assemble element stifness matrix and force array into global
!        left-hand-side matrix and right-hand side vector
!
      lmLocal(:)=reshape(estrutSistEqF_%lm(:,:,nel),(/nee/))
      lsym=.true.
      if (optSolver_=='GaussSkyline')   then
         call addlhs   (estrutSistEqF_%alhs, eleff, estrutSistEqF_%idiag, lmLocal, nee, diag, lsym) 
      endif
      if (optSolver_=='PardisoEsparso') then
        ! call addlhsCSR01 (estrutSistEqF_, eleff, nee) 
         call addlhsCSR   (estrutSistEqF_%alhs, eleff, estrutSistEqF_%Ap, estrutSistEqF_%Ai, lmLocal,  nee)
      endif
      if (optSolver_=='HYPREEsparso')   then
          call addnslHYPRE   (estrutSistEqF_%A_HYPRE, eleff, estrutSistEqF_%idiag, lmLocal, nee, diag, lsym)
      endif
      call addrhs     (estrutSistEqF_%brhs, elref, lmLocal, nee)

  500 continue
  
      if (optSolver_=='HYPREEsparso')   then
       do i = 1, estrutSistEqF_%neq
           estrutSistEqF_%rows(i) = i
       end do
      call adicionarValoresVetorHYPRE(estrutSistEqF_%b_HYPRE, 1, estrutSistEqF_%neq, estrutSistEqF_%rows, estrutSistEqF_%brhs)
         call fecharMatrizHYPRE            (estrutSistEqF_)  !%A_HYPRE)
         call fecharVetorHYPRE             (estrutSistEqF_) !%b_HYPRE)
         call fecharVetorHYPRE             (estrutSistEqF_)  !%u_HYPRE)

      endif



      return
      end subroutine

      end module

