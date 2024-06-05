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
       module mUtilSistemaEquacoes

       use mEstruturasDadosSistEq !, only : colht
       use mSolverHYPRE, only : criarMatrizHYPRE, criarVetorHYPRE
       use mSolverHYPRE, only : destruirMatrizHYPRE, destruirVetorHYPRE


!funcoes e subrotinas
        public :: load, dirichletConditions
        public :: kdbcF
        public :: montarEstrutDadosSistEqAlq

      contains

!**** new **********************************************************************
!
      subroutine montarEstrutDadosSistEqAlq(umaEstSistEq_)
!
      use mMalha,            only: x, conecNodaisElem, numnp, numel, nen, nsd
      use mMalha,            only: criarListaVizinhos, listaDosElemsPorNoCSR, numConexoesPorElem
      use mMalha,            only: criarListaVizinhosCRS
      use mMalha,            only: criarListaVizinhosCRS00
      use mLeituraEscrita,   only: iecho
      use mSolverHypre, only: mpi_comm
!
      implicit none
      type(estruturasArmazenamentoSistemaEq), intent(inout) :: umaEstSistEq_

      logical   :: simetria
      integer*4 :: meanbw
      real*8    :: t1, t2, t3
      integer*4 :: i, localSize

      write(*,'(a,i2,a)',advance='NO') " em montarEstrutDadosSistEqAlq, fev2019 "
      write(*,*) " "
    
      call timing(t1)
      allocate(umaEstSistEq_%lm(umaEstSistEq_%ndof,nen,numel))
      call formlm(umaEstSistEq_%id,conecNodaisElem,umaEstSistEq_%lm,umaEstSistEq_%ndof,umaEstSistEq_%ndof,nen,numel)
      allocate(umaEstSistEq_%idiag(umaEstSistEq_%neq+1)); umaEstSistEq_%idiag = 0  ! posicoes dos elementos da diagonal principal no vetor alhs    
      
  if(umaEstSistEq_%optSolver=="GaussSkyline") then
      call colht(umaEstSistEq_%idiag,umaEstSistEq_%lm,umaEstSistEq_%ndof,nen,numel,umaEstSistEq_%neq);
      call diag (umaEstSistEq_%idiag,umaEstSistEq_%neq,umaEstSistEq_%nalhs); 
      !stop
      write(*,*) " nalhs =", umaEstSistEq_%nalhs
      if(.not.associated(umaEstSistEq_%alhs)) then 
          allocate(umaEstSistEq_%alhs(umaEstSistEq_%nalhs));
          umaEstSistEq_%alhs=0.0 
      endif
      if(.not.umaEstSistEq_%simetria.and.(.not.associated(umaEstSistEq_%alhs)))then 
          allocate(umaEstSistEq_%clhs(umaEstSistEq_%nalhs));
          umaEstSistEq_%clhs=0.0 
      end if
	  umaEstSistEq_%Ap=> umaEstSistEq_%idiag; !posicoes dos elementos da diagonal principal no vetor alhs

      if(.not.associated (listaDosElemsPorNoCSR) ) then
        allocate(listaDosElemsPorNoCSR(nen,numnp)); listaDosElemsPorNoCSR=0

      call criarListaVizinhosCRS(nen,numnp,numel,conecNodaisElem,umaEstSistEq_%nVizinMax)
      !simetria=.true.
		 print*, "umaEstSistEq_%numCoefPorLinha = ", umaEstSistEq_%numCoefPorLinha
!     call criarGrafoEquacoesPorNo(umaEstSistEq_,  nsd, conecNodaisElem, listaDosElemsPorNoCSR, &
!                                     numnp, nen, numConexoesPorElem)
      end if
	  
   end if 
  
  if(umaEstSistEq_%optSolver=="PardisoEsparso") then

      if(.not.associated (listaDosElemsPorNoCSR) ) then
        allocate(listaDosElemsPorNoCSR(nen,numnp)); listaDosElemsPorNoCSR=0
      end if
      call criarListaVizinhosCRS(nen,numnp,numel,conecNodaisElem,umaEstSistEq_%nVizinMax)
      !simetria=.true.
      umaEstSistEq_%Ap=> umaEstSistEq_%idiag; !posicoes dos elementos da diagonal principal no vetor alhs
      call criarPonteirosMatEsparsa_CSR(umaEstSistEq_, nsd, conecNodaisElem, listaDosElemsPorNoCSR, & 
                                                             numnp, nen, numConexoesPorElem)       
      !write(*,    6000) 'Informacao do sistema de equacoes ',  &
      !   umaEstSistEq_%neq, umaEstSistEq_%nalhs, meanbw,       &
      !  (8.0*umaEstSistEq_%nalhs)/1000.0/1000.0,  umaEstSistEq_%optSolver
      if(.not.associated(umaEstSistEq_%alhs)) then
         allocate(umaEstSistEq_%alhs(umaEstSistEq_%nalhs));
         umaEstSistEq_%alhs=0.0
      end if
  endif 

  if(umaEstSistEq_%optSolver=="HYPREEsparso") then
     ! call destruirMatrizHYPRE       (umaEstSistEq_%A_HYPRE)
     ! call destruirVetorHYPRE        (umaEstSistEq_%b_HYPRE)
     ! call destruirVetorHYPRE        (umaEstSistEq_%u_HYPRE)
    ! este procedimento necessita de conhecer  Flower, Fupper
      umaEstSistEq_%Flower = 1 !- 1
      umaEstSistEq_%Fupper = umaEstSistEq_%neq!-1
      localSize=umaEstSistEq_%FUpper-umaEstSistEq_%Flower+1
      if(.not.associated(umaEstSistEq_%rows)) allocate(umaEstSistEq_%rows(localSize))
      do i = 1, localSize
      umaEstSistEq_%rows(i) = i !- 1
      end do
      call criarMatrizHYPRE (umaEstSistEq_%A_HYPRE, umaEstSistEq_%Flower, umaEstSistEq_%Fupper, mpi_comm )
      call criarVetorHYPRE  (umaEstSistEq_%b_HYPRE, umaEstSistEq_%Flower, umaEstSistEq_%Fupper, mpi_comm )
      call criarVetorHYPRE  (umaEstSistEq_%u_HYPRE, umaEstSistEq_%Flower, umaEstSistEq_%Fupper, mpi_comm )

       write(*,'(a,a,", ", i0,", ",i0)') "optSolver_  =", umaEstSistEq_%optSolver, umaEstSistEq_%Flower,umaEstSistEq_%Fupper
  endif 

      if(.not.associated(umaEstSistEq_%brhs))then
       allocate(umaEstSistEq_%brhs(umaEstSistEq_%neq));
       umaEstSistEq_%brhs=0.0
      end if
! 
      call timing(t2)
      write(*,*) "criar estruturas de dados,  tempo de parede = ", t2 - t1
      
      meanbw = umaEstSistEq_%nalhs/umaEstSistEq_%neq
      write(*,    6000) 'Informacao do sistema de equacoes ',  &
         umaEstSistEq_%neq, umaEstSistEq_%nalhs, meanbw,       &
        (8.0*umaEstSistEq_%nalhs)/1000.0/1000.0,  umaEstSistEq_%optSolver
      write(iecho,6000) 'Informacao do sistema de equacoes ',  &
         umaEstSistEq_%neq, umaEstSistEq_%nalhs, meanbw,       &
        (8.0*umaEstSistEq_%nalhs)/1000.0/1000.0,  umaEstSistEq_%optSolver
         
 6000 format(a///&
     ' e q u a t i o n    s y s t e m    d a t a              ',  //5x,&
     ' number of equations . . . . . . . . . . . . (neq    ) = ',i8//5x,&
     ' number of terms in left-hand-side matrix  . (nalhs  ) = ',i12//5x,&
     ' mean half bandwidth . . . . . . . . . . . . (meanbw ) = ',i8//5x,&
     ' memoria necessaria para a matriz do sistema (Mbytes)  = ',e10.2//5x, &
     ' Solver escolhido                                      = ',a)

      end subroutine montarEstrutDadosSistEqAlq
!
        
!**** new **********************************************************************
!
      subroutine load(id,f,brhs,ndof,numnp,nlvect)
!
!.... program to accumulate nodal forces and transfer into
!        right-hand-side vector
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4 :: ndof, numnp, nlvect
      integer*4 :: id(ndof,*)
      real*8    :: f(ndof,numnp,*),brhs(*)
!
      integer*4 :: nlv
      integer*4 :: i, j, k
!
      do 300 i=1,ndof
!
      do 200 j=1,numnp
      k = id(i,j)
      ! k = abs(id(i,j))
      if (k.gt.0) then
!
         do 100 nlv=1,nlvect
           brhs(k) = brhs(k) + f(i,j,nlv)
  100    continue
!
      endif
!
  200 continue
!
  300 continue
!
      do  i=1,ndof  ! BD
      do  j=1,numnp
        k = id(i,j)
        !!write(*,*) 'j=',j,', id=', id(i, j), ', brhs=',  brhs (k)
      end do
      end do
!
      return
      end subroutine
!
!**** new **********************************************************************
!
      subroutine dirichletConditions(id,d,f,ndof,numnp,nlvect)
!
!.... program to compute displacement boundary conditions
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: ndof, numnp, nlvect
      integer*4:: id(ndof,*)
      real*8  :: d(ndof,*),f(ndof,numnp,*)
!
      integer*4:: i, j, k, lv
      real*8  :: val
!
      write(*,*) " em subroutine dirichletConditions(id,d,f,ndof,numnp,nlvect)"
!
      do 300 i=1,ndof
!
            do 200 j=1,numnp
!
            k = id(i,j)
            if (k.gt.0) go to 200
       !     if (k.lt.0) go to 200
            val = 0.d0
                  do 100 lv=1,nlvect
                  val = val + f(i,j,lv)
100               continue
!
            d(i,j) = val
            !write(*,*) i, j, d(i,j)
!
  200       continue
!
 300  continue
!     stop
      write(*,*) " fim de dirichletConditions(id,d,f,ndof,numnp,nlvect)"
      return
      end subroutine
!

!**** new **********************************************************************
      subroutine btod(id,d,brhs,ndof,numnp)
!
!.... program to perform transfer from r.h.s. to displacement array
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: ndof, numnp
      integer*4:: id(ndof,*)
      real*8  :: d(ndof,*),brhs(*)
!
      integer*4:: i, j, k
!
      do 200 i=1,ndof
!
         do 100 j=1,numnp
         !k = id(i,j)
         k = abs(id(i,j))
         if (k.gt.0) then 
             d(i,j) = brhs(k)
         end if
  100    continue
!
  200    continue
!
      return
      end subroutine btod
!
!**** new **********************************************************************
!
      subroutine pivots(a,idiag,neq,nsq,iecho,*)
!
!.... program to determine the number of zero and negative terms in
!        array d of factorization a = u(transpose) * d * u
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                
      real*8  :: a(*)                                                       
      integer*4:: idiag(*)
      integer*4:: neq 
      integer*4:: nsq, iecho
!
      integer*4:: iz, in, n, i
!
      iz = 0
      in = 0
!
      do 100 n=1,neq
      i = idiag(n)
      if (a(i).eq.0.) iz = iz + 1
      if (a(i).lt.0.) in = in + 1
  100 continue
!
      write(iecho,1000) nsq,iz,in
!
      return 1
!
 1000 format(' ',&
     ' zero and/or negative pivots encountered                ', ///5x,&
     ' time sequence number   . . . . . . . . . . . (nsq  ) = ',i10//5x,&
     ' number of zeroes . . . . . . . . . . . . . . . . . . = ',i10//5x,&
     ' number of negatives  . . . . . . . . . . . . . . . . = ',i10//5x)
!
      end subroutine
!
!**** new **********************************************************************
!
      subroutine btdb(elstif,b,db,nee,nrowb,nstr)
!
!.... program to multiply b(transpose) * db taking account of symmetry
!        and accumulate into element stiffness matrix
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: nee,nrowb
      real*8  :: elstif(nee,*),b(nrowb,*),db(nrowb,*)
      integer*4:: nstr
!
      integer*4:: i,j
!
      do 200 j=1,nee
!
      do 100 i=1,j
      elstif(i,j) = elstif(i,j) + coldot(b(1,i),db(1,j),nstr)   
  100 continue
!
  200 continue
!
      return
      end subroutine

!**** new **********************************************************************
      function coldot(a,b,n)
!
!.... program to compute the dot product of vectors stored column-wise
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: n
      real*8  :: a(n),b(n)

      real*8  :: coldot
      integer*4:: i

      real*8  :: dot_product
      real*8  :: ddot
    
!
      coldot = 0.0d0
!
      do 100 i=1,n
        coldot = coldot + a(i)*b(i)
  100 continue

!     coldot = dot_product(a,b) ! intrinsec fortran funtion
!     coldot = ddot(n,a,1,b,1)    ! external blas function 

!
      return
      end function
!**** new **********************************************************************
      subroutine matadd00(a,b,c,ma,mb,mc,m,n,iopt)
!
!.... program to add rectangular matrices
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: ma,mb,mc,m,n,iopt
      real*8  :: a(ma,*),b(mb,*),c(mc,*)
!
      integer*4:: i,j
!
      go to (1000,2000,3000),iopt
!
!.... iopt = 1, add entire matrices
!
 1000 do 1200 j=1,n
!
      do 1100 i=1,m 
      c(i,j) = a(i,j) + b(i,j)
 1100 continue
!
 1200 continue
      return
!
!.... iopt = 2, add lower triangular and diagonal elements
!
 2000 do 2200 j=1,n
!
      do 2100 i=j,m 
      c(i,j) = a(i,j) + b(i,j)
 2100 continue
!
 2200 continue
      return
!
!.... iopt = 3, add upper triangular and diagonal elements
!
 3000 do 3200 j=1,n
!
      do 3100 i=1,j 
      c(i,j) = a(i,j) + b(i,j)
 3100 continue
!
 3200 continue
      return
!
      end subroutine matadd00

!**** new **********************************************************************
      subroutine kdbc(eleffm,elresf,dl,lm,nee)
      implicit none
      integer*4:: nee, lm(*)
      real*8  :: eleffm(nee,*),elresf(*),dl(*)
      integer*4:: i,j
      real*8  :: val
      !write(*,*) " em kdbc .............. "
      do j=1,nee
        val=dl(j)
        if(lm(j).gt.0) cycle !if(val.eq.0.0d0)  cycle
        do i=1,nee
          elresf(i)=elresf(i)-eleffm(i,j)*val
        end do
      end do
      return
      end subroutine kdbc
!**** new **********************************************************************
      subroutine kdbcF(eleffm,elresf,dl,lm,nee)
      implicit none
      integer*4:: nee, lm(*)
      real*8  :: eleffm(nee,*),elresf(*),dl(*)
      integer*4:: i,j
      real*8  :: val
      do j=1,nee
        val=dl(j)
        if(lm(j).gt.0) cycle !    if(val.eq.0.0d0) cycle 
        do i=1,j-1
          elresf(i)=elresf(i)-eleffm(i,j)*val
        end do
        do i=j+1,nee
          elresf(i)=elresf(i)-eleffm(i,j)*val
        end do
      end do
      do j=1,nee
        val=dl(j)
        if(lm(j).gt.0) cycle 
        elresf(j)=val
        eleffm(1:nee,  j)=0.0
        eleffm(j,  1:nee)=0.0
        eleffm(j,  j)=1.0
        !eleffm(1:j-1,  j)=0.0
        !eleffm(j+1:nee,j)=0.0
        !eleffm(j,  1:j-1)=0.0
        !eleffm(j,j+1:nee)=0.0
      end do
      return
      end subroutine kdbcF
 end module
