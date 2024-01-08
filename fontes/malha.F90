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

      module mMalha

      implicit none
      
      real*8,  allocatable :: x(:,:) 
      integer*4, allocatable :: conecNodaisElem(:,:)
      integer*4, allocatable :: listaDosElemsPorNo(:,:)
      integer*4, pointer  :: listaDosElemsPorNoCSR(:,:) => null()
      integer*4:: nsd, numel, numnp, nen
      integer*4:: numConexoesPorElem 


!
!     funcoes e subrotinas
      public :: nconec
      public :: genfl, genfl1
      public :: gensh, gensh1, gensh2, gensh3, igen

!
       contains
!
!**** new **********************************************************************
!
      subroutine nconec(nen,numnp,numel,conecElem,listaDosElemsPorNo)
!
      implicit none
!     
!     Objetivo: cria a matriz listaDosElemsPorNo: .... indices dos elementos que possuem o no
!
      integer*4 :: nen,numnp,numel
      integer*4, dimension(nen,numel)  :: conecElem
      integer*4, dimension(nen,numnp)  :: listaDosElemsPorNo
      integer*4 :: i,j,no,nel,ni,ns


      character(len=22) :: nomeArq="listaDosElemsPorNo.dat"
      integer*4:: ilistaDosElemsPorNo=30
      logical :: existe


      inquire (file=nomeArq,exist=existe)
      open(UNIT=ilistaDosElemsPorNo,FILE=nomeArq)
!      open(UNIT=ilistaDosElemsPorNo,FILE=nomeArq,FORM='UNFORMATTED')

      if(existe.eqv..false.) then
!     
!      nmax=max(nelx,nely)
!      nmin=min(nelx,nely)
!      nmax=nelx*nely
!
      listaDosElemsPorNo=0

!
!$OMP PARALLEL DO PRIVATE(J, NO, NEL)
      do i=1,numnp
!     
!     centra no noh
!
!         ni=i*nmax
!         if(ni.lt.0) ni=1
!     
         ni=1
!
!         ns=i+nmax
!         if(ns.gt.numel) ns=numel
!
         ns=numel
!
         do nel=ni,ns
!     
            do j=1,nen
!     
               no=conecElem(j,nel)
               if(no.eq.i) listaDosElemsPorNo(j,i) = nel
!     
            end do
!     
         end do
!     
      end do
!$OMP END PARALLEL DO

       do i=1, numnp
       write(ilistaDosElemsPorNo) listaDosElemsPorNo(1:nen,i)
       end do
!      write(ilistaDosElemsPorNo,*), listaDosElemsPorNo ! teste com arq.formatado
!      write(ilistaDosElemsPorNo), listaDosElemsPorNo ! teste com arq.nao-formatado


      else
          do i=1, numnp
            read(ilistaDosElemsPorNo,*) (listaDosElemsPorNo(j,i),j=1,nen)
          end do
!       read(ilistaDosElemsPorNo) listaDosElemsPorNo
      end if

      close(ilistaDosElemsPorNo)
    
      end subroutine
!
      subroutine genfl(a,nra,iin)   
      use mGlobaisEscalares  
!                                                                       
!.... program to read and generate floating-point nodal data            
!                                                                       
!         a       = input array                                         
!         nra     = number of rows in a (le.6)                          
!         n       = node number                                         
!         numgp   = number of generation points                         
!         ninc(i) = number of increments for direction i                
!         inc(i)  = increment for direction i                           
!                                                                       
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                        
      integer*4:: nra, iin                                               
      real*8  :: a(nra,*)
!
      integer*4:: n,numgp,ninc(3),inc(3)   
      integer*4:: i, j, m, mgen
      integer*4, parameter :: dim1 = 6, dim2 = 20
      real*8  :: temp(dim1,dim2)
      !real*8  :: temp(6,20)
!                                                                       
  100 continue                                                          
      read(iin,1000) n,numgp,(temp(i,1),i=1,nra)      
!       write(*,1000) n,numgp,(temp(i,1),i=1,nra)      

      if (n.eq.0) return                                                
 !     call move(a(1,n),temp,nra)           
      a(1:nra,n) = temp(1:nra,1)

      if (numgp.ne.0) then
         do 200 j=2,numgp
!                                                                       
         read(iin,1000) m,mgen,(temp(i,j),i=1,nra)    

!        if (mgen.ne.0) call move(temp(1,j),a(1,m),nra) 
         if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m) 


!                                                                       
  200    continue                               
         read(iin,2000) (ninc(i),inc(i),i=1,3)
         call genfl1(a,nra, temp, n, numgp, ninc, inc)                                                
      endif
      go to 100                                                         
!                                                                       
! 1000 format(2i10,6f10.0)                                                
 1000 format(2i10,6f10.0)                                                
 2000 format(16i10)                                                      
!                                                                       
      end  subroutine
!
!******************************************************************************
!
      subroutine genfl1(a,nra, temp, n, numgp, ninc, inc)     
      use mGlobaisEscalares                                      
!                                                                       
!.... program to generate floating-point nodal data 
!        via isoparametric interpolation         
!                                                                       
!         iopt = 1, generation along a line                             
!              = 2, generation over a surface                           
!              = 3, generation within a volume                            
!                                                                       
      implicit none                                          
!                                                                       
!.... remove above card for single-precision operation                  
!
      integer*4:: nra
      real*8  :: a(nra,*)
      !real*8  :: temp(6,20)
      integer*4, parameter :: dim1 = 6, dim2 = 20, um = 1
      real*8  :: temp(dim1,dim2)
      integer*4:: n,numgp,ninc(3),inc(3)                 
!
      real*8  :: sh(20), dr, ds, dt, r, s, t
      integer*4:: iopt, ni, nj, nk, ii, jj, kk, i, j, k
!
      iopt = 3                                                          
      if (ninc(3).eq.0) iopt = 2                                        
      if (ninc(2).eq.0) iopt = 1                                        
!                                                                       
      dr = zero                                                           
      ds = zero                                                           
      dt = zero                                                           
!                                                                       
      if (ninc(1).ne.0) dr = two/ninc(1)                                
      if (ninc(2).ne.0) ds = two/ninc(2)                                
      if (ninc(3).ne.0) dt = two/ninc(3)                                
!                                                                       
      ii = ninc(1)+1                                                    
      jj = ninc(2)+1                                                    
      kk = ninc(3)+1                                                    
!                                                                       
      ni = n                                                            
      nj = n                                                            
      nk = n                                                            
!                                                                       
      t = -one                                                          
      do 300 k=1,kk                                                     
!
      s = -one                                                          
      do 200 j=1,jj                                                     
!
      r = -one                                                          
      do 100 i=1,ii                                                     
!                                                                       
      call gensh(r,s,t,sh,numgp,iopt)                                   
      call multab(temp,sh,a(1,ni),dim1,dim2,nra,numgp,nra,um,um)               
      ni = ni + inc(1)                                                      
      r = r + dr                                                            
  100 continue                                                          
!                                                                       
      nj = nj + inc(2)                                                      
      ni = nj                                                             
      s = s + ds                                                            
  200 continue                                                          
!                                                                       
      nk = nk + inc(3)                                                      
      ni = nk                                                             
      t = t + dt                                                            
  300 continue                                                          
!                                                                       
      return                                                            
      end  subroutine
!
!****************************************************************************** 
!   
      subroutine gensh(r,s,t,sh,numgp,iopt)                             
!                                                                       
!.... program to call shape function routines         
!        for isoparametric generation         
!                                                                       
      implicit none                                      
!                                                                       
!.... modify above card for single-precision operation               
!               
      real*8  :: r, s, t, sh(*)                                                        
      integer*4:: numgp, iopt
!                                                                       
      go to (100,200,300),iopt                                                
!                                                                       
  100 call gensh1(r,sh,numgp)                                           
      return                                                            
!                                                                       
  200 call gensh2(r,s,sh,numgp)                                         
      return                                                            
!                                                                       
  300 call gensh3(r,s,t,sh,numgp)                                       
      return                                                            
!                                                                       
      end subroutine
!
!******************************************************************************
!
      subroutine gensh1(r,sh,n)  
      use mGlobaisEscalares                                       
!                                                                       
!.... program to compute 1d shape functions           
!        for isoparametric generation                     
!                                                                       
      implicit none                                          
!                                                                       
!.... modify above card(s) for single-precision operation               
!                                                                       
      real*8  :: r, sh(*)                                                   
      integer*4:: n
!                                                                       
      sh(2) = pt5*r                                                       
      sh(1) = pt5 - sh(2)                                                   
      sh(2) = pt5 + sh(2)                                                   
      if (n.eq.3) then
         sh(3) = one - r*r                                                     
         sh(1) = sh(1) - pt5*sh(3)                                             
         sh(2) = sh(2) - pt5*sh(3)                                             
      endif
!                                                                       
      return                                                            
      end subroutine
!
!******************************************************************************
!
      subroutine gensh2(r,s,sh,n)     
      use mGlobaisEscalares                                  
!
!.... program to compute 2d shape functions 
!        for isoparametric generation    
!                                                                       
      implicit none                                         
!
!.... modify above card for single-precision operation               
!                                                                       
      real*8  :: r, s, sh(*)                                                   
      integer*4:: n    
!
      real*8  :: r1, r2, r3, s1, s2, s3
!
      r2 = pt5*r                                                          
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s                                                          
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      sh(1) = r1*s1                                                       
      sh(2) = r2*s1                                                       
      sh(3) = r2*s2                                                       
      sh(4) = r1*s2                                                       
      if (n.eq.4) return                                                
!                                                                       
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      sh(5) = r3*s1                                                       
      sh(6) = s3*r2                                                       
      sh(7) = r3*s2                                                       
      sh(8) = s3*r1                                                       
      sh(1) = sh(1) - pt5*(sh(5) + sh(8))
      sh(2) = sh(2) - pt5*(sh(6) + sh(5))
      sh(3) = sh(3) - pt5*(sh(7) + sh(6))
      sh(4) = sh(4) - pt5*(sh(8) + sh(7))
!                                                                       
      return
      end subroutine
!
!******************************************************************************
!
      subroutine gensh3(r,s,t,sh,n) 
      use mGlobaisEscalares                                    
!                                                                       
!.... program to compute 3d shape functions            
!        for isoparametric generation   
!                                                                       
      implicit none                                         
!                                                                       
!.... modify above card for single-precision operation               
!                                                                       
      real*8  :: r, s, t, sh(*)                                                   
      integer*4:: n    
!
      real*8  :: r1, r2, r3, rs1, rs2, rs3, rs4
      real*8  :: s1, s2, s3, t1, t2, t3                                              
!                                                                       
      r2 = pt5*r
      r1 = pt5 - r2                                                         
      r2 = pt5 + r2                                                         
      s2 = pt5*s
      s1 = pt5 - s2                                                         
      s2 = pt5 + s2                                                         
      t2 = pt5*t
      t1 = pt5 - t2                                                         
      t2 = pt5 + t2                                                         
!                                                                       
      rs1 = r1*s1                                                         
      rs2 = r2*s1                                                         
      rs3 = r2*s2                                                         
      rs4 = r1*s2                                                         
      sh(1) = rs1*t1                                                      
      sh(2) = rs2*t1                                                      
      sh(3) = rs3*t1                                                      
      sh(4) = rs4*t1                                                      
      sh(5) = rs1*t2                                                      
      sh(6) = rs2*t2                                                      
      sh(7) = rs3*t2                                                      
      sh(8) = rs4*t2                                                      
      if (n.eq.8) return                                                 
!                                                                       
      r3 = one - r*r                                                        
      s3 = one - s*s                                                        
      t3 = one - t*t                                                        
      sh(17) = t3*rs1                                                     
      sh(18) = t3*rs2                                                     
      sh(19) = t3*rs3                                                     
      sh(20) = t3*rs4                                                     
      rs1 = r3*s1                                                         
      rs2 = s3*r2                                                         
      rs3 = r3*s2                                                         
      rs4 = s3*r1                                                         
      sh( 9) = rs1*t1                                                     
      sh(10) = rs2*t1                                                     
      sh(11) = rs3*t1                                                     
      sh(12) = rs4*t1                                                     
      sh(13) = rs1*t2                                                     
      sh(14) = rs2*t2                                                     
      sh(15) = rs3*t2                                                     
      sh(16) = rs4*t2                                                     
!                                                                       
      sh(1) = sh(1) - pt5*(sh( 9) + sh(12) + sh(17))
      sh(2) = sh(2) - pt5*(sh( 9) + sh(10) + sh(18))
      sh(3) = sh(3) - pt5*(sh(10) + sh(11) + sh(19))
      sh(4) = sh(4) - pt5*(sh(11) + sh(12) + sh(20))
      sh(5) = sh(5) - pt5*(sh(13) + sh(16) + sh(17))
      sh(6) = sh(6) - pt5*(sh(13) + sh(14) + sh(18))
      sh(7) = sh(7) - pt5*(sh(14) + sh(15) + sh(19))
      sh(8) = sh(8) - pt5*(sh(15) + sh(16) + sh(20))
!                                                                       
      return                                                            
      end subroutine                                                              
!**** new **********************************************************************
      subroutine igen(ia,m,iin)
      use mGlobaisEscalares
!
!.... program to read and generate integer*4, nodal data
!
!        ia = input array
!         m = number of rows in ia
!         n = node number
!        ne = end node in generation sequence
!        ng = generation increment
!    
      integer*4:: m, ia(m,*), iin
!
      integer*4:: ib(13)
      integer*4:: n, ne, ng
      integer*4:: i
!
  100 continue
      read(iin,1000) n,ne,ng,(ib(i),i=1,m)
!      write(*,1000) n,ne,ng,(ib(i),i=1,m)
      if (n.eq.0) return

      if (ng.eq.0) then
         ne = n
         ng = 1
      else
         ne = ne - mod(ne-n,ng)
      endif
!
      do 200 i=n,ne,ng
!     call imove(ia(1,i),ib,m)
      ia(1:m,i)=ib(1:m)
  200 continue
!
      go to 100
!
 1000 format(16i10)
      end subroutine
!     

!**** new **********************************************************************
      subroutine local(conectElem,x,xl,nen,nrowx,nrowxl)
!
!.... program to localize a global array
!
!        note: it is assumed nrowxl.le.nrowx
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: conectElem(*)
      integer*4:: nrowx, nrowxl, nen
      double precision :: x(nrowx,*),xl(nrowxl,*)
!
      integer*4:: i, j, node
!
      do 200 j=1,nen
      node = conectElem(j)
!
      do 100 i=1,nrowxl
      xl(i,j)= x(i,node)
  100 continue
!
  200 continue
!
      return
      end subroutine

!**** new **********************************************************************
      subroutine multab(a,b,c,ma,mb,mc,l,m,n,iopt)
!
!.... program to multiply two matrices
!
!        l = range of dot-product index
!        m = number of active rows in c
!        n = number of active columns in c
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                                                                       
      integer*4:: ma,mb,mc,l,m,n,iopt
      real*8  :: a(ma,*),b(mb,*),c(mc,*)
!
      integer*4:: i,j
!
      go to (1000,2000,3000,4000),iopt
!
!.... iopt = 1, c(i,j) = a(i,k)*b(k,j) , (c = a * b)
!
 1000 do 1200 i=1,m
!
      do 1100 j=1,n
      c(i,j) = rcdot(a(i,1),b(1,j),ma,l)
 1100 continue
!
 1200 continue
      return
!                                            t
!.... iopt = 2, c(i,j) = a(k,i)*b(k,j) (c = a  * b)
!
 2000 do 2200 i=1,m
!
      do 2100 j=1,n
      c(i,j) = dot_product (a(:,i),b(:,j)) !,l)
      !c(i,j) = coldot(a(1,i),b(1,j),l)
 2100 continue
!
 2200 continue
      return
!                                                t
!.... iopt = 3, c(i,j) = a(i,k)*b(j,k) (c = a * b )
!
 3000 do 3200 i=1,m
!
      do 3100 j=1,n
      c(i,j) = rowdot(a(i,1),b(j,1),ma,mb,l)
 3100 continue
!
 3200 continue
      return
!                                            t    t
!.... iopt = 4, c(i,j) = a(k,i)*b(j,k) (c = a  * b )
!
 4000 do 4200 i=1,m
!
      do 4100 j=1,n
      c(i,j) = rcdot(b(j,1),a(1,i),mb,l)
 4100 continue
!
 4200 continue
!
      return
      end subroutine

!**** new **********************************************************************
      function rcdot(a,b,ma,n)
!
!.... program to compute the dot product of a vector stored row-wise
!        with a vector stored column-wise
!
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                                                                       
      integer*4:: ma, n
      real*8  :: a(ma,*),b(*)
!
      real*8  :: rcdot
      integer*4:: i
!
      rcdot = 0.0d0
!
      do 100 i=1,n
      rcdot = rcdot + a(1,i)*b(i)
  100 continue
!
      return
      end function
!**** new **********************************************************************
      function rowdot(a,b,ma,mb,n)
!
!.... program to compute the dot product of vectors stored row-wise
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer*4:: ma, mb, n
      real*8  :: a(ma,*),b(mb,*)
!
      real*8  :: rowdot
      integer*4:: i
!
      rowdot = 0.0d0
!
      do 100 i=1,n
      rowdot = rowdot + a(1,i)*b(1,i)
  100 continue
!
     return
   end function
!
!=======================================================================
!     
      subroutine criarListaVizinhos(nen,numnp,numel,conecElem,listaDosElemsPorNo)
      implicit none
!     
!     Objetivo: cria a matriz listaDosElemsPorNo: .... indices dos elementos que possuem o no
!
      integer*4 :: nen,numnp,numel
      integer*4, dimension(nen,numel)  :: conecElem
      integer*4, dimension(nen,numnp)  :: listaDosElemsPorNo
      integer*4:: no,nel,l

      do nel=1, numel
         do l=1, nen
            no=conecElem(l,nel)
            listaDosElemsPorNo(l,no) = nel
         end do
      end do
    
      end subroutine

      end module
