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
      module mLeituraEscrita
!
      implicit none
!
      integer*4:: iin, iecho, icoords, iconects, iconectsL
      integer*4:: ignuplotPotencial, ignuplotFluxo, iparaview
!
!     funcoes e subrotinas
      public :: echo, leituraValoresCondContorno, leituraGeracaoCoordenadas
      public :: printf, printd, printp, prntel
      public :: printResultado, prtgnup, prtvB
      public :: escreverArqParaview, escreverPontosNodais, escreverConectividades
      public :: escreverTiposElementos, escreverEscalaresNodais
      public :: leituraGeracaoConectividades
      public :: genel1
!
      contains

!**** new **********************************************************************
      subroutine echo

      implicit none
!
!.... program to echo input data
!
      character*4 ia(50)
      integer*4:: iech, i
!
!     cabeçalho
      write(iecho,500)

      read(iin,1000) iech
      if (iech.eq.0) return
!
      write(iecho,2000) iech
      backspace iin
!
      do 100 i=1,100000
      read(iin,3000,end=200) ia
      if (mod(i,50).eq.1) write(iecho,4000)
      write(iecho,5000) ia

  100 continue
!
  200 continue
      rewind iin
      read(iin,1000) iech
!
      return
!
 500  format('programa de elementos finitos em fortran 90 baseado em:',// &
      'The Finite Element Method, Hughes, T. J. R., (2003)'//)
 1000 format(16i10)
 2000 format('1',' i n p u t   d a t a   f i l e               ',  //5x,&
     ' echo print code . . . . . . . . . . . . . . (iecho ) = ',i10//5x,&
     '    eq. 0, no echo of input data                        ',   /5x,&
     '    eq. 1, echo input data                              ',   ///)
 3000 format(20a4)
 4000 format(' ',8('123456789*'),//)
 5000 format(' ',20a4)

      end subroutine echo

!******************************************************************************
      subroutine leituraGeracaoConectividades(conecElem,mat,nen, iin)   
!
      use mGlobaisEscalares 
!
      implicit none
!                                                                       
!.... program to read and generate element node and material numbers    
!                                                                       
!         conecElem(nen,numel) = element node numbers                         
!         mat(numel)     = element material numbers                     
!         nen            = number of element nodes (le.27)              
!         n              = element number                               
!         ng             = generation parameter                         
!         nel(i)         = number of elements in direction i            
!         incel(i)       = element number increment for direction i     
!         inc(i)         = node number increment for direction i        
!   
      integer*4:: nen, iin
      integer*4:: conecElem(nen,*),mat(*)
!
      integer*4:: m,ng, i, itemp(27)
      integer*4:: n,nel(3),incel(3),inc(3)        
!                        
!                                                                       
  100 continue                                                          
      read(iin,1000) n,m,(itemp(i),i=1,nen),ng     
                
      if (n.eq.0) return                                                
      conecElem(1:nen,n)=itemp(1:nen)                                    
      mat(n)=m                                                          
      if (ng.ne.0) then
!                                                                       
!....... generate data                                                     
!                                                                       
         read(iin,1000) (nel(i),incel(i),inc(i),i=1,3)                     
         call genel1(conecElem,mat,nen,n,nel,incel,inc)                     
      endif
      go to 100                                                         
!                                                                       
 1000 format(16i10,10x,14i10)                                             
!                                                                       
      end subroutine                                                               
!******************************************************************************
      subroutine genel1(conecElem,mat,nen, n,nel,incel,inc) 
!                                                                       
!.... program to generate element node and material numbers             
!                                                                       
      integer*4 ::  nen, conecElem(nen,*),mat(*)                                       
      integer*4:: n,nel(3),incel(3),inc(3)                          

      integer*4:: i,j,k,ii,jj,kk,ie,le,je,ke
!                                                                       
!.... set defaults                                                      
!                                                                       
      call geneld                                                 
!                                                                       
!.... generation algorithm                                              
!                                                                       
      ie = n                                                            
      je = n                                                            
      ke = n                                                            
!                                                                      
      ii = nel(1)                                                       
      jj = nel(2)                                                       
      kk = nel(3)                                                       
!                                                                       
      do 300 k=1,kk                                                     
!
      do 200 j=1,jj                                                     
!
      do 100 i=1,ii                                                     
!                                                                       
      if (i.ne.ii) then
         le = ie                                                           
         ie = le + incel(1)                                                
         call geneli(conecElem(1,ie),conecElem(1,le),inc(1),nen)                       
         mat(ie) = mat(le)                                                 
      endif
  100 continue                                                          
!                                                                       
      if (j.ne.jj) then
         le = je                                                           
         je = le + incel(2)                                                
         call geneli(conecElem(1,je),conecElem(1,le),inc(2),nen)                       
         mat(je) = mat(le)                                                 
         ie = je                                                           
      endif
  200 continue                                                          
!                                                                       
      if (k.ne.kk) then
         le = ke                                                           
         ke = le + incel(3)                                                
         call geneli(conecElem(1,ke),conecElem(1,le),inc(3),nen)                       
         mat(ke) = mat(le)                                                 
         ie = ke
         je=ke                                                           
      endif
  300 continue                                                          
!                                                                      
      return                                                            
      contains 
!******************************************************************************
      subroutine geneld                                                 
!                                                                       
!.... program to set defaults for element node       
!        and material number generation                              
!                                                                       
      if (nel(1).eq.0) nel(1) = 1                                       
      if (nel(2).eq.0) nel(2) = 1                                       
      if (nel(3).eq.0) nel(3) = 1                                       
!                                                                       
      if (incel(1).eq.0) incel(1) = 1                                   
      if (incel(2).eq.0) incel(2) = nel(1)                              
      if (incel(3).eq.0) incel(3) = nel(1)*nel(2)                       
!                                                                       
      if (inc(1).eq.0) inc(1) = 1                                       
      if (inc(2).eq.0) inc(2) = (1+nel(1))*inc(1)                       
      if (inc(3).eq.0) inc(3) = (1+nel(2))*inc(2)                       
!                                                                       
      return                                                            
      end subroutine
!
!******************************************************************************
!
      subroutine geneli(conecElem2,conecElem1,inc,nen)                              
!                                                                       
!.... program to increment element node numbers                         
!                                                                       
      integer*4 ::  conecElem1(*),conecElem2(*)
      integer*4:: inc, nen

      integer*4:: i
!
      do 100 i=1,nen                                                    
      if (conecElem1(i).eq.0) then
         conecElem2(i) = 0
      else
         conecElem2(i) = conecElem1(i) + inc                         
      endif
  100 continue                                                          
!                                                                       
      return                                                            
      end  subroutine
      end subroutine genel1                                                              
!
!**** new **********************************************************************
      subroutine leituraGeracaoCoordenadas(x, nsd, numnp, iin, icoords, iprtin)
      use mMalha, only: genfl
!
!.... program to read, generate and write coordinate data
!
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4, intent(in)   :: nsd, numnp, iin, icoords, iprtin
      real*8, intent(inout) ::  x(nsd,*)
!
      integer*4:: i, n
!      
      call genfl(x,nsd,iin)
!
      if (iprtin.eq.1) return
!

         write(icoords,*) "# Coordenadas "
         do n=1,numnp
            write(icoords,2000) n,(x(i,n),i=1,nsd)   
         end do
!
      return
!
!  1000 format('1',' n o d a l   c o o r d i n a t e   d a t a '///5x,&
!      ' node no.',3(13x,' x',i1,' ',:)//)
 2000 format(6x,i12,10x,3(1pe15.8,2x))
      end subroutine
!**** new **********************************************************************
      subroutine leituraCodigosCondContorno(id, ndof, numnp, neq, iin, iecho,iprtin)
!
!.... program to read, generate and write boundary condition data
!        and establish equation numbers
!
      use mMalha, only: igen

      integer*4, intent(in) :: ndof, numnp, iin, iecho, iprtin
      integer*4, intent(inout) :: neq
      integer*4, intent(inout) :: id(ndof,numnp)
!
      integer*4:: nn, n, i
      logical pflag
!
      id(1:ndof,1:numnp) = 0
      call igen(id,ndof, iin)
!
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
         pflag = .false.
!
         do 100 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
    !       !print*, i, 'id =',  id(i,n)
  100    continue
!
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
            write(iecho,2000) n,(id(i,n),i=1,ndof)
         endif
  200    continue
      endif
!
!.... establish equation numbers
!
      neq = 0
!
      do 400 n=1,numnp
!
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif

!
  300 continue
!
  400 continue

!
      return
!
 1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n & 
              c o  d e s'/// &
      5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
!
      end subroutine
!**** new **********************************************************************
      subroutine leituraCodigosCondContornoBidu(id, ndof, numnp, neq, iin, iecho,iprtin)
!
!.... program to read, generate and write boundary condition data
!        and establish equation numbers
!
      use mMalha, only: igen

      integer*4, intent(in) :: ndof, numnp, iin, iecho, iprtin
      integer*4, intent(inout) :: neq
      integer*4, intent(inout) :: id(ndof,numnp)
!
      integer*4:: nn, n, i
      logical pflag
!
      write(*,*) "subroutine leituraCodigosCondContornoB(id, ndof, numnp, ..."
      id(1:ndof,1:numnp) = 0
      call igen(id,ndof, iin)
!
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
         pflag = .false.
!
         do 100 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
           !!print*, i, 'id =',  id(i,n)
  100    continue
!
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
            write(iecho,2000) n,(id(i,n),i=1,ndof)
         endif
  200    continue
      endif
!
!.... establish equation numbers
!
      neq = 0
!
      do 400 n=1,numnp
!
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif

!
  300 continue
!
  400 continue

       do i =1, numnp    !BD
           !print*,'i=', i, 'id =',  id(1,i)
       end do
    i = 3
    i = i+1; id(1,i) = 4 + 1
    i = i+1; id(1,i) = 6 + 1
    i = i+1; id(1,i) = 11 + 1
    i = i+1; id(1,i) = 1 + 1
    i = i+1; id(1,i) = 5 + 1
    i = i+1; id(1,i) = 9 + 1
    i = i+1; id(1,i) = 10 + 1
    i = i+1; id(1,i) = 2 + 1
    i = i+1; id(1,i) = 3 + 1
    i = i+1; id(1,i) = 8 + 1
    i = i+1; id(1,i) = 7 + 1
    i = i+1; id(1,i) = 0 + 1

       do i =1, numnp    !BD
           !print*,'i=', i, 'id =',  id(1,i)
       end do
!
      return
!
 1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n & 
              c o  d e s'/// &
      5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
!
      end subroutine
!**** new **********************************************************************
      subroutine leituraValoresCondContorno(f,ndof,numnp,j,nlvect,iprtin)
!
!.... program to read, generate and write nodal input data
!
!        f(ndof,numnp,nlvect) = prescribed forces/kinematic data (j=0)
!                             = nodal body forces(j=1)
!
      use mMalha, only : genfl
      implicit none
!
!.... remove above card for single-precision operation
!
      integer*4:: ndof, numnp, j, nlvect, iprtin
      real*8 f(ndof,numnp,nlvect)

      logical lzero
      integer*4 nlv
      character(len=30) :: rotulo

!
!     call clear(f,nlvect*numnp*ndof)
      f=0.0
      

      do 100 nlv=1,nlvect
      call genfl(f(1,1,nlv),ndof,iin)
      call ztest(f(1,1,nlv),ndof*numnp,lzero)
!
      if (iprtin.eq.0) then
!
         if (lzero) then
            if (j.eq.0) write(iecho,1000) nlv
            if (j.eq.1) write(iecho,2000)
         else
            if (j.eq.0) call printf(f,ndof,numnp,nlv)
!
            if (j.eq.1) then
               rotulo=" n o d a l  b o d y  f o r c e s  "
               call printd (rotulo, f,ndof,numnp,iecho)
            end if
!
         endif
      endif
!
  100 continue
!
      return
 1000 format('1'//,' there are no nonzero prescribed forces and ',&
         'kinematic boundary conditions for load vector number ',i10)
 2000 format('1'//,' there are no nonzero nodal body forces')
      end subroutine

!**** new **********************************************************************
      subroutine printf(f,ndof,numnp,nlv)
!
!.... program to print prescribed force and boundary condition data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer*4 :: ndof, numnp, nlv
      real*8 :: f(ndof,numnp,*)
!
      logical lzero
      integer*4:: nn, n, i
!
      nn = 0
!
      do 100 n=1,numnp
      call ztest(f(1,n,nlv),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) write(iecho,1000) nlv,(i,i=1,ndof)
         write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
      endif
  100 continue
!
      return
!
 1000 format('1',&
     ' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',&
     '  b o u n d a r y   c o n d i t i o n s'//5x,&
     ' load vector number = ',i10///5x,&
     ' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i10,10x,6(1pe15.8,2x))
      end subroutine
      
!**** new **********************************************************************
      subroutine printd(name,dva,ndof,numnp,icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer*4:: ndof,numnp, icode
      character (LEN=*) ::  name
      real*8 dva(ndof,*)
!
      logical lzero
      integer*4 nn, n, i
!
      nn = 0
!
      do 100 n=1,numnp
      call ztest(dva(1,n),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) &
           write(icode,1000) name,(i,i=1,ndof)
         write(icode,2000) n,(dva(i,n),i=1,ndof)
      endif
  100 continue
!
      return
!
 1000 format('1',11a4//1x,'node',6(11x,'dof',i1)/)
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end subroutine
!
!**** new **********************************************************************
!
      subroutine printp(a,idiag,neq,nsq,*)
      use mGlobaisEscalares
!
!.... program to print array d after Crout factorization 
!        a = u(transpose) * d * u
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer*4:: neq, nsq
      real*8 :: a(*)
      integer*4:: idiag(*)
!
      integer*4:: n, i
!
      do 100 n=1,neq
      if (mod(n,50).eq.1) write(iecho,1000) nsq
      write(iecho,1000)
      i = idiag(n)
      write(iecho,2000) n,a(i)
  100 continue
!
      return 1
!
 1000 format('1',' array d of factorization',/&
     ' a = u(transpose) * d * u ',                                //5x,&
     ' time sequence number   . . . . . . . . . . . . (nsq) = ',i10//5x)
 2000 format(1x,i10,4x,1pe20.8)
      end subroutine
!
!**** new **********************************************************************
!
      subroutine prntel(mat,conectElem,nen,numel,tipo)
      implicit none
!
!.... program to print data for element with "nen" nodes
!
!        note: presently the label formats are limited to
!              elements with one to nine nodes
!
      integer*4:: nen, numel
      integer*4:: mat(*),conectElem(nen,*)
      integer*4:: tipo
!
      integer*4 n, i
!
      if(tipo==1) then
	write(iconects,*) "# Conectividades nodais"
	do n=1,numel
	  write(iconects,2000) n,mat(n),(conectElem(i,n),i=1,nen)
	end do
      end if

      if(tipo==2) then
	write(iconectsL,*) "# Conectividades ladais"
	do  n=1,numel
	  write(iconectsL,3000) n,mat(n),(conectElem(i,n),i=1,nen)
	end do
      end if
!
      return
!
 1000 format(///10x,&
     ' e l e m e n t   d a t a',//1x,&
     ' element   material ',9('node ',i1,1x),/1x,&
     '  number    number',//)
 2000 format(1x,i10,9(2x,i10))
 3000 format(1x,i10,7(2x,i10))
      end subroutine

!**** new **********************************************************************
      subroutine printResultado(dva, ndof, numnp, inicio, fim, icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer*4:: ndof, numnp, inicio, fim, icode
      real*8  :: dva(ndof,numnp)
!
      integer*4:: n, i
!
      write(icode,*) "# Solucao"
      do 100 n=inicio,fim
         write(icode,2000) n,(dva(i,n),i=1,ndof)
         !write(*,*) n,(dva(i,n),i=1,ndof)
  100 continue
!
      return
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end subroutine


!**** new **********************************************************************
      subroutine prtgnup(name,x,dva,nsd,ndof,numnp,icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer*4:: nsd, ndof, numnp, icode
      character(len=*) :: name
      real*8 :: x(nsd,*),dva(ndof,*)
!
      integer*4:: n
      real*8  :: tmp, t1, t2
!
      call timing(t1)
      write(icode,'(a)') name
      tmp=0.0
      do 100 n=1,numnp
         write(icode,2000) x(1:nsd,n), dva(1:ndof,n)
         if(x(nsd,n+1)>x(nsd,n).and.n<numnp) write(icode,*) ""
  100 continue
!       do 100 n=1,numnp
!          write(icode,2000) (x(j,n),j=1,nsd), (dva(i,n),i=1,ndof)
!   100 continue
!
      call timing(t2)
      write(*,*) 'tempo de escrita em arquivo = ', t2 - t1
      return
!
 2000 format(6(1pe13.6,2x))
      end subroutine
!**** new **********************************************************************
    subroutine escreverArqParaview(campo, dim1, dim2, nen, conectElem, tipo, rotulo, tamRot)
    use mMalha, only: x, nsd, numel

    implicit none
    integer*4, intent(in) :: dim1, dim2
    double precision, intent(in) :: campo(dim1, dim2)
    integer*4:: nen
    integer*4:: conectElem(nen,numel)
    integer*4:: tipo  !tipo=1 para elemento, e tipo=2 para no
    integer*4:: tamRot

    character(len=tamRot) :: rotulo
    
    write(iparaview,'(a)')'# vtk DataFile Version 3.0'
    write(iparaview,'(a)')'vtk output'
    write(iparaview,'(a)')'ASCII'
    write(iparaview,'(a)')'DATASET UNSTRUCTURED_GRID'
    write(iparaview,'(a,i10,a)')'POINTS', dim2,' float '

    call escreverPontosNodais  (x, dim2, nsd)
! 
    write(iparaview,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
    call escreverConectividades(conectElem, numel, nen, nsd)
! 
    write(iparaview,'(a,i10)')'CELL_TYPES ', numel
    call escreverTiposElementos(numel, nsd)
! 

    if(tipo==1) write(iparaview,'(a,i10)')'CELL_DATA ', numel

    if(tipo==2) write(iparaview,'(a,i10)')'POINT_DATA',  dim1*dim2

    write(iparaview,'(3a)')'SCALARS ', trim(rotulo), ' float '
    write(iparaview,'(a)')'LOOKUP_TABLE default'

    call escreverEscalaresNodais(campo, dim1, dim2,rotulo,tamRot)

    end subroutine escreverArqParaview

!**** new **********************************************************************
      subroutine escreverPontosNodais  (coords, numnp, nsd)
      implicit none
      integer*4, intent(in) :: numnp, nsd
      real*8,  intent(in) :: coords(nsd,numnp)
!
      real*8  :: coordZ = 0.0 
      integer*4:: d, i
!
      if(nsd==2) then
	  do i=1,numnp
	      write(iparaview,'(3(1x, 1pe15.8))') (coords(d,i),d=1,nsd), coordZ 
	  end do
      end if

      if(nsd==3) then
	  do i=1,numnp
	      write(iparaview,'(3(1x, 1pe15.8))') (coords(d,i),d=1,nsd)
	  end do
      end if
      end subroutine escreverPontosNodais


!**** new **********************************************************************
      subroutine escreverConectividades(conectElem, numel, nen, nsd)
      implicit none
      integer*4, intent(in)  :: numel, nen, nsd
      integer*4, intent(in)  :: conectElem(nen,numel)
!
      integer*4 n, i
!
      if(nsd==2) then
	do  n=1,numel
	  write(iparaview,'(i10,9(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
	end do
      end if

      if(nsd==3) then
	do  n=1,numel
	  write(iparaview,'(i10,18(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
	end do
      end if

 end subroutine escreverConectividades

!**** new **********************************************************************
      subroutine escreverTiposElementos(numel, nsd)
      implicit none
      integer*4, intent(in)   :: numel, nsd
!
      integer*4:: i
!
      if(nsd==2) then
	do  i =1,numel
	  write(iparaview,'(a)') '9'!trim(adjustl(tipo))
	end do
      end if 
!
      if(nsd==3) then
	do  i =1,numel
	  write(iparaview,'(a)') '12'!trim(adjustl(tipo))
	end do
      end if 
      end subroutine escreverTiposElementos

!**** new **********************************************************************
      subroutine escreverEscalaresNodais(v, tam1, tam2, rotulo, tamRot)
      implicit none
      integer*4, intent(in)  :: tam1,tam2
      real*8, intent(in)   :: v(tam1,tam2)
      integer*4:: tamRot
      character(len=tamRot) :: rotulo
!
      character(len=tamRot+5) ::  rotuloN
      integer*4:: i,j
      character(len=5):: eixo
      real*8 :: limite

      limite=1.e-20
      do i=1,tam1

!        if(i>1) then
           write(eixo,'(i0)') i
           if(rotulo.ne.'potencial') then
           rotuloN=trim(rotulo)//'Dir'//trim(eixo)
           write(iparaview,'(3a)')'SCALARS ', trim(rotuloN), ' float '
           write(iparaview,'(a)')'LOOKUP_TABLE default'
           endif
!        endif

        do j=1, tam2
             write(iparaview,*) v(i,j)
        end do
       
      end do

      end subroutine escreverEscalaresNodais

!**** new **********************************************************************

    subroutine escreverArqParaviewIntermed(campo, dim1, dim2, rotulo, tamRot)
    use mMalha, only: x, nsd, numel, numnp

    implicit none
    integer*4, intent(in) :: dim1, dim2
    double precision, intent(in) :: campo(dim1, dim2)

    integer*4:: tamRot
    character(len=tamRot) :: rotulo

     if(rotulo.ne.'velocidade') then
     write(iparaview,'(3a)')'SCALARS ', trim(rotulo), ' float '
     write(iparaview,'(a)')'LOOKUP_TABLE default'
     endif

     call escreverEscalaresNodais(campo, dim1, dim2, rotulo, tamRot)

     end subroutine escreverArqParaviewIntermed


!**** new **********************************************************************
! *** BEGIN -- ROTINA PARA ESCREVER MALHA PARA O PARAVIEW (AUTOR: KADU/PROMEC 06/2011)

    subroutine escreverArqParaview_geo()

    use mMalha, only: x, conecNodaisElem, nsd, numel, numnp

    implicit none
    integer*4:: i , iparaview_geo  

      iparaview_geo  = 31
      open(unit=iparaview_geo , file= 'modelo.geo')

      write(iparaview_geo,'(a)')'Title1'
      write(iparaview_geo,'(a)')'Title2'
      write(iparaview_geo,'(a)')'node id given'
      write(iparaview_geo,'(a)')'element id given'
      write(iparaview_geo,'(a)')'coordinates'
      write(iparaview_geo,'(i8)')  numnp

      do i = 1, numnp
      WRITE (iparaview_geo,'(I8,3E12.5)') I,x(1,i),x(2,i),0.d0
      enddo

      WRITE (iparaview_geo,'(A,/,A,/,A,/,I8)')                     &
                                'part 1'           ,    &
                                'malha'            ,    &
                                'quad4'            ,    &
                                 numel


      WRITE (iparaview_geo,'(5I8)')  (I,conecNodaisElem(1,i),conecNodaisElem(2,i),conecNodaisElem(3,i), &
                                                      conecNodaisElem(4,i),i=1, numel ) 

     end subroutine escreverArqParaview_geo


!**** new **********************************************************************
! *** BEGIN -- ROTINA PARA ESCREVER ESCALAR POR ELEMENTO PARA O PARAVIEW (AUTOR: KADU/PROMEC 06/2011)

    subroutine escreverArqParaview_res(campo, ndim , npasso, var)

    implicit none
    integer*4, intent(in) :: ndim, npasso
    double precision, intent(in) :: campo(ndim)
    character*1, intent(in) :: var

    integer*4:: i
    character*30 :: filename
    integer*4, parameter :: iparaview_res=41

    write (filename(1:3),'(i3.3)') npasso
    if (var=='p') write (filename,'(a)')   'potencial'//'.'//filename(1:3)
    if (var=='s') write (filename,'(a)') 'saturation'//'.'//filename(1:3)

    open  (unit=iparaview_res,file=filename,form='formatted')

       write (iparaview_res,'(a,i5,1x,a)') 'Ensight Escalar passo ',npasso
       write (iparaview_res,'(a/a)') 'part 1','quad4' 
       write (iparaview_res,'(1p,6e12.5)') (campo(i),i=1,ndim)

    close  (unit=iparaview_res)


    end subroutine escreverArqParaview_res


!
!=======================================================================
!     
      subroutine prt(nen,nsd,numel,conectElem,t0,u,iunit)
!      
      use mMalha, only: x
      use mMalha, only: local
!      
      implicit none
!     
!     imprime campos escalares para o gnuplot ou para o matlab
!
      integer*4                  :: nen,nsd,numel
      integer*4                  :: conectElem(nen,numel)
      real(8), dimension(*)     :: u
      real(8)                   :: t0
!     
      integer*4:: nel
      real(8) :: xx,yy,uu
!     
      integer*4:: iunit
!
      real*8 :: xg, yg
      real*8 :: xl(nsd,nen)
!     
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
      write(iunit,*)
!     
      do nel=1,numel

         call local(conectElem(1,nel),x,xl,nen,nsd,nsd)
         xg = sum(xl(1,1:nen))/nen
         yg = sum(xl(2,1:nen))/nen
         xx=xg !xc(1,nel)
         yy=yg !xc(2,nel)
!
         uu=u(nel)
!
         write(iunit,"(5(f25.15,2x))") xx,yy,uu 
      end do
!     
      write(iunit,*)
!     
      end subroutine
!     
!=======================================================================
!
      subroutine prtvB(nen,nsd,numel,conecNodaisElem,t0,velocLadal,ndofV, numLados, conecLadaisElem, numLadosElem, iunit)
!
      use mMalha, only: x
      use mMalha, only: local
!
      implicit none
!
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      integer*4                  :: nen,numel,nsd, ndofV ,numLados,numLadosElem
      real(8), dimension(ndofV,numLados) :: velocLadal
      real(8)                   :: t0
      integer*4                  :: conecNodaisElem(nen,numel),conecLadaisElem(numLadosElem,numel)
!
      integer*4:: nel
      real*8 :: xl(nsd,nen)
      real(8) :: xg,yg, vcx,vcy
       integer*4:: lado1, lado2, lado3, lado4
!
      integer*4:: iunit
!
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!
      do nel=1,numel
!
        call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
        xg = sum(xl(1,1:nen))/nen
        yg = sum(xl(2,1:nen))/nen

        lado1 = conecLadaisElem(1,nel);  lado2 = conecLadaisElem(2,nel);
        lado3 = conecLadaisElem(3,nel);  lado4 = conecLadaisElem(4,nel);
        vcx = (velocLadal(1,lado2)+velocLadal(1,lado4))/2.0
        vcy = (velocLadal(1,lado1)+velocLadal(1,lado3))/2.0
        write(iunit,"(5(f25.15,2x))") xg,yg,vcx,vcy ! velocLadal(1,nel),velocLadal(2,nel)

      end do
!
      write(iunit,*)
!
      end subroutine

!     
!=======================================================================
!     
      subroutine prtv(nen,nsd,numel,ndofV,numLados,conecElem,t0,u,iunit)
!      
      use mMalha, only: x
      use mMalha, only: local
!
      implicit none
!     
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      integer*4                  :: nen,numel,nsd, ndofV, numLados
      real(8), dimension(ndofV,numLados) :: u
      real(8)                   :: t0
      integer*4                  :: conecElem(nen,numel)
!
      integer*4:: nel
      real*8 :: xl(nsd,nen)
!     
      real(8) :: xx,yy
!     
      integer*4:: iunit
!     
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!     
      do nel=1,numel
!     
         call local(conecElem(1,nel),x,xl,nen,nsd,nsd)
         xx = sum(xl(1,1:nen))/nen
         yy = sum(xl(2,1:nen))/nen
!
         write(iunit,"(5(f25.15,2x))") xx,yy,u(1,nel),u(2,nel)
!          write(*,"(5(f25.15,2x))") xx,yy,u(1,nel),u(2,nel)
      end do
!     
      write(iunit,*)
!     
      end subroutine

!
!=================================================================================
!
    subroutine abrirArquivos(comDS_)

!        iin    = input unit number
!        iecho  = output unit of input data
!        iouter  = output unit of error norms
!
    logical, intent(in) :: comDS_

    character(len=20) :: nomeIn, nomeEcho
!
      iin        = 15
      iecho      = 16 
      icoords    = 18
      iconects   = 19
!
      ignuplotPotencial = 30
      ignuplotFluxo   = 31
      iparaview       = 32
!
      nomeIn='input.dat'
      if(comDS_.eqv..true.) nomeIn='inputDS.dat'
      write(*,*) " +++", nomeIn
      nomeEcho='echo.dat'
!
      open(unit=iin,    file=nomeIn, status='old', err=100) 
      open(unit=iecho , file=nomeEcho)
!
      open(unit=icoords   ,file= 'coordenadas.dat')
      open(unit=iconects  ,file= 'conectsNodais.dat')
!
      open(unit=ignuplotPotencial ,file= 'resultadoP.dat')
      open(unit=ignuplotFluxo     ,file= 'resultadoF.dat')
      open(unit=iparaview ,file= 'resultado.vtk')

      return 

      100 continue
      write(*,*) ' arquivo: ', nomeIn, ' NAO encontrado'
      stop 1
       
   end subroutine abrirArquivos
!
!=================================================================================
!
    subroutine fecharArquivos()

      close(iin      )
      close(iecho    )
      close(icoords  )
      close(iconects )
      close(ignuplotPotencial )
      close(ignuplotFluxo   )
      close(iparaview)
    end subroutine fecharArquivos

      end module
