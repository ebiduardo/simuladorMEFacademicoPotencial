   module mGlobais
   
   integer ::  nsd
   integer :: ndofF, ndofP
   integer, parameter :: zero=0
   
   
   contains
   
   !
!**** new **********************************************************************
!
      subroutine imprimirInfChar(label, informacao, iout)
      implicit none
      
      character(len=80), intent(in) :: label, informacao
      integer, intent(in) :: iout
      
      write(iout,'(a)') label
      write(iout,'(a)') informacao
      write(iout,'(a)') "}"
   
      end subroutine
!
!**** new **********************************************************************
!
      subroutine imprimirInfInt(label, informacao, iout)
      implicit none
      
      character(len=80), intent(in) :: label
      integer ::  informacao
      integer, intent(in) :: iout
      
      write(iout,'(a)') label
      write(iout,'(i8)') informacao
      write(iout,'(a)') "}"
   
      end subroutine
!
!**** new **********************************************************************
!
      subroutine imprimirInfReal(label, informacao, iout)
      implicit none
      
      character(len=80), intent(in) :: label
      real*8,  intent(in) ::  informacao
      integer, intent(in) :: iout
      
      write(iout,'(a)') label
      write(iout,'(f10.4)') informacao
      write(iout,'(a)') "}"
   
      end subroutine
      
!
!**** new **********************************************************************
!
      subroutine readstat(flag,ifile,x)
!
      implicit none
      integer :: ifile
      character(len=80) :: x
      character(len=80) :: flag,flag1
!
      read(ifile,"(a)") flag1
      if(trim(flag).eq.trim(flag1)) then
      read(ifile,*) x
      else
      write(*,*) "Erro na leitura de ", flag
      stop
      end if
!
      end subroutine
      
!
!**** new **********************************************************************
!
      subroutine imprimir(label, informacao, iout)
      implicit none
      
      character(len=80), intent(in) :: label, informacao
      integer, intent(in) :: iout
      
      write(iout,'(a)') label
      write(iout,'(a)') informacao
   
      end subroutine

      
    end module
    
    

   program main

      implicit none
      
      integer :: inData, inInput, inRand, inNovoInput
      character(len=28) :: input, datain, randfiles, novoInput

      
      integer :: i, numLinhasDataIn, cont

      input    ="./input.dat";     inInput=20     
      novoInput="inputDS.dat"; inNovoInput =23
      
      OPEN(UNIT=inInput,FILE=input)     
      OPEN(UNIT=inNovoInput ,FILE=novoInput)

      call alterarFormatoInput(inInput, inNovoInput)
     

      CLOSE(inInput)

   end program main
!
!**** new **********************************************************************
!
      subroutine alterarFormatoInput(iin, iout)
      
      use mGlobais
   
      implicit none
      
      integer :: iin, iout
!      
      character(len=80) :: flag, informacao
!
      integer :: iech,exec, iprtin
      integer :: numnp, numel, nen, npint
      integer :: nlvectP, nlvectF
      character(len=80):: title
      integer :: nummat, i, m, numLadosElem
      real*8:: c(3), grav(3)

   
      read(iin,'(20a4)') iech
      if (iech.eq.0) return
   
      read(iin,'(a)') title
      read(iin,'(3i10)') exec, iprtin, nsd
!             
      read(iin,'(5i10)') numnp, numel, nen, npint
      read(iin,'(3i10)') nlvectP, nlvectF
   
      flag="*title_pc{" 
      call imprimirInfChar(flag, title, iout)
      
      flag="*exec_pc{"
      call imprimirInfInt(flag, exec, iout)
      
      flag="*iprtin_pc{"
      call imprimirInfInt(flag, iprtin, iout)
      
      flag="*nsd_pc{"
      call imprimirInfInt(flag, nsd, iout)
   
      flag="*numnp_pc{"
      call imprimirInfInt(flag, numnp, iout)

      flag="*numel_pc{"
      call imprimirInfInt(flag, numel, iout)

      flag="*nen_pc{"
      call imprimirInfInt(flag, nen, iout)
      
      flag="*npint_pc{"
      call imprimirInfInt(flag, npint, iout)
      
      flag="*nlvectP_pc{"
      call imprimirInfInt(flag, nlvectP, iout) 

      flag="*nlvectF_pc{"
      call imprimirInfInt(flag, nlvectF, iout) 
          
      call lerImprimirCoordenadas(nsd, iin, iout)
      
      call lerEscreverConectsNodais(nen,nsd,iin,iout)
      
      write(iout,'(a)') "*codigos_cond_contorno_potencial_pc{"
      ndofP=1
      call lerImprimirCodigosCondContorno(ndofP,iin,iout)

      if(nlvectP.gt.0) then
         write(iout,'(a)') "*valores_cond_contorno_potencial_pc{"
         call lerImprimirValoresCondContorno(ndofP, iin, iout)
      endif

      write(iout,'(a)') "*codigos_cond_contorno_fluxo_pc{"
      ndofF=nsd
      call lerImprimirCodigosCondContorno(ndofF,iin,iout)
      if(nlvectF.gt.0) then
         write(iout,'(a)') "*valores_cond_contorno_fluxo_pc{"
         call lerImprimirValoresCondContorno(ndofF, iin, iout)
      endif
      
      read(iin,'(i10)') nummat
      flag="*nummat_pc{"
      call imprimirInfInt(flag, nummat, iout)   
      
      write(iout,'(a)') "*prop_fisica_meio_pc{"
      read (iin,  5000) m,(c(i),i=1,3)
      write(iout,5000) m,(c(i),i=1,3)
      write(iout,'(a)') "}"
      
      write(iout,'(a)') "*grav_pc{"
      read  (iin,  '(3e10.2)') (grav(i),i=1,3)
      write (iout,'(3e10.2)') (grav(i),i=1,3)
      write(iout,'(a)') "}"
      
      
 5000  format(i10,5x,3e10.2)
   contains
!
!**** new **********************************************************************
!
      subroutine lerEscreverConectsLadais3D(nen, nsd, iin, inout)   
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
      integer :: nen, nsd, iin, inout
      integer :: itemp(27) 
!
      integer :: n,m,ng,i,nel(3),incel(3),inc(3)
!
      write(inout,'(a)') "*conectividades_ladais_pc{"
  100 continue                                                          
      read(iin,1000) n,m,(itemp(i),i=1,nen),ng     
      if(n==0) then
         write(inout,'(i10)') zero
         write(inout,'(a)') "}"
      else
         if(ng==0) then
            write(inout,1000) n,m,(itemp(i),i=1,nen)
         else
            write(inout,1000) n,m,(itemp(i),i=1,nen), ng
         endif   
      endif
                
      if (n.eq.0) return                                                
      if (ng.ne.0) then
! !                                                                       
! !....... generate data                                                     
! !                                                                       
!          read(iin,1000) (nel(i),incel(i),inc(i),i=1,nsd)     
!           write(inout,1000) (nel(i),incel(i),inc(i),i=1,nsd)                  
      endif
      go to 100                                                         
!                                                                       
 1000 format(9i10)
!                                                                       
      end subroutine   
!
!**** new **********************************************************************
!
      subroutine lerEscreverConectsLadais2D(nen, nsd, iin, inout)   
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
      integer :: nen, nsd, iin, inout
      integer :: itemp(27) 
!
      integer :: n,m,ng,i,nel(3),incel(3),inc(3)
!
      
      write(inout,'(a)') "*conectividades_ladais_pc{"
  100 continue                                                          
      read(iin,1000) n,m,(itemp(i),i=1,nen),ng     
      if(n==0) then
         write(inout,'(i10)') zero
         write(inout,'(a)') "}"
      else
         if(ng==0) then
            write(inout,1000) n,m,(itemp(i),i=1,nen)
         else
            write(inout,1000) n,m,(itemp(i),i=1,nen), ng
         endif
         
      endif
                
      if (n.eq.0) return                                                
      if (ng.ne.0) then
!                                                                       
!....... generate data                                                     
!                                                                       
         read(iin,1000) (nel(i),incel(i),inc(i),i=1,nsd)     
          write(inout,1000) (nel(i),incel(i),inc(i),i=1,nsd)                  
      endif
      go to 100                                                         
!                                                                       
 1000 format(16i10,10x,14i10)                                             
!                                                                       
      end subroutine   
!
!**** new **********************************************************************
!
      subroutine lerEscreverConectsNodais(nen, nsd, iin, inout)   
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
      integer :: nen, nsd, iin, inout
      integer :: itemp(27) 
!
      integer :: n,m,ng,i,nel(3),incel(3),inc(3)
!
      
      write(inout,'(a)') "*conectividades_nodais_pc{"
  100 continue                                                          
      read(iin,1000) n,m,(itemp(i),i=1,nen),ng     
      if(n==0) then
         write(inout,'(i10)') zero
         write(inout,'(a)') "}"
      else
         if(ng==0) then
            write(inout,1000) n,m,(itemp(i),i=1,nen)
         else
            write(inout,1000) n,m,(itemp(i),i=1,nen), ng
         endif 
      endif
                
      if (n.eq.0) return                                                
      if (ng.ne.0) then
!                                                                       
!....... generate data                                                     
!                                                                       
         read(iin,1000) (nel(i),incel(i),inc(i),i=1,nsd)     
          write(inout,1000) (nel(i),incel(i),inc(i),i=1,nsd)                  
      endif
      go to 100                                                         
!                                                                       
 1000 format(16i10,10x,14i10)                                             
!                                                                       
      end subroutine   
!
!**** new **********************************************************************
!
      subroutine lerImprimirValoresCondContorno(nra,iin,inout)   
!                                                                       
!.... program to read floating-point nodal data            
!                                                                                                   
!         nra     = number of rows
!         n       = node number
!         numgp   = number of generation points
!         ninc(i) = number of increments for direction i
!         inc(i)  = increment for direction i
!                                                                       
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                        
      integer :: nra, iin, inout
!
      real*8  :: temp(6,20)
      integer :: n,numgp,ninc,inc
      integer :: i, j, m, mgen    

      
  100 continue                                                          
      read(iin,1000) n,numgp,(temp(i,1),i=1,nra) 
      if(n==0) then
         write(inout,'(i10)') zero
         write(inout,'(a)') "}"
      endif
     
      if (n.eq.0) return                                                
      write(inout,1000) n,numgp,(temp(i,1),i=1,nra)    

      if (numgp.ne.0) then
         do 200 j=2,numgp
!                                                                       
         read(iin,1000) m,mgen,(temp(i,j),i=1,nra)
         write(inout,1500) (temp(i,j),i=1,nra)        
                 
  200    continue                               
         read(iin,2000) ninc,inc
         write(inout,2000) ninc,inc

      endif
      
      go to 100  
      write(inout,2000) zero
      write(inout,'(a)') "}"

!                                                                       
 1000 format(2i10,3e10.2)
 1500 format(20x,3e10.2)                                                
 2000 format(16i10)                                                      
!                                                                       
      end  subroutine   
!
!**** new **********************************************************************
!
      subroutine lerImprimirCodigosCondContorno(m,iin,inout)
!
!.... program to read integer nodal data
!
      implicit none
!
      integer :: m, iin, inout
!
      integer :: ib(m)
      integer :: n, ne, ng
      integer :: i
!
  100 continue
      read(iin,1000) n,ne,ng,(ib(i),i=1,m)
      if(n==0) then
         write(inout,'(i10)') zero
         write(inout,'(a)') "}"
      else
         write(inout,1000) n,ne,ng,(ib(i),i=1,m)
      endif

      if (n.eq.0) return
!
      go to 100
!
 1000 format(16i10)
      end subroutine
!
!**** new **********************************************************************
!
      subroutine lerImprimirCoordenadas(nra,iin,inout)   
!                                                                       
!.... program to read floating-point nodal data            
!                                                                                                   
!         nra     = number of rows
!         n       = node number
!         numgp   = number of generation points
!         ninc(i) = number of increments for direction i
!         inc(i)  = increment for direction i
!                                                                       
      implicit none
!                                                                       
!.... remove above card for single-precision operation               
!                        
      integer :: nra, iin, inout
!
      real*8  :: temp(6,20)
      integer :: n,numgp,ninc(3),inc(3)   
      integer :: i, j, m, mgen    
!                                                                       
      write(inout,'(a)') "*coordenadas_nodais_pc{"
      
  100 continue                                                          
      read(iin,1000) n,numgp,(temp(i,1),i=1,nra) 
      if(n==0) then
         write(inout,'(i10)') zero
         write(inout,'(a)') "}"
      endif
     
      if (n.eq.0) return                                                
      write(inout,1000) n,numgp,(temp(i,1),i=1,nra)    

      if (numgp.ne.0) then
         do 200 j=2,numgp
!                                                                       
         read(iin,1000) m,mgen,(temp(i,j),i=1,nra)
         write(inout,1500) (temp(i,j),i=1,nra)        

  200    continue                               
         read(iin,2000) (ninc(i),inc(i),i=1,nra)
         write(inout,2000) (ninc(i),inc(i),i=1,nra)

      endif
      
      go to 100  
      write(inout,2000) zero
      write(inout,'(a)') "}"

!                                                                       
 1000 format(2i10,6e10.3)
 1500 format(20x, 6e10.3)                                                
 2000 format(16i10)                                                      
!                                                                       
      end  subroutine   

   
   end subroutine alterarFormatoInput
