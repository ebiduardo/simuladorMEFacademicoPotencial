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
!
!      subroutine ztest(a,n,lzero)
!      use mGlobaisEscalares, only : zero
!      implicit none
!      integer*4, intent(in)    :: n 
!      real*8, intent(in)     :: a(n)
!      logical, intent(inout) :: lzero
!      integer*4:: i
!      lzero = .true.
!      do 100 i=1,n
!      if (a(i).ne.zero) then
!         lzero = .false.
!         return
!      endif
!  100 continue
!      end subroutine ztest
!
!**** new **********************************************************************
!
      subroutine timing(tempoDeParede)
      implicit none     
!
!.... program to determine elapsed cpu time
!
      real*8, intent(inout) :: tempoDeParede
!
      character(LEN=8)  :: date
      character(LEN=10) :: minhaHora
      character(LEN=5)  :: zone
      integer*4,dimension(8) :: values
      integer*4::  horas, minutos, segundos, milesimosSeg
!
      call date_and_time(date,minhaHora,zone,values);
!
      horas=values(5);      minutos=values(6); 
      segundos=values(7); milesimosSeg=values(8);    
!
      tempoDeParede = (60*horas+minutos)*60.0+segundos+milesimosSeg/1000.00
!
      return
      end subroutine timing
!
