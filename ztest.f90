

      subroutine ztest(a,n,lzero)
!
!.... program to determine if an array contains only zero entries
!
      implicit none
!
!.... remove above card for single precision operation
!
      real*8 :: zero = 0.0
      integer*4, intent(in)    :: n
      real*8, intent(in)     :: a(n)
      logical, intent(inout) :: lzero
!
      integer*4:: i
!
      lzero = .true.
!
      do 100 i=1,n
      if (a(i).ne.zero) then
         lzero = .false.
         return
      endif
  100 continue
!
      end subroutine ztest

      program ztestClient

      real*8 :: zero = 0.0
      real*8 :: f(2,4,3)
      logical :: lzero
      f=zero
      f(:,:,2)=zero+1.

      write(*,'(a,*(f3.1,1x))') "f = ", f
      write(*,*) "lzero = ", sum(f(:,:,1)) == zero
      write(*,*) "lzero = ", sum(f(:,:,2)) == zero
      write(*,*) "lzero = ", sum(f(:,:,3)) == zero
      lzero = sum(f(:,:,1)) == zero
      write(*,*) "lzero = ", lzero
      lzero = sum(f(:,:,2)) == zero
      write(*,*) "lzero = ", lzero

      f(1,2,3)=1
      write(*,'(a,*(f3.1,1x))') "f(:,:,3) = ", f(:,:,3) 
      write(*,*) "any non zero?  ", (f(:,:,3) /= zero) 
      write(*,*) "lzero = ", .not. any(f(:,:,3) /= zero) 

      end program ztestClient
