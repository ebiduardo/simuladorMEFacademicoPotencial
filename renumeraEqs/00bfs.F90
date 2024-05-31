module grafo

implicit none

type, public :: vertice
    integer num, peso 
  contains
    procedure :: acessarNum 
    procedure :: atribuirNum
    procedure :: mostrarConteudoV
end type vertice

type, public :: verticeL
    type(vertice) , pointer :: V    => null()
    type(verticeL), pointer :: next => null()
    type(verticeL), pointer :: pai  => null()
  contains
    procedure :: mostrarConteudoLV
    !procedure :: incluirVerticeL
end type verticeL

   type(verticeL), allocatable :: adjArray(:)
   type(verticeL), allocatable :: adjArrayC(:)
   type(verticeL), pointer :: toVisit=> null() , visitedL=>null()
   logical, allocatable :: visited(:)

contains

integer pure function acessarNum(this)
    class(vertice), intent(in) :: this
    acessarNum = this%num
end function acessarNum

subroutine atribuirNum(this,num_)
    class(vertice), intent(out) :: this
    integer, intent(in) :: num_
    this%num=num_
end subroutine atribuirNum

subroutine mostrarConteudoV(this)
    class(vertice), intent(in) :: this
    write(*,'(i3,", ")',advance='No') acessarNum(this)
end subroutine mostrarConteudoV

subroutine mostrarConteudoLV(this)
    class(verticeL), intent(in) :: this
    call  mostrarConteudoV(this%V)
end subroutine mostrarConteudoLV

subroutine mostrarConteudoL(lV_)
 type(verticeL), pointer :: lV_

 type(verticeL), pointer :: aux
 write(*,'(a)',advance='no') ': ';
 aux=>lV_
   do while(associated(aux))
     call mostrarConteudoLV(aux)
     aux=>aux%next
   end do
 print*, '...';
end subroutine mostrarConteudoL

subroutine mostrarConteudoG(adjArray_)
 type(verticeL), intent(in) :: adjArray_(:)
 integer :: eq
 print*, "em  mostrarConteudoG, num vertices = ", size(adjArray_)
 do eq = 1, size(adjArray_)
  write(*,'(a, i2, ", ")', advance='no') 'adjs a ', eq; call mostrarConteudoL(adjArray_(eq)%next )
 ! write(*,'(a, i2, ", ")', advance='no') 'bandaLocal= ', maiorValor(adjArray_(eq)%next)-menorValor(adjArray_(eq)%next)
 end do
 print'(/a,i5)', ':::..., banda maxima=', bandaMax(adjArray_)
end subroutine mostrarConteudoG

integer function maiorValor(lV_)
 type(verticeL), pointer :: lV_
 type(verticeL), pointer :: aux
 integer::maiorN
 aux=>lV_
   maiorValor=acessarNum(aux%V) 
   do while(associated(aux))
     if(acessarNum(aux%V)>maiorValor) maiorValor=acessarNum(aux%V)
     aux=>aux%next
   end do
end function maiorValor

integer function menorValor(lV_)
 type(verticeL), pointer :: lV_
 type(verticeL), pointer :: aux
 integer::maiorN
 aux=>lV_
   menorValor=acessarNum(aux%V) 
   do while(associated(aux))
     if(acessarNum(aux%V)<menorValor) menorValor=acessarNum(aux%V)
     aux=>aux%next
   end do
end function menorValor

integer function bandaMax(adjArray_)
 type(verticeL), intent(in) :: adjArray_(:)
 integer :: eq, banda
 eq=1
 banda=maiorValor(adjArray_(eq)%next)-menorValor(adjArray_(eq)%next)+1
 do eq = 2, size(adjArray_)
  banda=maiorValor(adjArray_(eq)%next)-menorValor(adjArray_(eq)%next)+1
  if(banda>bandaMax) bandaMax=banda
 end do
end function bandaMax

subroutine incluirVerticeL(this, v_)
   type(verticeL), pointer, intent(inout) :: this
   type(vertice), target, intent(in) :: v_

   type(verticeL), pointer  :: novo
   !call mostrarConteudoV(v_ )
   allocate(novo)
   novo%V=>v_ 
   novo%next=>this;
   this=>novo;
   return
end subroutine incluirVerticeL

function excluirVerticeMaisAntigo(c_)
       type(verticeL), pointer :: excluirVerticeMaisAntigo
       type(verticeL), intent(inout), pointer:: c_

       type(verticeL), pointer :: p
       type(verticeL), pointer :: ant

       !allocate(ant,p)!,ant)
       print*, " em excluirVerticeMaisAntigo(VertList** c_) "

       if(.not.associated(c_)) then
            print*, " ... lista vazia! "
            excluirVerticeMaisAntigo=>null()
            return
       end if
       ant=>c_;
       if(.not.associated(c_%next)) then
           print*, " QUASE vazia, "
           write(*,'(a)',advance='no') "A, excluindo:"; call mostrarConteudoV(c_%V); print*;
           c_=>null()
           excluirVerticeMaisAntigo=>ant
           return
       end if

       p=>c_
       do 
          if(.not.associated(p%next)) exit
          ant=>p
          p=>p%next;
       end do 
       ant%next=>NULL(); 
       excluirVerticeMaisAntigo=>p;
       return
end function  excluirVerticeMaisAntigo

function bfs (adjArray_, neq_, origem_, destino_)
   type(verticeL), pointer :: bfs 
   type(verticeL), intent(in) :: adjArray_(:)
   type(verticeL), pointer :: current, adj, inicio
   integer, intent(in) :: neq_, origem_, destino_ 
   logical, allocatable :: included(:)

   integer ::  i 
   print*," em bfs, origem_,", origem_, ", destino_ ", destino_

   allocate(visited(neq_));  visited =.false.
   allocate(included(neq_)); included=.false.

   toVisit=>null()
   inicio=>adjArray_(origem_)%next
   do i = 1, 3
     inicio=>inicio%next
   end do
   write(*,'(a)'   , advance='no' ) "incluindo :" 
   write(*,'(a,i5)', advance='yes') ',', acessarNum(inicio%V);
   call incluirVerticeL(toVisit, inicio%V) 
   included(acessarNum(inicio%V) )=.true.

   write(*,'(a)',advance='no') '0, toVisit, '; call mostrarConteudoL(toVisit ); print*
   allocate(adj)  
   do while(associated(toVisit))
      current=>excluirVerticeMaisAntigo(toVisit)
      adj => adjArray_(acessarNum(current%V))%next 
      write(*,'(a)',   advance='no') "adjacents from ="; call mostrarConteudoV(current%V);
      call mostrarConteudoL(adjArray_(acessarNum(current%V))%next);
      write(*,'(a)',   advance='no') "incluindo :" 
      do while(associated(adj))
         if(.not.included(acessarNum(adj%V)) .and. acessarNum(adj%V)/=acessarNum(current%V)) then
            write(*,'(a,i5)',advance='no')',', acessarNum(adj%V); 
            call incluirVerticeL(toVisit, adj%V)
            included(acessarNum(adj%V) )=.true.
         endif
         adj=>adj%next
      end do 
      print*
      call incluirVerticeL(visitedL, current%V)
      visited(acessarNum(current%V))=.true.
      print*, "included = ", included
      print*, "visited  = ", visited
      write(*,'(a)',advance='no') '2, toVisit : '; call mostrarConteudoL(toVisit );
      write(*,'(a)',advance='no') '   visitedL: '; call mostrarConteudoL(visitedL );
   end do
      !call mostrarConteudoL(visitedL)
      bfs => visitedL
end function bfs 

 subroutine montarAdjArray(  adjArray_ ,LMstencilEq_, listaVertices_, neq_, numMaxVizEq_ )
 type(verticeL), intent(in) :: adjArray_(:)
 integer :: LMstencilEq_(neq_,numMaxVizEq_)
 type(vertice) :: listaVertices_(0:) 
 integer, intent(in) :: neq_, numMaxVizEq_
 integer :: eq, i
  do eq = 1, neq_
    allocate(adjArray(eq)%next)
    adjArray(eq)%next=>null()
    i=1
    do while(i<=numMaxVizEq_)
      if(LMstencilEq_(eq,i)>0) then
         call incluirVerticeL(adjArray(eq)%next, listaVertices_(LMstencilEq_(eq,i)))
      endif
      i=i+1
    end do
  end do 
 end subroutine montarAdjArray

end module grafo

program renum
 use grafo
 implicit none
 integer, parameter ::numnp=21,neq = 15, numMaxVizEq=9, ndof=1

 integer :: LMstencilEq(neq,numMaxVizEq)
 integer :: id(ndof,numnp)
 integer :: numNova(neq)
 integer :: eq, i, n

 !https://stackoverflow.com/questions/8900336/arrays-of-pointers
 !does not define an array of pointers, as you might think, but a pointer to an array.
 !type(verticeL), pointer :: adjArray(:)

 type(vertice), allocatable, target :: listaVertices(:) 
 type(verticeL), pointer :: caminho, aux
 allocate(listaVertices(0:neq))
 allocate(adjArray(1:neq)); !adjArray=>null()
 allocate(adjArrayC(1:neq)); !adjArray=>null()

 call setId()
 call setLMstencil()

 call montarAdjArray(adjArray, LMstencilEq, listaVertices, neq, numMaxVizEq)
 call mostrarConteudoG(adjArray)

 caminho => bfs (adjArray, neq, 1, 1)
 call mostrarConteudoL(caminho)
 call mostrarConteudoG(adjArray)
 call renumerar()
 call mostrarConteudoG(adjArray)

contains 

subroutine renumerar()
 write(*, '(a)', advance='no') "numeracao 1 = " 
 do i = 1, neq;
    write(*, '(i0,", ")', advance="no") acessarNum(listaVertices(i))
 end do; print*

 i = neq
 do while (associated(caminho)) 
    adjArrayC(i)=adjArray(acessarNum(caminho%V))
    print'(i0,a,i0)', i,' vai para posição ', acessarNum(caminho%V);
    call atribuirNum(caminho%V,i) ! altera os num dos vertices da listaVertices tambem 
    caminho=>caminho%next
    i = i-1
 end do; 

 adjArray=adjArrayC

 write(*, '(a)', advance='no') "numeracao nova = " 
 do i = 1, neq;
    write(*, '(i0,", ")', advance="no") acessarNum(listaVertices(i))
 end do; print*
end subroutine renumerar

 subroutine setLMstencil()
 integer :: eq, i, n
LMstencilEq =0
n=4; eq=1; LMstencilEq(eq, 1:n) = (/1,2,6,7 /)
print'(a, 30i3)', "LMstencilEq(1, 1:4) =", LMstencilEq(eq, 1:n)
n=6; eq=2; LMstencilEq(eq, 1:n) = (/1,2,3,6,7,8 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=6; eq=3; LMstencilEq(eq, 1:n) = (/2,3,4,7,8,9 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
!return
n=6; eq=4; LMstencilEq(eq, 1:n) = (/3,4,5,8,9,10 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=4; eq=5; LMstencilEq(eq, 1:n) = (/4,5,9,10 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=6; eq=6; LMstencilEq(eq, 1:n) = (/1,2,6,7,11,12 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=9; eq=7; LMstencilEq(eq, 1:n) = (/1,2,3,6,7,8,11,12,13 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=9; eq=8; LMstencilEq(eq, 1:n) = (/2,3,4,7,8,9,12,13,14 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=9; eq=9; LMstencilEq(eq, 1:n) = (/3,4,5,8,9,10,13,14,15 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=6; eq=10; LMstencilEq(eq, 1:n) = (/4,5,9,10,14,15 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=4; eq=11; LMstencilEq(eq, 1:n) = (/6,7,11,12 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=6; eq=12; LMstencilEq(eq, 1:n) = (/6,7,8,11,12,13 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=6; eq=13; LMstencilEq(eq, 1:n) = (/7,8,9,12,13,14 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=6; eq=14; LMstencilEq(eq, 1:n) = (/8,9,10,13,14,15 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
n=4; eq=15; LMstencilEq(eq, 1:n) = (/9,10,14,15 /)
print'(i3,a,30i3)',eq,',', LMstencilEq(eq, 1:n)
end subroutine setLMstencil

subroutine setId()
 integer :: eq, i, n
n=1; id(1:ndof, n)=(/0/); n=2; id(1:ndof, n)=(/1/);
n=3; id(1:ndof, n)=(/2/); n=4; id(1:ndof, n)=(/3/);
n=5; id(1:ndof, n)=(/4/); n=6; id(1:ndof, n)=(/5/);
n=7; id(1:ndof, n)=(/0/); n=8; id(1:ndof, n)=(/0/);
n=9; id(1:ndof, n)=(/6/); n=10; id(1:ndof, n)=(/7/);
n=11; id(1:ndof, n)=(/8/); n=12; id(1:ndof, n)=(/9/);
n=13; id(1:ndof, n)=(/10/); n=14; id(1:ndof, n)=(/0/);
n=15; id(1:ndof, n)=(/0/); n=16; id(1:ndof, n)=(/11/);
n=17; id(1:ndof, n)=(/12/); n=18; id(1:ndof, n)=(/13/);
n=19; id(1:ndof, n)=(/14/); n=20; id(1:ndof, n)=(/15/);
n=21; id(1:ndof, n)=(/0/);
write(*, '(a)', advance='no') "id = " 
write(*, '(100(i0,", "))') id(:,1:numnp)
 do i = 1, numnp
   call atribuirNum(listaVertices(id(1,i)),id(1,i))
 end do
 do i = 1, neq; call mostrarConteudoV(listaVertices(i)); end do; print*
 end subroutine setId
end program renum
