  module mEstruturasDadosSistEq
  
  type estruturasArmazenamentoSistemaEq 
        integer*4           :: neq, nalhs   ! skyline + CSR + HYPRE
        integer*4           :: ndof, nlvect ! skyline + CSR + HYPRE
        integer*4, pointer  :: lm(:,:,:)=> null() ! skyline + CSR + HYPRE
        integer*4, pointer  :: id(:,:)  => null() ! skyline + CSR + HYPRE   

        integer*4, pointer  :: idiag(:) => null() ! Skyline + CSR

        real*8, pointer :: u(:,:), f(:,:,:) ! skyline + CSR + HYPRE   
        real*8, pointer :: uTempoAnt(:,:), uIterAnt(:,:)
        real*8, pointer :: brhs(:)=>null()  ! skyline + CSR + HYPRE
        real*8, pointer :: alhs(:)=>null()  ! skyline + CSR 
        real*8, pointer :: clhs(:)=>null()  ! skyline 

        integer*4, pointer  :: rows(:)          => NULL()                                    ! HYPRE
        real*8   , pointer  :: initialGuess(:)  => null()                       ! HYPRE

        integer*4,  pointer :: Ap(:)=>null()               ! CSR
        integer*4,  pointer :: Ai(:)=>null()               ! CSR
        integer*4,  pointer :: LMstencilEq(:,:) => null()  ! CSR
        integer*4           :: numCoefPorLinha, nVizinMax  ! CSR
        integer*4           :: pt(64), iparm(64)
        REAL*8              :: dparm(64)
             
        integer*8           :: A_HYPRE, parcsr_A, b_HYPRE, par_b, u_HYPRE, par_u ! HYPRE
        integer*8           :: solver                                            ! HYPRE
        integer*4           :: solver_id, precond_id                             ! HYPRE
        integer*4           :: Clower, Cupper                                    ! HYPRE
        integer*4           :: Flower, Fupper                                    ! HYPRE
        integer*4           :: num_iterations
        real*8              :: tol, final_res_norm
  
        real*8 :: uInicial
        logical           :: simetria
        character(len=80) :: optSolver=""
        integer*4, pointer :: listaElemPorNo(:,:)
        character(LEN=15) :: label

        character(len=4)  :: eliminate=""
        
        
   end type

      contains
!
!=============================================================================
!
      subroutine formlm (id,conecElem,lm,ndof,ned,nen,numel)
!
      implicit none
!
!.... program to form lm array
!
      integer*4:: ndof,ned,nen,numel
      integer*4:: id(ndof,*),conecElem(nen,*),lm(ned,nen,*)
!
      integer*4:: i,j,k,node
!
      do 300 k=1,numel
!
      do 200 j=1,nen
      node=conecElem(j,k)
!
      do 100 i=1,ndof

      lm(i,j,k) = id(i,node)

  100 continue
!
  200 continue
!
  300 continue

!
      return
      end subroutine
      !**** new **********************************************************************
!
      subroutine diag(idiag,neq,n)
!
      implicit none
!
!.... program to compute diagonal addresses of left-hand-side matrix
!
      integer*4:: neq
      integer*4:: n
      integer*4:: idiag(neq)
!
      integer*4:: i

      n = 1
      idiag(1) = 1 
      if (neq.eq.1) return
!
      do 100 i=2,neq
      idiag(i) = idiag(i) + idiag(i-1) + 1
  100 continue
      n = idiag(neq)
    !write(*,'(a, 20i4)') ", idiag=", idiag(1:10)
    !write(*,'(a, 20i9)') ", idiag=", idiag(neq-10:neq)
!
      return
      end subroutine
!**** new **********************************************************************
      subroutine colht(idiag,lm,ned,nen,numel,neq)
!
!.... program to compute column heights in global left-hand-side matrix
!
      implicit none
      integer*4:: ned, nen, numel
      integer*4:: neq
      integer*4:: idiag(*),lm(ned,nen,*)
!
      integer*4:: i, j, k
      integer*4:: m, minimo, num
!
!     write(*,*) "lm=", lm(:,:,1:numel)

      do 500 k=1,numel
      minimo = neq
      do 200 j=1,nen
      do 100 i=1,ned
      num = abs(lm(i,j,k))
      if (num.gt.0) minimo = min(minimo,num)  ! min= min0(min,num) 
  100 continue
  200 continue

!     write(*,'(a,4i3)', advance='NO') ' lm=', lm(:,:,k)
!     write(*,'(a,i2,a,i2)', advance='NO') ' nel=', k, ', minimo = ', minimo
!
      do 400 j=1,nen
      do 300 i=1,ned
      num = abs(lm(i,j,k))
      if (num.gt.0) then
         m = num - minimo
         if (m.gt.idiag(num)) then
              idiag(num) = m
         endif
      endif
  300 continue
  400 continue
  500 continue
    !write(*,'(a, 20i2)') ", idiag=", idiag(1:10)
    !write(*,'(a, 20i9)') ", idiag=", idiag(neq-10:neq)

  !  write(*,'(a, i5)') ", neq=", neq
  !  write(*,'(a, 200i5)') ", idiag=", idiag(1:200); stop
!
      return
      end subroutine
      
!******************************************************************************
      subroutine criarGrafoEquacoesPorNo(estrutSistEq_,  nsd, conectsElem, listaDosElems, &
                                                 numConexoes, nen, numConexoesPorElem)
      use mMalha, only : numnp
      implicit none 
!
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
        
      integer*4,  intent(in) :: nsd, numConexoes, nen, numConexoesPorElem
      integer*4,  intent(in) :: conectsElem(nen,*)
      integer*4 :: listaDosElems(numConexoesPorElem,numnp)
      integer*4 :: i
	  
      write(*,*) " em subroutine criarGrafoEquacoesPorNo(nsd, ndof, neq, num"
!

         estrutSistEq_%simetria = .false.
         !if(.not.associated(LMstencilEq_))
         print*, " simetria =", estrutSistEq_%simetria

         !estrutSistEq_%numCoefPorLinha = (estrutSistEq_%nVizinMax)*estrutSistEq_%ndof
		 		 print*, "estrutSistEq_%nVizinMax = ",estrutSistEq_%nVizinMax
         estrutSistEq_%numCoefPorLinha = (estrutSistEq_%nVizinMax+1)*estrutSistEq_%ndof*nen;
		 		 print*, "estrutSistEq_%numCoefPorLinha = ", estrutSistEq_%numCoefPorLinha
!         write(*,*) " estrutSistEq_%nVizinMax       =", estrutSistEq_%nVizinMax 
!         write(*,*) " estrutSistEq_%neq             =", estrutSistEq_%neq 
!         write(*,*) " estrutSistEq_%numCoefPorLinha =", estrutSistEq_%numCoefPorLinha 
         allocate(estrutSistEq_%LMstencilEq(estrutSistEq_%neq,estrutSistEq_%numCoefPorLinha)); 
         estrutSistEq_%LMstencilEq=0

!         print*, "shape(estrutSistEq_%LMstencilEq)=",shape(estrutSistEq_%LMstencilEq)
!		 do i = 1, numnp
!          print*, listaDosElems(:,i)
!         end do
	 	 !stop 10
         print*, "+++numCoefPorLinha=", estrutSistEq_%numCoefPorLinha
         call montarLmStencilNodal_CSR (estrutSistEq_%LMstencilEq, listaDosElems, estrutSistEq_%id, &
                   estrutSistEq_%lm, &
                   conectsElem, estrutSistEq_%numCoefPorLinha,  &
                   estrutSistEq_%ndof, numConexoes, nen, estrutSistEq_%nVizinMax, &
                   estrutSistEq_%neq, estrutSistEq_%simetria)
         do i = 1, estrutSistEq_%neq
          print*, i, estrutSistEq_%LMstencilEq(i,:)
        end do
      stop 22
      return
      end subroutine criarGrafoEquacoesPorNo !(estrutSistEq_ ,  nsd, conectsElem, listaDosElems, &
!
!******************************************************************************
      subroutine criarPonteirosMatEsparsa_CSR(estrutSistEq_,  nsd, conectsElem, listaDosElems, &
                                                 numConexoes, nen, numConexoesPorElem)
!                                        
      use mMalha, only : numnp
      implicit none 
!
      type (estruturasArmazenamentoSistemaEq) :: estrutSistEq_
        
      integer*4,  intent(in) :: nsd, numConexoes, nen, numConexoesPorElem
      integer*4 :: conectsElem(nen,*), listaDosElems(numConexoesPorElem,numnp)
      integer*4 :: i

      write(*,*) " em subroutine criarPonteirosMatEsparsa_CSR(nsd, ndof, neq, num"
!
         estrutSistEq_%numCoefPorLinha = (estrutSistEq_%nVizinMax+1)*estrutSistEq_%ndof;
         allocate(estrutSistEq_%LMstencilEq(estrutSistEq_%neq,estrutSistEq_%numCoefPorLinha)); 

         call montarLmStencilNodal_CSR (estrutSistEq_%LMstencilEq, listaDosElems, estrutSistEq_%id, &
                   estrutSistEq_%lm, &
                   conectsElem, estrutSistEq_%numCoefPorLinha,  &
                   estrutSistEq_%ndof, numConexoes, nen, estrutSistEq_%nVizinMax, &
                   estrutSistEq_%neq, estrutSistEq_%simetria)

         call montarPonteiroAp_CSR(estrutSistEq_%Ap, estrutSistEq_%LMstencilEq, estrutSistEq_%numCoefPorLinha,&
                                                                        estrutSistEq_%neq, estrutSistEq_%nalhs)

         allocate(estrutSistEq_%Ai(estrutSistEq_%nalhs)); estrutSistEq_%Ai=0
         call montarPonteiroAi_CSR(estrutSistEq_%Ai, estrutSistEq_%LMstencilEq, estrutSistEq_%numCoefPorLinha,&
                                                                        estrutSistEq_%neq, estrutSistEq_%nalhs)

!        call escreverEstruturaEsparsa( AlhsP, brhsP, ApPotencial, AiPotencial, neq, nalhs
! 
              return
              end subroutine

      !
!**** new *************************************************************
      subroutine montarLmStencilNodal_CSR(LMstencilEq, listaDosElems, id, lm, conectsElem,&
                        numCoefPorLinha, ndof, numConexoes, nen, nVizinMax, neq, simetria) 
      use mMalha, only : numnp, numel
      implicit none
      integer*4,  intent(out) :: LMstencilEq(neq,numCoefPorLinha)
      integer*4,  intent(in)  :: listaDosElems(nVizinMax,*)
      integer*4,  intent(in)  :: id(ndof,numnp), conectsElem(nen,*)
      integer*4,  intent(in)  :: lm(nen,ndof,numel)
      integer*4,  intent(in)  :: numCoefPorLinha, ndof, numConexoes
      integer*4,  intent(in)  :: nen, nVizinMax, neq
      logical,    intent(in)  :: simetria
!
      integer*4 :: nc, eq, pId, no, elemViz, i, j, k, nel
      integer*4 :: numVizAtualizado, numZerosID, numEqs 
      integer*4 :: LMstencilEqId(nVizinMax*(nen)*ndof) 
      integer*4 :: copiaIdNode(nen*ndof)
      integer*4 :: LMLocal(nen*ndof)
      logical   :: simetriaB

      write(*,*) " em subroutine montarLmStencilNodal_CSR(...., "
      print*, " simetria =", simetria, ", numCoefPorLinha=", numCoefPorLinha
      !stop 11111
    !call escreverId_Lm()

      LMstencilEq=0
      do nc=1, numConexoes !nc = numnp
       do j = 1, ndof
          eq = abs(id(j,nc))
          if (eq > 0) then
            numVizAtualizado=nVizinMax
!			print*," eq=", eq," numVizAtualizado =" ,numVizAtualizado
               LMstencilEqId=0
           ! construindo uma lista de numeros de equacoes dos nos dos elementos
           ! vizinhos ao nó que possui a equacao eq
           ! incluindo: numero de eq. repetidos e zeros do vetor original com tamanho exagerado 
            do k=1, nVizinMax ! percorre as equacoes dos vizinhos do nó nc 
               elemViz=listaDosElems(k,nc)
               if(elemViz>0) then
                  LMLocal(:)=reshape(lm(:,:,elemViz),(/nen*ndof/))
                  pId = 1 + nen * ndof * (k - 1)  
     ! print*, " pId =", pId, ", nen = ", nen, ", ndof = ", ndof
     ! print*, " LMstencilEqId(pId:pId+nen*ndof  =",  LMstencilEqId(pId:pId+nen*ndof-1)
!      print*, " abs(LMLocal(:) =", abs(LMLocal(:))
                  LMstencilEqId(pId:pId+nen*ndof-1) = abs(LMLocal(:))
     ! stop 22
               else
                  numVizAtualizado=numVizAtualizado-1
               end if
            end do ! k=1, nVizinMax 
			
    !B1   write(*,'(a,i0,a,81(1x,i0))'), "1, noh=", nc, ", LmStencilEqId, ->", LMstencilEqId(1:numVizAtualizado*nen*ndof)
           ! ordenando a lista de numeros de equacoes para retirar os numeros repetidos que
           ! ficarão um ao lado do outro
            call ordenarLMstencil(LMstencilEqId,numVizAtualizado*nen*ndof)

    !B1  write(*,'(a,i0,a,81(1x,i0))'), "1.1, noh=", nc, ", LmStencilEqId, ->", LMstencilEqId(1:numVizAtualizado*nen*ndof)
           ! eliminando os numeros de equacoes repetidos, atribuindo zero no lugar
    !       if(simetria) then !elimina num. eq. repetidas
!	print*,eq,"b, LMstencilEqId=", LMstencilEqId
            do i=1, numVizAtualizado*nen*ndof-1 
             if(LMstencilEqId(i)==LMstencilEqId(i+1))LMstencilEqId(i)=0!elimina num. eq. repetidas
            end do
    !       end if
           ! ordenando a lista de numeros de equacoes SEM repeticao 
           !     para agrupar os zeros aa esquerda  e uma lista correta de num.
           !     de equacoes aa direita
            call ordenarLMstencil(LMstencilEqId,numVizAtualizado*nen*ndof)
!			print*,eq,"c, LMstencilEqId=", LMstencilEqId
           ! eliminando os numeros das equacoes menores de o valor de eq
           ! contando o numero de zeros para separar o joio do trigo  
            numZerosId=0
            simetriaB=simetria
      !      simetriaB=.true.
            do i=1, numVizAtualizado*nen*ndof 
               if(simetria .and. LMstencilEqId(i)<eq) then
                LMstencilEqId(i)=0;
               end if
                if ( LMstencilEqId(i) == 0 ) numZerosId=numZerosId+1
           end do
!		   print*, " numZerosID       =" ,numZerosID  
!		   print*, " numVizAtualizado*nen*ndof-numZerosID =" ,numVizAtualizado*nen*ndof-numZerosID
		   
       !B1 write(*,'(a,i0,a,81(1x,i0))'), "1.3, noh=", nc, ", LmStencilEqId, ->", LMstencilEqId(1:numVizAtualizado*nen*ndof)
       !if(eq > 2) stop
            ! guardando somente os valores diferentes de zero 
!	print'(i3,a,30i3)' ,eq, "c, LMstencilEqId=", LMstencilEqId(1:)
             !  write(*,'(4(i4 ))',  advance='NO') numVizAtualizado,nen,ndof,numZerosID; 
    !set2022           write(*,'(2(a,i0,a ))',  advance='NO') ", n=",numVizAtualizado*nen*ndof-numZerosID, ";"
               !write(*,'(2(a,i0))',  advance='YES') "numZerosId+1=",numZerosId+1, ', fim=',numVizAtualizado*nen*ndof!-numZerosID
!	print'(i3,a,30i3)' ,eq, "d, LMstencilEqId=", LMstencilEqId(numZerosId+1:numVizAtualizado*nen*ndof)
            !if(eq == 2) stop
            LMstencilEq(eq,1:numVizAtualizado*nen*ndof-numZerosID)=(LMstencilEqId(numZerosId+1:numVizAtualizado*nen*ndof))
	!print'(i3,a,30i3)' ,eq, "e, LMstencilEq  =", LMstencilEq (eq,1:)!numVizAtualizado*nen*ndof-numZerosID)         

        !call escreverLMstencilEq()
         end if !if (eq > 0) then
      end do
     end do
   return
    contains
    subroutine escreverId_Lm()
      do i=1, numnp
         write(*,'(a,i0,a,i0,a)') "noh =",i, "; id(1:ndof, n)=(/", id(1:ndof, i ),'/)' 
      end do !  
      do nel = 1, numel
          LMLocal(:)=reshape(lm(:,:,nel),(/nen*ndof/))
	  write(*,'(a,i0,a)', advance='NO') "nel=", nel, "; LM(1:ndof*nen,nel) = (/"
          do i =1, nen*ndof-1
               write(*, '(i0,a)', advance='NO') LMLocal (i), ','         
          end do 
          write(*, '(i0,a)', advance='yes') LMLocal (i), '/)'         
      end do
!      print*, " nVizinMax =", nVizinMax
!      print*, "shape(LMstencilEq)=",shape(LMstencilEq)
!      print*, "shape(LMstencilEqId)=",shape(LMstencilEqId)
     end subroutine escreverId_Lm
    subroutine escreverLMstencilEq()
            write(*,'(a,i0,a)', advance='NO') "eq=", eq, "; LMstencilEq(eq, 1:n) = (/"
            do i =1, numVizAtualizado*nen*ndof-numZerosID-1
               write(*, '(i0,a)', advance='NO') LMstencilEq (eq,i), ','         
            end do 
            write(*, '(i0,a)', advance='yES') LMstencilEq (eq,i),'/)'
		   !stop 33
       !B  write(*,'(a,i0,a,81(1x,i0))'), "em 2, eq=", eq, ", LmStencilEq, ->", LMstencilEq(eq,1 :numVizAtualizado*nen*ndof-numZerosID)
         end subroutine escreverLMstencilEq
end subroutine montarLmStencilNodal_CSR     !
    subroutine ordenarLMstencil(LMstencilEq,numCoefPorLinha)
      implicit none
      integer*4,  intent(in)    :: numCoefPorLinha      
      integer*4,  intent(inout) :: LMstencilEq(numCoefPorLinha)
      integer*8 :: menorEq, n, nn , tmp
      do n = 1, numCoefPorLinha
        menorEq=n
        do nn = n+1, numCoefPorLinha
           if((LMstencilEq(nn))<(LMstencilEq(menorEq))) menorEq = nn
           !if(abs(LMstencilEq(nn))<abs(LMstencilEq(menorEq))) menorEq = nn
        end do
        if(n == menorEq) cycle
        tmp                  = LMstencilEq(n)
        LMstencilEq(n)       = LMstencilEq(menorEq)
        LMstencilEq(menorEq) = tmp
      enddo
    end subroutine ordenarLMstencil
      !
      !**** new *************************************************************
      !
    subroutine montarPonteiroAp_CSR (Ap, LMstencilEq, numCoefPorLinha, neq, nonzeros)
      implicit none
      integer*4, intent(inout) :: Ap(:)
      integer*4, intent(in)    :: LMstencilEq(neq,numCoefPorLinha)
      integer*4, intent(in)    :: numCoefPorLinha,  neq
      integer*4, intent(inout) :: nonzeros
      integer*8 :: l,j
      ! Montando Ap
      write(*,'(a)') "em subroutine montarPonteiroAp_CSR (Ap, LMstencilEq, numCoefPorLinha, neq, nonzeros)"
      !write(*,*) neq, numCoefPorLinha, nonzeros 
      call montarListaPonteiros(Ap, LMstencilEq,neq,numCoefPorLinha)
      !B1 write(*,*) "Ap =", Ap(1:neq+1)
      ! Contando os valores nao nulos
      nonzeros=0
      do l=2, neq+1
            nonzeros=nonzeros+(Ap(l)-Ap(l-1))
      end do
      write(*,*) "neq, nonZeros, numCoefPorLinha, nonZeros/neq"
      write(*,*) neq, nonZeros, numCoefPorLinha, nonZeros/neq
    end subroutine montarPonteiroAp_CSR
      !
      !**** new *************************************************************
      !
    subroutine montarPonteiroAi_CSR (Ai, LMstencilEq, numCoefPorLinha, neq, nonzeros)
      implicit none
      integer*4,  intent(out) :: Ai(:)
      integer*4 :: nonzeros, neq
      integer*4 :: numCoefPorLinha
      !integer*4 :: LMstencilEq(neq,numCoefPorLinha)
      integer*4 :: LMstencilEq(9,numCoefPorLinha)
        call montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)
    end subroutine montarPonteiroAi_CSR
      !
      !**** new *************************************************************
      !
      subroutine montarListaPonteiros(Ap, LMstencilEq, neq, numCoefPorLinha)
      implicit none
      integer*4,  intent(in) :: neq
      integer*4,  intent(out) :: Ap(neq+1)
      integer*4,  intent(in) :: numCoefPorLinha
      integer*4,  intent(in) :: LMstencilEq(neq,numCoefPorLinha)
      integer*4 :: LMstencilEqTemp(0:numCoefPorLinha)
      integer*8 :: n, l
      integer*4 :: posPonteiro, numCoefAlhs
      posPonteiro=0
      numCoefAlhs=0
      write(*,*) " em subroutine montarListaPonteiros(Ap, LMstencilEq, neq)"
      do l=1,  neq
     !B1    write(*,'(a, 59i3)') "lmstencilEQ =", lmstencilEQ (l,:); 
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=LMstencilEq(l,:)
           posPonteiro=posPonteiro+1
           Ap(posPonteiro)=numCoefAlhs+1
         !B1  write(*,'(a, 19i3)') "lmStencilEqTemp =", LmStencilEqTemp (:) 
           do n = 1, numCoefPorLinha
            if(LMstencilEqTemp(n)==0.or.LMstencilEqTemp(n) == LMstencilEqTemp(n-1)) cycle
             numCoefAlhs = numCoefAlhs+1
           enddo    
         if(l==neq) then
            posPonteiro=posPonteiro+1
            Ap(posPonteiro)=numCoefAlhs+1
         end if
       !B1  write(*,'(a,13i3)')"parcial, Ap =", Ap (1:neq+1); 
      end do
       !B1 write(*,'(a,13i3)')"final, Ap =", Ap (1:neq+1); 
      end subroutine montarListaPonteiros
      !
      !**** new *************************************************************
      !
      subroutine montarListaIndices(Ai, LMstencilEq, neq, nonzeros, numCoefPorLinha)
      implicit none
      integer*4,  intent(in)    :: nonzeros, neq
      integer*4,  intent(inout) :: Ai(nonzeros)
      integer*4,  intent(in)    :: numCoefPorLinha
      integer*4,  intent(in)    :: LMstencilEq(neq,numCoefPorLinha)
!
      integer*4 :: LMstencilEqTemp(0:numCoefPorLinha)
      integer*4 :: i, n, posColunas, k
!
      write(*,*) "em subroutine montarListaIndices(Ai, LMstencilEq, neq, nonzeros, ..."
      posColunas=0
      do i=1, neq
         LMstencilEqTemp=0
         LMstencilEqTemp(1:numCoefPorLinha)=abs(LMstencilEq(i,:))
         if(sum(lmstencilEqTemp)>0) then 
         !B1 write(*,'(a,16(1x,i0))') "LMstencilEqTemp =", LMstencilEqTemp
         do n = 1, numCoefPorLinha
             if(LMstencilEqTemp(n).ne.LMstencilEqTemp(n-1).and.LMstencilEqTemp(n).ne.0 ) then
                 posColunas=posColunas+1
                 Ai(posColunas)= LMstencilEqTemp(n)
             end if
         enddo
         end if !if(sum(lmstencilEqTemp)>0)
      end do
     !B1 write(*,*) "Ai = ",   ai
     end subroutine montarListaIndices
!
!   end subroutine criarPonteirosMatEsparsa_CSR
!=======================================================================

      subroutine escreverSistemaAlgCSRemMTX(alhs, brhs, Ap, Ai,  nAlhs, neq, nomeArq)
      implicit none 
      real*8,  intent(in)   :: Alhs(:), brhs(:)
      integer*4,  intent(in)  :: Ap(:), Ai(:)
      integer*4,  intent(in)  :: neq, nAlhs
      character (len=*)  :: nomeArq
!234567
      integer*8 :: i, j, k
      integer*8 :: luSist = 1836 
      character(len=40), parameter :: formatoEscritaA='(1x,2(i0,1x),e23.16)'
      character(len=40), parameter :: formatoEscritaAL='(2(i0,1x),e15.8)'
      character(len=40), parameter :: formatoEscritaALB='(3(i0,1x),e15.8)'
      character(len=40), parameter :: formatoEscritaB='(e15.8)'
      character(len=40), parameter :: formatoEscritaC='(e15.8,a)'
      real*8 :: t1, t2, t3, t4
      integer*8 :: nalhsSemZeros 
!
      write(*,*) " em subroutine escreverSistemaAlgCSRemMTX(alhs, brhs, Ap, Ai, ... "
      write(*,*)" file=", trim(nomeArq)
      call timing(t1)
      open(file=nomeArq, unit=luSist) 
!
      nalhsSemZeros = nAlhs
      k = 1
      do i = 1, neq
           !write(*, *) Ap(i) ,  Ap(i+1) 
           do j = Ap(i) ,  Ap(i+1) - 1  
                if(alhs(k) .eq. 0.0d0) then
                        nalhsSemZeros = nalhsSemZeros - 1
                 !B1       write(*, formatoEscritaALB  ) k,  i, Ai(j), alhs(k)
                endif
                k = k + 1
           end do
      end do
      !do i=1,Ap(neq+1)-1
      !   !do j = Ap(i) ,  Ap(i+1) - 1  
         ! if(alhs(i) .eq. 0.0d0)  then
      !    if(alhs(i) .ne. 0.0d0) then
      !        k = k;
      !    else
      !        nalhsSemZeros = nalhsSemZeros - 1
      !        !write(*,'(1x,i0,",")', ADVANCE="NO") i
      !        write(*,'(2(1x,i0,","))') k, i, Ai(Ap(i)),   
      !        k = k+1
      !    endif
         !end do
      !end do

!      write(luSist,'(a)')'%% matriz A de coeficientes reais simetrica positiva definida ' 
      write(luSist,'(a)')'%%MatrixMarket matrix coordinate real symmetric' 
      write(luSist,'(a)')'% produzida pelo metodo classico de galerkin para o metodo de elementos finitos '
      write(luSist,'(a)')'% armazenamento esparso CSR '
      write(luSist,'(a,3(i0,a))' )  '%  matriz ',neq, 'X',neq, ' com ', nAlhs, ' coefs diferentes de zero' 
      write(luSist,'(  3(i0,1x))' )  neq, neq, nalhsSemZeros 
      write(*,'(  3(a,i0,1x))' )  "neq = ", neq, ", nAlhs=", nAlhs, ", menos alguns zeros ", nalhsSemZeros 
      !write(*,*) size(alhs), size(brhs), size(Ap), size(Ai)  
      k = 1
      do i = 1, neq
           !write(*, *) Ap(i) ,  Ap(i+1) 
           do j = Ap(i) ,  Ap(i+1) - 1  
                if(alhs(k) .ne. 0.0d0) then
                        write(luSist, formatoEscritaAL  ) i, Ai(j), alhs(k)
                        if(mod(j,nAlhs/10)==1.or.j==nAlhs) write(*, formatoEscritaALB  ) k,  i, Ai(j), alhs(k)
                endif
              !  write(*, formatoEscritaA   ) i, Ai(j), alhs(k)
                !write(luA+luAux, '( 3(i10,2x), e20.10, 2x,i5)'     ) i, j,  Ai(j), alhs(k), k
                k = k + 1
           end do
      end do
!      write(luSist,'( (a,i0,a))') '% lado direito com ',neq,' elementos' 
!      write(luSist,'(  3(i10))' )  neq 
      i = 1
      write(luSist, formatoEscritaC )  brhs(i), "  BRHS"
      do i = 2, neq
          write(luSist, formatoEscritaB )  brhs(i)
      end do
      close(luSist)
      call timing(t2)
      write(*,*) " tempo de escrita =", t2 - t1

     end subroutine escreverSistemaAlgCSRemMTX

     subroutine lerSistemaAlgMTXemCSR(alhs, brhs, Ap, Ai, neq, nonzeros,  nomeArq) 
     implicit none 
     real*8   ,intent(inout) :: alhs(:) ! matriz do sistema 
     real*8   ,intent(inout) :: brhs(:) ! lado direito do sistema
     integer*4,intent(inout) :: Ap(:) ! vetor apontadores para  inicio de uma linha
     integer*4,intent(inout) :: Ai(:) ! vetor com as colunas dos coefcientes de alhs
     integer*4,intent(inout) :: nonzeros
     integer*4,intent(inout) :: neq 
     character(len=*), intent(in) :: nomeArq
!
     integer   :: luSist = 1000, luSistLido = 2000 
     integer*8 :: i, j, k, ieq, ieqAnterior
     character(len=100) :: label 
     character(len=40), parameter :: formatoEscritaALB='(3(i0,1x),e15.8)'
     character(len=40), parameter :: formatoEscritaA='(2(i0,1x),e15.8)'
     character(len=40), parameter :: formatoEscritaB='(1(i0,1x),e15.8)'
!
     write(*,*) " em subroutine lerSistemaMTXemCSR(alhs, brhs, ..."
     write(*,*) " FUNCIONANDO PARA MATRIZES SIMETRICAS COM ELEMENTOS FORNECIDOS POR LINHAS"
     open(file=nomeArq, unit=luSist) 

     read(luSist,'(a)' ) label; write(*,* ) label 
     read(luSist,'(a)' ) label; write(*,* ) label 
     read(luSist,'(a)' ) label; write(*,* ) label 
     read(luSist,'(a)' ) label; write(*,* ) label 
     read(luSist,* )     neq, neq, nonzeros
     write(*,'(a,3(i10))' ) "+++", neq, neq, nonzeros

! o arquivo que será lido estah escrito com os coeficientes de uma linha 
!      considerando matrix simetrica
! as linhas estao em ordem crescente de 1 ateh neq (num de equacoes) 

     Ap(:) = 1
     k = 0; ieqAnterior = 0
     do i = 1, nonzeros
          k = k + 1 ! contador de nao zeros lidos
          read(luSist, * ) ieq,  Ai(k), alhs(k)
          if(mod(k,nonZeros/10)==1.or.k==nonZeros) write(*, formatoEscritaALB  ) k,  ieq, Ai(k), alhs(k)
          if(ieq > ieqAnterior) Ap(ieq) = k  
          ieqAnterior = ieq
     end do
     k = k + 1
     Ap(ieq+1) = k
   
     write(*,'( (a,i0,a))') 'lado direito com ',neq,' elementos' 
     do i = 1, neq
         read(luSist, * ) brhs(i)
     end do

      write(luSistLido,'(  3(i0,1x))' )  neq, neq, nonzeros 
      k = 1
      do i = 1, neq
           do j = Ap(i) ,  Ap(i+1) - 1  
                write(luSistLido, formatoEscritaA   ) i, Ai(j), alhs(k)
                k = k + 1
           end do
      end do
      do i = 1, neq
          write(luSistLido, formatoEscritaB ) i, brhs(i)
      end do

     close(luSist); 
     end subroutine lerSistemaAlgMTXemCSR
     
     end module mEstruturasDadosSistEq

