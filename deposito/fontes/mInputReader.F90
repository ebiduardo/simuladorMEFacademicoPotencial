!!>
!!         programa de elementos finitos em fortran 90
!!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!!
!!         Pilato Jr , Vinicius Pessoa  {pilato@deepsoft.com.br, pessoa@deepsoft.com.br}
!!
!!         Desenvolvido por DeepSoft para LNCC/MCT
!!         Rio de Janeiro, 11.2013-04.2014
!!
!!=================================================================================

!> Modulo responsavel por reunir subrotinas para leitura do arquivo de entrada.
module mInputReader


    implicit none
    !> Armazena as linhas do arquivo de input.
    character(len=200), allocatable :: file_lines(:)
    !> Armazena o numero de linhas no arquivo.
    integer*4 number_of_lines

    contains

    !> Leitura e geracao de valores de condicoes de contorno.
    !!
    !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
    !! @param f             f.
    !! @param ndof          ndof
    !! @param numnp         numnp
    !! @param j             j
    !! @param nlvect        nlvect
    !! @param iprtin        iprtin
    subroutine leituraValoresCondContornoDS(keyword_name,f,ndof,numnp,j, nlvect, iprtin, iecho)

        use mleituraEscrita, only: printf, printd
        use mGlobaisEscalares, only : zero

        implicit none

        integer*4 :: ndof, numnp, j, nlvect, iprtin, iecho, keyword_line
        character(len=50) :: keyword_name
        real*8 f(ndof,numnp,nlvect)

        logical lzero
        integer*4 nlv
        character(len=35) :: rotulo

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.-1) return

        f=0.0

        do 100 nlv=1,nlvect
            call genflDS(f(1,1,nlv),ndof,keyword_line)
            !call ztest(f(1,1,nlv),ndof*numnp,lzero)
            lzero = sum(f(:,:,nlv)) == zero
!
            if (iprtin.eq.0) then
                 if (lzero) then
                    if (j.eq.0) write(iecho,1000) nlv
                    if (j.eq.1) write(iecho,2000)
                 else
                    if (j.eq.0) call printf(f,ndof,numnp,nlv,iecho)
                    if (j.eq.1) then
                        rotulo=" n o d a l  b o d y  f o r c e s "
                        call printd (rotulo, f,ndof,numnp,iecho)
                    end if
                endif
            endif
        100 continue
        return
 1000 format('1'//,' there are no nonzero prescribed forces and ',&
         'kinematic boundary conditions for load vector number ',i10)
 2000 format('1'//,' there are no nonzero nodal body forces')
    end subroutine leituraValoresCondContornoDS !**********************************************************************


    !> Leitura e geracao de codigos de condicoes de contorno.
    !!
    !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
    !! @param id            id.
    !! @param ndof          numnp
    !! @param numnp         numnp
    !! @param neq           neq
    !! @param iecho         iecho
    !! @param iprtin        iprtin
    subroutine leituraCodigosCondContornoDS(keyword_name, id, ndof, numnp, neq, iecho, iprtin)

        integer*4, intent(in) :: ndof, numnp, iecho, iprtin
        integer*4, intent(out) ::  neq
        integer*4, intent(inout) :: id(ndof,numnp)
        character(len=50) keyword_name

        integer*4 :: keyword_line
        integer*4 :: nn, n, i
        logical pflag

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.-1) return
        id = 0
        call igenDS(id,ndof, keyword_line)
        if (iprtin.eq.0) then
            nn=0
            do n=1,numnp
                pflag = .false.
                do i=1,ndof
                    if (id(i,n).ne.0) pflag = .true.
                end do
                if (pflag) then
                    nn = nn + 1
                    if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
                    write(iecho,2000) n,(id(i,n),i=1,ndof)
                endif
            end do
        endif
!.... establish equation numbers
      neq = 0
      do  n=1,numnp
      do  i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif
      end do
      end do

    return

      print *, "id =", id
      do  n=1,numnp
      do  i=1,ndof
      if (id(i,n).ne.0) then
         id(i,n) = neq - id(i,n) + 1
      endif
      end do
      end do
      print *, "id =", id
!
    return
!
1000 format(' n o d a l   b o u n d a r y   c o n d i t i o n &
     &         c o  d e s'/// &
      5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
!
    end subroutine leituraCodigosCondContornoDS !********************************************************************




    !> Leitura e geracao de codigos de condicoes de contorno.
    !!
    !! @param keyword_name  Palavra-chave que indica o tipo da condicao de contorno.
    !! @param id            id.
    !! @param ndof          numnp
    !! @param numnp         numnp
    !! @param neq           neq
    !! @param iecho         iecho
    !! @param iprtin        iprtin
    subroutine leituraCodigosCondContornoDSF(keyword_name, id, ndof, numnp, neq, iecho, iprtin)
!        use mInputReader,     only: findKeyword, igenDS


        integer*4, intent(in) :: ndof, numnp, iecho, iprtin
        integer*4 ::  neq, keyword_line
        integer*4, intent(inout) :: id(ndof,numnp)
        character(len=50) keyword_name

        integer*4:: nn, n, i
        logical pflag

        write(*,*) " em subroutine leituraCodigosCondContornoDSF(keyword_name,"

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.-1) return

        id = 0
        call igenDS(id,ndof, keyword_line)
!
        if (iprtin.eq.0) then
            nn=0
            do 200 n=1,numnp
                pflag = .false.
!
                do 100 i=1,ndof
                    if (id(i,n).ne.0) pflag = .true.
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
                 neq = neq + 1
                if (id(i,n).eq.0) then
                !    neq = neq + 1
                    id(i,n) = neq
                else
                    !id(i,n) = 1 - id(i,n)
                    id(i,n) = -neq !1 - id(i,n)
                endif
!
            300 continue
!
               !     write(*,2000) n,(id(i,n),i=1,ndof)
        400 continue
!
    return
!
1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n &
     &         c o  d e s'/// &
      5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
!
    end subroutine leituraCodigosCondContornoDSF !********************************************************************

    !> Leitura e geracao de coordenadas
    !!
    !! @param x          Matriz onde serao armazeadas as coordenadas.
    !! @param nsd        Corresponde ao nao de linhas da matriz x.
    !! @param numnp      Numero de xxx
    !! @param iprtin     iprtin
    subroutine leituraGeracaoCoordenadasDS(x, nsd, numnp, iprtin, icoords)
            
      implicit none

      integer*4, intent(in)   :: nsd, numnp, iprtin, icoords
      real*8, intent(inout) ::  x(nsd,*)
      integer*4:: i, n, keyword_line
      character(len=50) keyword_name

      keyword_name = 'coordenadas_nodais'
      keyword_line = findKeyword(keyword_name)
      call genflDS(x,nsd,keyword_line)

      if (iprtin.eq.1) return 
! 
      open(unit=icoords    , file= 'coordenadas.dat')
      do n=1,numnp 
         if (mod(n,50).eq.1) write(icoords,1000) (i,i=1,nsd) 
         write(icoords,2000) n,(x(i,n),i=1,nsd)   
      enddo 
! 
      return 
! 
 1000 format(///,' n o d a l   c o o r d i n a t e   d a t a '///5x, &
     &' node no.',3(13x,' x',i1,' ',:)//) 
 2000 format(6x,i10,10x,3(1pe15.8,2x)) 
    end subroutine leituraGeracaoCoordenadasDS!***************************************************************************************************


    !---------------------------------------------------------------------------------
    !> Le arquivo de input e armazena seu conteudo em um array.
    !! @param file_name Nome do arquivo a ser lido.
    subroutine readInputFileDS(iin)
        implicit none
        integer :: iin

        character(len=200) file_line
        character(len=200), allocatable :: main_file_lines(:)
        character(len=200), allocatable :: include_files(:)
        character(len=200)              :: include_file

        integer*4 :: success, lines_count
        integer*4 :: main_number_of_lines, main_number_of_includes
        integer*4 :: include_index, inc_nlines, inc_inc, merge_lines, i
        integer*4, allocatable :: include_indexes(:), include_number_of_lines(:)

        main_number_of_lines    = 0
        main_number_of_includes = 0

        call analyzeFileInput(main_number_of_lines, main_number_of_includes, iin)

        if (main_number_of_includes.eq.0) then
            call createSimpleInputFile(iin)
            return
        end if

        allocate(main_file_lines(main_number_of_lines))
        allocate(include_indexes(main_number_of_includes))
        allocate(include_number_of_lines(main_number_of_includes))
        allocate(include_files(main_number_of_includes))

        lines_count = 1
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            main_file_lines(lines_count) = file_line
            lines_count = lines_count + 1
        end do
        rewind(iin)

        number_of_lines = main_number_of_lines

        !Number of lines
        do i=1, main_number_of_includes
            include_index = findInclude(i, main_file_lines, main_number_of_lines)
            read(main_file_lines(include_index), '(A)') include_file

            include_file = adJustl(include_file)
            call analyzeFile(include_file, inc_nlines, inc_inc)
            number_of_lines = number_of_lines + inc_nlines
            include_indexes(i) = include_index
            include_number_of_lines(i) = inc_nlines
            include_files(i) = include_file
        end do

        !Prepare final struct.
        call prepareFileLines(include_indexes, include_number_of_lines, main_number_of_includes, main_file_lines)

        !Merge contensts.
        merge_lines = 0
        do i=1, main_number_of_includes
            call mergeIncludeContents(include_files(i), include_indexes(i) + merge_lines)
            merge_lines = merge_lines + include_number_of_lines(i)
        end do

        deallocate(main_file_lines)
        deallocate(include_indexes)
        deallocate(include_number_of_lines)
        return
    end subroutine readInputFileDS !***********************************************************************************


    !> Cria a estretura de input usando um arquivo de entrada sem includes
    !! @param file_name Nome do arquivo a ser lido.
    subroutine createSimpleInputFile(iin)
        implicit none
        integer :: iin

        integer*4          :: success, lines_count
        character(len=200) :: file_line
        number_of_lines = 0
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
        end do
        rewind(iin)

        allocate(file_lines(number_of_lines))
        !TO-DO avoid two-times read
        lines_count = 1
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            file_lines(lines_count) = file_line
            lines_count = lines_count + 1
        end do
        rewind(iin)
    end subroutine createSimpleInputFile !*****************************************************************************

    !> Le o conteudo do arquivo de include e armazena no array principal.
    !! @param   include_index   O index do include.
    !! @param   include_files   Array com includes.
    !! @param   include_line    A linha do include.
    subroutine mergeIncludeContents(include_file, include_line)
        implicit none
        integer*4          ::  include_line
        character(len=200) ::  include_file

        character(len=200) file_line
        integer*4 file_channel, success, current_index

        file_channel = 1
        current_index = include_line

!         open(unit=file_channel, file=include_file)       
        open(unit=file_channel, file=include_file, status='old', err=100)             
        do
            read(file_channel, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            file_lines(current_index) = file_line
            current_index = current_index + 1
        end do
        close(file_channel)
        return
        100 print*, "falta o arquivo: ", include_file, ", mencionado em inputDS.dat"
        
    end subroutine mergeIncludeContents !******************************************************************************

    !> Efetua a aloca��o da estrutura definitiva, preparando a linha dos arquivos originais para receber os includes
    !! @param   include_indexes             Array os indices de ocorrencias dos includes.
    !! @param   include_number_of_lines     Array com o numero de linhas de cada include
    !! @param   number_of_includes          Numero de includes.
    !! @param   original_file_lines         Linhas do arquivo de entrada original.
    subroutine prepareFileLines(include_indexes, include_number_of_lines, number_of_includes, original_file_lines)
        integer*4 number_of_includes, number_of_original_lines, line_index, shift_lines
        integer*4 include_indexes(:), include_number_of_lines(:)
        character(len=200) original_file_lines(:)

        integer*4 current_include_index, original_index

        allocate(file_lines(number_of_lines))

        current_include_index = 1
        original_index = 1
        line_index = 1
        shift_lines = 0
        do while ( line_index <= number_of_lines)
            if (original_index.eq.(include_indexes(current_include_index))) then
                line_index = line_index + include_number_of_lines(current_include_index)
                current_include_index = current_include_index + 1
            end if
            file_lines(line_index) = original_file_lines(original_index)
            line_index = line_index + 1
            original_index = original_index + 1
        end do


    end subroutine prepareFileLines !**********************************************************************************

    !> Efetua algumas an�lises no arquivo recebido.
    !! @param   number_of_lines     N�mero de linhas.
    !! @param   number_of_include   N�mero de ocorr�ncias da palavra include.
    subroutine analyzeFileInput(number_of_lines, number_of_includes, iin)

        integer*4 number_of_lines, number_of_includes
        integer :: iin

        character(len=50) include_keyword, formated_keyword
        character(len=200) file_line
        integer*4 keyword_len, success

        include_keyword  = "include"
        keyword_len      = len(trim(include_keyword)) + 2
        formated_keyword = trim('*' // trim(include_keyword) // '{')

        number_of_lines    = 0
        number_of_includes = 0
!        lunitInicial = 15
        do
            read(iin, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                number_of_includes = number_of_includes + 1
            end if
        end do
        rewind(iin)

    end subroutine analyzeFileInput !***************************************************************************************

    !> Efetua algumas an�lises no arquivo recebido.
    !! @param   file_name           O nome do arquivo.
    !! @param   number_of_lines     N�mero de linhas.
    !! @param   number_of_include   N�mero de ocorr�ncias da palavra include.
    subroutine analyzeFile(file_name, number_of_lines, number_of_includes)
        character(len=200) file_name, file_line
        integer*4 number_of_lines, number_of_includes

        character(len=50) include_keyword, formated_keyword
        integer*4 keyword_len, file_channel, success, lunitInicial

        include_keyword = "include"
        keyword_len = len(trim(include_keyword)) + 2
        formated_keyword = trim('*' // trim(include_keyword) // '{')

        number_of_lines = 0
        number_of_includes = 0

        file_channel = 2
        lunitInicial = 15
        file_channel = lunitInicial

!         open(unit=file_channel, file=file_name, status='old', err=100)        
        open(unit=file_channel, file=file_name, status='old', err=100)   
        do
            read(file_channel, "(A)", iostat=success) file_line
            if (success.ne.0) exit
            number_of_lines = number_of_lines + 1
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                number_of_includes = number_of_includes + 1
            end if
        end do
        close(file_channel)
        return
        100 print*, "file", file_name, " - incluido no arquivo inputDS.dat - não existe"
   
    end subroutine analyzeFile !***************************************************************************************

    !> Procura a n-esima palavra-chave include.
    !! @param  position         Corresponde a posicao desejada.
    !! @param  file_lines       Linhas do arquivo.
    !! @param  number_of_lines  Numero de linhas atuais.
    !! @return O indice da palavra-chave no array que contem as linhas do arquivo de entrada.
    integer*4 function findInclude(position, file_lines, number_of_lines)
        implicit none
        integer*4 position, number_of_lines, current_position
        character(len=200) file_lines(:)
        character(50) keyword, formated_keyword
        character(len=120) file_line
        integer*4 i, keyword_len

        keyword = "include"
        keyword_len = len(trim(keyword)) + 2
        formated_keyword = trim('*' // trim(keyword) // '{')
        current_position = 0

        do i=1, number_of_lines
            file_line = file_lines(i)
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                current_position = current_position + 1
                if (current_position.eq.position) then
                    findInclude = i + 1
                    return
                end if
            end if
        end do
        findInclude = 0
        return
    end function findInclude !*****************************************************************************************

    !> Procura uma palavra-chave.
    !! @param  keyword A palavra-chave.
    !! @return O indice da palavra-chave no array que contem as linhas do arquivo de entrada.
    integer*4 function findKeyword(keyword)
        implicit none
        character(50) keyword, formated_keyword
        character(len=120) file_line
        integer*4 i, keyword_len
        findKeyword = 0
        do i=1, number_of_lines, 1
            file_line = file_lines(i)
            keyword_len = len(trim(keyword)) + 2
            formated_keyword = trim('*' // trim(keyword) // '{')
            if (formated_keyword.eq.file_line(1:keyword_len)) then
                findKeyword = i + 1
                return
            end if
        end do
        return
    end function findKeyword !*****************************************************************************************

    !> Efetua a leitura de uma palavra-chave to tipo inteiro. Se nao encontrado, associa o valor defualt fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde o valor inteiro sera atribuido.
    !! @param default_value Valor default.
    subroutine readIntegerKeywordValue(keyword, target, default_value, ierr)
        implicit none
        character(50) :: keyword
        integer*4     :: target, default_value
        integer       :: ierr
!
        character(120) file_line
        integer*4 :: keyword_line = 0
        character(20) :: origem  
        target = default_value
        origem = 'valor default, ' ;ierr=1
        keyword_line = findKeyword(keyword)
        if (keyword_line.ne.0) then
            file_line = adjustL(trim(file_lines(keyword_line)))
            read(file_line, *) target
            origem='valor lido, ';    ierr=0
        end if
        write(*,'(a, a, a, i0)') origem, trim(keyword), '=', target 
        return
    end subroutine readIntegerKeywordValue !***************************************************************************

    !> Efetua a leitura de uma palavra-chave to tipo string. Se nao encontrado, associa o valor defualt fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde a string sera atribuido.
    !! @param default_value Valor default.
    subroutine readStringKeywordValue(keyword, target, defaultValue, ierr)
        implicit none
        character(50) keyword
        character(120) file_line
        character(len=*) :: target, defaultValue
        integer :: ierr
        !
        integer*4 keyword_line
        character(20) :: origem  
        
        origem = 'valor default, ' 
!        print*, "em readStringKeywordValue, para ler: ", keyword
        target = defaultValue
!        print*, value
        ierr=1
        keyword_line = findKeyword(keyword)
        !print*, "origem ............................= ", origem
        !print*, " ..... keyword_line=",keyword_line
        if (keyword_line.ne.0) then
            file_line = adjustL(trim(file_lines(keyword_line)))
            read(file_line, '(a)') target
            origem='valor lido, ';    ierr=0
        end if
        !read(file_lines(keyword_line), '(a)') value
        write(*,'(a, a, a, a)') origem, trim(keyword), '=', target 
       ! if(defaultValue=="YES") stop
        return
    end subroutine readStringKeywordValue 
!****************************************************************************
    !> Efetua a leitura de uma palavra-chave to tipo real. Se nao encontrado, associa o valor defualt fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde o real sera atribuido.
    !! @param default_value Valor default.
    subroutine readRealKeywordValue(keyword, target, default_value, ierr)
      implicit none
      character(50) keyword
      real(8) target
      real(8) default_value
      integer*4 keyword_line
      integer :: ierr

      character(120) file_line
      character(20) :: origem 
      target = default_value
      origem = 'valor default, ' 
      ierr=1
      keyword_line = findKeyword(keyword)
      if (keyword_line.ne.0) then
          file_line = adjustL(trim(file_lines(keyword_line)))
          read(file_line, *) target
          origem='valor lido, ';    ierr=0
        end if
      write(*,'(a, a, a, e15.7)') origem, trim(keyword), '=', target 
      return
    end subroutine readRealKeywordValue
!****************************************************************************
    !> Efetua a leitura de valores definidos para propriedades de sa�da. Na pr�tica s�o lidas 3 vari�veis.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        A variavel destino.
    !! @param var_out       O valor da variavel out.
    !! @param var_n         Valor n.
    subroutine readOutFlagKeyword(keyword, targetb, var_out, var_n, dominio, ierr)
        implicit none
        
        character(50) keyword
        integer*4 targetb, var_n
        character(len=128) :: var_out
        character(len=*) ::  dominio

        integer*4 keyword_line
        character(len=120) :: file_line
        integer :: ierr
            character*1 tab
            tab = char(11)
!
        keyword_line = findKeyword(keyword)
        if (keyword_line.eq.0) then
            ierr=1
            return
        end if
        read(file_lines(keyword_line), *) targetb
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line), "(a)") var_out
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line), *) var_n
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line), "(a)") dominio
        write(*,'(a, a, a, 1x, 2(1x,i5,1x, a) )') 'valores lidos, ', trim(keyword),'=', &
                 targetb, trim(var_out), var_n, trim(dominio)
        
    end subroutine readOutFlagKeyword

    ! ALTERADO PAT
    ! NOVO !    
    !! Efetua a leitura de uma palavra-chave to tipo logico. 
    !! Se nao encontrado, associa o valor defualt fornecido.
    !! @param keyword       A palavra-chave a ser encontrada.
    !! @param target        Variavel onde o real sera atribuido.
    !! @param default_value Valor default.
    subroutine readLogicalKeywordValue(keyword, target, default_value, ierr)
       implicit none
       character(50) keyword
       character(120) file_line
       logical target, default_value
       integer :: ierr

       integer*4 keyword_line

       character(20) :: origem

       target = default_value
       origem = 'valor default, '
       ierr=1

       keyword_line = findKeyword(keyword)
       if (keyword_line.ne.0) then
       file_line = adjustL(trim(file_lines(keyword_line)))
       read(file_line, *) target
       origem='valor lido, ';    ierr=0
       end if
       write(*,'(a, a, a, L)') origem, trim(keyword), '=', target
       return
    end subroutine readLogicalKeywordValue !****************************************************************************
    
    !> Efetua a geracao de coordeadas, de acordo com parametros.
    !!
    !! @param a      Matriz onde serao armazenados os dados.
    !! @param nra    Inteiro indicando nra
    !! @param nLinhaArqInput    Indice da linha onde as coordenadas estao posicionadas no array linhas no arquivo de entrada.

    subroutine genflDS(a,nra, nLinhaArqInput )
      use mMalha, only: genfl1   
!      use mGlobaisEscalares

      implicit none

      integer*4:: nra, nLinhaArqInput
      real*8  :: a(nra,*)
      real*8  :: temp(6,20)
      integer*4:: n,numgp,ninc(3),inc(3)
      integer*4:: i, j, m, mgen
      character(80) :: novaLinha
    

      write(*,*) " em subroutine genflDS(a,nra, nLinhaArqInput )"
      !write(*,*)"nLinhaArqInput=", nLinhaArqInput, ", file_lines =", trim(file_lines(nLinhaArqInput) )
      100 continue
      novaLinha=trim(file_lines(nLinhaArqInput))
      !write(*,*) novaLinha, "..."
      !read(novaLinha,1002) n,numgp,(temp(i,1),i=1,nra)
      read(file_lines(nLinhaArqInput),1002) n,numgp,(temp(i,1),i=1,nra)
      1002 format(2i10,6f10.7)
      !write(*,*) n,numgp,(temp(i,1),i=1,nra)
      nLinhaArqInput = nLinhaArqInput + 1
      !write(*,*)"nLinhaArqInput=", nLinhaArqInput, ", file_lines =", trim(file_lines(nLinhaArqInput) )

      if (n.eq.0) return
     ! call move(a(1,n),temp,nra)
      a(1:nra,n) = temp(1:nra,1)

      if (numgp.ne.0) then
      do 200 j=2,numgp
        read(file_lines(nLinhaArqInput),1002) m,mgen,(temp(i,j),i=1,nra)
         nLinhaArqInput = nLinhaArqInput + 1
      !   write(*,*)"nLinhaArqInput=", nLinhaArqInput, ", file_lines =", trim(file_lines(nLinhaArqInput) )
        if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m) 
        ! temp(1:nra,j)=a(1:nra,m) 
      200 continue

         read(file_lines(nLinhaArqInput),2000) (ninc(i),inc(i),i=1,3)
         !read(file_lines(nLinhaArqInput),2000) (ninc(i),inc(i),i=1,3)
         nLinhaArqInput = nLinhaArqInput + 1

         call genfl1(a,nra, temp, n, numgp, ninc, inc)
      endif
      go to 100

      1000 format(2i10,6e15.7)
      1001 format(20x,6e15.7)
      2000 format(16i10)

    end  subroutine !**************************************************************************************************

    subroutine genflDS_ptj(a,nra, nLinhaArqInput )
      use mMalha, only: genfl1   
      use mGlobaisEscalares

      implicit none

      integer*4:: nra, nLinhaArqInput
      real*8  :: a(nra,*)
      real*8  :: temp(6,20)
      integer*4:: n,numgp,ninc(3),inc(3)
      integer*4:: i, j, m, mgen

      100 continue
      read(file_lines(nLinhaArqInput),1000) n,numgp,(temp(i,1),i=1,nra)
      nLinhaArqInput = nLinhaArqInput + 1

      if (n.eq.0) return
     ! call move(a(1,n),temp,nra)
      a(1:nra,n) = temp(1:nra,1)

      if (numgp.ne.0) then
         do 200 j=2,numgp
         read(file_lines(nLinhaArqInput),1000) m,mgen,(temp(i,j),i=1,nra)
         nLinhaArqInput = nLinhaArqInput + 1

        if (mgen.ne.0) temp(1:nra,j)=a(1:nra,m)         !B
         200    continue
         read(file_lines(nLinhaArqInput),2000) (ninc(i),inc(i),i=1,3)
         nLinhaArqInput = nLinhaArqInput + 1

         call genfl1(a,nra, temp, n, numgp, ninc, inc)
      endif
      go to 100

      1000 format(2i10,6f15.0)
      2000 format(16i10)

    end  subroutine !**************************************************************************************************

    !> Subrotina para ler e gerar dados nodais inteiros.
    !! @param ia              Array de entrada.
    !! @param m               Numero de linhas na matriz de entrada.
    !! @param nLinhaArqInput  Indice da linha no array de linhas do arquivo de entrada.
    subroutine igenDS(ia, m, nLinhaArqInput)
!        use mGlobaisEscalares

        integer*4:: m, ia(m,*), nLinhaArqInput
        integer*4:: ib(m)
        integer*4:: n, ne, ng
        integer*4:: i
        !
        100 continue
        read(file_lines(nLinhaArqInput),1000) n,ne,ng,(ib(i),i=1,m)
        nLinhaArqInput = nLinhaArqInput + 1

        if (n.eq.0) return

        if (ng.eq.0) then
            ne = n
            ng = 1
        else
            ne = ne - mod(ne-n,ng)
        endif
        !
        do 200 i=n,ne,ng
        ia(:,i)=ib
        200 continue
        !
        go to 100
        !
        1000 format(16i10)
    end subroutine igenDS !********************************************************************************************

    !> Subrotina responsavel por ler e gerar conectividades nodais e ladais.
    !> @param keyword_name  O nome da keyword associada.
    !> @param conecElem     C�digo do elemento
    !> @param mat           C�digo do material
    !> @param nen           N�mero de elementos.
    subroutine leituraGeracaoConectividadesDs(keyword_name, conecElem, mat, nen, ierr)
        use mleituraEscrita , only: genel1
        integer*4:: n,nel(3),incel(3),inc(3)

        integer*4:: nen
        integer*4:: conecElem(nen,*),mat(*) 
        integer*4:: itemp(27) !B
        character(len=50) keyword_name
        !
        integer*4:: m,ng,i, keyword_line
        integer :: ierr

!        write(*,'(a)') " em subroutine leituraGeracaoConectividadesDS(keyword_name, conecElem, mat, nen)"
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
           ierr=1
           return
        endif

        100 continue
        read(file_lines(keyword_line),1000) n,m,(itemp(i),i=1,nen),ng
        keyword_line = keyword_line + 1

        if (n.eq.0) return
        !call imove(conecElem(1,n),itemp,nen)
        conecElem(1:nen,n)=itemp(1:nen)
        mat(n)=m
        if (ng.ne.0) then
            !Generate data
            read(file_lines(keyword_line),1000) (nel(i),incel(i),inc(i),i=1,3)
            keyword_line = keyword_line + 1
            call genel1(conecElem,mat,nen,n,nel,incel,inc)
        else

        endif
        go to 100
        
        1000 format(16i10,10x,14i10)

    end subroutine !***************************************************************************************************

    !> Subrotina respons�vel por ler e gerar elementos de face.
    !> @param keyword_name     A plavra-chave associada.
    !> @param conecElem        C�digo dos n�s dos elementos.
    !> @param nen              N�mero de lementos.
    !> @param nelx             N�mero de elmentos em x.
    !> @param nely             N�mero de elmentos em y.
    !> @param nelz             N�mero de elmentos em z.
    subroutine genelFacesDS(keyword_name, conecElem, nen, nelx, nely, nelz, ierr)
        use mGlobaisEscalares
        use mMalha, only: numel

        implicit none
        integer*4:: nen, nelx, nely, nelz
        integer*4:: conecElem(nen,*)
        character(len=50) keyword_name
        integer*4:: ierr

        integer*4:: ng, n, m, nel, i
        integer*4:: condicao, condicao2, keyword_line
        
        print*, "nelx", nelx
        print*, "nely", nely

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
           ierr=1
           return
        endif

        read(file_lines(keyword_line),1000) n,m,(conecElem(i,1),i=1,nen),ng
        keyword_line = keyword_line + 1
!
        condicao=0
        condicao2=0

        do nel=2, numel
            if(condicao==0.and.condicao2==0) then
                do i=1, nen
                    conecElem(i,nel)=conecElem(i,nel-1)+1
                end do
            else
                if(condicao==1.and.condicao2==0) then
                    do i=1, nen
                        if(i<=4) then
                            conecElem(i,nel)=conecElem(i,nel-1)+(nelx+1)+1
                        else
                            conecElem(i,nel)=conecElem(i,nel-1)+1
                        endif
                    end do
                else
                    if(condicao==1.and.condicao2==1) then
                        do i=1, nen
                            if(i<=4) then
                                conecElem(i,nel)=conecElem(i,nel-1)+(nelx+1)+nelx+(nelx*nely)+1
                            else
                                conecElem(i,nel)=conecElem(i,nel-1)+(nelx*(nely+1))+(nely*(nelx+1))+1
                            endif
                        enddo
                    end if
                end if
            end if
            if(mod(nel, nelx)==0) then
                condicao=1
            else
                condicao=0
            end if
            if(mod(nel, nelx*nely)==0) then
                condicao2=1
            else
                condicao2=0
            end if
        end do
        1000 format(16i10,10x,14i10)
    end subroutine


    subroutine readNodeElementsDS
        use mGlobaisEscalares
        use mMalha,          only: nen
        implicit none
        character(len=50) keyword_name
        integer :: ierr

        keyword_name = "ntype"
        call readIntegerKeywordValue(keyword_name, ntype, 1, ierr)
        keyword_name = "numat"
        call readIntegerKeywordValue(keyword_name, numat, 1, ierr)
        keyword_name = "nen"
        call readIntegerKeywordValue(keyword_name, nen, 4, ierr)
        keyword_name = "nicode"
        call readIntegerKeywordValue(keyword_name, nicode, nen, ierr)
    end subroutine !***************************************************************************************************

    
 !***************************************************************************************************
    !> Efetua a leitura de propriedades de materiais.
    !> @param keyword_name  Keyword especifica das  propriedades de materiais.
    subroutine readMaterialPropertiesDS(keyword_name, iecho, ierr)
        use mGlobaisEscalares
        use mGlobaisArranjos
        
        implicit none
        character(len=50) keyword_name
        integer*4 n, m, i, keyword_line
        integer :: iecho, ierr

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
           ierr=1
           return
        endif

        write(*,*) "lido, ", keyword_name

        do 400 n=1,numat
        if (mod(n,50).eq.1) write(iecho,4000) numat
        read (file_lines(keyword_line),  5000) m,(c(i,m),i=1,3)
        keyword_line = keyword_line + 1
        write(iecho,6000) m,(c(i,m),i=1,3)
        write(*,6000) m,(c(i,m),i=1,3)
        400 continue
        5000  format(i10,5x,5f10.0)
        6000  format(2x,i3,1x,5(1x,1pe11.4))
        4000  format(///,&
                ' m a t e r i a l   s e t   d a t a                      ',  //5x,&
                ' number of material sets . . . . . . . . . . (numat ) = ',i10///,2x,'set',4x,'Kx ',&
                10x,'Ky',10x,'Kz')


    end subroutine readMaterialPropertiesDS !**************************************************************************

    !> Faz a leitura de constant body forces.
    !> @param keyword_name  Keyword especifica para constant body forces.
    subroutine readConstantBodyForcesDS(keyword_name, iecho, ierr)
        use mGlobaisArranjos, only: grav

        implicit none
        character(len=50) keyword_name
        integer*4 i, keyword_line
        integer :: iecho, ierr

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
           ierr=1
           return
        endif

        write(*,*) "lido, ", keyword_name
        read  (file_lines(keyword_line),  7000) (grav(i),i=1,3)
        write (iecho,8000) (grav(i),i=1,3)
        write (*,8000) (grav(i),i=1,3)

        7000 format(8 f10.0)
        8000 format(///,&
         ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
         ' direcao 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
         ' direcao 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
         ' direcao 3............................ = ',      1pe15.8,//)
    end subroutine readConstantBodyForcesDS 
    
    !**************************************************************************
    subroutine leituraRegiaoDs(keyword_name, geoform, numel, ierr)
        character(len=50) keyword_name
        integer :: numel
        CHARACTER*12 :: GEOFORM(numel)
!
        integer*4:: m,ng,i, keyword_line
        integer :: ierr

        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
           ierr=1
           write(*,'(a, a)') 'nao encontrado, ', keyword_name
           return
        endif

        !read(file_lines(keyword_line),1000) (GEOFORM(i),i=1,numel)
        do i=1, numel
           read(file_lines(keyword_line),1000) GEOFORM(i)
           write(*,'(a,a)')   "em leituraRegiaoDs, valor lido =",  GEOFORM(i)
           keyword_line = keyword_line + 1
        end do
        1000 format(a12)
    end subroutine    
    

end module mInputReader

