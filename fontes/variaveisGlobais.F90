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
      module mGlobaisArranjos
!
        integer*4, allocatable :: npar(:)
        real*8               :: etime(6)
        character(len=80)    :: title
        integer*4, allocatable :: mat(:)
        real*8,  allocatable :: grav(:), bf(:,:), c(:,:)
        logical :: listaSolverDisponivel(3)

      end module ! mGlobaisArranjos

!
      module mGlobaisEscalares
!
        integer*4:: exec,iprtin
        integer*4:: numParElem=15
        integer*4:: ntype,numat,nrowsh,nicode,npint

        real*8, parameter  :: zero=0.0d0, pt25=0.25d0, pt5=0.5d0
        real*8, parameter  :: four=4.0d0, five=5.0d0, six=6.0d0
        real*8, parameter  :: pt1667=0.1666666666666667d0
        real*8, parameter  :: one=1.0d0, two=2.0d0, three=3.0d0
        real*8, parameter  :: pt8=0.8d0
        real*8, parameter  :: pt45=0.45d0
        real*8, parameter  :: pi=3.14159265359

        integer*4:: id0
        integer*4:: optSolver

      end module mGlobaisEscalares


