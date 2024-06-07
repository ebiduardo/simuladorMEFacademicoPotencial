dirFontes=${1:-${dirFontesD}}
opcaoA=${2:-"1"}
opcaoB=${3:-"1"}
opcaoC=${4:-"1"}
dirFontesD=fontes
    #dirFontes=${4:-"fontes"}
exeSufixExtra=${5:-""}

     listaSolvers=(zero Gauss Pardiso HYPRE PH)
listaCompiladores=(zero gfortran ifort) # pgf90)
     listaOPTIMIZ=(zero "-g -O0" "-g -O0" "-fast")
     listaOPTIMIZ=(zero "-g -O0" "-O4" "-fast")

 solver=${listaSolvers[opcaoA]}
     FC=${listaCompiladores[opcaoB]}
OPTIMIZ=${listaOPTIMIZ[opcaoC]}
maquina=desconhecida
maquina="bidu-debian"
maquina="DESKTOP-PRUDTOO"
maquina="pc012137"
maquina=$(hostname)

#HYPRE_DIR="/mnt/c/Users/EduardoLMGarcia/lib/hypre-2.11.2"
#HYPRE_DIR="/cygdrive/c/Users/EduardoLMGarcia/"
#HYPRE_DIR="/cygdrive/c/Users/EduardoLMGarcia/lib/hypre-2.11.2/appWindows"
#HYPRE_DIR="/cygdrive/c/Users/EduardoLMGarcia/lib/hypre-2.11.2/cygWin"
#HYPRE_DIR="/mnt/c/Users/$USER/OneDrive/aLncc/lib/hypre-2.11.2/appWindows"
#HYPRE_DIR="/mnt/c/Users/$USER/OneDrive/aLncc/lib/hypre-2.11.2B/src/"
USER=EduardoLMGarcia
USER=$LOGNAME
HYPRE_DIR="/mnt/c/Users/$USER/OneDrive/aLncc/lib/hypre-2.11.2/src/"
HYPRE_DIR="/mnt/c/Users/$USER/OneDrive/aLncc/lib/hypre-2.11.2B/src/"
#LIBS="$LIBS -L/mnt/c/Users/bidu/OneDrive/aLncc/lib/hypre-2.11.2B/src/lib/ -lHYPRE "

     lS=(z G P H PH)
     lC=(z G I)  # P)
     lO=(z 0 4 f)
sufixoExec="${lS[opcaoA]}${lC[opcaoB]}${lO[opcaoC]}"

LIBPARDISO=MKL

echo "escolhas: $opcaoA	 $opcaoB	$opcaoC"
echo "escolhas: $solver	 $FC 		$OPTIMIZ ...  $maquina"
echo "sufixo  : $sufixoExec"
echo "sufixo Extra  : $exeSufixExtra" ...

dirBin=bin

OUTROS="-DmostrarTempos -Ddebug"
OUTROSF="-DmostrarTempos"
unset listaFontes
unset listaObjetos

     listaFontes=(${dirFontes}/variaveisGlobais.F90  \
        ${dirFontes}/malha.F90 \
	${dirFontes}/estruturasDadosSistEq.F90 \
	${dirFontes}/leituraEscrita.F90 \
        ${dirFontes}/funcoesDeForma.F90 \
 	${dirFontes}/utilitarios.F90 \
	${dirFontes}/mInputReader.F90 \
	${dirFontes}/solverGaussSkyline.F90 \
        ${dirFontes}/solverHypre.F90 \
        ${dirFontes}/solverPardisoCSR.F90 \
 	${dirFontes}/utilSistemaEquacoes.F90 \
	${dirFontes}/potencial.F90 ${dirFontes}/fluxo.F90 \
        ${dirFontes}/driver.F90 )

     echo ${listaFontes[*]} 
     i=0
     for f in ${listaFontes[*]} 
       do
           listaObjetos[$i]=${dirBin}/$(basename ${f/F90/o/})
#	   echo "${listaObjetos[$i]}"
           ((i=i+1))
       done

LIBS=""
if [ "$FC" = "ifort" ]; then 
   ARGINC="-w -module include";
   COMP="-fopenmp  -DwithOMP"
   LOMP="-fopenmp"
   INTEL_DIR=/opt/intel
   INTEL_DIR=/opt/intel/oneapi/
   INTEL_MPI=${INTEL_DIR}/oneapi/mpi/
   INTEL_MPI_LIB=$HYPRE_DIR/2021.9.0/lib/release # libmpi.so
   INTEL_MPI_LIB=$HYPRE_DIR/2021.10.0/lib/release # libmpi.so
   INTEL_MPI_INC="-I $INTEL_MPI/2021.9.0/include"  
   INTEL_MPI_INC="-I $INTEL_MPI/2021.10.0/include"  
   LIBS="$LIBS  -qmkl"; 
   INCS="$INCS $INTEL_MPI_INC"
fi
   comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH::${INTEL_MPI}"; eval $comando

if [ "$FC" = "gfortran" ]; then
   ARGINC="-I include";
   ARGINC="-M include";
   ARGINC="-fpic -ffree-line-length-none -J include"; # GNU Fortran (Debian 4.7.2-5) 4.7.2
   OUTROSF="-DmostrarTempos  -fbounds-check   "
   COMP="-fopenmp  -DwithOMP"
   LOMP="-fopenmp"
   INCS="$INCS $ARGINC"
fi

# /hpc/pardiso5.0
#Architecture X86-64, 64-bit, icc/ifort 13.01 libpardiso500-INTEL1301-X86-64.so
#Architecture X86-64, 64-bit, gcc/gfortran 4.7.2 libpardiso500-GNU472-X86-64.so

#Architecture X86-64, 64-bit, icc/ifort 13.01     Linux MPI libpardiso500-MPI-INTEL1301-X86-64.so
#Architecture X86-64, 64-bit, mpicc/mpif77 4.7.2     Linux MPI libpardiso500-MPI-GNU472-X86-64.so

   LIBS="" ;  

case $solver in
"Gauss" )
   ppSolver="-DSkyline";
;;

"Pardiso")
  if [ "$FC" = "gfortran" ]; then
    FC=gfortran
    HYPRE_DIR="/mnt/c/Users/bidu/OneDrive/aLncc/lib/hypre-2.11.2/src"
    HYPRE_LIB="$HYPRE_DIR/lib/appWindows/"
    HYPRE_INC="${HYPRE_DIR}/include/"
    comando="export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/hpc/pardiso5.0/"; eval $comando
    PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-GNU472-X86-64 -lblas -llapack -fopenmp -lpthread -lm" 
    PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-GNU472-X86-64 -L${LAPACK_DIR}     -lblas -llapack -fopenmp -lpthread -lm" 
  fi
  if [ "$FC" = "ifort" ]; then
    FC=ifort
    if [ "$LIBPARDISO" = "MKL" ]; then
       PARDISOLNK="-qmkl"; 
    else
       PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-INTEL1301-X86-64 -lblas -llapack -fopenmp -lpthread -lm" 
       PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-INTEL1301-X86-64  -L${LAPACK_DIR}    -lblas -llapack -fopenmp -lpthread -lm" 
    fi
  fi
   comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PARDISO_DIR}"; eval $comando
   ppSolver="-DwithPardiso -Dwithcrs -DSkyline  $COMP ";  
   ppSolver="-DwithPardiso -Dwithcrs            $COMP ";  
   LIBS="${PARDISOLNK}"; 
;;

"HYPRE")
  ppSolver="-DwithHYPRE -DwithPardiso -Dwithcrs -DSkyline $COMP ";
  ppSolver="-DwithHYPRE -Dwithcrs $COMP"; 
  if [ "$FC" = "gfortran" ]; then
    FC=mpif90
    HYPRE_DIR="/mnt/c/Users/bidu/OneDrive/aLncc/lib/hypre-2.11.2/src"
    HYPRE_LIB="${HYPRE_DIR}/lib/appWindows" 
    HYPRE_LNK="-L${HYPRE_LIB} -lHYPRE " 
    HYPRE_INC="${HYPRE_DIR}/include/"
  fi
  if [ "$FC" = "ifort" ]; then
    FC=mpiifort
    HYPRE_DIR="/mnt/c/Users/bidu/OneDrive/aLncc/lib/hypre-2.11.2_intel"
    HYPRE_LIB="$HYPRE_DIR/src/lib"
    HYPRE_LNK="-L${HYPRE_LIB} -lHYPRE " 
    HYPRE_INC="${HYPRE_DIR}/src/include/"
    INC="-I /opt/intel/oneapi/mpi/2021.6.0/include"
    INC="-I /opt/intel/oneapi/mpi/2021.10.0/include"
    INC="-I $INTEL_MPI/include"  
    INC="-I /opt/intel/oneapi/mpi/2021.9.0/include"  
    INC="-I /opt/intel/oneapi/mpi/2021.10.0/include"  
    FC="mpiifort $INC"
  fi
  LIBS="$LIBS  ${HYPRE_LNK} "; 
  comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH::${HYPRE_LIB}"; eval $comando
  echo ::: $LIBS
;;

"PH")
  ppSolver="-DwithHYPRE -DwithPardiso -Dwithcrs -DSkyline $COMP ";
  ppSolver="-DwithHYPRE -DwithPardiso -Dwithcrs           $COMP ";
  HYPRELNK="-L${HYPRE_DIR}/lib -lHYPRE " 
  HYPRE_INC="${HYPRE_DIR}/include/"
  if [ "$FC" = "gfortran" ]; then
    FC=mpif90
    PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-GNU472-X86-64 -L${LAPACK_DIR} -lblas -llapack -fopenmp -lpthread -lm"
    LIBS="${HYPRE_LNK} ${PARDISOLNK} "; 
  fi
  if [ "$FC" = "ifort" ]; then
    FC=mpiifort
    PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-INTEL1301-X86-64     -lblas -llapack -fopenmp -lpthread -lm"
    LIBS="${HYPRE_LNK} -qmkl"; 
  fi

  if [ "$LIBPARDISO" = "MKL" ]; then
       PARDISOLNK="-mkl";
  else
       PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-INTEL1301-X86-64 -lblas -llapack -fopenmp -lpthread -lm"
  fi
   comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PARDISO_DIR}:${HYPRE_DIR}"; eval $comando
;;
*)
  echo "...erro na escolha do solver....em compilar.sh" 
  exit;;
esac

#echo +++  maquina FC solver ppSolver ARGINC COMP LOMP LIBS
comando="($FC --version)"; echo $comando; eval $comando

if [ $# == 0 ]; then
nomeExecutavel=simulador${exeSufixExtra}.exe
else 
nomeExecutavel=simulador${sufixoExec}${exeSufixExtra}.exe
fi
#echo "+++ digite Enter para gerar o executavel: $nomeExecutavel "
#read
comando="rm -rf $dirBin/$nomeExecutavel"
echo $comando
eval $comando



for i in $(seq 1 ${#listaFontes[*]})
 do
   FFLAGS="${OPTIMIZ}  ${OUTROSF}"
   FFLAGS="${OPTIMIZ}  ${OUTROSF} ${COMP}"
   comando="${FC} -c ${ARGINC} ${FFLAGS} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}  " ; 
   comando="${FC} -c ${ARGINC} ${ppSolver} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}" ; 

if [ "${listaFontes[i-1]}" = "${dirFontes}/solverHypre.F90" ]; then
   comando="${FC} -c ${ARGINC} ${ppSolver} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}" ; 
fi

if [ ".${listaFontes[i-1]}." = ".${dirFontes}/solverPardisoCSR.F90." ]; then
   comando="${FC} -c ${ARGINC} ${ppSolver} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}" ; 
fi

if [ "${listaFontes[i-1]}" = "${dirFontes}/driver.F90" ]; then
   comando="${FC} -c ${ARGINC} ${ppSolver}  ${FFLAGS} ${listaFontes[i-1]} -o ${listaObjetos[i-1]}" ; 
fi
   echo $comando;  eval $comando
done

LFLAGS="${LOMP} ${OPTIMIZ}"
comando="${FC} ${LFLAGS} -o ${dirBin}/${nomeExecutavel} ${listaObjetos[*]} ${LIBS} ";
echo $comando; eval $comando
#mpif90 -fopenmp -g -O0 -o bin/prototipoHYPRE.exe bin/*.o -L"/mnt/c/Users/bidu/OneDrive/aLncc/lib/hypre-2.11.2B/src/lib/" -lHYPRE 


echo +++
echo +++ executavel criado
ls -ltr ${dirBin}/${nomeExecutavel}
echo +++

#rm $dirBin/*.o include/*
#echo +++ aguardando um sinal
expDir="exp05x023D"
expDir="exp07"
expDir="exp05x02"
numThreads="1"
comando="./rodarExperimento.sh $expDir $numThreads $opcaoA $opcaoB $opcaoC "
#echo $comando; 
#read
echo +++
#eval $comando
 
