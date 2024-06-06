#    opcaoA=${2:-"1"}
#    opcaoB=${3:-"1"}
#    opcaoC=${4:-"1"}
#sufixoTela=${4:-""};

if [ $# == 0 ]; then
  echo "informe o expDIR obrigatoriamente  " 
  echo "ex: $0 exp05x02 " 
  comandoRUN="(export OMP_NUM_THREADS=1; cd exp05x02 ; time  $PWD/bin/simulador.exe ; |tee -a telaSimulador.txt)";
  echo $comandoRUN 
  echo "                 e se quiser inclua outros argumentos opcionais " 
  echo "ex: $0 expDir nomeExecutavel pathExecutavel numThreads marcaNomeArqTela" 
  comandoRUN="(export OMP_NUM_THREADS=numThreads; cd expDir ; time  /pathExecutavel/nomeExecutavel; |tee -a telaMarcaNomeArqTela.txt)";
  echo $comandoRUN 
  exit
fi

expDir=${1}
nomeExecutavel=${2:-"simulador.exe"}
binDir=${3:-"$(pwd)/bin"}
numThreads=${4:-"1"}
EXEC=$binDir/$nomeExecutavel
sufixoExec=${5:-"$(basename $nomeExecutavel |cut -d"." -f1)"}

arqTela="tela_${sufixoExec}_${numThreads}Thr${sufixoTela}";

comandoRUN="(export OMP_NUM_THREADS=$numThreads; cd $expDir ; time ${EXEC}   |tee -a ${arqTela}.txt)";
echo $comandoRUN 
#read -p "+++ digite ENTER para executar o simulador "  

#LIBPARDISO=""
#LIBPARDISO=MKL; sufixoExec="${sufixoExec}${LIBPARDISO}"
#PARDISO_DIR="/usr/local/lib/pardiso5.0"
#HYPRE_DIR="/usr/local/hypre-2.10.1/"
#export PARDISO_LIC_PATH=$PARDISO_DIR
#export  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PARDISO_DIR
#PARDISOLNK="-L${PARDISO_DIR} -lpardiso500-MPI-INTEL1301-X86-64 -lblas -llapack -fopenmp -lpthread -lm" 
#PARDISOLNK="-mkl"; 
#comando="LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${PARDISO_DIR}:${HYPRE_DIR}";echo $comando; eval $comando

echo +++ |tee $expDir/${arqTela}.txt

echo +++ |tee -a $expDir/${arqTela}.txt 
date     |tee -a $expDir/${arqTela}.txt 
echo +++ |tee -a $expDir/${arqTela}.txt 
echo +++  RODAR EXPERIMENTO COM AS SEGUINTES ESCOLHAS: |tee -a $expDir/${arqTela}.txt
echo +++ |tee -a $expDir/${arqTela}.txt 
echo +++  maquina:    $maquina                  |tee -a $expDir/${arqTela}.txt
echo +++  compilador: $FC                       |tee -a $expDir/${arqTela}.txt
echo +++  solver:     $solver                   |tee -a $expDir/${arqTela}.txt
echo +++  executavel: ${binDir}/$nomeExecutavel |tee -a $expDir/${arqTela}.txt
echo +++  diretorio:  ${expDir}                 |tee -a $expDir/${arqTela}.txt

echo +++ |tee -a $expDir/${arqTela}.txt
echo +++ |tee -a $expDir/${arqTela}.txt
echo +++ |tee -a $expDir/${arqTela}.txt

echo $comandoRUN |tee -a $expDir/${arqTela}.txt
formato="\"\t%E real,\t%U user,\t%S sys\" "
#echo $formato
#time -f $formato eval $comando 2> tempoMedido.txt
dataInicio=$(date)
  start_time="$(date -u +%s)"
    eval $comandoRUN 2> tempoMedido.txt
  end_time="$(date -u +%s)"
  elapsed=$((end_time-start_time))
dataFinal=$(date)
echo " +++ inicio da simulacao: $dataInicio" |tee -a $expDir/${arqTela}.txt
echo " +++ fim .. da simulacao: $dataFinal " |tee -a $expDir/${arqTela}.txt
echo "Total of $elapsed seconds elapsed"
cat tempoMedido.txt                          |tee -a $expDir/${arqTela}.txt; 
rm tempoMedido.txt

