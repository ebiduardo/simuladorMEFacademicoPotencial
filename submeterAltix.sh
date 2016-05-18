#$ -S /bin/bash
#  $ -q prjmurad 
#$ -q fila 
# $ -q k20 
# $ -pe  threads 12 
# $ -pe  mpit6  2 
#$ -cwd
#$ -o saida.txt
#$ -l h_vmem=24G

 
# Rsrv675x027Caps6kx6k/

#Declare array with 4 elements
numThreads=(  8 4 2 )
numThreads=(  6  3  1 )
numThreads=( '12'  6  3  1 )
numThreads=(  '1'  )
# get number of elements in the array
ELEMENTS=${#numThreads[@]}

#for (( i=0;i<$ELEMENTS;i++)); do
#
#  declare -x OMP_NUM_THREADS=${numThreads[${i}]}
#
#  if [ ${OMP_NUM_THREADS}  -lt 10 ]; then
#     declare -x arqTela="Tela_0${OMP_NUM_THREADS}THRK20A_Hypre.txt"
#  else
#     declare -x arqTela="Tela_${OMP_NUM_THREADS}THRK20.txt"
#  fi 
#
  declare -x arqTela="Tela_filaAltix-xe.txt"
  expDir=/prj/prjedlg/bidu/Experimentos/ExperimentosImCoopLncc/2D/exp0500x0500/
  expDir=/prj/prjedlg/bidu/Experimentos/ExperimentosImCoopLncc/2D/exp1000x1000/
  comando="./rodarExperimento.sh 2 1 2 ${expDir} >> ${arqTela}" 
#
  echo 'data     :' $(date)                     >  ${arqTela}
  echo 'local    :' $(hostname)                 >> ${arqTela}
  echo 'diretorio:' $(pwd)                      >> ${arqTela}
  echo 'comando  :' ${comando}                  >> ${arqTela}

  echo em submeterAltix.sh: ${comando}

  comando="./rodarExperimento.sh 1 1 2 ${expDir}" 
  #eval ${comando}
  comando="./rodarExperimento.sh 1 2 2 ${expDir}" 
  #eval ${comando}
  comando="./rodarExperimento.sh 2 1 2 ${expDir}" 
  eval ${comando}
  comando="./rodarExperimento.sh 2 2 2 ${expDir}" 
  eval ${comando}
  comando="./rodarExperimento.sh 3 1 2 ${expDir}" 
  eval ${comando}
  comando="./rodarExperimento.sh 3 2 2 ${expDir}" 
  eval ${comando}
  comando="./rodarExperimento.sh 4 1 2 ${expDir}" 
  eval ${comando}
  comando="./rodarExperimento.sh 4 2 2 ${expDir}" 
  eval ${comando}

#done 
