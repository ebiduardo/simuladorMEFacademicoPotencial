#$ -S /bin/bash
# $ -q prjmurad 
# $ -q fila 
#$ -q k20 
# $ -pe  threads 11 
# $ -pe  mpit6  1 
#$ -cwd
# $ -o saida0250x02503D_H.txt
#$ -l h_vmem=32G    ###onde 'x' é a quantidade de memória necessária para cada processo.
numThreads=4
comando="(export OMP_NUM_THREADS=${numThreads}; cd $expDir; /usr/bin/time $simulador  &> $arqTela) "

#binDir=/media/sf_LNCC/simuladorPotencial/bin/
binDir=/prj/prjedlg/bidu/simuladorPotencial/bin/
binDir=/prj/prjedlg/bidu/simuladorReservatorioGeomecB/bin

simulador=simuladorPotHYPRETOLMAIOR;   arqTela00=telaADeformedElemHYPRE21_
simulador=simuladorPotPardiso; arqTela00=telaADeformedElemPardiso_
simulador=simuladorHYPRE.exe;   arqTela00=telaHYPRE
#./compilar.sh 2 2 2
#source ./confCompilador.sh 2 2 2
source /hpc/modulos/bash/intel-cluster_studio_xe_2013.sh

simuladorH=simuladorHI4.exe;   arqTela00H=telaHI4_A
simuladorP=simuladorPI4.exe;   arqTela00P=telaPI4_A

expDir=../imCoopLncc/3d/exp0250x0250x0016/
expDir=../imCoopLncc/3d/exp0250x0250x0001/
expDirBase=/prj/prjedlg/bidu/simuladorReservatorioGeomecB/malhaGrossaTempo30dias
expDirBase=/prj/prjedlg/bidu/Experimentos/ExperimentosImCoopLncc/2D/


#listaExpDir=($expDirBase/exp0250x0250x0001/ $expDirBase/exp0250x0250x0002/ $expDirBase/exp0250x0250x0004/ $expDirBase/exp0250x0250x0008/ $expDirBase/exp0250x0250x0016/ $expDirBase/exp0250x0250x0032/)

#listaExpDir=($expDirBase/exp0500x0500x0001/ $expDirBase/exp0500x0500x0002/ $expDirBase/exp0500x0500x0004/ $expDirBase/exp0500x0500x0008/ $expDirBase/exp0500x0500x0016) # $expDirBase/exp0500x0500x0032 )
#listaExpDir=( exp100x21x20_young09_2way exp100x21x20_young09_1way)
listaEspessuras=("4.0e+02" "4.0e+01" "4.0e+00"           "4.0e-01" "4.0e-02")
listaEspessuras=("4.0e+02" "4.0e+01" "4.0e+00" "2.0e+00" "4.0e-01" "4.0e-02")
listaEspessuras=("2.0e-02" "4.0e-02" "8.0e-02" "1.6e-01" "3.2e-01")
listaEspessuras=("2.0e+00" "4.0e+00" "8.0e+00" "1.6e+01" "3.2e+01")
listaEspessuras=("2.0e-02" "2.0e-02" "2.0e-02" "2.0e-02" "2.0e-02")

listaExpDir=$(ls $expDirBase/exp0*)
listaExpDir=$(find $expDirBase -name "exp*000" |sort)

binDir=/prj/prjedlg/bidu/BTsimuladorMEFacademico/bin

#echo ${listaExpDir[*]}; exit
#for novaEspessura in  $( echo ${listaEspessuras[*]} );
#for j in  $( seq 0 0 ); #$( seq 6 -1 0 );
for expDir in  $( echo ${listaExpDir[*]} );
 do
 #expDir=${expDirBase}/${listaExpDir[j]} ;
 #expDir=${expDirBase}/$expDir;
   #for i in $(seq 0 4); #  4 -1 0);
   for i in $(seq 0 0); #  4 -1 0);
    do
     (
#      novaEspessura=${listaEspessuras[j]}; 
      arqTela="${arqTela00}.txt"
      rm -r $arqTela
      #cd $expDir;
#      sed "s/espessura/$novaEspessura/g" inputB00.dat > input.dat;
      hostname             |tee  $arqTela ;
      date                 |tee -a  $arqTela ;
      pwd                  |tee -a  $arqTela ;
      echo  ... $expDir    |tee -a $arqTela ;
#      echo "espessura = ",  $novaEspessura ", copia da tela em: " $arqTela |tee -a $arqTela ;
#      head -n 15 input.in |tee -a $arqTela ;
      export OMP_NUM_THREADS=${numThreads}
      comando="( export OMP_NUM_THREADS=${numThreads}; cd $expDir; time  ${binDir}/$simulador |tee -a $arqTela)"
      comando="(time  ${binDir}/$simulador |tee -a $arqTela)"

      comandoH="(./rodarExperimento.sh 3 2 2  ${expDir})"
      comandoP="(./rodarExperimento.sh 2 2 2  ${expDir})"

      echo $comandoH        |tee -a $arqTela;
      eval $comandoH
      echo $comandoP        |tee -a $arqTela;
      eval $comandoP

   #   ls -ltr ${arqTela00}*
    )
    done
 done
#exit
