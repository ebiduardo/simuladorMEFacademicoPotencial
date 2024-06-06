#rm bin/simuladorTransiente.exe
rm telaVarios.txt 

nomeExec="simuladorTransiente.exe"
function compilarB(){
comando="./compilar.sh 3 2 1 src/;  "
  echo $comando |tee -a telaVarios.txt
  eval $comando |tee -a telaVarios.txt
  echo $comando; eval $comando 

comando="mpiifort -c -w -module include -DwithPardiso -Dwithcrs -DwithHYPRE -fopenmp -DwithOMP src/driver.F90 -o bin/driver.o "
  echo $comando |tee -a telaVarios.txt
  eval $comando |tee -a telaVarios.txt

comando="mpiifort -c -w -module include -DwithPardiso -Dwithcrs             -fopenmp -DwithOMP src/solverPardisoCSR.F90 -o bin/solverPardisoCSR.o "
  echo $comando |tee -a telaVarios.txt
  eval $comando |tee -a telaVarios.txt

comando="mpiifort -fopenmp -g -O0 -o bin/simuladorHI0.exe bin/variaveisGlobais.o bin/malha.o bin/estruturasDadosSistEq.o bin/leituraEscrita.o bin/funcoesDeForma.o bin/utilitarios.o bin/mInputReader.o bin/solverGaussSkyline.o bin/solverHypre.o bin/solverPardisoCSR.o bin/utilSistemaEquacoes.o bin/potencial.o bin/fluxo.o bin/driver.o -L/mnt/c/Users/bidu/OneDrive/aLncc/lib/hypre-2.11.2_intel/src/lib -lHYPRE -qmkl "
  echo $comando |tee -a telaVarios.txt
  eval $comando |tee -a telaVarios.txt

exe="$(ls -tr bin/simulador* |tail -n1)"; 
echo $exe
comando="mv $exe bin/$nomeExec"
  echo $comando |tee -a telaVarios.txt
  eval $comando |tee -a telaVarios.txt

}
compilarB

echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
comando="source configHypreIntel.sh"
echo $comando |tee -a telaVarios.txt
#eval $comando |tee -a telaVarios.txt
echo LD_LIBRARY_PATH=$LD_LIBRARY_PATH
comando="ldd ./bin/${nomeExec} |tee telaVarios.txt"
echo $comando |tee -a telaVarios.txt
eval $comando |tee -a telaVarios.txt

#exp06x02F_Dirichlet_Gauss - malha 6x2, CC dirihclet, solver direto, eliminação de gauss, armazenamento skyline, F ull (SEM eliminacao equacoes nos Graus de liberdade prescritos)
#exp06x02E_Dirichlet_Gauss - malha 6x2, CC dirihclet, solver direto, eliminação de gauss, armazenamento skyline, E     (COM eliminacao equacoes nos Graus de liberdade prescritos)
#exp06x02E_Dirichlet_HYPRE - malha 6x2, CC dirihclet, solver iterativo, Gradiente conjugado, F ull (SEM eliminacao equacoes) 
#exp06x02F_Dirichlet_HYPRE - malha 6x2, CC dirihclet, solver iterativo, Gradiente conjugado, E     (COM eliminacao equacoes)
#
#exp06x02E_Fonte_Gauss - malha 6x2, força aplicada, solver iterativo, Gradiente conjugado, F ull (SEM eliminacao equacoes)
#exp06x02F_Fonte_Gauss - malha 6x2, força aplicada, solver iterativo, Gradiente conjugado, F ull (SEM eliminacao equacoes)
#exp06x02E_Fonte_HYPRE - malha 6x2, força aplicada, solver iterativo, Gradiente conjugado, F ull (SEM eliminacao equacoes)
#exp06x02F_Fonte_HYPRE - malha 6x2, força aplicada, solver iterativo, Gradiente conjugado, E     (COM eliminacao equacoes)

ls ./exp06x02[EF]_[DF]*_[GH]* -d

listaDirG=(exp06x02F_Dirichlet_Gauss) 
listaDirG=(\
exp06x02F_Dirichlet_Gauss exp06x02E_Dirichlet_Gauss \
exp06x02F_Fonte_Gauss     exp06x02E_Fonte_Gauss \
)

listaDirH=(exp06x02F_Fonte_HYPRE)
listaDirH=(\
exp06x02F_Dirichlet_HYPRE exp06x02E_Dirichlet_HYPRE \
exp06x02F_Fonte_HYPRE     exp06x02E_Fonte_HYPRE \
)

listaDirM=(exp06x02F_Dirichlet_MKL)
listaDirM=(\
exp06x02F_Dirichlet_MKL exp06x02E_Dirichlet_MKL \
)


newList=${listaDir[@]:4:4} 
newList=${listaDir[@]:0:4} 
newList=${listaDir[@]:8:2} 
#newList=(exp06x02F_Dirichlet_Gauss exp06x02F_Dirichlet_HYPRE)
#newList=(exp06x02E_Dirichlet_Gauss exp06x02E_Dirichlet_HYPRE)

newList=${listaDirG[@]} 
newList=${listaDirH[@]} 

listaOptSolver=(Gauss HYPRE MKL) 


echo --- ${listaOptSolver[@]}
#echo --- ${newList[@]}

for s in ${listaOptSolver[@]} 
   do
   if [ "$s" = "Gauss" ]; then newList=${listaDirG[@]} ; fi 
   if [ "$s" = "HYPRE" ]; then newList=${listaDirH[@]} ; fi 
   if [ "$s" = "MKL"   ]; then newList=${listaDirM[@]} ; fi 
   echo ${newList[@]}; 

   for d in  ${newList[@]} 
      do
      echo -----TRANSIENTE-------------- $d  >  telaT.txt
      sed -i -e '/estacionario/{n;d}'  $d/input.dat; sed -i -e '/estacionario/a .FALSE.' $d/input.dat
      comando="(./rodarExperimento.sh  $d ${nomeExec} )|grep \":\\|lador\\|pleto,\\|em sol\\|Euc\\|neq\" -a -A 3|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\\|Espar\\|subr\\|coef\\|criar\\|Info\\|Beg\"   >> telaT.txt "
      echo $comando; eval $comando 
      comando="head -44  telaT.txt |tee -a telaVarios.txt "
      echo $comando; eval $comando 
      comando="tail -16   telaT.txt |tee -a telaVarios.txt "
      echo $comando; eval $comando 
#  read
      echo -----ESTACIONARIO------------ $d |tee  telaE.txt
      sed -i -e '/estacionario/{n;d}'  $d/input.dat; sed -i -e '/estacionario/a .TRUE.' $d/input.dat
      comando="(cd $d; pwd; ../bin/${nomeExec}; pwd)|grep \"lador\\|pleto,\\|em sol\\|Euc\\|neq\" -a -A 1|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\"  >> telaE.txt "
      comando="(./rodarExperimento.sh  $d ${nomeExec} )|grep \":\\|lador\\|pleto,\\|em sol\\|Euc\\|neq\" -a -A 3|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\"  >> telaE.txt "
      echo $comando; eval $comando 
      comando="tail -24  telaE.txt |tee -a telaVarios.txt "
      echo $comando; eval $comando 
      #  read
      done # for d in  ${newList[@]} 
done # for s in  ${listaOptSolver[@]} 
#rm tela.txt

