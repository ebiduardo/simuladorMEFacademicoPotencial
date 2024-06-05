#rm bin/simuladorTransiente.exe

#./compilar.sh 3 2 1 src/; 
#exe="$(ls -tr bin/simulador* |tail -n1)"; echo $exe
#mv $exe bin/simuladorTransiente.exe

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


newList=${listaDir[@]:4:4} 
newList=${listaDir[@]:0:4} 
newList=${listaDir[@]:8:2} 
#newList=(exp06x02F_Dirichlet_Gauss exp06x02F_Dirichlet_HYPRE)
#newList=(exp06x02E_Dirichlet_Gauss exp06x02E_Dirichlet_HYPRE)

newList=${listaDirG[@]} 
newList=${listaDirH[@]} 

#. compilar.sh 3 1 2 srcTransiente/
EXEC=simuladorTransienteHG4.exe
EXEC=simuladorTransiente.exe

rm telaVarios.txt 
echo --- ${newList[@]}
for d in  ${newList[@]} 
  do
  echo -----TRANSIENTE-------------- $d |tee  telaT.txt
  sed -i -e '/estacionario/{n;d}'  $d/input.dat; sed -i -e '/estacionario/a .FALSE.' $d/input.dat
  comando="(cd $d; pwd; ../bin/${EXEC}; pwd)|grep \"lador\\|pleto,\\|em sol\\|Euc\\|neq\" -a -A 1|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\"  >> telaT.txt "
  comando="(./rodarExperimento.sh  $d 1 1 1 ${EXEC})|grep \":\\|lador\\|pleto,\\|em sol\\|Euc\\|neq\" -a -A 1|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\"  >> telaT.txt "
  comando="(./rodarExperimento.sh  $d 1 1 1 ${EXEC})|grep \":\\|lador\\|pleto,\\|em sol\\|Euc\\|neq\" -a -A 3|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\"  >> telaT.txt "
  echo $comando; eval $comando 
  comando="head -27  telaT.txt |tee -a telaVarios.txt "
  echo $comando; eval $comando 
  comando="tail -10   telaT.txt |tee -a telaVarios.txt "
  echo $comando; eval $comando 
  echo -----ESTACIONARIO------------ $d |tee  telaE.txt
  sed -i -e '/estacionario/{n;d}'  $d/input.dat; sed -i -e '/estacionario/a .TRUE.' $d/input.dat
  comando="(cd $d; pwd; ../bin/${EXEC}; pwd)|grep \"lador\\|pleto,\\|em sol\\|Euc\\|neq\" -a -A 1|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\"  >> telaE.txt "
  comando="(./rodarExperimento.sh  $d 1 1 1 ${EXEC})|grep \":\\|lador\\|pleto,\\|em sol\\|Euc\\|neq\" -a -A 3|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\"  >> telaE.txt "
  echo $comando; eval $comando 
  comando="tail -15  telaE.txt |tee -a telaVarios.txt "
  echo $comando; eval $comando 
done
#rm tela.txt

