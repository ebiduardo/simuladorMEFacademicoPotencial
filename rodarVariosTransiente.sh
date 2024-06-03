#exp06x02F_Dirichlet_Gauss - malha 6x2, CC dirihclet, solver direto, eliminação de gauss, armazenamento skyline, F ull (SEM eliminacao equacoes)
#exp06x02E_Dirichlet_Gauss - malha 6x2, CC dirihclet, solver direto, eliminação de gauss, armazenamento skyline, E (COM eliminacao equacoes)
#exp06x02E_Dirichlet_HYPRE - malha 6x2, CC dirihclet, solver iterativo, Gradiente conjugado, F ull (SEM eliminacao equacoes) 
#exp06x02F_Dirichlet_HYPRE - malha 6x2, CC dirihclet, solver iterativo, Gradiente conjugado, E     (COM eliminacao equacoes)
#
#exp06x02E_Fonte_Gauss - malha 6x2, força aplicada, solver iterativo, Gradiente conjugado, F ull (SEM eliminacao equacoes)
#exp06x02F_Fonte_Gauss - malha 6x2, força aplicada, solver iterativo, Gradiente conjugado, F ull (SEM eliminacao equacoes)
#exp06x02E_Fonte_HYPRE - malha 6x2, força aplicada, solver iterativo, Gradiente conjugado, F ull (SEM eliminacao equacoes)
#exp06x02F_Fonte_HYPRE - malha 6x2, força aplicada, solver iterativo, Gradiente conjugado, E     (COM eliminacao equacoes)

ls ./exp06x02[EF]_[DF]*_[GH]* -d
listaDir=(\
exp06x02F_Dirichlet_Gauss exp06x02E_Dirichlet_Gauss \
exp06x02E_Dirichlet_HYPRE exp06x02F_Dirichlet_HYPRE \
exp06x02E_Fonte_Gauss     exp06x02F_Fonte_Gauss \
exp06x02E_Fonte_HYPRE     exp06x02F_Fonte_HYPRE \
)


newList=${listaDir[@]:4:4} 
newList=${listaDir[@]:0:4} 
newList=${listaDir[@]:8:2} 
#newList=(exp06x02F_Dirichlet_Gauss exp06x02F_Dirichlet_HYPRE)
#newList=(exp06x02E_Dirichlet_Gauss exp06x02E_Dirichlet_HYPRE)

newList=${listaDir[@]} 

#. compilar.sh 3 1 2 srcTransiente/
EXEC=simuladorTransienteHG4.exe
mv bin/simuladorHG4.exe bin/$EXEC
EXEC=simuladorTransiente.exe

rm tela.txt
echo --- ${newList[@]}
for d in  ${newList[@]} 
  do
  echo ----------------------------- $d
  comando=". rodarExperimento.sh $d "
  comando=". rodarExperimento.sh $d 3 1 2 | grep \"lador\\|pleto,\\|em sol\\|Euc\" -a -A 1|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\" |tail -20"
  comando="(cd $d; pwd; ../bin/${EXEC}; pwd )|grep \"lador\\|pleto,\\|em sol\\|Euc\" -a -A 1|grep -v \"alloc\\|fact\\|passo\\|seg\\|--\" |tail -20"
  echo $comando
  eval $comando |tee -a tela.txt
  echo $?
done

