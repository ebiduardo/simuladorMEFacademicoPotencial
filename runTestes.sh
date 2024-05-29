simDir=/mnt/c/Users/bidu/OneDrive/aLncc/simuladorMEFacademicoPotencial2023Fev
simuladorN=$simDir/bin/simuladorI.exe

expDir=../prototipos/matrizesPotencial2D/exp01000x01000/

listaExpDir=(exp06x02F_Fonte exp06x02F_Dirichlet exp06x02E_Fonte exp06x02E_Dirichlet) 

listaS=(                            GaussSkyline)
listaS=(               HYPREEsparso GaussSkyline)
listaS=(PardisoEsparso HYPREEsparso GaussSkyline)


echo --- ${listaS[@]}
echo --- ${listaExpDir[@]}
for expDir in  ${listaExpDir[@]}
  do 
  echo --- ${expDir}
  for solverN in  ${listaS[@]}
    do
    echo - - -  ${expDir} +  ${solverN}
    #read # continue
    lin1=$(grep solver\_pressao $expDir/inputDS.dat -n|cut -d ":" -f 1) 
    #echo lin1=$lin1;
    comando="awk -F\":\" 'BEGIN{print \"sed \\\" \", $lin1+1, \"s/.*/$solverN/\\\" $expDir/inputDS.dat -i;\"}'"
    #echo :::  $comando
    comandoB=$(eval $comando)
    #echo :::  $comandoB
    bash -c "$comandoB"
    lin2=$(grep solver_fluxo $expDir/inputDS.dat -n|cut -d ":" -f 1) 
    comando="awk -F\":\" 'BEGIN{print \"sed \\\" \", $lin2+1, \"s/.*/$solverN/\\\" $expDir/inputDS.dat -i;\"}'"
    #echo :::  $comando
    comandoB=$(eval $comando)
    #echo :::  $comandoB
    bash -c "$comandoB"
    grep solver $expDir/inputDS.dat -A 1
    comando=" (cd $expDir/; $simuladorN)"
    echo $comando
    bash -c "$comando"

  done
done

#exit
