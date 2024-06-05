rm bin/simuladorEstacionario.exe
rm bin/simuladorTransiente.exe
#./compilar.sh 1 1 1 srcEstacionario/ ; mv bin/simuladorGG0.exe bin/simuladorEstacionario.exe
#./compilar.sh 1 1 1 srcTransiente/   ; mv bin/simuladorGG0.exe bin/simuladorTransiente.exe
./compilar.sh 1 1 1 src/   ; mv bin/simuladorGG0.exe bin/simuladorTransiente.exe
expDir=exp06x02F_Dirichlet_Gauss
echo TRANSIENTE, $expDir
(cd $expDir/; ../bin/simuladorTransiente.exe   ) |head -17 ; read
(cd $expDir/; ../bin/simuladorTransiente.exe   ) |grep -a  ", potencial" -A 6 |head -7 
(cd $expDir/; ../bin/simuladorTransiente.exe   ) |grep -a  ", potencial" -A 6 |tail -7 
echo ESTACIONARIO, $expDir
(cd $expDir/; ../bin/simuladorEstacionario.exe ) |grep -a  ", potencial" -A 6
expDir=exp06x02F_Fonte_Gauss
echo TRANSIENTE, $expDir
(cd $expDir/; ../bin/simuladorTransiente.exe   ) |grep -a  ", potencial" -A 6 |head -7 
(cd $expDir/; ../bin/simuladorTransiente.exe   ) |grep -a  ", potencial" -A 6 |tail -7 
echo ESTACIONARIO, $expDir
(cd $expDir/; ../bin/simuladorEstacionario.exe ) |grep -a  ", potencial" -A 6
expDir=exp06x02FDG
