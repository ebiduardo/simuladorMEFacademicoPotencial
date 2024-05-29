# bidu@lncc.br, abril 2023, proj monan
# para escolha de parametros para o fpm
# uso do arquivo projectOpt.rsp 
# ou por parametro na linha de comando


HYPRE_DIR=/mnt/c/Users/bidu/OneDrive/aLncc/lib/hypre-2.11.2/
HYPRE_LIB=$HYPRE_DIR/src/lib
export    LIBRARY_PATH=$LIBRARY_PATH:$HYPRE_LIB

FPM="../monan/fpm/fpm"
EXE="simuladorP"
EXEF="simulador.exe"
echo Criando o executavel 
comando="echo y | fpm  clean; rm build/cache.*  "
#echo $comando; eval $comando
echo +++ $? +++ 
comando="fpm  build  @projectOpt"
comando="fpm  build  -flag=\"-fPIC -O0 -g -fcheck=bounds -fbacktrace -ffree-line-length-200\""
comando="$FPM  build  --compiler \"mpif90\"\
       	-flag=\"-DwithHYPRE -fPIC -O0 -g -fcheck=bounds -fbacktrace -ffree-line-length-200\" \
       	--link-flag \"-lHYPRE\" -verbose"
echo $comando; eval $comando
if [ "$?" != "0" ]; then
	echo erro! interrompendo o sript FPM.sh
       	exit
fi 

comando="find build -name $EXE |xargs -I{} cp \"{}\" bin/$EXEF"
#EXE="$(find build -name $EXE|head -n 1)"
echo EXE=$EXE
EXE="$(find build -name $EXE|grep app)"
echo EXE=$EXE; # exit
echo copiando $EXE para bin/$EXEF 
comando="cp $EXE bin/$EXEF"
echo $comando; eval $comando
comando="ls -ltr bin/$EXEF"
echo $comando; eval $comando
