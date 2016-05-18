#script para geracao de arquivo de entrada 2D com numeracao com maior largura de banda
#!/bin/bash

nsd=2 # numero de dimensoes do modelo: 1D, 2D ou 3D -> 1,  2 ou 3 respectivamente
printCoord=0 # 1 nao 0 sim

nelDir1p=8  # numero na direcao X padrao
nelDir1=${1:-${nelDir1p}} # atribuir valor padrao aa variavel nelDir1

nelDir2p=2  # numero na direcao Y padrao
export nelDir2=${2:-${nelDir2p}} # atribuir valor padrao aa variavel nelDir2

#echo "  gerando ${arqSaida} de entrada com malha: ${nelDir1} X ${nelDir2}"
#echo "  a partir do arquivo:  ${arqEntrada}"
#echo " Confirma?"
#read 

####### definicao de parametros (numero de nos e numero de elementos
numnp=$(( (nelDir1+1) * (nelDir2+1) )) # numero total de nos
numEl=$(( nelDir1 * nelDir2)) # numero total de elementos  
nen=4
npint=4

domI_X=0.0
domF_X=1.0
domI_Y=0.0
domF_Y=1.0


printf "%5s \n" 0 
printf "  experimento inicial com programa modular nodal+ladal, CC dirichlet \n"
printf "%10s %9s %9s\n" 1 ${printCoord} ${nsd}
printf "%10s %9s %9s %9s \n" ${numnp} ${numEl} ${nen} ${npint}
printf "%10s %9s \n" 1 0 

#Coordenadas
printf "%10s %9s %9s %9s\n" 1 4   ${domI_X} ${domI_Y} 
printf "%10s %9s %9s %9s\n" "" "" ${domF_X} ${domI_Y} 
printf "%10s %9s %9s %9s\n" "" "" ${domF_X} ${domF_Y} 
printf "%10s %9s %9s %9s\n" "" "" ${domI_X} ${domF_Y} 
printf "%10s %9s %9s %9s\n" ${nelDir1} "1" ${nelDir2} $(( nelDir1+1 ))
printf "%10s \n" 0 


#Conectividades
n1=1
n2=$(( n1 + 1 ))
n3=$(( nelDir1+1+1+1    ))
n4=$(( n3-1    ))
printf "%10s %9s %9s %9s %9s %9s %9s\n" "1" "1" ${n1} ${n2} ${n3} ${n4} "1"

incElemX=1
incNoX=1
incElemY=${nelDir1}
incNoY=$(( nelDir1+1 ))
printf "%10s %9s %9s %9s %9s %9s\n" ${nelDir1} ${incElemX} ${incNoX} ${nelDir2} ${incElemY} ${incNoY}
printf "%10s \n" 0 


#Condicoes de Contorno (ID)
printf "%10s %9s %9s %9s\n" "1" $(( numnp-nelDir1 ))  $(( nelDir1+1 ))  "1"
printf "%10s %9s %9s %9s\n"  $(( nelDir1+1 )) ${numnp}  $(( nelDir1+1 )) "1" 
printf "%10s \n" 0 


#Condicoes de Contorno (F)
valI=" 1.0d-0"
valF="-1.0d-0"
#esquerda
printf "%10s %9s %9s\n" "1" "2"  ${valI}  
printf "%10s %9s %9s\n" "" ""    ${valI}
printf "%10s %9s \n" ${nelDir2}  $(( nelDir1+1 ))    
#direita
printf "%10s %9s %9s\n" $(( nelDir1+1 )) "2"  ${valF}  
printf "%10s %9s %9s\n" "" ""    ${valF}
printf "%10s %9s \n" ${nelDir2}  $(( nelDir1+1 ))    
printf "%10s \n" 0 


#parametros diversos
printf "%10s \n" 0 
printf "%10s \n" 1
alfa="1.00e+00"
beta="1.00e+00"
gama="1.00e+00"
printf "%10s %14s %9s %9s\n"  "1"  ${alfa} ${beta} ${gama} 
printf "%10s %9s %9s\n"  "0.0"  "0.0" "0.0" 


printf "*end \n"
