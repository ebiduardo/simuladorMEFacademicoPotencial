arquivoOriginal=coefSistEqAlgEsparsoP.mtx
arquivoDestino=matASistEqAlgEsparsoP.mtx

cp $arquivoOriginal $arquivoDestino
sed -i 's/%%/% %/' $arquivoDestino 
sed -i '/% %/i %%MatrixMarket matrix coordinate real symmetric' $arquivoDestino 
sed -i '/BRHS/,/$!/d'  $arquivoDestino 

arquivoOriginal=coefSistEqAlgEsparsoF.mtx
arquivoDestino=matASistEqAlgEsparsoF.mtx

cp $arquivoOriginal $arquivoDestino
sed -i 's/%%/% %/' $arquivoDestino 
sed -i '/% %/i %%MatrixMarket matrix coordinate real symmetric' $arquivoDestino 
sed -i '/BRHS/,/$!/d'  $arquivoDestino 

octave  < criarImagemMatrizA.m
