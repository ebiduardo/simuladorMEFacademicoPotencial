Bp = load    ('vetBSistEqAlgEsparsoP.mtx');
Ap = mmreadE ("matASistEqAlgEsparsoP.mtx");
up = Ap\Bp

permut = symrcm (Ap);
ApP   = Ap(permut,permut);
BpP    = Bp(permut);
upP=ApP\BpP
[nr, npBretorno]=size(permut)
for i=1:npBretorno
   pUpRetorno(permut(i)) = i; 
endfor
upOrigOrdering=upP(pUpRetorno)
up-upOrigOrdering

