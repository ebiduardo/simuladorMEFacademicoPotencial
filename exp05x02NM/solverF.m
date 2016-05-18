 Bf = load    ('vetBSistEqAlgEsparsoF.mtx');
Af = mmreadE ("matASistEqAlgEsparsoF.mtx");
uf = Af\Bf

permut = symrcm (Af);
AfP   = Af(permut,permut);
BfP   = Bf(permut);
ufP=AfP\BfP
[nr, npBretorno]=size(permut)
for i=1:npBretorno
  pUfRetorno(permut(i)) = i;
endfor
ufOrigOrdering=ufP(pUfRetorno)
uf-ufOrigOrdering

