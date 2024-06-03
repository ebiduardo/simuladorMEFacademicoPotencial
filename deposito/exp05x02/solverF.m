 Bf = load    ('vetBSistEqAlgEsparsoF.mtx');
Af = mmreadE ("matASistEqAlgEsparsoF.mtx");
uf = Af\Bf
# disp(u)


Bf = load    ('vetBSistEqAlgEsparsoF.mtx');
Af = mmreadE ("matASistEqAlgEsparsoF.mtx");
permut = symrcm (Af);

AfP   = Af(permut,permut);
BfP    = Bf(permut);
ufP=AfP\BfP
#for i=1:12 pBretorno(permut(i)) = i; endfor
[nr, npBretorno]=size(permut)
for i=1:npBretorno pBfRetorno(permut(i)) = i; endfor
uOrig=ufP(pBfRetorno)
uf-uOrig

