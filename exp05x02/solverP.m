 Bp = load    ('vetBSistEqAlgEsparsoP.mtx');
Ap = mmreadE ("matASistEqAlgEsparsoP.mtx");
u = Ap\Bp
# disp(u)


Bp = load    ('vetBSistEqAlgEsparsoP.mtx');
Ap = mmreadE ("matASistEqAlgEsparsoP.mtx");
permut = symrcm (Ap);

ApP   = Ap(permut,permut);
BpP    = Bp(permut);
uP=ApP\BpP
# disp(uP)
for i=1:12 pBretorno(permut(i)) = i; endfor
for i=1:size(permut) pBretorno(permut(i)) = i; endfor
uOrig=uP(pBretorno)
#u-uP(permV)
u-uOrig

