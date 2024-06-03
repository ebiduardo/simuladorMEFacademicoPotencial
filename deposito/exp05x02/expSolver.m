
Ny=3;
Nx=3;
ex = ones(Nx,1);
Axx = spdiags([ex -2*ex ex], -1:1, Nx, Nx);
Ix = speye(Nx);
ey = ones(Ny,1);
Ayy = spdiags([ey -2*ey ey], -1:1, Ny, Ny);
Iy = speye(Ny);
A = kron(Axx,Iy)+kron(Ix,Ayy);
C=symrcm (A); %this is the actual RCM
% just to visualize things 
spy (A)
spy (A(C,C))
b=ones(9,1);
% store the permutated matrix
Amut= (A(C,C));
% solve and compare
x=A\b
xmut=Amut\b
x==xmut % results in avector full of zeros, e.g. False.
