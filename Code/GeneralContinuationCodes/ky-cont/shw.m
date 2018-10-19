function [F, J] = shw(y0)
global dx dy y x N Nx Ny n lu dxu z1 tang  
global chi mu mu0 k0 a III Ga tau BDx Mx 


[F, Vphi, Vy, Jin]  = shwf(y0);


Vphi = Vphi(1:N);
Vy = Vy(1:N);

w = y0(1:N); kx = y0(N+1); ky = y0(N+2);


%Jww = Lw(kx,ky) - spdiags(-2*tau*(w + chi*Vphi)+3*(w + chi*Vphi).^2,0,N,N)+Jin;
Jww = Lw(kx,ky) - spdiags(3*(w + chi*Vphi).^2,0,N,N)+Jin;
%Jphw = Vy(III)';
Jkall = tang';

J = speye(N+2,N+2);
J(1:N,1:N) = Jww;
%J(N+1,N-1*n+1-10:N-10) = Jphw; 
J(N+1,III) =dx*dy*Vy(III)';%(III)';
%J(N+1,1:N) =dx*dy*Vy(:)';%(III)';
J(N+2,:) = Jkall;
%J(N+1,1:N) = Jphw;
del = 1e-7;del2 = del*2;
dkx = zeros(N+2,1);dkx(N+1,1) = del;
dky = zeros(N+2,1);dky(N+2,1) = del;
J(1:N+2,N+1) = (shwf(y0-2*dkx)/12 - (2/3)*shwf(y0-dkx) + (2/3)*shwf(y0+dkx) - shwf(y0+2*dkx)/12)/del;%
%J(1:N+2,N+1)=(shwf(y0+dkx) - shwf(y0-dkx))/del2;
J(1:N+2,N+2) = (shwf(y0-2*dky)/12 - (2/3)*shwf(y0-dky) + (2/3)*shwf(y0+dky) - shwf(y0+2*dky)/12)/del; %
%J(1:N+2,N+2)=(shwf(y0+dky) - shwf(y0-dky))/del2;
