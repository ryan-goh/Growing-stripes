function [F, Vphi, Vy,Jinh] = shwf(y0)
global dx dy dxu Lapu Dxu xu x y lu xuext BDx BDy N Nx Ny n Xi Yi tang z1
global mu0 mu BDx chi a III Ga tau Mx znew ky



w = y0(1:N); kx = y0(N+1); c = y0(N+2);  ;
ur = getroll(sqrt(kx^2+ky^2));
%urext = [ur; ur; ur];
urext = repmat(ur,Mx+1,1);
Upol = griddedInterpolant(xuext, urext,'spline');
Uphi = Upol(kx*Xi-Yi);
Vphi = reshape(Uphi,Nx*Ny,1);
Vy = BDy*Vphi;
%Jinh = 0*speye(N+2,N+2);
Jinh =  Ga;

UU = (w+chi*Vphi);
%Fw = (Lw(kx,ky)+Ga)*UU+ tau*(UU).^2 - UU.^3
%-chi*((Lw(kx,ky)+Ga)*Vphi + tau*(Vphi).^2 - Vphi.^3);
LL = Lw(kx,c)+Ga;
Fw = LL*UU - UU.^3 - chi*(LL*Vphi - Vphi.^3);
Fpc = dx*dy*Vy(III)'*w(III);%dx*dy*Vy'*w;%dx*dy*Vy(III)'*w(III);
Fk = (y0-z1)'*tang;

F = [Fw; Fpc; Fk];
Vphi=[Vphi;0;0];
Vy=[Vy;0;0];


end
