function [UU] = sprof(x,y,kky,znew)
global dxu Lapu xu n mu0
%%%input a core solution and other system parameters, use far-field solver
%%%to build full solution
mu0;
Lx = x(end);
Ly = y(end);
dx = x(2) - x(1);
dy = y(2) - y(1);
Nx = floor(Lx/dx);
Ny = floor(Ly/dy);
N = Nx*Ny;



%%Need this to build full solution from core solution an kx, ky,c
igL = 1;[0:dy:Ly]';
igR = floor(Nx*8/10);
igW = floor(Nx*7/10);
hp = -mu0;
hm = mu0;

l = x(igW);
chi1=heaviside(x-l-4);
chi1=chi1+1/64*(3*(x-l).^5-(x-l).^6).*(1-heaviside(x-l-2)).*heaviside(x-l);
chi1=chi1+(1+1/64*(3*(x-l-4).^5+(x-l-4).^6)).*(1-heaviside(x-l-4)).*heaviside(x-l-2);
chi1 = chi1(end:-1:1);
chi2 = spdiags(chi1,0,Nx,Nx);
chi = kron(speye(Ny),chi2);


%%%Build stuff for inital roll 
%%%%Create initial guess and interpolation to get rolls on whole domain

%roll domain parameters
%nn = 200;
  lu=2*pi;  dxu=dx;
  xu =[dx:dx:2*pi];%[dx/2:dx:2*pi-dx/2];[0:dxu:2*pi];
  n = length(xu);

%  %generate L_mu^per
  e = ones(n,1);
  Lapu = spdiags([e -2*e e], -1:1, n, n);%Dirichlet
%Lapu(1,1) =-1;Lapu(n,n)=-1; %change toNeumann
Lapu(1,n) = 1;Lapu(n,1)=1; %change to periodic
Lu=-(dxu^(-2)*Lapu+speye(n))^2+mu0*speye(n);

%%%Pseudo-spectral far-field laplacian, if desired=
	lu = 2*pi; dxu = dy;
	xu = [dxu:dxu:2*pi]; n=length(y); n2 = Ny;
	Lapu= toeplitz([-pi^2/(3*dxu^2)-1/6 0.5*(-1).^(2:n)./sin(dxu*(1:n-1)/2).^2]);
	Lapu = ((dxu)^2)*Lapu;



size(xu)





%znew = PROF(:);
UI = reshape(znew(1:end-2),Nx,Ny)';
%figure(10)
%
%n
%size(xu)
%size(Lapu)



tic
kkx = znew(end-1);
k0 = sqrt(kkx^2+kky^2);
u = getroll(k0);


toc

%add a grid point left and right to improve interpolation 
Mx = floor(Lx/2/pi);
%xuext=[-lu+dxu:dxu:Mx*lu];[xu-lu xu xu+lu];
xuext = -lu+xu;
for ii=0:Mx-1
xuext = [xuext,xu+ii*lu];
end
uext=repmat(u,Mx+1,1);%[u; u; u];
Upol=griddedInterpolant(xuext,uext,'spline');
xi = x;mod(x,lu);
[Xi,Yi] = ndgrid(xi,y);
Ufull = Upol(kkx*Xi-Yi);
% figure(1)
% imagesc(x,y,Ufull');

%%Cutoff far-field
      vI = (1*chi)*reshape(Ufull,Nx*Ny,1);
      VI = reshape(vI(1:end),Nx,Ny)';
      
      
 %figure(2)
 %imagesc(x,y,UI+VI)
 
 
UU = UI+VI;



end

