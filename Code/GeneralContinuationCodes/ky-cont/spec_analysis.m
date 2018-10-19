%%Script to analyze spectral data from one continuation profile

global mu0 mu ky Lapu dy dx dxu xu xuext n Nx Ny N Xi Yi chi k0 z1 tang  
global BDx BD2x BDy BD2y BI a III Ga tau Dyu Mx znew

%%%number of eigenvalues to track
pp = 20;






%%Pick a certain ky-slice of data
%load('ContData1c.mat')
%load('ContDataZZ9.825658e-01.mat')
%load('ContData9.980035e-01.mat')
%load('ContDataPerp9.980035e-01.mat')
load('ContDataPerp1.120102e+00.mat')
PROF1 = PROF;

[cmm,iimx] = min(PROF(end,:));
PROF = [PROF(:,iimx:-1:1),PROF(:,iimx+1:end)];


x_higher_order = false;


%Need to formula numerous things to get the jacobian back out
Lx = x(end);
Ly = y(end);
dx = x(2) - x(1);
dy = y(2) - y(1);
Nx = floor(Lx/dx);
Ny = floor(Ly/dy);
N = Nx*Ny


%%Far-field cut-off copied from cont.m
%Build quenching interface/trigger
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

        % Solve for a roll
tic
kkx = KDAT(1,end-1);
k0 = sqrt(kkx^2+ky^2);
u = getroll(k0);
toc
size(u)
%add a grid point left and right to improve interpolation 
Mx = floor(Lx/2/pi)
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
figure(100)
imagesc(x,y,Ufull');
%%%%%%

%%%%%
%%%%Build Finite difference matrices for full solution -> need neumann

Ix = speye(Nx);
Iy = speye(Ny);
E = ones(Nx,1);
Dx = spdiags([-E/2 zeros(Nx,1) E/2],-1:1,Nx,Nx);
%Neumann left
%Neumann left
Dx(1,2) = 0; Dx(end,end-1) = 0;
%Dx(1,1) = Dx(1,1) + -1/2; Dx(end,end) = Dx(end,end) + 1/2; %Neumann
if x_higher_order == true %%1/12	−2/3	0	2/3	−1/12
    clear Dx
    Dx = spdiags([(1/12)*E (-2/3)*E 0*E (2/3)*E (-1/12)*E],-2:2,Nx,Nx); ...
     Dx(1,1) = 0; Dx(1,2) = 0; Dx(1,3) = 0;
    Dx(2,1) = -1;Dx(2,2) = 0; Dx(2,3) = 1;Dx(2,4) = 0;
    %Dx(end-1,end-2) = 0;Dx(end,end) = 0;
end
Dx = Dx/dx;


D2x = spdiags([E -2*E E],-1:1,Nx,Nx);
D2x(1,1) = -1;D2x(end,end) = -1; %Neumann
if x_higher_order == true  %%−1/12	4/3	−5/2	4/3  −1/12
                           clear D2x
    D2x = spdiags([(-1/12)*E (4/3)*E (-5/2)*E (4/3)*E (-1/12)*E],- ...
                  2:2,Nx,Nx);
    D2x(1,1) = -1;D2x(1,2) = 1;D2x(1,3) = 0;
    D2x(2,1) = 1; D2x(2,2) = -2; D2x(2,3) = 1;
    D2x(end-1,end-2) = 1; D2x(end-1,end-1) = -2;D2x(end-1,end) = 1;
    D2x(end,end-1) = 0;D2x(end,end-1) = 1;D2x(end,end) = -1;
end

D2x = D2x/(dx^2);
BDx = kron(Iy,Dx);  %x-Derivative
BD2x = kron(Iy,D2x); %xx-derivative

%%y pseudo differentiation matrices
 dsp = 2*pi/Ny; sp = dsp*(1:Ny)'; Ny2 = Ny/2;
 D2t = toeplitz([-pi^2/(3*dsp^2)-1/6 0.5*(-1).^(2:Ny)./sin(dsp*(1:Ny- ...
                                                   1)/2).^2]);
 D2y =( (2*pi/Ly)^2)*D2t;

 colum  = [0 0.5*(-1).^(1:Ny-1).*cot((1:Ny-1)*dsp/2)];
 Dt = toeplitz(colum,colum([1 Ny:-1:2]));
 D1y = (2*pi/Ly)*Dt;
 
 BDy = kron(D1y,Ix); %y-derivative
BD2y = kron(D2y,Ix); %yy-derivative

%BL = -(kron(Iy,Ix) + kron(Iy,D2x)+eps*kron(D2y,Ix))^2;
BI = kron(Ix,Iy);

%%%%%

%%%%%Loop through speed c values finding spectrum
CC = KDAT(:,end);
KX = KDAT(:,end-1);
L2N = [];
SUM = [];

EFCN = cell(length(CC),1);
EVALS = cell(length(CC),1);

%%loop over all the values
for jj=1:length(PROF(1,:))-1
INP =jj;Ind_prof(jj);%(jj-1)*6+1;

%%get far-field for KX(i),ky values
kx = KX(INP)
c = CC(INP)

test = c - PROF(end,jj)

ur = getroll(sqrt(kx.^2+ky.^2));
urext = repmat(ur,Mx+1,1);
Upol = griddedInterpolant(xuext, urext,'spline');
Uphi = Upol(kx*Xi-Yi);
Vphi = reshape(Uphi,Nx*Ny,1);
Vy = BDy*Vphi;

%combine core and ff to get solution
UU = PROF1(1:end-2,jj)+chi*Vphi; %%solution
DxUU = BDx*(PROF1(1:N,jj)+chi*Vphi);%x-derivative
DyUU = BDy*(PROF1(1:N,jj)+chi*Vphi);%y-derivative
VV = reshape(UU,Nx,Ny)';
DxVV=reshape(DxUU,Nx,Ny)';
DyVV=reshape(DyUU,Nx,Ny)';
%%plot solution and its derivatives
figure(101)
imagesc(x,y,VV)
xlabel('x')
ylabel('y')
title(['speed c = ' num2str(c)])
% figure(102)
% subplot(2,1,1)
% imagesc(x,y,DxVV)
% caxis([-.1 .1])
% title(['\partial_x u and \partial_y u'])
% subplot(2,1,2)
% imagesc(x,y,DyVV)
% drawnow

%note the l2 norm
l2 = norm(UU,2);
sss = sum(UU);
SUM(jj) = sss;
L2N(jj) = l2;


%%%Plot where we are on the speed-wavenumber curve the solution 
figure(3)
plot(CC,KX,'.-',c,kx,'o')

%%Form linearization
J = Lw(kx,c) + (mu-mu0*speye(N)) - spdiags(3*UU.^2,0,N,N);

[VEC,VAL] =eigs(J,pp,0.4,'MaxIteration',300,'Tolerance',1e-11);%,'Display',1);
%[VEC,VAL] = eig(full(J));
VALD = diag(VAL);
figure(20)
plot(real(VALD),imag(VALD),'+')
title(['Spectrum of linearization for k_x = ', num2str(kx),' c = ' ...
       num2str(c)])
drawnow

[VALS,INS] = sort(real(VALD));
                 INS = INS(end:-1:1);
 SVALS = VALD(INS);
 SVEC = VEC(:,INS);
 EFCN{jj} = SVEC;
 EVALS{jj} = SVALS;
 
 
 %%Plot the eigenvectors
      
                 EVEC1= reshape(VEC(1:end,INS(1)),Nx,Ny)';
	         EVEC2= reshape(VEC(1:end,INS(2)),Nx,Ny)';
                 EVEC3= reshape(VEC(1:end,INS(3)),Nx,Ny)';
                 EVEC4 = reshape(VEC(1:end,INS(4)),Nx,Ny)';
    figure(21)
	         subplot(3,1,1)
           	imagesc(x,y,real(EVEC1))
            xlabel(['$y$'],'Interpreter','latext','FontSize',20)
            ylabel(['$y$'],'Interpreter','latext','FontSize',20)
 	       %  title(['Lam_1 = ',num2str(VALD(INS(1))),...
              %          ' Lam_2 = ',num2str(VALD(INS(2))),'Lam_3 = ',...
             %           num2str(VALD(INS(3))), ' Lam_4 = ',num2str(VALD(INS(4)))])
	         subplot(3,1,2)
	         imagesc(x,y,real(EVEC2))
            subplot(3,1,3)
	          imagesc(x,y,real(EVEC3))
                  % subplot(3,1,4)
	          %imagesc(x,y,real(EVEC4))

end


%%%Analyze first 4 evalues and eigenvectors
MAIN = [];
for jj = 1:length(CC)-1
    DUM = EVALS{jj};
    MAIN(:,jj) = DUM(1:3);
end


%%%Linear interpolate to approx. fold bifurcation location:
for jj=1:length(CC)-1
    if  MAIN(1,jj)>0.001
        jfold = jj;
        break
    end 
end

mf = (MAIN(1,jfold) - MAIN(2,jfold-1))/(CC(jfold) - CC(jfold-1));
cf = -(MAIN(2,jfold-1)/mf) + CC(jfold-1)
EE = EFCN{jfold};
Ef = EE(:,1);
Ef1 = reshape(Ef,Nx,Ny)';
figure(17)
imagesc(x,y,real(Ef1))
pbaspect([2.5 10*Ly/Lx 2])
xlabel(['x'])
ylabel(['y'])
drawnow





%%Same to find pitchfork location:
for jj=1:length(CC)-1
    if  MAIN(2,jj)>0.001
        jpf = jj;
        break
    end 
end
mp = (MAIN(2,jpf) - MAIN(3,jpf-1))/(CC(jpf) - CC(jpf-1));
cp = -(MAIN(2,jpf)/mf) + CC(jpf)
EE = EFCN{jpf+1};
Ef = EE(:,2);
Ef1 = reshape(Ef,Nx,Ny)';
figure(118)
imagesc(x,y,real(Ef1))
pbaspect([2.5 10*Ly/Lx 2])
xlabel(['x'])
ylabel(['y'])
drawnow


%%Plot eigenvalue curves
figure(4) 
plot(CC(1:end-1),[real(MAIN(1:3,:))],'x-')
xlabel(['c'])
ylabel(['Re(\lambda)'])

IND_p = [1 20 40 68];
hp=figure(50)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(CC(1:end-1),L2N,'LineWidth',2)
plot(CC(IND_p),L2N(IND_p),'o','LineWidth',2)
plot(CC(jfold),L2N(jfold),'k.','MarkerSize',20)
plot(CC(jpf),L2N(jpf),'k.','MarkerSize',20)
text(CC(IND_p(1)),L2N(IND_p(1))+0.02,'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}$', 'Interpreter', 'latex','FontSize',18)
text(CC(IND_p(2)),L2N(IND_p(2))+0.02,'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}$', 'Interpreter', 'latex','FontSize',18)
text(CC(IND_p(3)),L2N(IND_p(3))+0.02,['$\raisebox{.5pt}{\' ...
                    'textcircled{\raisebox{-.9pt} {3}}}$'], 'Interpreter', 'latex','FontSize',18)
text(CC(IND_p(4)),L2N(IND_p(4))+0.02,['$\raisebox{.5pt}{\' ...
                    'textcircled{\raisebox{-.9pt} {4}}}$'], 'Interpreter', 'latex','FontSize',18)
text(CC(jfold)+0.006,L2N(jfold)+0.02,['$\raisebox{.5pt}{\' ...
                    'raisebox{-.9pt} {S.N.}}$'], 'Interpreter', 'latex','FontSize',18)
text(CC(jpf)-0.06,L2N(jpf),['$\raisebox{.5pt}{\' ...
                    'raisebox{-.9pt} {Pitchfork}}$'], 'Interpreter', 'latex','FontSize',18)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$\| u\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
hold off
kky = ky;

hs=figure(51)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii  = 1:length(IND_p)
    %ufull = full{jj}(:,Nx*IND_s(ii)+1:Nx*(IND_s(ii)+1));
    ufull = sprof(x,y,kky,PROF(:,IND_p(ii)));%moduli{jj}(:,IND_s(ii)));
    subplot(4,1,ii)
    imagesc(x,y,ufull)
    xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    title([num2str(ii)])
end
hold off
drawnow










%load('ContData9.980035e-01.mat')
load('ContDataPerp9.980035e-01.mat')
PROF1 = PROF;




x_higher_order = false;


%Need to formula numerous things to get the jacobian back out
Lx = x(end);
Ly = y(end);
dx = x(2) - x(1);
dy = y(2) - y(1);
Nx = floor(Lx/dx);
Ny = floor(Ly/dy);
N = Nx*Ny




%%%%%%Plot some perp solutions before and after the fold curve
IND_s = [1 30 65];
hb = figure(200)

hold on 
plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
plot(PROF(end,IND_s),PROF(end-1,IND_s),'o','LineWidth',2)
text(PROF(end,IND_s(1)),PROF(end-1,IND_s(1))+0.02,'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}$', 'Interpreter', 'latex','FontSize',18)
text(PROF(end,IND_s(2))+0.02,PROF(end-1,IND_s(2))+0.02,'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}$', 'Interpreter', 'latex','FontSize',18)
text(PROF(end,IND_s(3)),PROF(end-1,IND_s(3))+0.02,'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}$', 'Interpreter', 'latex','FontSize',18)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
hold off

kky = ky;


hp=figure(201)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii  = 1:length(IND_s)
    %ufull = full{jj}(:,Nx*IND_s(ii)+1:Nx*(IND_s(ii)+1));
    ufull = sprof(x,y,kky,PROF(:,IND_s(ii)));%moduli{jj}(:,IND_s(ii)));
    subplot(4,1,ii)
    imagesc(x,y,ufull)
    xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    title([num2str(ii)])
end
hold off
drawnow





%%%%%%%%%%%load the oblique stripesload('ContData9.980035e-01.mat') and
%%%%add to the plot
load('ContData9.980035e-01.mat')
Lx = x(end);
Ly = y(end);
dx = x(2) - x(1);
dy = y(2) - y(1);
Nx = floor(Lx/dx);
Ny = floor(Ly/dy);
N = Nx*Ny
IND_s1 = 140;

[kxmm,iimx] = min(PROF(end-1,:));
PROF = [PROF(:,iimx:-1:1),PROF(:,iimx+1:end)];


figure(200)
hold on
plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
plot(PROF(end,IND_s1),PROF(end-1,IND_s1),'o','LineWidth',2)
text(PROF(end,IND_s1),PROF(end-1,IND_s1)-0.02,['$\raisebox{.5pt}{\' ...
                    'textcircled{\raisebox{-.9pt} {4}}}$'], ...
     'Interpreter', 'latex','FontSize',18)
xlim([-0.0001 max(PROF(end,:))+0.01])
ylim([-0.015 0.4092])
hold off



kky = ky
figure(201)
hold on
ufull = sprof(x,y,kky,PROF(:,IND_s1));%moduli{jj}(:,IND_s(ii)));
    subplot(4,1,4)
    imagesc(x,y,ufull)
    xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    title([num2str(4)])
    hold off
    drawnow














