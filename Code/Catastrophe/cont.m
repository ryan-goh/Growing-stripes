function [KDAT,PROF,FULL_PROF,dxx,dyy,Lx,Ly] = cont(dtcont,cc,PROF_I)


global mu0 mu ky Lapu dy dx dxu xu xuext n Nx Ny N Xi Yi chi k0 z1 tang  
global BDx BD2x BDy BD2y BI a III Ga tau Dyu Mx znew


%%Discretization parameters
dt = 0.001;
tf = 300;
Lx = 40*pi;
Ly = 2*pi;
dxx = 0.25; 
Nx = floor(Lx/dxx)+1
dx = Lx/Nx;
dyy = 0.25; 
Ny = floor(Ly/dyy)+1
dy = Ly/Ny;
a = 0;
mu0 = 0.25;
tau=0;
plotbif = true;false;true; true;
plotprof =  true;false;true;true;
TIM =true;false;false;false; %time individual steps of continuation alg.
STEPS = 5000;
SAVE = true;true;
save_int = 2;
psuedoff = true; %%Use pseudo spectral instead of finite diff for
                  %%far-field solve
x_higher_order =true; %%higher order differences in x
ydisc_fd = false;  %true = finite diff. in y, false  = psedo spec in y
spectrum = false
tau = 0;
%dir = 1;%%1 for right, -1 for left

x =[dx:dx:Lx]';[dx/2:dx:Lx-dx/2]';  Nx = length(x);[0:dx:Lx]';
y = dy*[1:Ny]';[dy:dy:Ly]';; Ny = length(y);[0:dy:Ly]';
N = Nx*Ny

%Initial wavenumbers
%mu0 = 0.25;
%ky_init= -0.001;0.895;-0.03;.87;-0.03;0.9375;-0.03;895;-0.02;0.9025;0.9093;-0.05;0.9093;0.8967;-0.01; -0.05;%.2%1.082;  %%start off finding a 1-D modulated traveling wave
%kx_init =1;0.02725;1;.04;1;0.03892;1.0;0.2014;1;0.158;0.1498;1;0.1498;0.2;1;1.0;%0.04218;sqrt(0.978^2-ky_init^2);

%c = 2;
c_init = cc;
ky_init = PROF_I(end);
kx_init = PROF_I(end-1);
ky = ky_init;

k0 = sqrt(ky_init^2+kx_init^2);


%Build quenching interface/trigger
igL = 1;[0:dy:Ly]';
igR = floor(Nx*8/10);
igW = floor(Nx*7/10);
hp = -mu0;
hm = mu0;

mu0
x(igR)

mu1 = mu0*2*(heaviside(-(x-x(igR)+dx/2))-0.5);    
%mu1= (hp)*tanh(4*(x - x(igR)-dx/2))+hp+mu0; %(-hp/2)*(tanh(10.5*(x-x(igL)))-tanh(10.5*(x-x(igR)))) + hp+mu0;
mu2 = spdiags(mu1,0,Nx,Nx);
mu = kron(speye(Ny),mu2);

%%Build far-field cutoff function
%l=Lx-40;
l = x(igW);
chi1=heaviside(x-l-4);
chi1=chi1+1/64*(3*(x-l).^5-(x-l).^6).*(1-heaviside(x-l-2)).*heaviside(x-l);
chi1=chi1+(1+1/64*(3*(x-l-4).^5+(x-l-4).^6)).*(1-heaviside(x-l-4)).*heaviside(x-l-2);
chi1 = chi1(end:-1:1);
chi2 = spdiags(chi1,0,Nx,Nx);
chi = kron(speye(Ny),chi2);


cut1 = (2-tanh(100*(x - x(igR))))/2 - 0.5;
cut = kron(speye(Ny),spdiags(cut1,0,Nx,Nx));

figure(30)
plot(x,diag(chi2),x,diag(mu2))


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

%%%Pseudo-spectral far-field laplacian, if desired
if psuedoff
	lu = 2*pi; dxu = dy/2;
	xu = [dxu:dxu:2*pi]; n=2*length(y); n2 = Ny/2;
	Lapu= toeplitz([-pi^2/(3*dxu^2)-1/6 0.5*(-1).^(2:n)./sin(dxu*(1:n-1)/2).^2]);
	Lapu = ((dxu)^2)*Lapu;
        
end


% Solve for a roll
tic
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
Ufull = Upol(kx_init*Xi-Yi);
%figure(1)
%imagesc(x,y,Ufull');

vI = cut*(1*speye(N)-1*chi)*reshape(Ufull,Nx*Ny,1);%%reshape grid of pattern to vector each spatial slice one col at a time

z1 = [PROF_I(1:end-1);c_init];%[vI; kx_init;c_init];


%test
UI = reshape(vI,Nx,Ny)';
figure(6)
imagesc(x,y,UI);
drawnow



%%%%Build Finite difference matrices for large domain

Ix = speye(Nx);
Iy = speye(Ny);
E = ones(Nx,1);
Dx = spdiags([-E/2 zeros(Nx,1) E/2],-1:1,Nx,Nx);
 %Dx(1,1) = Dx(1,1) + -1/2; Dx(end,end) = Dx(end,end) + 1/2; %Neumann
if x_higher_order == true %%1/12	−2/3	0	2/3	−1/12
    clear Dx
    Dx = spdiags([(1/12)*E (-2/3)*E 0*E (2/3)*E (-1/12)*E],-2:2,Nx,Nx);
end
Dx = Dx/dx;


D2x = spdiags([E -2*E E],-1:1,Nx,Nx);
%D2x(1,1) = -1;D2x(end,end) = -1; %Neumann
if x_higher_order == true  %%−1/12	4/3	−5/2	4/3  −1/12
                           clear D2x
    D2x = spdiags([(-1/12)*E (4/3)*E (-5/2)*E (4/3)*E (-1/12)*E],-2:2,Nx,Nx);
end

D2x = D2x/(dx^2);
BDx = kron(Iy,Dx);  %x-Derivative
BD2x = kron(Iy,D2x); %xx-derivative





if ydisc_fd 
     E = ones(Ny,1);
     D2y = spdiags([E -2*E E],-1:1,Ny,Ny);  D2y(1,end) = 1; D2y(end,1) = 1;
     D2y = D2y/(dy^2); 
         
     lowerDiag = -1/2*ones(Ny,1);
     centerDiag = -0*ones(Ny,1);
     upperDiag = 1/2*ones(Ny,1);


     D1y = spdiags([lowerDiag centerDiag upperDiag],-1:1,Ny,Ny);
     D1y(1,end) = D1y(1,end)+lowerDiag(1);
     D1y(end,1) = D1y(end,1)+upperDiag(end);
     D1y = D1y/(dy);
     
          if x_higher_order
             D1y= spdiags([(1/12)*E (-2/3)*E 0*E (2/3)*E (-1/12)* ...
                           E],-2:2,Ny,Ny);
             D1y(1,end) = -2/3;
             D1y(1,end-1) = 1/12;
             D1y(2,end) = 1/12;
             D1y(end,1) = 2/3;
             D1y(end,2) = -1/12;
             D1y(end-1,1)=-1/12;
             D1y = D1y/dy;
             
             
             D2x = spdiags([(-1/12)*E (4/3)*E (-5/2)*E (4/3)*E (-1/12)*E],-2:2,Ny,Ny);
             D2x(1,end) = 4/3;
             D2x(1,end-1) = -1/12;
             D2x(2,end) =-1/12;
             D2x(end-1,1) = -1/12;
             D2x(end,1) = 4/3;
             D2x(end,2) = -1/12;
             end
             
             
              
else
    
 dsp = 2*pi/Ny; sp = dsp*(1:Ny)'; Ny2 = Ny/2;
 D2t = toeplitz([-pi^2/(3*dsp^2)-1/6 0.5*(-1).^(2:Ny)./sin(dsp*(1:Ny- ...
                                                   1)/2).^2]);
 D2y =( (2*pi/Ly)^2)*D2t;

 colum  = [0 0.5*(-1).^(1:Ny-1).*cot((1:Ny-1)*dsp/2)];
 Dt = toeplitz(colum,colum([1 Ny:-1:2]));
 D1y = (2*pi/Ly)*Dt;
 end
 
 
BDy = kron(D1y,Ix); %y-derivative
BD2y = kron(D2y,Ix); %yy-derivative

%BL = -(kron(Iy,Ix) + kron(Iy,D2x)+eps*kron(D2y,Ix))^2;
BI = kron(Ix,Iy);



%%%Index for phase condition, basically one x,y period, located
%%%near the left boundary
%PHT =floor((Lx/2/pi-1)/2)-2;
PHT=0;%floor(Lx/2/pi)-5;%floor(Lx/4/pi)*1.0-3;

III = []; NXper = floor(2*pi/dx)
for jj = 1:Ny
    III = [III,NXper*PHT+[floor(1):NXper]+Nx*(jj-1)];
end
III=III';


%%%Test to make sure it's doing what it's supposed to
TT= reshape(Ufull,Nx*Ny,1);
TT(III) = 0;
TI = reshape(TT,Nx,Ny)';
figure(600)
imagesc(x,y,TI);
drawnow







%%%Initial newton run to get starting point without trigger
a = 1; Ga = a*(mu-mu0*speye(N));
znew = z1;
tang = zeros(Nx*Ny+2,1); tang(Nx*Ny+2:N+2) = 1;

options = optimset('Jacobian','on','Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxIter',50,'Algorithm','trust-region-reflective');
[znew,DUM,exitflag] = fsolve(@shw,z1,options);
exitflag;
%znew(end)
%znew(end-1)
%test




%%
UI = reshape(znew(1:end-2),Nx,Ny)';
u = getroll(sqrt(znew(end-1)^2+ky^2));
%xuext=[xu-lu xu xu+lu];
uext=repmat(u,Mx+1,1);%[u; u; u];
Upol=griddedInterpolant(xuext,uext,'spline');
Ufull = Upol(znew(end-1)*Xi-Yi);
%figure(1)
%imagesc(x,y,Ufull');

vI = (1*chi)*reshape(Ufull,Nx*Ny,1);
VI = reshape(vI(1:end),Nx,Ny)';



figure(1)
subplot(3,1,1)
imagesc(x,y, VI);
colorbar
subplot(3,1,2)
imagesc(x,y,UI);
colorbar
drawnow
subplot(3,1,3)
imagesc(x,y,UI+VI);
colorbar
drawnow

figure(2)
plot(shw(znew))


%%%%Code to homotope a term tau*u^2 into the nonlinearity,
%%%%hypothesized that it would break a symmetry
% %%%%%
% dtau = 0.01;


% for jj = 1:50
% 	tau = tau+dtau
% 	tang = zeros(Nx*Ny+2,1); tang(Nx*Ny+1:N+2) = 1;

% 	options = optimset	('Jacobian','on','Display','iter','TolFun',1e-6,'TolX',1e-6,'MaxIter',50,'Algorithm','trust-region-reflective');
% 	[z1,DUM,exitflag] = fsolve(@shw,z1,options);
% 	exitflag

% 	z1(end)
% 	z1(end-1)

% 	UI = reshape(z1(1:end-2),Nx,Ny)';
% 	u = getroll(sqrt(z1(end-1)^2+z1(end)^2));
% 	%add a grid point left and right to improve interpolation 
% 	xuext=[xu-lu xu xu+lu];
% 	uext=[u; u; u];
% 	Upol=griddedInterpolant(xuext,uext,'spline');
% 	xi = mod(x,lu);
% 	[Xi,Yi] = ndgrid(xi,y);
% 	Ufull = Upol(Xi-Yi);
% 	%figure(1)
% 	imagesc(x,y,Ufull');

% 	vI = (1*chi)*reshape(Ufull,Nx*Ny,1);
% 	VI = reshape(vI(1:end),Nx,Ny)';

% 	figure(1)
% 	subplot(3,1,1)
% 	imagesc(x,y, VI);
% 	colorbar
% 	subplot(3,1,2)
% 	imagesc(x,y,UI);
% 	colorbar
% 	drawnow
% 	subplot(3,1,3)
% 	imagesc(x,y,UI+VI);
% 	colorbar
% 	drawnow

% 	figure(2)
% 	plot(shw(z1))	
% end
%%%%


   

%%%Continuation algorithm; step in k_y, solve for kx and ky
%KDAT = [z1(1:Nx)',z1(N+1), z1(N+2)];
%PROF = z1;
Ind_prof = 1;
z10 = znew;
z100 = znew;
KDAT = [];
PROF = [];
FULL_PROF = [];
for direct = [1,-1]%dir%[-1 1]
dt = dtcont;


em2=zeros(N+2,1);em2(N+2)=1; 
tic
[FF,JJ]=shw(znew);
tangn=JJ\em2;
LinSolveTIme=toc
tang=1*tangn/norm(tangn,2); 



options = optimset('Jacobian','on','Display','iter','TolFun',1e-8,'TolX',1e-8,'MaxIter',50,'Algorithm','trust-region-reflective');
z1=z10+dt*0.001*tang;

tic
[znew,fdum,exitflag]=fsolve(@shw,z1,options);%solution
[FF,JJ]=shw(znew);
exitflag



tangn=JJ\em2;  
tang=tangn/norm(tangn,2);    
timer =toc  
z1 = znew;
KDAT = [KDAT;znew(1:Nx)',znew(N+1), znew(N+2)];
PROF = [PROF,znew];


UI = reshape(znew(1:end-2),Nx,Ny)';
      u = getroll(sqrt(znew(end-1)^2+ky^2));
      uext=repmat(u,Mx+1,1);
      Upol=griddedInterpolant(xuext,uext,'spline');
      Ufull = Upol(znew(end-1)*Xi-Yi);
      vI = (1*chi)*reshape(Ufull,Nx*Ny,1);
      VI = reshape(vI(1:end),Nx,Ny)';
 FULL_PROF = [FULL_PROF,UI+VI];

dt = dtcont;
if direct*tang(N+2)<0
    tang = -tang;
    display('switch')
end

z1=znew+dt*0.01*tang;
z2 = znew;
dt=dtcont;

options = optimset('Jacobian','on','Display','off','TolFun',1e-8, ...
                   'TolX',1e-8,'MaxIter',10,'Algorithm','trust-region-reflective');

%%Rough estimate of time for full number of continuation steps to run
s_co = 0;
TIME = STEPS*timer/(3600)
%%%Start the loop
counter = 1;
ConditionNum = condest(JJ);
dt = 0.0001;
for j = 1:STEPS

      

      
      
      if TIM    
	 tic
      end
 
      [znew,fdum,exitflag]=fsolve(@shw,z1,options);

      if exitflag < 1
          for iii=1:4
              dt = dt/5;
              iii;
              z1 = z2 + dt*tang;
              [znew,fdum,exitflag]=fsolve(@shw,z1,options);
		
              if exitflag>0
                % display(['reduced step size to ', num2str(dt)])
                  dt = dt*3*iii;
                  break
              end
          end
      end
      
              

      [FF,JJ]=shw(znew);tangn=JJ\em2; 

      
     if (tangn(N+1:N+2)'*tang(N+1:N+2)<0)
         tangn=-tangn;
     end
 % if (tang(N+2)<0)
 %   tang=-tang;
 % end


        % if znew(N+1)<0.015 && tang(N+2)<0.0001 && s_co<1
        %   tang(N+1:N+2) = [-tang(N+2);tang(N+1)];
        %   s_co = s_co+1;
        %   display('switch')
        %   end
       
dt =dtcont;
          
      if abs(znew(N+2))<0.025 && dt>0.9
        dt = 0.3;
	%disp('reduce step')
      end
                            
                            %if abs(znew(N+2))<0.1
                            %dt = 1;
                            %%disp('reduce step')
                            %%end
                            
                           % if abs(znew(N+2))<0.02
                          %  dt = 0.5;
                            %disp('reduce step')
                          %  end

      %if znew(N+2)>0.95
      %  dt = 0.02;
      %else
      %   dt = dtcont;
      %end

 
     if abs(znew(N+1))<0.1
          dt = 0.2;
      end
	if abs(znew(N+1))<0.02
        dt = 0.1;
	end

	%disp('reduced')

      tang=tangn/norm(tangn,2);  

      %ConditionNumNew = condest(JJ)
      % if ConditionNumNew > 5*ConditionNum 
      %     tang(N+1:N+2) = [-tang(N+2),tang(N+1)];
      %     % tang = -tang;
      %     disp('branch ??')
      % end
      % ConditionNum = ConditionNumNew;
      
      % if znew(N+1)<0.1 && znew(N+1)>0.025 && tang(N+1)>0
      %     tang = -tang;
      %     disp('switch')
      %     %tang(N+1:N+2) = [-tang(N+2),tang(N+1)];
      % end
      
      
      z1=znew;
      z2 = znew;
      z1 = z1 + dt*tang;
      KDAT = [KDAT;znew(1:Nx)', znew(N+1), znew(N+2)];
      if TIM    
          timer = toc
	end
      
      if plotbif
          % %%find any Eigenvalues and eigenvectors
              if spectrum
              if j>1 &&  mod(j,4) ==0
      	          %LinCond = cond(JJ(1:N+2,1:N+2))    
                  pp = 10; %%number of eigenvalues
       	         [VEC,VAL] = eigs(JJ(1:N+2,1:N+2),pp,0);%,'SubspaceDimension',3*pp+5);
                  	    VALD = diag(VAL);
     	           figure(20)
       	         plot(real(VALD),imag(VALD),'+')
 	         title(['k_x = ', num2str(z1(end-1)),' c = ' num2str(z1(end))])
 	         [VALS,INS] = sort(abs(VALD));
	         EVEC1= reshape(VEC(1:end-2,1),Nx,Ny)';
	         EVEC2= reshape(VEC(1:end-2,2),Nx,Ny)';
 	         figure(21)
	         subplot(1,2,1)
           	imagesc(x,y,real(EVEC1))
 	         title(['Lam = ',num2str(VAL(1)),'VecEnd= ' num2str(EVEC1...
	         	(end-1,1)),'VecEnd2= ' num2str(EVEC1(end,1))])
	         subplot(1,2,2)
	         imagesc(x,y,real(EVEC2))
      	         end
              end
              
      %%Plot bif
      figure(4)
      subplot(1,2,1)
      plot(KDAT(:,end),KDAT(:,end-1),'.-')
      title(['\mu =',num2str(mu0),',  c =', num2str(znew(end)), ', ky = ', ...
             num2str(ky),  ])
      xlabel('c')
      ylabel('k_x')
      
       subplot(1,2,2)
      plot(KDAT(:,end),sqrt(KDAT(:,end-1).^2+ky.^2))
      title(['|k| = ' num2str(sqrt(znew(end-1)^2+ky^2))])
      xlabel('c')
      ylabel('|k|')
      drawnow
      
      end
      
      if plotprof &&  mod(j,5) ==0
      %%Plotting
      
          if ConditionNum>1e18
            disp('Condition Number too large')
            break
            end
    
      UI = reshape(znew(1:end-2),Nx,Ny)';
      u = getroll(sqrt(znew(end-1)^2+ky^2));
      %add a grid point left and right to improve interpolation 
      % xuext=[xu-lu xu xu+lu];
      %uext=[u; u; u];
      uext=repmat(u,Mx+1,1);
      Upol=griddedInterpolant(xuext,uext,'spline');
      %xi = mod(x,lu);
      % [Xi,Yi] = ndgrid(xi,y);
      Ufull = Upol(znew(end-1)*Xi-Yi);
      %figure(1)
      %imagesc(x,y,Ufull');

      vI = (1*chi)*reshape(Ufull,Nx*Ny,1);
      VI = reshape(vI(1:end),Nx,Ny)';

      figure(1)
      subplot(3,1,1)
      imagesc(x,y,VI);
      %surface(Xi,Yi,VI,'EdgeColor','None');
      colorbar
      ylabel(['\chi u_p'])
      subplot(3,1,2)
      imagesc(x,y,UI)
      %surface(Xi,Yi,UI,'EdgeColor','None');
      colorbar
      ylabel(['w'])
      subplot(3,1,3)
      imagesc(x,y,UI+VI)
      %surface(Xi',Yi',UI+VI,'EdgeColor','None');
      colorbar
      ylabel(['u = w+\chi u_p'])
      drawnow
      
      end
      
      PROF = [PROF,znew];
      FULL_PROF = [FULL_PROF,UI+VI];

      
       if SAVE == true && mod(j,save_int)==0
                   PROF = [PROF,znew];
                   Ind_prof = [Ind_prof,j+1];
          
 
            fname = sprintf('ContDataZZ%d.mat',ky)
 
           save(fname, 'ky','mu0', 'KDAT','PROF','Ind_prof','x','y','mu','igL','igR','igW')
 
       end
      
      
      if exitflag<1
      exitflag
      break
      end
      
    if znew(end) < 1e-5 && counter>30
        display('c< 1e-5')
        break
    end
    
        if znew(end-1) <-0.05
        display('kx<-0.05')
        break
    end
    

    if max(znew(floor(4*Nx/8):Nx:N))<0.05            %%%%track when front detaches
          LinfNorm=max(znew(floor(3.5*Nx/8):Ny:N))
          disp('spreading')
          break
    end
    
        if max(znew(floor(0.5*Nx/8):Nx:N))>0.2 %%%%detect boundary defects
          LinfNorm=max(znew(floor(4.25*Nx/8):Ny:N))
		disp('boundary')
          break
    end
      
    counter = j;
end

counter
end









%%%


















