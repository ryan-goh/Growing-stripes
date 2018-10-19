format long
clear all
close all

%%%Script to continue perpendicular stripes, if ky< absolute spectrum triple point
%%% continues nearby oblique stripes at reattachment point

%%%Defining variables to store moduli surface in:  Will be a curve
%%%of data for each k_x,ky point on KDAT curve.
%%USES MATLAB STRUCTURES for kx and c
I = 50;  %%%number of ky values to run

KY_I = linspace(0.78,1.22,I); %%KY -values
KX_I = 0*KY_I;  %%initial kx values, equal to zero for perpendicular stripes
c = 0;  %Initial speed guess

%%%or load values from kx-ky run
  %load('ContData_c-3.000000e-01.mat','PROF','c')
  %C = c;
  %KX_II = PROF(end-1,:);
  %KY_II = PROF(end,:);

  %clear PROF
  %IND = [1 20 40 60 80 100 120 130 140 150 160 200 300];

  %KX_I = KX_II(IND);
  %KY_I = KY_II(IND);

  % figure(20) %%Plot to see which values are ones starting continuation
  % plot(KY_II,KX_II,KY_I,KX_I,'o-')


%kx_i = 1;
%ky_i = -0.03;
%KX_I = [.99;.9842;.9685;.9406;.9027;.846];
%KY_I = [0;.1;.2;.3042;.4004;.5];
%I = length(KX_I);
moduli = cell(1,I);


dtcont = 1.5;

%%%Example parallel poopl configurations, (have to define in MATLAB pool preferences) to make multiple runs at once. Beware of memory restraints if too many parallel workers run at once
%p = gcp;
%delete(p)
%parpool('3x3',3)  %%three workers, each has access to two computational threads (on an intel core i7)
%parpool('4x2',4)   %%four workers, two computational threads each
%parpool('local')  %%default parallel pool
parfor i = 1:I  %%parallel, set
%for i = 1:I  %%serial, run one continuation at a time.
    i
    s = struct;
	cadd = 0.1;  %%set c start point slightly different for each ky value
	if KY_I(i)<0.85
		cadd = 0.02;
	end
	if KY_I(i)<0.8
		cadd = 0.03;
	end
    [KDAT,dxx,dyy,Lx,Ly] = cont(dtcont,c+cadd,KX_I(i),KY_I(i),[1,-1]);  %%makes a continuation run in c, perp stripes
    %%Perp stripes

    [c_det,i_det] = max(KDAT(:,end)); %%Find fold point in c, pitchfork will be nearby.

   [KDATob,dxx,dyy,Lx,Ly] = cont(dtcont,c_det+0.01,KX_I(i)+sqrt(0.01),KY_I(i),[-1,1]);  %%continue in oblique stripes, just make an initial guess with kx non-zero.
   %%Obliquestripes

   KDAT = [KDAT;KDATob];

   moduli{i} = [KDAT]; %%Combine data and save in moduli, to get
                       %%solution profiles, see data files saved via
                       %%each run in cont.m
end





l = [];
for jj = 1:I
l = [l;length(moduli{jj}(:,end))];
end

    KX = [];%zeros(sum(l),1);
    KY = KX;
    CC = KX;
for jj = 1:I
KX = [KX;moduli{jj}(:,end-1)];
CC = [CC; moduli{jj}(:,end)];
KY = [KY;KY_I(jj)*ones(l(jj),1)];
end

    KM = sqrt(KX.^2 + KY.^2);

%%Scatter plot of moduli Data
    figure(1)
    scatter3(KY,CC,KM,'.')
    xlabel('ky')
    ylabel('c')
    ylim([0 max(CC)])
    zlabel('k_x')
    drawnow
%%Compilation of line:
figure(11)
hold on
for ii = 1:I
  plot3(KY_I{ii}(l(jj,1)),moduli{ii}(:,end),moduli{ii}(:,end-1),'LineWidth',2,'Color',[0 .44 .72])
  xlabel('ky')
  ylabel('c')
  ylim([0 max(CC)])
  zlabel('k_x')
  drawnow
end
  hold off




%%Calculate linear predictions for detachment curve, from pinched double-root analysis,
    %%perp transition speed:
    KYI = linspace(min(KY_I),max(KY_I),100);
C0 = (3+(-8).*KYI.^2+4.*KYI.^4).^(-1).*((6/7).^(1/2).*((-6)+6.*KYI.^2+ ...
  3.^(1/2).*((-9)+32.*KYI.^2+(-16).*KYI.^4).^(1/2)).^(1/2)+(-11).*( ...
  2/21).^(1/2).*KYI.^2.*((-6)+6.*KYI.^2+3.^(1/2).*((-9)+32.*KYI.^2+( ...
  -16).*KYI.^4).^(1/2)).^(1/2)+4.*(6/7).^(1/2).*KYI.^4.*((-6)+6.* ...
  KYI.^2+3.^(1/2).*((-9)+32.*KYI.^2+(-16).*KYI.^4).^(1/2)).^(1/2)+( ...
  -4).*(2/21).^(1/2).*KYI.^6.*((-6)+6.*KYI.^2+3.^(1/2).*((-9)+32.* ...
  KYI.^2+(-16).*KYI.^4).^(1/2)).^(1/2)+(1/14).*42.^(-1/2).*((-6)+6.* ...
  KYI.^2+3.^(1/2).*((-9)+32.*KYI.^2+(-16).*KYI.^4).^(1/2)).^(3/2)+( ...
  -10/21).*(2/21).^(1/2).*KYI.^2.*((-6)+6.*KYI.^2+3.^(1/2).*((-9)+ ...
  32.*KYI.^2+(-16).*KYI.^4).^(1/2)).^(3/2)+(5/21).*(2/21).^(1/2).* ...
  KYI.^4.*((-6)+6.*KYI.^2+3.^(1/2).*((-9)+32.*KYI.^2+(-16).*KYI.^4) ...
  .^(1/2)).^(3/2)+(-10/147).*(2/21).^(1/2).*((-6)+6.*KYI.^2+3.^(1/2) ...
  .*((-9)+32.*KYI.^2+(-16).*KYI.^4).^(1/2)).^(5/2)+(10/147).*(2/21) ...
  .^(1/2).*KYI.^2.*((-6)+6.*KYI.^2+3.^(1/2).*((-9)+32.*KYI.^2+(-16) ...
  .*KYI.^4).^(1/2)).^(5/2)+(-1/294).*42.^(-1/2).*((-6)+6.*KYI.^2+ ...
  3.^(1/2).*((-9)+32.*KYI.^2+(-16).*KYI.^4).^(1/2)).^(7/2));


  %%Lin spreading speed
 C1= 2.*(2/3).^(1/2).*((-2)+2.*KYI.^2+2.^(1/2).*(5+(-4).*KYI.^2+2.* ...
  KYI.^4).^(1/2)).^(1/2)+(-2).*(2/3).^(1/2).*KYI.^2.*((-2)+2.* ...
  KYI.^2+2.^(1/2).*(5+(-4).*KYI.^2+2.*KYI.^4).^(1/2)).^(1/2)+(1/3).* ...
  (2/3).^(1/2).*((-2)+2.*KYI.^2+2.^(1/2).*(5+(-4).*KYI.^2+2.*KYI.^4) ...
  .^(1/2)).^(3/2);

  KX1=(2.*(2/3).^(1/2).*((-2)+2.*KYI.^2+2.^(1/2).*(5+(-4).*KYI.^2+2.* ...
  KYI.^4).^(1/2)).^(1/2)+(-2).*(2/3).^(1/2).*KYI.^2.*((-2)+2.* ...
  KYI.^2+2.^(1/2).*(5+(-4).*KYI.^2+2.*KYI.^4).^(1/2)).^(1/2)+(1/3).* ...
  (2/3).^(1/2).*((-2)+2.*KYI.^2+2.^(1/2).*(5+(-4).*KYI.^2+2.*KYI.^4) ...
  .^(1/2)).^(3/2)).^(-1).*(3.^(-1/2).*(6+(-6).*KYI.^2+2.^(1/2).*(5+( ...
  -4).*KYI.^2+2.*KYI.^4).^(1/2)).^(1/2).*((-2)+2.*KYI.^2+2.^(1/2).*( ...
  5+(-4).*KYI.^2+2.*KYI.^4).^(1/2)).^(1/2)+(-1).*3.^(-1/2).*KYI.^2.* ...
  (6+(-6).*KYI.^2+2.^(1/2).*(5+(-4).*KYI.^2+2.*KYI.^4).^(1/2)).^( ...
  1/2).*((-2)+2.*KYI.^2+2.^(1/2).*(5+(-4).*KYI.^2+2.*KYI.^4).^(1/2)) ...
  .^(1/2)+(1/8).*3.^(-1/2).*(6+(-6).*KYI.^2+2.^(1/2).*(5+(-4).* ...
  KYI.^2+2.*KYI.^4).^(1/2)).^(1/2).*((-2)+2.*KYI.^2+2.^(1/2).*(5+( ...
  -4).*KYI.^2+2.*KYI.^4).^(1/2)).^(3/2));


 C2=(-2).*(2/3).^(1/2).*((-2)+2.*KYI.^2+(-1).*(13+(-32).*KYI.^2+16.* ...
  KYI.^4).^(1/2)).^(1/2)+2.*(2/3).^(1/2).*KYI.^2.*((-2)+2.*KYI.^2+( ...
  -1).*(13+(-32).*KYI.^2+16.*KYI.^4).^(1/2)).^(1/2)+(-1/3).*(2/3).^( ...
  1/2).*((-2)+2.*KYI.^2+(-1).*(13+(-32).*KYI.^2+16.*KYI.^4).^(1/2)) ...
  .^(3/2);


    %%Overlay of scatter plot data with all-stripe detachment predictions
        figure(31)
    scatter3(KY,CC,KX,'.')
    hold on
    plot3(KYI,C0,KYI*0)
    plot3(KYI,C1,KX1)
    plot3(KYI,C2,KYI*0)
    hold off
    xlabel('ky')
    ylabel('c')
    %ylim([0 .6])
    zlabel('kx')
    legend(['Data'],[ 'C0 -> Abs. Spec. trans'],['C1 -> obl. lin. sp.'],['C2->Perp lin. sp.'])
    drawnow


%%To plot interpolated surfaces, see file /Code/PerpDetachment/plot-fig.m

%%5Save data
    mu0 = 0.25;
	Lx = 40*pi; Ly = 2*pi;dx = 0.25; dy = 0.25; %%have to fill these out everytime  :()
      fname = sprintf('moduli_dat_perp.mat')  %%make sure to label the mu value
      save(fname,'moduli','KX','KY','CC','mu0','c','I','l','dx','dy','Lx','Ly')
   % save(fname, 'moduli','x','y')
%%%
