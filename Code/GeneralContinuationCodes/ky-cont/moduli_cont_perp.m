format long
clear all
close all

%%%Script to continue perp stripes in ky, with range of c fixed

%%%Defining variables to store moduli surface in:  Will be a curve
%%%of data for each k_x,ky point on KDAT curve.
%%USES MATLAB STRUCTURES for kx and c
II = 50;
Cold = linspace(0.05,2.05,II);  %%same speed discretization as for
                                %%oblique stripes, want to match it for when perpendicular stripes exist
icc = find(Cold<0.548);
C = Cold(icc);

C = [0,0.01,0.02,0.03,0.04,C];
I = length(C)
moduli = cell(1,I);
kdat = cell(1,I);

kx_i = 0;
ky_i = 1.13%1.19;


dtcont = 2;

%p = gcp;
%delete(p)
%parpool('3x3',3)
parpool('4x2',4)
%parpool('local')
parfor i = 1:I  %%parallel
%for i=3;%i = 1:I  %%serial
    i
    c = C(i);
    s = struct;
    [KDAT,SPROF,Ind_prof,dxx,dyy,Lx,Ly] = cont(dtcont,c,kx_i,ky_i,[1,-1]);  %Continue up then down in $c$
	kdat{i} = KDAT;  %store moduli data
	moduli{i} = SPROF; %% store solution profiles, beware this variable can get very, very large


    end



l = [];
for jj = 1:I
l = [l;length(kdat{jj}(:,end))];
end

    KX = [];%zeros(sum(l),1);
    KY = KX;
    CC = KX;
for jj = 1:I
KX = [KX;kdat{jj}(:,end-1)];
KY = [KY; kdat{jj}(:,end)];
CC = [CC;C(jj)*ones(l(jj),1)];
end

    KM = sqrt(KX.^2 + KY.^2);


    figure(1)
    scatter3(KY,CC,KM,'.')
    xlabel('ky')
    ylabel('c')
    ylim([0 max(C)])
    zlabel('|k|')
    drawnow

    figure(2)
    hold on
    for ii = 1:I
      plot3(kdat{ii}(:,end),C(ii)*ones(l(ii),1), kdat{ii}(:,end-1))
      xlabel('ky')
      ylabel('c')
      ylim([0 max(C)])
      zlabel('kx')
      drawnow

    end
    hold off

    %%Define linear spreading speed detatchment curve, for kx\neq0
    KY_det = [0:0.01:1.2];
    QU1 = 1-KY_det.^2; QU1 = 2*QU1;
    QU2 = -3/4 + 2*KY_det.^2-KY_det.^4;
    KX_det = 3*(3*QU1+sqrt(7*QU1.^2+24*QU2)).^(3/2);
    KX_det = KX_det./(8*sqrt(2)*(2*QU1+sqrt(7*QU1.^2+24*QU2)));
    C_det = sqrt(2/27)*sqrt(-QU1 + sqrt(7*QU1.^2 + 24*QU2)).*(2*QU1+sqrt(7*QU1.^2+24*QU2));

        figure(3)
    scatter3(KY,CC,KX,'.')
    xlabel('ky')
    ylabel('c')
    %ylim([0 .6])
    zlabel('kx')
    drawnow

    mu0 = 0.25;
	Lx = 40*pi; Ly = 2*pi;dx = 0.25; dy = 0.25; %%have to fill these out everytime
      fname = sprintf('moduli_dat_full_perp.mat')  %%make sure to label the mu value
      save(fname,'moduli','kdat','KX','KY','CC','mu0','C','I','l','dx','dy','Lx','Ly')%,'x','y')
   % save(fname, 'moduli','x','y')  %%To save on file size, do not save moduli
%%%
