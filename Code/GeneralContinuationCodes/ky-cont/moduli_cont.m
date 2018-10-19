format long
clear all
close all



%%%Defining variables to store moduli surface in:  Will be a curve
%%%of data for each k_x,ky point on KDAT curve.
%%USES MATLAB STRUCTURES for kx and c

%%Set of speed values to continue on
I = 50;  %%number of c values
moduli = cell(1,I);
C = linspace(0.05,2.05,I); %%define c values


%%Initial wavenumber guess, start ar the ky = 0 slice of the moduli space
kx_i = 0.995;
ky_i = 0.0001;


%%Continuation step-size
dtcont = 2;


%%%%%parallel parameters and setup, change to parpool('local') to
%%%%%use defaul.  parpool('4x2',4)  -> assigns 4 workers, each
%%%%%getting 2 threads (for intel processor with hyperthreading (i7)) to use

%%Can also do in serial, changing to a for loop
%p = gcp;
%delete(p)
%parpool('3x3',3)
parpool('4x2',4)
%parpool('local')
parfor i = 1:I  %%parallel
%for i = 1:I  %%serial
    i
    c = C(i);
    s = struct;
    [KDAT,SPROF,Ind_prof,dxx,dyy,Lx,Ly] = cont(dtcont,c,kx_i,ky_i,1);  %% Continuation step
	kdat{i} = KDAT;    %%%Wavenumbers (kx,ky)
	moduli{i} = SPROF; %%%%%full solution profile (w,kx,ky)
    end



l = [];
for jj = 1:I
l = [l;length(kdat{jj}(:,end))];
end
    KX = [];
    KY = KX;
    CC = KX;
for jj = 1:I-3
KX = [KX;kdat{jj}(:,end-1)];
KY = [KY; kdat{jj}(:,end)];
CC = [CC;C(jj)*ones(l(jj),1)];
end

    KM = sqrt(KX.^2 + KY.^2);

%%Plot moduli space using bulk wavenumber
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


%%Plot scatter plot of moduli space
        figure(3)
    scatter3(KY,CC,KX,'.')
    xlabel('ky')
    ylabel('c')
    %ylim([0 .6])
    zlabel('kx')
    drawnow

    mu0 = 0.25;
	Lx = 40*pi; Ly = 2*pi;dx = 0.25; dy = 0.25; %%have to fill these out everytime
      fname = sprintf('moduli_dat_full.mat')  %%make sure to label the mu value
      save(fname,'moduli','kdat','KX','KY','CC','mu0','C','I','l','dx','dy','Lx','Ly','x','y')  %%To save on file size, do not save moduli, note solution profiles are saved in data files for each run
   % save(fname, 'moduli','x','y')
%%%
