format long
clear all
close all



%%%Defining variables to store moduli surface in:  Will be a curve
%%%of data for each k_x,ky point on KDAT curve.
%%USES MATLAB STRUCTURES for kx and c
%I = 40;
%C  = linspace(0.005,0.04,I); %full
%C = linspace(0.4,1.95,I); %avoid small
%C = linspace(0.05,2.0,I);

%%Read values off of moduli_full_kx.dat
%C = 0.4808;

%%or load values from kx-ky run
load('ContData0.mat','PROF','c')  %%%%This loads a set of solutions
                                  %%%%from the ky continuation code
                                  %%%%(i.e. run
                                  %%%% "cont(1,0,0.001,1.02,-1)"
                                  %%%% in the GeneralModuliSpace directory

C = c;
KX_II = PROF(end-1,:);
KY_II = PROF(end,:);

IK = length(KY_II);

%IND = [1 20 40 60 80 100 120 130 140 150 160 200 300];
%IND =[140:20:600];
%IND = [IND,500:50:1400];%[100:7:156];
%IND = [1:33:450];
%IND = [IND,450:20:600,650,700,800];


IND = 400:20:700;
IND = [IND,750]; %%%%Pick set of ky values on zig-zag c = 0
                 %%%%critical curve to continue in c>0,


%IND = [200:33:450];
%IND = [IND,470:20:600,650];

KX_I = KX_II(IND);
KY_I = KY_II(IND);
PROF_I = PROF(:,IND);


%%See which values you are starting from
figure(20)
plot(KY_II,KX_II,KY_I,KX_I,'o-')

clear PROF

%kx_i = 1;
%ky_i = -0.03;
%KX_I = [.99;.9842;.9685;.9406;.9027;.846];
%KY_I = [0;.1;.2;.3042;.4004;.5];
I = length(KX_I);
moduli = cell(1,I);
kdat = cell(1,I);
full = cell(1,I);


%8.714651404000000e+03 sec for I = 10
%parpool('local',4)  %%if no parallel profile saved (need at least 4 open compute threads)

dtcont =0.5;
%mu0 = 0.25  %%-> entered in cont.m

%moduli.kx = [];
%moduli.c = [];

%p = gcp;
%delete(p)
%parpool('3x3',3)
parpool('4x2',4)
%parpool('local')
%for i = 1:I  %%parallel
parfor i = 1:I  %%serial
%for i = 4
    i
    c =C+0.0001; %C(i);
    s = struct;
    [KDAT,SPROF,FPROF,dxx,dyy,Lx,Ly] = cont(dtcont,c,PROF_I(:,i));
    moduli{i} = SPROF;
    kdat{i} = KDAT;
    full{i} = FPROF;


    end



l = [];
I_sn = [];
KY_sn = [];
KX_sn = [];
C_sn = [];


%%Locate the saddle-node curve on the kink-drag buble.
KY_sn = PROF_I(end,:);
for jj = 1:I
l = [l;length(kdat{jj}(:,end))];
[ff,gg] = min(kdat{jj}(:,end-1));
I_sn = [I_sn;gg];
KX_sn = [KX_sn;kdat{jj}(gg,end-1)];
C_sn = [C_sn; kdat{jj}(gg,end)];

end

    KX = [];%zeros(sum(l),1);
    KY = KX;
    CC = KX;
for jj = 1:I
KX = [KX;kdat{jj}(:,end-1)];
CC = [CC; kdat{jj}(:,end)];
KY = [KY;KY_I(jj)*ones(l(jj),1)];
ax = gca;
ax.ColorOrderIndex = 1;
figure(15)
hold on
plot3(KY_I(jj)*ones(l(jj),1),kdat{jj}(:,end),kdat{jj}(:,end-1),'.-','LineWidth',1)

end
ax = gca;
ax.ColorOrderIndex = 3;
plot3(KY_II,KY_II*0,KX_II,'LineWidth',2)

    plot3(KY_sn,C_sn,KX_sn,'g','LineWidth',2)
    hold off
    xlabel('ky')
    ylabel('c')
    xlim([0.989 1])
    ylim([0 3e-3])
    zlim([0 .15])
    zlabel('kx')
    drawnow

hold off

    KM = sqrt(KX.^2 + KY.^2);
    KM_II = sqrt(KX_II.^2 + KY.^2);


    figure(1)
    scatter3(KY,CC,KM,'.')
    hold on
    plot3(KY_II,KY_II*0,KM_II)
    hold off
    xlabel('ky')
    ylabel('c')
    % ylim([0 max(C)])
    zlabel('|k|')
    drawnow


figure(5)
 scatter3(KY,CC,KX,'.','LineWidth',0.1)
    hold off
    xlabel('k_y')
    ylabel('c')
    %ylim([0 2])
    zlabel('k_x')
    ax = gca;
    ax.FontSize = 16;



        figure(3)
    scatter3(KY,CC,KX,'.')
    hold on
    plot3(KY_II,KY_II*0,KX_II,'LineWidth',2)
    plot3(KY_sn,C_sn,KX_sn,'g','LineWidth',2)
    hold off
    xlabel('ky')
    ylabel('c')
   xlim([0.89 1])
    zlabel('kx')
    drawnow




    mu0 = 0.25;
	Lx = 40*pi; Ly = 2*pi;dx = 0.25; dy = 0.25; %%have to fill these out everytime
      fname = sprintf('moduli_dat_zz.mat')  %%make sure to label the mu value
      save(fname,'kdat','moduli','full','KX_sn','KY_sn','C_sn','KX','KY','CC','KY_I','mu0','C','I','l','dx','dy','Lx','Ly')  %%Note: omit 'moduli' and 'full' variables to reduce file size, though plotting scripts plot_fig.m needs these.
   % save(fname, 'moduli','x','y')
%%%
