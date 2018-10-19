%%%%Script to process data files aftre the fact:

rerun =true; %%%filled in some data
touch_up = false;  %%Filled in one "all-detachment" slice



%%or load values from kx-ky run
load('ContData0.mat','PROF','c')
%load('ContData0Lx28.mat','PROF','c')
C = c;
KX_II = PROF(end-1,:);
KY_II = PROF(end,:);
C_II = 0*KX_II;

IK = length(KY_II);

%IND = [1 20 40 60 80 100 120 130 140 150 160 200 300];
%IND =[140:20:600];
%IND = [IND,500:50:1400];%[100:7:156];
%IND = [1:33:450];
%IND = [IND,450:20:600,650,700,800];
%IND = [200:33:450];
%IND = [IND,470:20:600,650];

%IND = 400:20:700;%%%latest
%IND = [IND,750];


%%Get ky values close to those in detachment moduli space,
KY_d = linspace(.78,1.22,50);
for ii = 1:10
     [dm,im] = min(abs(KY_II -KY_d(ii))) ;
IND_cl(ii) = im;
end

IND_cl = IND_cl(end:-1:1);
%IND_c = zeros(1,2*length(IND_cl));
%IND_c(2:2:end) = IND_cl;
%IND_c(1:2:end-1) = IND_cl-20;


IND = [500:100:1400];
IND = [IND,IND_cl,1900:100:2600];

%IND = [800:50:2000];
%IND = [IND,1700:150:IK];
%IND=find(abs(KY_II-0.842)<1e-4);


KX_I = KX_II(IND);
KY_I = KY_II(IND);
PROF_I = PROF(:,IND);
I = length(KY_I);
KY_else = [];
IND_else = [];
IND_yas = [];

count = 0;
for ii = 1:I
	kky = KY_I(ii);
	fname = sprintf('ContDataZZ%d.mat',kky);
if exist(fname,'file')
	count = count+1;
	IND_yas = [IND_yas,ii];
		display('yes')
	else
	KY_else = [KY_else,kky];
	IND_else = [IND_else,ii];
end
end









ID = 1000;
JD = count;
CG = [];
KYG = KY_I(IND_yas)'*ones(1,ID);
KXG = [];

CGt = [];
KYGt = KY_I(1:15)'*ones(1,ID);
KXGt = [];
CGb = [];
KYGb = KY_I(1:15)'*ones(1,ID);
KXGb = [];


C_sn1 =[];
KX_sn1 =[];
KY_sn1 = KY_I(IND_yas);



KX = [];
KY = [];
CC = [];

figure(20)
plot(KY_II,KX_II,KY_I,KX_I,'-o')

clear PROF



for ii = IND_yas;
	kky = KY_I(ii)
	fname = sprintf('ContDataZZ%d.mat',kky);
		load(fname)
		kxd = PROF(end-1,:);
		cd = PROF(end,:);
		    [id] = find(cd<0);
    		cd(id) = [];
    		kxd(id) = [];
		if kky<0.845
			kxd = kxd(1:end-9);
			cd = cd(1:end-9);
		end

   		 %%Find Saddle-node point
                 if ii==8
                   [csn,isn] = max(cd(1:2000));
                   kxsn = kxd(isn);
                else
                    [csn,isn] = max(cd);
                    kxsn = kxd(isn);
                end


    		C_sn1 = [C_sn1,csn];
    		KX_sn1 = [KX_sn1,kxsn];

    		tt = linspace(0,2*pi,length(cd));
    		Y = [cd;kxd];
    		sl0 = [cd(2) - cd(1);kxd(2) - kxd(1)];
    		sle = [cd(end) - cd(end-1);kxd(end) - kxd(end-1)];
    		Y = [sl0,Y,sle];
    		pp = spline(tt,Y);
    		yy = ppval(pp,linspace(0,2*pi,ID));
    		figure(1)
    		plot(yy(1,:),yy(2,:))
    		drawnow
    		CG = [CG;yy(1,:)];
    		KXG = [KXG; yy(2,:)];
		CC = [CC;cd'];
		KX = [KX;kxd'];
		KY = [KY;ky*ones(length(cd),1)];
		%%%%end of existence of file loop;


                if ii <16;

                if ii==8
                   [csn,isn] = max(cd(1:2000));
                   kxsn = kxd(isn);
                    ct = cd(1:isn);
                    kxt = kxd(1:isn);
                    cb = cd(isn+1:end);
                    kxb = kxd(isn+1:end);
                else
                  ct = cd(1:isn);
                cb = cd(isn+1:end);
                kxt = kxd(1:isn);
                kxb = kxd(isn+1:end);
                end

                %interpolate surface for top
    		tt = linspace(0,2*pi,length(ct));
    		Y = [ct;kxt];
    		sl0 = [ct(2) - ct(1);kxt(2) - kxt(1)];
    		sle = [ct(end) - ct(end-1);kxt(end) - kxt(end-1)];
    		%Y = [sl0,Y,sle];
    		pp = pchip(tt,Y);%spline(tt,Y);
    		yy = ppval(pp,linspace(0,2*pi,ID));
    		figure(1)
    		plot(yy(1,:),yy(2,:))
    		drawnow
    		CGt = [CGt;yy(1,:)];
    		KXGt = [KXGt; yy(2,:)];

                 %interpolate surface for bottom
    		tt = linspace(0,2*pi,length(cb));
    		Y = [cb;kxb];
    		sl0 = [cb(2) - cb(1);kxb(2) - kxb(1)];
    		sle = [cb(end) - cb(end-1);kxb(end) - kxb(end-1)];
    		%Y = [sl0,Y,sle];
    		pp = pchip(tt,Y);%spline(tt,Y);
    		yy = ppval(pp,linspace(0,2*pi,ID));
    		figure(2)
    		plot(yy(1,:),yy(2,:))
    		drawnow
    		CGb = [CGb;yy(1,:)];
    		KXGb = [KXGb; yy(2,:)];
                end

end






%%%Add one more additional data point(from run below), insert
%%%inbetween 15th and 16th ky slice (found by hand)
	kky = KY_II(IND(15)+3); %read off from below
	fname = sprintf('ContDataZZ%d.mat',kky);
		load(fname)
		kxd = PROF(end-1,:);
		cd = PROF(end,:);
		    [id] = find(cd<0);
    		cd(id) = [];
    		kxd(id) = [];
		if kky<0.845
			kxd = kxd(1:end-9);
			cd = cd(1:end-9);
		end

   		 %%Find Saddle-node point
    		[csn,isn] = max(cd);
    		kxsn = kxd(isn);


    		C_sn1 = [C_sn1(1:15),csn,C_sn1(16:end)];%[C_sn1,csn];
    		KX_sn1 =[KX_sn1(1:15),kxsn,KX_sn1(16:end)];% [KX_sn1,kxsn];

    		tt = linspace(0,2*pi,length(cd));
    		Y = [cd;kxd];
    		sl0 = [cd(2) - cd(1);kxd(2) - kxd(1)];
    		sle = [cd(end) - cd(end-1);kxd(end) - kxd(end-1)];
    		Y = [sl0,Y,sle];
    		pp = spline(tt,Y);
    		yy = ppval(pp,linspace(0,2*pi,ID));
    		figure(1)
    		plot(yy(1,:),yy(2,:))
    		drawnow
    		CG = [CG(1:15,:);yy(1,:);CG(16:end,:)];%[CG;yy(1,:)];
    		KXG = [KXG(1:15,:);yy(2,:);KXG(16:end,:)];%[KXG;
                                                      %yy(2,:)];
                KYG = [KYG(1:15,:);ky*ones(1,ID);KYG(16:end,:)];

		CC = [CC;cd'];
		KX = [KX;kxd'];
		KY = [KY;ky*ones(length(cd),1)];






%%now start a plotting

h = figure(1)
box on
set(gcf, 'Position', [400, 400, 700, 450]);
set(gcf,'PaperPositionMode','auto')
scatter3(KY,CC,KX,'.','LineWidth',2)



figure(2)
surf(KYG,CG,KXG,'FaceAlpha',0.75)
shading interp
hold on
plot3(KY_sn1(1:15),C_sn1(1:15),KX_sn1(1:15),'r','LineWidth',3)
hold off
view([142 33])
xlim([0.78 1.02])
ylim([0 0.2])


figure(3)
box on
set(gcf, 'Position', [400, 400, 700, 450]);
set(gcf,'PaperPositionMode','auto')
scatter3(KY,CC,KX,'.','LineWidth',2)
hold on
plot3(KY_II,C_II,KX_II,'k','LineWidth',2)
surf(KYG,CG,KXG,'FaceAlpha',0.75)
shading interp
hold off
view([142 33])


figure(4)
scatter3(KY,CC,KX,'.','LineWidth',2)
hold on
surf(KYGt,CGt,KXGt,'FaceAlpha',0.75)
shading interp
surf(KYGb,CGb,KXGb,'FaceAlpha',0.75)
shading interp
xlim([0.85 1])
hold off

%%Save just the raw data, and surf data:
save('moduli_dat_ZZBubble.mat','KYG','CG','KXG','KYGt','KXGt','CGt','KXGb','KYGb','CGb','KY','CC','KX','KY_sn1','C_sn1','KX_sn1','KX_II','KY_II','C_II')




figure(55)
surf(KYG(1:15,:),CG(1:15,:),KXG(1:15,:),'FaceAlpha',0.75)
shading interp
hold on
surf(KYG(16:end,:),CG(16:end,:),KXG(16:end,:),'FaceAlpha',0.75)
shading interp
hold off


%%%%One little extra run (comment out after getting the right data set, to
%%%%fill in a missing slice
%%%
%[KDATad,dxx,dyy,Lx,Ly] = cont(0.5,.14,PROF_III(:,i));
if touch_up
load('ContData0.mat','PROF','c')
C = c;
KX_II = PROF(end-1,:);
KY_II = PROF(end,:);

Iad = IND(15)+3;
KX_III = KX_II(Iad);
KY_III = KY_II(Iad);
PROF_III = PROF(:,Iad);
III = length(KY_III);

[KDAT,dxx,dyy,Lx,Ly] = cont(0.5,0.0001,PROF_III);
end

%%%%
%%%%
%%%% rerun the KY_II's that were missed on first run, shouldn't have to do
%%%% this any more
IND_else;
IND_else = [IND_else,9,13];
if rerun


load('ContData0.mat','PROF','c')
C = c;
KX_II = PROF(end-1,:);
KY_II = PROF(end,:);


KX_III = KX_II(IND(IND_else));
KY_III = KY_II(IND(IND_else));
PROF_III = PROF(:,IND(IND_else));
III = length(KY_III);


figure(20)
plot(KY_II,KX_II,KY_III,KX_III,'-o')

clear PROF

%kx_i = 1;
%ky_i = -0.03;
%KX_I = [.99;.9842;.9685;.9406;.9027;.846];
%KY_I = [0;.1;.2;.3042;.4004;.5];
kdat = cell(1,III);



%8.714651404000000e+03 sec for I = 10
%parpool('local',4)

dtcont =0.5;
%mu0 = 0.25  %%-> entered in cont.m

%moduli.kx = [];
%moduli.c = [];

%p = gcp;
%delete(p)
%parpool('3x3',3)
parpool('4x2',4)
%parpool('6x1',6)
%parpool('local')
%for i = 1:I  %%parallel
close all
parfor i = 1:III  %%serial
%for i = 1:III
    i
    c =C+0.0001; %C(i);
    s = struct;
    %[KDAT,SPROF,FPROF,dxx,dyy,Lx,Ly] = cont(dtcont,c,PROF_I(:,i));
	[KDAT,dxx,dyy,Lx,Ly] = cont(dtcont,c,PROF_III(:,i));
    %moduli{i} = SPROF;
    kdat{i} = KDAT;
    %full{i} = FPROF;


    end







end
