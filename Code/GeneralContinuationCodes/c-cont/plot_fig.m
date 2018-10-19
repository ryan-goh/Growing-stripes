%%%%Plotting figures for reattachment curves
global mu0
%%loads the ky continuation (with c=0) which the bubble emanates from, then
%%clears other data
makevideo = 0;
flipkyc = 0;  %%plots versions of most figures with ky and c axes flipped.
makesolutionvideo = 0;

load('moduli_dat_perp.mat')
dx = 0.25; dy = 0.25; Lx = 40*pi;
II = length(moduli);
KY_II = linspace(1/sqrt(2)+0.08,1.22,II);
I = II-6;
KY_I=KY_II(1:I);

moduli_ob = cell(1,I);
moduli_perp = cell(1,I);

%%Sorts moduli data into perp & oblique
for ii = 1:I
kxs = moduli{ii}(:,end-1);
cs = moduli{ii}(:,end);
kys = KY_I(ii);
[ikbig] = find(abs(kxs)>1e-5 & cs>=0);
[iksm] = find(abs(kxs)<=1e-5 & cs>=0);
kxso = kxs(ikbig);
cso = cs(ikbig);
kxsp = kxs(iksm);
csp = cs(iksm);

moduli_ob{ii} = [kxso,cso];
moduli_perp{ii} = [kxsp,csp];



% figure(111)
% plot(cso,kxso,'o')
% hold on
% plot(csp,kxsp,'x')
% plot(moduli{ii}(:,end),moduli{ii}(:,end-1))
% hold off

end
%%%%
%%%%



%%%FULL moduli data vectors
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

    %%%%%%perp surface
    %%%%%
    CGp = [];
    KYGp = [];
    KXGp = [];
    ID = I;
    JD = 1000;
    KYGp = KY_I(1:ID)'*ones(1,JD);


    C_sn = [];
    C_pf = [];
    KY_sn = KY_I;


    for ii = 1:I
        cd = moduli_perp{ii}(:,end);
        kxd = moduli_perp{ii}(:,end-1);
        md = length(cd);


        [cf,iif] = max(cd);
        kx_f = kxd(iif)
        C_sn = [C_sn,cf];

%         %%find pitchfork point
%         [kx_pf,ipf] = min(abs(kxd(jd:end)));
%         cpf = cd(ipf+jd-1);
%         C_pf = [C_pf;cpf];

        %[id] = find(kxd<0);
        %kxd(id) = [];
        %cd(id) = [];
        tt = linspace(0,2*pi, length(cd));
        Y = [cd';kxd'];
        sl0 = [1,0];
        sle = [cd(end) - cd(end-1);kxd(end) - kxd(end-1)];
        pp = spline(tt,Y);
        yy = ppval(pp,linspace(0,2*pi,JD));
        figure(1)
%         plot(yy(1,:),yy(2,:),'.',cd,kxd,'x')
%
%         figure(3)
%         plot(cd,kxd,'.',cf,kx_f,'x')
%
        CGp = [CGp;yy(1,:)];
        KXGp = [KXGp;yy(2,:)];

    end

    %%%%
    %%%%oblique surface
        CGo = [];
    KYGo = [];
    KXGo = [];
    ID = I;
    JD = 1000;
    KYGo = KY_I(1:ID)'*ones(1,JD);



    C_pf = [];
    KY_sn = KY_I;


    for ii = 1:I
        cd = moduli_ob{ii}(:,end);
        kxd = moduli_ob{ii}(:,end-1);
        md = length(cd);
        [kdum,ikm]  = min(kxd);
        kxd = [kxd(ikm:-1:1);kxd(ikm+1:end)];
        cd = [cd(ikm:-1:1);cd(ikm+1:end)];
         if ii ==30 || ii==31
            cd = cd(1:end-1);
            kxd = kxd(1:end-1);
         end

        %%find pitchfork point
        [kx_pf,ipf] = min(abs(kxd));
        cpf = cd(ipf);
        C_pf = [C_pf;cpf];

        %[id] = find(kxd<0);
        %kxd(id) = [];
        %cd(id) = [];
        tt = linspace(0,2*pi, length(cd));
        Y = [cd';kxd'];
        sl0 = [1,0];
        sle = [cd(end) - cd(end-1);kxd(end) - kxd(end-1)];
        pp = spline(tt,Y);
        yy = ppval(pp,linspace(0,2*pi,JD));
%         figure(1)
%         plot(yy(1,:),yy(2,:),'.')
%
%         figure(3)
%         plot(cd,kxd,'.',cpf,kx_pf,'x')
%
        CGo = [CGo;yy(1,:)];
        KXGo = [KXGo;yy(2,:)];

    end







    figure(1)
    surf(KYGo,CGo,KXGo)
    shading interp


















    %%%%%load in add on data to fill in some holes
    %%%%%
    load('moduli_dat_perp1.mat','moduli_ob','moduli_perp','minus')
dx = 0.25; dy = 0.25; Lx = 40*pi;
Ia = length(moduli_perp)
KY_Ia = linspace(1.14,1.23,Ia);%linspace(1/sqrt(2)+0.08,1.22,I);

l = [];
l_ob = [];
l_perp = [];

for jj = 1:Ia
l_ob = [l_ob;length(moduli_ob{jj}(:,end))];
l_perp = [l_perp;length(moduli_perp{jj}(:,end))];
end
l = l_ob+l_perp;
    KXa = [];%zeros(sum(l),1);
    KYa = KXa;
    CCa = KXa;
for jj = 1:Ia
KX = [KX;moduli_perp{jj}(:,end-1);moduli_ob{jj}(1:end+minus(jj),end-1)];
CC = [CC; moduli_perp{jj}(:,end);moduli_ob{jj}(1:end+minus(jj),end)];
KY = [KY;KY_Ia(jj)*ones(l(jj)+minus(jj),1)];
end
    KM = sqrt(KX.^2 + KY.^2);


    CGob = [];
    KYGob = [];
    KXGob= [];

        CGperp = [];
    KYGperp = [];
    KXGperp = [];

    C_sna = [];
    C_pfa = [];
    KY_sna = KY_I;
    IDa = Ia;
    JDa = 1000;
    KYGa = KY_Ia(1:IDa)'*ones(1,JD);

    skipcounter = 0;

    for ii = 1:Ia

        %%find fold

        kxs = moduli_perp{ii}(:,end-1);
    cs = moduli_perp{ii}(:,end);
    kys = KY_I(ii);
    [ikbig] = find(abs(kxs)>1e-5 & cs>=0);
    [iksm] = find(abs(kxs)<=1e-5 & cs>=0);
    kxso = kxs(ikbig);
    cso = cs(ikbig);
    kxsp = kxs(iksm);
    csp = cs(iksm);



%         cd = moduli_perp{ii}(:,end);
%         kxd = moduli_perp{ii}(:,end-1);
%         md = length(cd);

        [cf,iif] = max(csp);
        kx_f = kxd(iif);
        C_sn = [C_sn,cf];

                tt = linspace(0,2*pi, length(kxsp));
        Y = [csp';kxsp'];
        sl0 = [1,0];
        sle = [csp(end) - csp(end-1);kxd(end) - kxd(end-1)];
        pp = spline(tt,Y);
        yy = ppval(pp,linspace(0,2*pi,JD));
%         figure(1)
%         plot(yy(1,:),yy(2,:),'.',csp(iif),kxsp(iif),'x')
           CGperp = [CGperp;yy(1,:)];
        KXGperp = [KXGperp;yy(2,:)];

        %%find pitchfork point
        minu = minus(ii);
        cdfs = moduli_ob{ii}(1:end+minu,end);
        kxf = moduli_ob{ii}(1:end+minu,end-1);
        [kdum,ikm]  = min(kxf);
        kxf = [kxf(ikm:-1:1);kxf(ikm+1:end)];
        cdfs = [cdfs(ikm:-1:1);cdfs(ikm+1:end)];


        [iiss] = find(cdfs>=0);
        cdf = cdfs(iiss);

        [kx_pf,ipf] = min(abs(kxf));
        cpf = cdf(ipf);
        C_pf = [C_pf;cpf];







         tt = linspace(0,2*pi, length(cdf));
        Y = [cdf';kxf'];
        sl0 = [1,0];
        sle = [cdf(end) - cdf(end-1);kxf(end) - kxf(end-1)];
        pp = spline(tt,Y);
        yy = ppval(pp,linspace(0,2*pi,JD));
%         figure(2)
%         plot(yy(1,:),yy(2,:),'.',cdf(ipf),kxf(ipf),'x')
%
%                 figure(3)
%         plot(cd,kxd,'.',cdf,kxf,'.',cpf,kx_pf,'x')

        if KY_Ia(ii)<sqrt(4+sqrt(3))/2;
            skipcounter = skipcounter+1;
        CGob = [CGob;yy(1,:)];
        KXGob = [KXGob;yy(2,:)];
        end

    end
    %%%%%%%%



  %%perp transition speed:
    KYI = linspace(0.0,1.23,2500);
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



[iic] = find(KYI>=sqrt(4+sqrt(3))/2);%find(abs(C2)==real(C2))

  icull = find(KY>1.189 & KX>4e-3 & CC>=0);
 KY(icull) = [];%%%removes extraneous points where Oblique-stripes don't exist, but one data point was saved
    CC(icull) = [];
    KX(icull) = [];

figure(22)
scatter3(KY,CC,KX,'.')
hold on
plot3([KY_I,KY_Ia],C_sn,0*C_sn,'LineWidth',2)
plot3([KY_I,KY_Ia],C_pf,0*C_pf,'LineWidth',2)
surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',0.7)
shading interp
surf(KYGa,CGperp,KXGperp,'FaceAlpha',0.7)
shading interp
hold off


    %     KXGob(:,6:10) = [];
 figure(34)
  box on
    set(gcf, 'Position', [100, 100, 900, 600]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',0.6)
 shading interp
    hold on
    shading interp
    surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
    shading interp
      area([KY_sn,KY_Ia],C_sn,'FaceColor',[0 0.0 0.5],'FaceAlpha',0.6)
      plot3(0,0,0)
      plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',2)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',2)
    plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',20)
    xlim([.76,1.225])
    hold off


    %%load data from Full Bubble  computations
    ZZs= load('moduli_dat_ZZBubble.mat');
    %%Continuation data points
    KXZZ =ZZs.KX;
    KYZZ = ZZs.KY;
    CCZZ = ZZs.CC;

     %%Approximate Saddle-Node Data
    KXZZ_sn = ZZs.KX_sn1;
    KYZZ_sn = ZZs.KY_sn1;
    CCZZ_sn = ZZs.C_sn1;

    %%Surface data
    KXZZG = ZZs.KXG;
    KYZZG = ZZs.KYG;
    CCZZG = ZZs.CG;

        %%Surface data
    KXZZGt = ZZs.KXGt;
    KYZZGt = ZZs.KYGt;
    CCZZGt = ZZs.CGt;

     KXZZGb = ZZs.KXGb;
    KYZZGb = ZZs.KYGb;
    CCZZGb = ZZs.CGb;

    %%zig-zag Critical curve
    KXZZr = ZZs.KX_II;
    KYZZr = ZZs.KY_II;
    CCZZr = ZZs.C_II






    figure(333)
      box on
    set(gcf, 'Position', [100, 100, 600, 400]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',0.6)
 shading interp
    hold on
    shading interp
    surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
    shading interp
      area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',0.6)
       scatter3(KY,CC,KX,'.')
       plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3)
    plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',20)
    hold off
    xlim([0.76, 1.225])
    ylim([0 , 1.4])
    zlim([0, 0.7])
    view([-28 54])
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
    colormap(parula(256))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)




    %%Different viewpoint
      figure(334)
      box on
    set(gcf, 'Position', [100, 100, 600, 400]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',0.6)
 shading interp
    hold on
    shading interp
    surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
    shading interp
      area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',0.6)
       scatter3(KY,CC,KX,'.')
       plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3)
    plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',20)
    hold off
    xlim([0.76, 1.225])
    ylim([0 , 1.4])
    zlim([0, 0.7])
    view([0 90])
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
    colormap(parula(256))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)





    if flipkyc

        %%Different viewpoint
      figure(9334)
          box on
    set(gcf, 'Position', [100, 100, 600, 400]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',0.6)
 shading interp
    hold on
    shading interp
    surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
    shading interp
      area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',0.6)
       scatter3(KY,CC,KX,'.')
       plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3)
    plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',20)
    hold off
    xlim([0.76, 1.225])
    ylim([0 , 1.4])
    zlim([0, 0.7])
    view([90 90])
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
    colormap(parula(256))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)


    end


%%With Zig-Zag buuble
jjj = find(   abs(KXGo(6,:) - .39)<2e-4);
jjf = find(abs(KXGo(7,:) - .4042)<2e-4);
jdf = jjf - jjj;
  JJ = 1:jjf;
JJ([2:jdf+1]*6) = [];



lll = find(   abs(KXGo(7,:) - .3396)<2e-4);
lll = lll(1);
llf = find(   abs(KXZZG(16,:) - .4119)<1e-4);
llf = llf(1)
ldf = (1000-llf) - (1000-lll);
LL = [llf:1000];
LL([1:ldf-10]*2) = [];
LL([1:10]*7) = [];
JF = [lll:1000];


MMM = [1:1000];
MMM(2:2:end) = [];
imm = 1000-llf-length(MMM);
MMM([1:imm]*3) = [];








    %%With Zig-Zag buuble-different view, just interpolated surface
    fa =1;
    figure(336)
         box on
    set(gcf, 'Position', [100, 100, 900, 600]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',fa)
 shading interp
    hold on
    plot3(-1,-1,-1)
        plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',40)
    plot3(KYZZ_sn(15),CCZZ_sn(15),KXZZ_sn(15),'g.','MarkerSize',40)
    surf(KYGo(1:6,:),CGo(1:6,:),KXGo(1:6,:),'FaceAlpha',fa)
    shading interp
    surf(KYGo(7:end-5,:),CGo(7:end-5,:),KXGo(7:end-5,:),'FaceAlpha',fa)
        shading interp
    surf([KYGa(1,:);KYGo(end-5,:)],[CGob(1,:);CGo(end-5,:)],[KXGob(1,:);KXGo(end-5,:)],'FaceAlpha',fa)
        shading interp
     shading interp
    surf(KYZZGt(1:15,:),CCZZGt(1:15,:),KXZZGt(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZGb(1:15,:),CCZZGb(1:15,:),KXZZGb(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZG(16:end,:),CCZZG(16:end,:),KXZZG(16:end,:),'FaceAlpha',fa)
        shading interp
    surf([KYGo(6,1:jjj);KYGo(7,JJ)],[CGo(6,1:jjj);CGo(7,JJ)],[KXGo(6,1:jjj);KXGo(7,JJ)],'FaceAlpha',fa)
        shading interp
    surf([KYGo(7,JF);KYZZG(16,LL)],[CGo(7,JF);CCZZG(16,LL)],[KXGo(7,JF);KXZZG(16,LL)],'FaceAlpha',fa)
        shading interp
    surf([KYZZGt(15,MMM);KYZZG(16,1:llf)],[CCZZGt(15,MMM);CCZZG(16,1:llf)],[KXZZGt(15,MMM);KXZZG(16,1:llf)],'FaceAlpha',0.6)
        shading interp
    area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',fa)
    %scatter3(KY,CC,KX,'.')
    plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3,'LineSmoothing','on')
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3,'LineSmoothing','on')
    plot3(KYZZ_sn(1:15),CCZZ_sn(1:15),KXZZ_sn(1:15),'r','LineWidth',3,'LineSmoothing','on')
   %scatter3(KYZZ,CCZZ,KXZZ,'.','MarkerEdgeColor', [0    0.4470    0.7410])
    plot3(KYZZr(9:end),CCZZr(9:end),KXZZr(9:end),'k','LineWidth',3,'LineSmoothing','on')
    hold off
    xlim([0.28, 1.225])
    ylim([0 , 2.11])
    zlim([0, 1])
    view([0 90])
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 28)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 28)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 28)
    colorbar
    colormap(parula(256))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',28)
    set(yaxis,'FontSize',28)



    if flipkyc
       figure(9336)
         box on
    set(gcf, 'Position', [100, 100, 900, 600]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',fa)
 shading interp
    hold on
    plot3(-1,-1,-1)
        plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',40)
    plot3(KYZZ_sn(15),CCZZ_sn(15),KXZZ_sn(15),'g.','MarkerSize',40)
    surf(KYGo(1:6,:),CGo(1:6,:),KXGo(1:6,:),'FaceAlpha',fa)
    shading interp
    surf(KYGo(7:end-5,:),CGo(7:end-5,:),KXGo(7:end-5,:),'FaceAlpha',fa)
        shading interp
    surf([KYGa(1,:);KYGo(end-5,:)],[CGob(1,:);CGo(end-5,:)],[KXGob(1,:);KXGo(end-5,:)],'FaceAlpha',fa)
        shading interp
     shading interp
    surf(KYZZGt(1:15,:),CCZZGt(1:15,:),KXZZGt(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZGb(1:15,:),CCZZGb(1:15,:),KXZZGb(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZG(16:end,:),CCZZG(16:end,:),KXZZG(16:end,:),'FaceAlpha',fa)
        shading interp
    surf([KYGo(6,1:jjj);KYGo(7,JJ)],[CGo(6,1:jjj);CGo(7,JJ)],[KXGo(6,1:jjj);KXGo(7,JJ)],'FaceAlpha',fa)
        shading interp
    surf([KYGo(7,JF);KYZZG(16,LL)],[CGo(7,JF);CCZZG(16,LL)],[KXGo(7,JF);KXZZG(16,LL)],'FaceAlpha',fa)
        shading interp
    surf([KYZZGt(15,MMM);KYZZG(16,1:llf)],[CCZZGt(15,MMM);CCZZG(16,1:llf)],[KXZZGt(15,MMM);KXZZG(16,1:llf)],'FaceAlpha',0.6)
        shading interp
    area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',fa)
    %scatter3(KY,CC,KX,'.')
    plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3,'LineSmoothing','on')
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3,'LineSmoothing','on')
    plot3(KYZZ_sn(1:15),CCZZ_sn(1:15),KXZZ_sn(1:15),'r','LineWidth',3,'LineSmoothing','on')
   %scatter3(KYZZ,CCZZ,KXZZ,'.','MarkerEdgeColor', [0    0.4470    0.7410])
    plot3(KYZZr(9:end),CCZZr(9:end),KXZZr(9:end),'k','LineWidth',3,'LineSmoothing','on')
    hold off
    xlim([0.28, 1.225])
    ylim([0 , 2.11])
    zlim([0, 1])
    view([90 90])
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 28)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 28)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 28)
    colorbar
    colormap(parula(256))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',28)
    set(yaxis,'FontSize',28)







    end


%%With Zig-Zag buuble--Just the bubble part
    fa = 0.9;
    figure(337)
           box on
    set(gcf, 'Position', [100, 100, 900, 600]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',fa)
 shading interp
    hold on
%          surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
%     shading interp
    plot3(-1,-1,-1,'.')
%         surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
%     shading interp
    surf([KYGo(6,1:jjj);KYGo(7,JJ)],[CGo(6,1:jjj);CGo(7,JJ)],[KXGo(6,1:jjj);KXGo(7,JJ)],'FaceAlpha',fa)
    shading interp
        surf([KYGo(7,JF);KYZZG(16,LL)],[CGo(7,JF);CCZZG(16,LL)],[KXGo(7,JF);KXZZG(16,LL)],'FaceAlpha',fa)
    shading interp
    surf(KYGo(1:6,:),CGo(1:6,:),KXGo(1:6,:),'FaceAlpha',fa)
    shading interp
    surf(KYGo(7:end,:),CGo(7:end,:),KXGo(7:end,:),'FaceAlpha',fa)
     shading interp
%       surf(KYZZG(1:15,:),CCZZG(1:15,:),KXZZG(1:15,:),'FaceAlpha',fa)
%       shading interp
    surf(KYZZGt(1:15,:),CCZZGt(1:15,:),KXZZGt(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZGb(1:15,:),CCZZGb(1:15,:),KXZZGb(1:15,:),'FaceAlpha',fa)
      shading interp
     surf(KYZZG(16:end,:),CCZZG(16:end,:),KXZZG(16:end,:),'FaceAlpha',fa)
    shading interp
    surf([KYZZGt(15,MMM);KYZZG(16,1:llf)],[CCZZGt(15,MMM);CCZZG(16,1:llf)],[KXZZGt(15,MMM);KXZZG(16,1:llf)],'FaceAlpha',fa)
        shading interp
      area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',fa)
%       scatter3(KY,CC,KX,'.')
       plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3)
    plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',25)
       plot3(KYZZ_sn(15),CCZZ_sn(15),KXZZ_sn(15),'g.','MarkerSize',25)
    plot3(KYZZ_sn(1:15),CCZZ_sn(1:15),KXZZ_sn(1:15),'r','LineWidth',3)
  %  scatter3(KYZZ,CCZZ,KXZZ,'.','MarkerEdgeColor', [0    0.4470    0.7410])
    plot3(KYZZr,CCZZr,KXZZr,'k','LineWidth',3)
    hold off
    xlim([0.7, 1.225])
    ylim([0 , 1.6])
    zlim([0, 0.8])
     view([54 63])
     colormap(parula(256))
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 28)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 28)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 28)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',30 )
    set(yaxis,'FontSize',30)
    set(zaxis,'FontSize',30)




%%With Zig-Zag buuble--Just the bubble part
    fa = 0.9
    figure(338)
           box on
    set(gcf, 'Position', [100, 100, 900, 600]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',fa)
 shading interp
    hold on
%          surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
%     shading interp
    plot3(-1,-1,-1,'.')
%         surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
%     shading interp
    surf([KYGo(6,1:jjj);KYGo(7,JJ)],[CGo(6,1:jjj);CGo(7,JJ)],[KXGo(6,1:jjj);KXGo(7,JJ)],'FaceAlpha',fa)
    shading interp
        surf([KYGo(7,JF);KYZZG(16,LL)],[CGo(7,JF);CCZZG(16,LL)],[KXGo(7,JF);KXZZG(16,LL)],'FaceAlpha',fa)
    shading interp
    surf(KYGo(1:6,:),CGo(1:6,:),KXGo(1:6,:),'FaceAlpha',fa)
    shading interp
    surf(KYGo(7:end,:),CGo(7:end,:),KXGo(7:end,:),'FaceAlpha',fa)
     shading interp
%       surf(KYZZG(1:15,:),CCZZG(1:15,:),KXZZG(1:15,:),'FaceAlpha',fa)
%       shading interp
    surf(KYZZGt(1:15,:),CCZZGt(1:15,:),KXZZGt(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZGb(1:15,:),CCZZGb(1:15,:),KXZZGb(1:15,:),'FaceAlpha',fa)
      shading interp
     surf(KYZZG(16:end,:),CCZZG(16:end,:),KXZZG(16:end,:),'FaceAlpha',fa)
    shading interp
    surf([KYZZGt(15,MMM);KYZZG(16,1:llf)],[CCZZGt(15,MMM);CCZZG(16,1:llf)],[KXZZGt(15,MMM);KXZZG(16,1:llf)],'FaceAlpha',fa)
        shading interp
      area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',fa)
%       scatter3(KY,CC,KX,'.')
       plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3)
    plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',25)
       plot3(KYZZ_sn(15),CCZZ_sn(15),KXZZ_sn(15),'g.','MarkerSize',25)
    plot3(KYZZ_sn(1:15),CCZZ_sn(1:15),KXZZ_sn(1:15),'r','LineWidth',3)
  %  scatter3(KYZZ,CCZZ,KXZZ,'.','MarkerEdgeColor', [0    0.4470    0.7410])
    plot3(KYZZr,CCZZr,KXZZr,'k','LineWidth',3)
    hold off
    xlim([0.7, 1.225])
    ylim([0 , 1.6])
    zlim([0, 0.8])
     view([112 60])
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 30)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 30)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 30)
colormap(parula(256))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',30)
    set(yaxis,'FontSize',30)
    set(zaxis,'FontSize',30)


%%With Zig-Zag buuble--Just the bubble part
    fa = 0.9
    figure(339)
           box on
    set(gcf, 'Position', [100, 100, 900, 600]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',fa)
 shading interp
    hold on
%          surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
%     shading interp
    plot3(-1,-1,-1,'.')
%         surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
%     shading interp
    surf([KYGo(6,1:jjj);KYGo(7,JJ)],[CGo(6,1:jjj);CGo(7,JJ)],[KXGo(6,1:jjj);KXGo(7,JJ)],'FaceAlpha',fa)
    shading interp
        surf([KYGo(7,JF);KYZZG(16,LL)],[CGo(7,JF);CCZZG(16,LL)],[KXGo(7,JF);KXZZG(16,LL)],'FaceAlpha',fa)
    shading interp
    surf(KYGo(1:6,:),CGo(1:6,:),KXGo(1:6,:),'FaceAlpha',fa)
    shading interp
    surf(KYGo(7:end,:),CGo(7:end,:),KXGo(7:end,:),'FaceAlpha',fa)
     shading interp
%       surf(KYZZG(1:15,:),CCZZG(1:15,:),KXZZG(1:15,:),'FaceAlpha',fa)
%       shading interp
    surf(KYZZGt(1:15,:),CCZZGt(1:15,:),KXZZGt(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZGb(1:15,:),CCZZGb(1:15,:),KXZZGb(1:15,:),'FaceAlpha',fa)
      shading interp
     surf(KYZZG(16:end,:),CCZZG(16:end,:),KXZZG(16:end,:),'FaceAlpha',fa)
    shading interp
    surf([KYZZGt(15,MMM);KYZZG(16,1:llf)],[CCZZGt(15,MMM);CCZZG(16,1:llf)],[KXZZGt(15,MMM);KXZZG(16,1:llf)],'FaceAlpha',fa)
        shading interp
      area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',fa)
%       scatter3(KY,CC,KX,'.')
       plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3)
    plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',25)
       plot3(KYZZ_sn(15),CCZZ_sn(15),KXZZ_sn(15),'g.','MarkerSize',25)
    plot3(KYZZ_sn(1:15),CCZZ_sn(1:15),KXZZ_sn(1:15),'r','LineWidth',3)
  %  scatter3(KYZZ,CCZZ,KXZZ,'.','MarkerEdgeColor', [0    0.4470    0.7410])
    plot3(KYZZr,CCZZr,KXZZr,'k','LineWidth',3)
    hold off
    xlim([0.7, 1.225])
    ylim([0 , 1.6])
    zlim([0, 0.8])
     view([-45 36])
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 30)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 30)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 30)
colormap(parula(256))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',30)
    set(yaxis,'FontSize',30)
        set(zaxis,'FontSize',30)



       fa =0.7;
    figure(400)
              box on
        set(gcf, 'Position', [100, 100, 900, 600]);
        set(gcf,'PaperPositionMode','auto')
        surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',0.6)
        shading interp
        hold on
        %        surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
        %   shading interp
        surf(KYGo(1:6,:),CGo(1:6,:),KXGo(1:6,:),'FaceAlpha',0.6)
        shading interp
    surf(KYGo(7:end-5,:),CGo(7:end-5,:),KXGo(7:end-5,:),'FaceAlpha',fa)
        shading interp
    surf([KYGa(1,:);KYGo(end-5,:)],[CGob(1,:);CGo(end-5,:)],[KXGob(1,:);KXGo(end-5,:)],'FaceAlpha',fa)
        shading interp
        surf(KYZZGt(1:15,:),CCZZGt(1:15,:),KXZZGt(1:15,:),'FaceAlpha',0.6)
        shading interp
        surf(KYZZGb(1:15,:),CCZZGb(1:15,:),KXZZGb(1:15,:),'FaceAlpha',0.6)
        shading interp
        surf(KYZZG(16:end,:),CCZZG(16:end,:),KXZZG(16:end,:),'FaceAlpha',0.6)
        shading interp
        surf([KYGo(6,1:jjj);KYGo(7,JJ)],[CGo(6,1:jjj);CGo(7,JJ)],[KXGo(6,1:jjj);KXGo(7,JJ)],'FaceAlpha',0.6)
        shading interp
        surf([KYGo(7,JF);KYZZG(16,LL)],[CGo(7,JF);CCZZG(16,LL)],[KXGo(7,JF);KXZZG(16,LL)],'FaceAlpha',0.6)
        shading interp
        surf([KYZZGt(15,MMM);KYZZG(16,1:llf)],[CCZZGt(15,MMM);CCZZG(16,1:llf)],[KXZZGt(15,MMM);KXZZG(16,1:llf)],'FaceAlpha',0.6)
        shading interp
        %surf([KYGo(7,:);KYZZG(16,:)],[CGo(7,:);CCZZG(16,:)],[KXGo(7,:);KXZZG(16,:)],'FaceAlpha',0.6)
        %shading interp

        area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',0.6)
        scatter3(KY,CC,KX,'.')
        plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3)
        plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3)
        plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',40)
        plot3(KYZZ_sn(15),CCZZ_sn(15),KXZZ_sn(15),'g.','MarkerSize',40)
        plot3(KYZZ_sn(1:15),CCZZ_sn(1:15),KXZZ_sn(1:15),'r','LineWidth',3)
        scatter3(KYZZ,CCZZ,KXZZ,'.','MarkerEdgeColor', [0    0.4470    0.7410])
        plot3(KYZZr(9:end),CCZZr(9:end),KXZZr(9:end),'k','LineWidth',3)
        hold off
        xlim([0.25, 1.225])
        ylim([0 , 2.1])
        zlim([0, 1])
        view([54 63])
        %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
        xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 28)
        yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 28)
        zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 28)
    colormap(parula(256))
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16);
        set(xaxis,'FontSize',28)
        set(yaxis,'FontSize',28)
        set(zaxis,'FontSize',28)



  %%%%%make an animation

  if makevideo

  v = VideoWriter('moduli.avi');
  v.Quality = 100;
open(v);

  for jj = 1:360
              fa =0.75;
    fig = figure(401)
              box on
        set(gcf, 'Position', [100, 100, 900, 600]);
        set(gcf,'PaperPositionMode','auto')
        surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',fa)
        shading interp
        hold on
        %        surf(KYGo,CGo,KXGo,'FaceAlpha',0.6)
        %   shading interp
        surf(KYGo(1:6,:),CGo(1:6,:),KXGo(1:6,:),'FaceAlpha',fa)
        shading interp
    surf(KYGo(7:end-5,:),CGo(7:end-5,:),KXGo(7:end-5,:),'FaceAlpha',fa)
        shading interp
    surf([KYGa(1,:);KYGo(end-5,:)],[CGob(1,:);CGo(end-5,:)],[KXGob(1,:);KXGo(end-5,:)],'FaceAlpha',fa)
        shading interp
        surf(KYZZGt(1:15,:),CCZZGt(1:15,:),KXZZGt(1:15,:),'FaceAlpha',fa)
        shading interp
        surf(KYZZGb(1:15,:),CCZZGb(1:15,:),KXZZGb(1:15,:),'FaceAlpha',fa)
        shading interp
        surf(KYZZG(16:end,:),CCZZG(16:end,:),KXZZG(16:end,:),'FaceAlpha',fa)
        shading interp
        surf([KYGo(6,1:jjj);KYGo(7,JJ)],[CGo(6,1:jjj);CGo(7,JJ)],[KXGo(6,1:jjj);KXGo(7,JJ)],'FaceAlpha',fa)
        shading interp
        surf([KYGo(7,JF);KYZZG(16,LL)],[CGo(7,JF);CCZZG(16,LL)],[KXGo(7,JF);KXZZG(16,LL)],'FaceAlpha',fa)
        shading interp
        surf([KYZZGt(15,MMM);KYZZG(16,1:llf)],[CCZZGt(15,MMM);CCZZG(16,1:llf)],[KXZZGt(15,MMM);KXZZG(16,1:llf)],'FaceAlpha',fa)
        shading interp
        %surf([KYGo(7,:);KYZZG(16,:)],[CGo(7,:);CCZZG(16,:)],[KXGo(7,:);KXZZG(16,:)],'FaceAlpha',0.6)
        %shading interp

        area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',fa)
        scatter3(KY,CC,KX,'.')
        plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3)
        plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3)
        plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',40)
        plot3(KYZZ_sn(15),CCZZ_sn(15),KXZZ_sn(15),'g.','MarkerSize',40)
        plot3(KYZZ_sn(1:15),CCZZ_sn(1:15),KXZZ_sn(1:15),'r','LineWidth',3)
        scatter3(KYZZ,CCZZ,KXZZ,'.','MarkerEdgeColor', [0    0.4470    0.7410])
        plot3(KYZZr(9:end),CCZZr(9:end),KXZZr(9:end),'k','LineWidth',3)
        hold off
        xlim([0.78, 1.225])
        ylim([0 , 1.5])
        zlim([0, 0.7])
        view([jj+30 30])
        lighting gouraud
        %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
        xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
        yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
        zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
        colormap(parula(256))
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16);
        set(xaxis,'FontSize',20)
        set(yaxis,'FontSize',20)
        frame = getframe(fig);
        writeVideo(v,frame);
  end

close(v)

end

% filename = 'testAnimated.gif'; % Specify the output file name
% for idx = 1:length(im)
%     [A,map] = rgb2ind(im{idx},256);
%     if idx == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.3);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.3);
%     end
% end




  h=figure(100)
    box on
    set(gcf, 'Position', [100, 100, 900, 450]);
    set(gcf,'PaperPositionMode','auto')
    scatter3(KY,CC,KX,'.')
    hold on
    %plot3(KYI(1:iic-1),C0(1:iic-1),KYI(1:iic-1)*0,'LineWidth',2)
    plot3(KYI(1:iic-1),C1(1:iic-1),KX1(1:iic-1),'LineWidth',2)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',2)
    hold off
    legend(['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)

    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)


    drawnow


    %%%Find and plot reattachment curve (ky,c)
    %%%Find and plot fold curve in (ky,c)
    figure(32)
plot([KY_I';KY_Ia'],[C_sn'],'x','LineWidth',2)
hold on
plot(KY_I(2:34),C_pf(2:34),'x','LineWidth',2)
plot([KY_Ia],C_pf(end-9:end),'x','LineWidth',2)
hold off

%%%Load NWS data from Arnd
% load location of pitchfork
load('pf_data.mat')
%  kygrid=[kygrid1';kygrid2'];
%
B=sortrows([kygrid0,ccrit0,cbranch0]);
kygrid0=B(:,1);
ccrit0=B(:,2);
cbranch0=(B(:,3)+1)/20;

h=figure(100)
box on
set(gcf, 'Position', [100, 100, 1800, 450]);
set(gcf,'PaperPositionMode','auto')
hold on
box on
plot(kygrid0,ccrit0,'-m', 'LineWidth',2)

load('sn_cont_4th_0.25.mat')
kygrid=[kygrid1';kygrid2'];

B=sortrows([kygrid,ccrit]);
kygridB=B(:,1);
ccritB=B(:,2);
plot(kygrid,ccrit,'-r','linewidth',2)


hold on
load('sn_cont_small_ky.mat')
kygrid=kys;

B=sortrows([ccrit/2,kygrid]);
B=B(1:170,:);
kygrid=B(:,2);
ccrit=B(:,1);
plot(kygrid,ccrit,'-r','linewidth',2)
plot(sqrt(4 + sqrt(3))/2,0.438691,'k.','MarkerSize',20)

hold off
axis([0.75 1.25 -0.0001 0.7])
xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)

xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)





h=figure(35)
box on
set(gcf, 'Position', [100, 100, 600, 400]);
set(gcf,'PaperPositionMode','auto')

plot([KY_I';KY_Ia'],[C_sn'],'o','LineWidth',2)
hold on
plot(KYI(1:iic(1)-1),C1(1:iic(1)-1),'LineWidth',2)
plot(KYI(iic),C2(iic),'LineWidth',2)
plot([KY_I(1:34),KY_Ia],[C_pf],'+','LineWidth',2,'Color',[0.,0.4,0.1])
plot(sqrt(4 + sqrt(3))/2,0.438691,'k.','MarkerSize',40)
plot([kygrid;kygridB],[ccrit;ccritB],'-r','linewidth',2)
plot(kygrid0,ccrit0,'-m','linewidth',2)
xlim([0.77 1.225])
ylim([0, 0.7])
hold off
xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)

xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
% lgd=legend({'Saddle-node curve', '$c_{lin}(k_y)$',...
%     '$c_{lin,0}(k_y), k_y>$','Oblique stripe reattachment','Abs. Spec. triple point','Saddle-node NWS','Pitchfork-NWS'},'Interpreter',...
%        'latex')
% lgd.Location = 'Northwest'



if flipkyc

h=figure(935)
box on
set(gcf, 'Position', [100, 100, 600, 400]);
set(gcf,'PaperPositionMode','auto')

plot([C_sn'],[KY_I';KY_Ia'],'o','LineWidth',2)
hold on
plot(C1(1:iic(1)-1),KYI(1:iic(1)-1),'LineWidth',2)
plot(C2(iic),KYI(iic),'LineWidth',2)
plot([C_pf],[KY_I(1:34),KY_Ia],'+','LineWidth',2,'Color',[0.,0.4,0.1])
plot(0.438691,sqrt(4 + sqrt(3))/2,'k.','MarkerSize',40)
plot([ccrit;ccritB],[kygrid;kygridB],'-r','linewidth',2)
plot(ccrit0,kygrid0,'-m','linewidth',2)
ylim([0.77 1.225])
xlim([0, 0.7])
hold off
yaxis=ylabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
ax = gca;
ax.YDir = 'reverse';
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)

xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)

% lgd=legend({'Saddle-node curve', '$c_{lin}(k_y)$',...
%     '$c_{lin,0}(k_y), k_y>$','Oblique stripe reattachment','Abs. Spec. triple point','Saddle-node NWS','Pitchfork-NWS'},'Interpreter',...
%        'latex')
% lgd.Location = 'Northwest'


end


  v = VideoWriter('test.avi');
  v.Quality = 100;
open(v);

if makesolutionvideo

    %%With Zig-Zag buuble-different view, just interpolated surface
    fa =1;
   fig =  figure(501)
         box on
    set(gcf, 'Position', [100, 100, 900, 600]);
    set(gcf,'PaperPositionMode','auto')
 surf(KYGa(1:skipcounter,:),CGob,KXGob,'FaceAlpha',fa)
 shading interp
    hold on
    plot3(-1,-1,-1)
        plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',40)
    plot3(KYZZ_sn(15),CCZZ_sn(15),KXZZ_sn(15),'g.','MarkerSize',40)
    surf(KYGo(1:6,:),CGo(1:6,:),KXGo(1:6,:),'FaceAlpha',fa)
    shading interp
    surf(KYGo(7:end-5,:),CGo(7:end-5,:),KXGo(7:end-5,:),'FaceAlpha',fa)
        shading interp
    surf([KYGa(1,:);KYGo(end-5,:)],[CGob(1,:);CGo(end-5,:)],[KXGob(1,:);KXGo(end-5,:)],'FaceAlpha',fa)
        shading interp
     shading interp
    surf(KYZZGt(1:15,:),CCZZGt(1:15,:),KXZZGt(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZGb(1:15,:),CCZZGb(1:15,:),KXZZGb(1:15,:),'FaceAlpha',fa)
      shading interp
    surf(KYZZG(16:end,:),CCZZG(16:end,:),KXZZG(16:end,:),'FaceAlpha',fa)
        shading interp
    surf([KYGo(6,1:jjj);KYGo(7,JJ)],[CGo(6,1:jjj);CGo(7,JJ)],[KXGo(6,1:jjj);KXGo(7,JJ)],'FaceAlpha',fa)
        shading interp
    surf([KYGo(7,JF);KYZZG(16,LL)],[CGo(7,JF);CCZZG(16,LL)],[KXGo(7,JF);KXZZG(16,LL)],'FaceAlpha',fa)
        shading interp
    surf([KYZZGt(15,MMM);KYZZG(16,1:llf)],[CCZZGt(15,MMM);CCZZG(16,1:llf)],[KXZZGt(15,MMM);KXZZG(16,1:llf)],'FaceAlpha',0.6)
        shading interp
    area([KY_sn,KY_Ia(1:end-1)],C_sn(1:end-1),'FaceColor',[0 0.0 0.5],'FaceAlpha',fa)
    %scatter3(KY,CC,KX,'.')
    plot3(KYI(1:iic(1)+1),C1(1:iic(1)+1),KX1(1:iic(1)+1),'LineWidth',3,'LineSmoothing','on')
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',3,'LineSmoothing','on')
    plot3(KYZZ_sn(1:15),CCZZ_sn(1:15),KXZZ_sn(1:15),'r','LineWidth',3,'LineSmoothing','on')
   %scatter3(KYZZ,CCZZ,KXZZ,'.','MarkerEdgeColor', [0    0.4470    0.7410])
    plot3(KYZZr(9:end),CCZZr(9:end),KXZZr(9:end),'k','LineWidth',3,'LineSmoothing','on')
    hold off
    xlim([0.28, 1.225])
    ylim([0 , 2.11])
    zlim([0, 1])
    view([0 90])
    %    legend([],[],[],['Data'],[' oblique lin. sp.'],['perp lin. sp.'])%,['C2->Perp lin. sp.'])
    xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 28)
    yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 28)
    zaxis = zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 28)
    colorbar
    colormap(parula(256))
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',28)
    set(yaxis,'FontSize',28)

        frame = getframe(fig);
        writeVideo(v,frame);
     close(v)
end



%%%%%Zig-Zag reattachment profiles
%%%%Plot a specific solution (unfortunately only have the oblique stripe
%%%%data for this set of continuations
jj = 2;
kky = KY_I(jj);
%%%for a perp solution

fname = sprintf('ContData%d.mat',kky)
load(fname)

%%finds negative end point of curve and rearranges data
[dum,iif] = min(PROF(end-1,:));
%kxplot = [PROF(end-1,iif:-1:1),PROF(end-1,iif+1:end)];
%cxplot =  [PROF(end,iif:-1:1),PROF(end,iif+1:end)];
PROFsh = [PROF(:,iif:-1:1),PROF(:,iif+1:end)];


IND_s = [20 175 200 240];
hb = figure(200)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(PROFsh(end,:),PROFsh(end-1,:),'LineWidth',2)
plot(PROFsh(end,IND_s),PROFsh(end-1,IND_s),'o','LineWidth',2)
text(PROFsh(end,IND_s(1))+0.001,PROFsh(end-1,IND_s(1))+0.015,'(1)','Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}$', 'Interpreter', 'latex','FontSize',18)
text(PROFsh(end,IND_s(2))+0.002,PROFsh(end-1,IND_s(2)),'(2)','Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}$', 'Interpreter', 'latex','FontSize',18)
text(PROFsh(end,IND_s(3))-0.005,PROFsh(end-1,IND_s(3))-0.01,'(3)','Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}$', 'Interpreter', 'latex','FontSize',18)
text(PROFsh(end,IND_s(4)),PROFsh(end-1,IND_s(4))-0.02,'(4)','Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {4}}}$', 'Interpreter', 'latex','FontSize',18)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
hold off
ylim([0, max(PROFsh(end-1,:))])
kky = ky;

hp=figure(201)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii = 1:length(IND_s)
kkx = PROFsh(end-1,IND_s(ii)+1);%KDAT(Ind_PROFsh(ii),end-1);
cc = PROFsh(end,IND_s(ii)+1);%KDAT(Ind_PROFsh(ii),end);
znew = PROFsh(:,IND_s(ii)+1);
solu = sprof(x,y,kky,znew);
%%%Plot some solutions
 subplot(4,1,ii)
    imagesc(x,y,solu)

    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);

    set(yaxis,'FontSize',20)
    title(['(',num2str(ii),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])
end
hold off
 xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
 set(xaxis,'FontSize',20)
drawnow
