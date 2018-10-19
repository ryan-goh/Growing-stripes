%%%%Plotting figures for reattachment curves
global mu0
%%loads the ky continuation (with c=0) which the bubble emanates from, then
%%clears other data

load('moduli_dat_full.mat')
%dx = 0.25; dy = 0.25; Lx = 40*pi;
I = length(moduli);
Io = I;
KYo = KY;
KXo = KX;
CCo = CC;
Co = C;

flipkyc = 1;

%%loads(''moduli','kdat','KX','KY','CC','mu0','C','I','l','dx','dy','Lx','Ly','x','y')
%KY_I = linspace(1/sqrt(2)+0.08,1.22,I);

% l = [];
% for jj = 1:I
% l = [l;length(moduli{jj}(:,end))];
% end

%     KX = [];%zeros(sum(l),1);
%     KY = KX;
%     CC = KX;
% for jj = 1:I
% KX = [KX;moduli{jj}(:,end-1)];
% CC = [CC; moduli{jj}(:,end)];
% KY = [KY;KY_I(jj)*ones(l(jj),1)];
% end
    KMo = sqrt(KX.^2 + KY.^2);


    CGo = [];
    KYGo = [];
    KXGo = [];

    %KY_sn = KY_I;
    IDo = I-3;
    JD = 1000;
    CGo = C(1:IDo)'*ones(1,JD);
    %KYG = KY_I(1:ID)'*ones(1,JD);

    KXmax = [];
    KYmax = [];
    Cmax = Co(1:IDo);

    for ii = 1:I-3
        kyd = moduli{ii}(end,:);
        kxd = moduli{ii}(end-1,:);
        md = length(kyd);
        for jj = 1:md
            if abs(kxd(jj))>1e-5;%%end of perp stripes
                jd = jj;
                break
            end
        end

        [kym,iim] = max(kyd);
        kxm = kxd(iim);
        KXmax = [KXmax;kxm];
        KYmax = [KYmax;kym];


        %%%%
        %%%
        tt = linspace(0,2*pi, length(kyd));
        Y = [kyd;kxd];
        sl0 = [1,0];
        sle = [kyd(end) - kyd(end-1);kxd(end) - kxd(end-1)];
        pp = spline(tt,Y);
        yy = ppval(pp,linspace(0,2*pi,JD));
%         figure(1)
%         plot(yy(1,:),yy(2,:),'.')
%
%         figure(3)
%         plot(kyd,kxd,'.')

        KYGo = [KYGo;yy(1,:)];
        KXGo = [KXGo;yy(2,:)];

    end


% figure(2)
% scatter3(KY,CC,KX,'.')
% hold on
% surf(KYGo,CGo,KXGo,'FaceAlpha',0.75)
% shading interp
% xlabel('k_y')
% ylabel('c')
% zlabel('k_x')
%  %   xlim([min(KY_I) max(KY_I)])
%  %   ylim([0 max(CC)*(1.2)])
%  %   zlim([0 max(KX)])
%     drawnow
%
% figure(22)
% scatter3(KYo,CCo,KXo,'.')



%%%Plot a specific solution
%ii = 100;
%jj = 3;
%Nx= floor(Lx/dx)+1;
%Ny = floor(Ly/dy)+1;
%dx = Lx/Nx;
%dy = Ly/Ny;
%y = dy*[1:Ny]';
%x =[dx:dx:Lx]';
%kkx = kdat{jj}(ii,end-1);
%kky = kdat{jj}();
%znew = moduli{jj}(:,ii);
%cc = C(jj);
%ufull = sprof(x,y,cc,znew);%full{jj}(:,Nx*ii+1:Nx*(ii+1));
%figure(5)
%imagesc(x,y,ufull)
%xlabel('x')
%ylabel('y')
%title([' c = ', num2str(cc),' ky = ', num2str(KY_I(jj)), ', kx = ', num2str(kkx)])



%%%%Run through one continuation run (can put a movie wrapper around this
%%%%loop to record, also plots the corresponding point on a cross-section
%%%%of the moduli surface
%jj = 10;
%for ii = 100;%1:3:length(moduli{jj}(1,:))
%    %%Plot a specific solution
%     Nx= floor(Lx/dx)+1;
%     Ny = floor(Ly/dy)+1;
%     y = dy*[1:Ny]';
%     x =[dx:dx:Lx]';
%    KXdum = moduli{jj}(end-1,:);
%    KYdum = moduli{jj}(end,:);
%    kkx =KXdum(ii);
%    cc = C(jj);
%    kky = KYdum(jj);
%    %%function to output
%    w = moduli{jj}(:,ii);
%    ufull = sprof(x,y,cc,w);

   % ufull = full{jj}(:,Nx*ii+1:Nx*(ii+1));
%    figure(3)
%    imagesc(x,y,ufull)
%    xlabel('x')
%    ylabel('y')
%    title([' c = ', num2str(cc),' ky = ', num2str(kky), ', kx = ', num2str(kkx)])


%    figure(4)
%    plot(KYdum,KXdum,'-',kky,kkx,'x')
%    xlabel('c')
%    ylabel('k_x')
%    title(['ky = ', num2str(kky)])
%    drawnow
%end








%%%Now load up the perp stripes
load('moduli_dat_full_perp.mat','kdat','moduli','I','KX','KY','CC','C')
Ip = I;
Cp = C;
KXp = KX;
KYp = KY;
CCp = CC;
KYpmax = [];
KYpmin = [];
Cpmax = Cp;
    CGp = [];
    KYGp = [];
    KXGp = [];

    %KY_sn = KY_I;
    IDp = Ip;
    JDp = 1000;
    CGp = Cp(1:IDp)'*ones(1,JDp);
    %KYG = KY_I(1:ID)'*ones(1,JD);



 for ii = 1:I
        kyd = kdat{ii}(:,end);
        kxd = kdat{ii}(:,end-1);
        md = length(kyd);
        for jj = 1:md
            if abs(kxd(jj))>1e-5;%%end of perp stripes
                jd = jj;
                break
            end
        end



%        [kym,iim] = max(kyd);
%        kxm = kxm(iim);
%        KXmax = [KXmax;kxm];
	KYpmin = [KYpmin;min(kyd)];
        KYpmax = [KYpmax;max(kyd)];


        %%%%
        %%%
        tt = linspace(0,2*pi, length(kyd));
        Y = [kyd';kxd'];
        sl0 = [1,0];
        sle = [kyd(end) - kyd(end-1);kxd(end) - kxd(end-1)];
        pp = spline(tt,Y);
        yy = ppval(pp,linspace(0,2*pi,JD));
%         figure(1)
%         plot(yy(1,:),yy(2,:),'.')
%
%         figure(3)
%         plot(kyd,kxd,'.')

        KYGp = [KYGp;yy(1,:)];
        KXGp = [KXGp;yy(2,:)];

    end

KXf = [KXp;KXo];
KYf = [KYp;KYo];
CCf = [CCp;CCo];


% %%%%%%%%%%
% %%%%%%%%%Now load up a few more perp stripes
% %%%Now load up the perp stripes,,,,
% %load('moduli_dat_full_perpadd.mat','kdat','moduli','I','KX','KY','CC','C')
% %Ipa = I;
% %Cp = C;
%
%
% KXpa = [];
% KYpa = [];
% CCpa = [];
%
%     CGpa = [];
%     KYGpa = [];
%     KXGpa = [];
%
% Cpa = [0.455405405405405   0.482432432432432   0.509459459459459 0.53 0.54]; %%Don't delete these
% Ipa = length(Cpa);
% Cpmax = [Cpmax,Cpa];
%
%
%     %KY_sn = KY_I;
%     IDpa = Ipa;
%     JDpa = 1000;
%     CGpa = Cpa(1:IDpa)'*ones(1,JDpa);
%
%
% for ii =1:Ipa
% fname = sprintf('ContDataPerp%d.mat',Cpa(ii));
% load(fname,'PROF');
% kyd = PROF(end,:);
% kxd = PROF(end-1,:);
%
%
%
% iijd = find(abs(kxd)<1e-6);
% CCpa = [CCpa;Cpa(ii)*ones(length(kxd(iijd)),1)];
% KXpa = [KXpa;kxd(iijd)'];
% KYpa = [KYpa;kyd(iijd)'];
% KYpmax = [KYpmax;max(kyd(iijd))];
% KYpmin = [KYpmin;min(kyd(iijd))];
% %%%%
%         %%%
%         tt = linspace(0,2*pi, length(kyd(iijd)));
%         Y = [kyd(iijd);kxd(iijd)];
%         sl0 = [1,0];
%         sle = [kyd(end) - kyd(end-1);kxd(end) - kxd(end-1)];
%         pp = spline(tt,Y);
%         yy = ppval(pp,linspace(0,2*pi,JDpa));
%         figure(1)
% %         plot(yy(1,:),yy(2,:),'.')
% %
% %         figure(3)
% %         plot(kyd,kxd,'.')
%
%         KYGpa = [KYGpa;yy(1,:)];
%         KXGpa = [KXGpa;yy(2,:)];
%
% end
%
%
%
% KXf = [KXf;KXpa];
% KYf = [KYf;KYpa];
% CCf = [CCf;CCpa];



  %%linear predictions for all-stripe detachment speed, perpendicular for ky\gtrsim1.19, oblique for ky\lesssim1.19
    KYI = linspace(0,1.23,3000);
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
    [iic] = find(abs(C2)==real(C2))


%%Plot full scatter plot data of moduli space, with surface interpolation
figure(32)
scatter3(KYf,CCf,KXf,'.')
hold on
surf([KYGp,[CGp,[KXGp,'FaceAlpha',0.75)
shading interp
%surf(KYGpa,CGpa,KXGpa,'FaceAlpha',0.75)
%shading interp
 xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
  yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    zaxis=zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
 %   xlim([min(KY_I) max(KY_I)])
    ylim([0.05 2])
 %   zlim([0 max(KX)])
     xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    drawnow






%only scatter plot data
figure(22)
scatter3(KY,CC,KX,'.')


%%Scatter plot data with linear predictions for detachment
        figure(31)
    scatter3(KYf,CCf,KXf,'.')
    hold on
    %plot3(KYI,C0,KYI*0,'LineWidth',2)
    plot3(KYI(1:iic(1)),C1(1:iic(1)),KX1(1:iic(1)),'LineWidth',2)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',2)
    hold off
    xlabel('ky')
    ylabel('c')
    ylim([0.05 2])
    zlabel('kx')
    %legend(['Data'],[ 'C1 -> obl. lin. sp.'],['C2 -> perp. lin. sp.'])%,['C2->Perp lin. sp.'])
    view([45 40])
    drawnow



    %%%Find and plot full detachment curve
    %%%Find and plot fold curve in (ky,c)
    figure(302)
    box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
plot(KYmax(13:end),Cmax(13:end),'x','LineWidth',2)
hold on
plot(KYI(1:end-50),C1(1:end-50),'LineWidth',2)
plot(KYI(iic),C2(iic),'LineWidth',2)
plot(KYpmax(1:end-3),Cpmax(1:end-3),'x','LineWidth',2)
ylim([0.05 2.11])
hold off
 xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
 yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)


    if flipkyc  %%flipkyc axes

    figure(9302)
    box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
plot(Cmax(13:end),KYmax(13:end),'x','LineWidth',2)
hold on
plot(C1(1:end-50),KYI(1:end-50),'LineWidth',2)
plot(C2(iic),KYI(iic),'LineWidth',2)
plot(Cpmax(1:end-3),KYpmax(1:end-3),'x','LineWidth',2)
xlim([0.05 2.11])
hold off
 yaxis=ylabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
 xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    ax = gca;
    ax.YDir = 'reverse'
    end






%%Plots scatter data with interpolation for oblique stripe surface, and filled in perpendicular stripe surface
figure(23)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
hold on
scatter3(KYf,CCf,KXf,'.')
surf(KYGo,CGo,KXGo,'FaceAlpha',0.75)
shading interp
%surf(KYGp,CGp,KXGp,'FaceAlpha',0.75)
%shading interp
%surf(KYGpa,CGpa,KXGpa,'FaceAlpha',0.75)
%shading interp
%surf([KYGp;KYGpa],[CGp;CGpa],[KXGp;KXGpa],'FaceAlpha',0.75)
%shading interp
area([KYpmin(3:end);KYpmax(end:-1:3)],[Cpmax(3:end),Cpmax(end:-1:3)]','FaceAlpha',0.5,'FaceColor',[0 0 0.1])
 xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
 yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    zaxis=zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    set(zaxis,'FontSize',20)
   xlim([0 1.23])
   ylim([0.05 2.1])
 %   zlim([0 max(KX)])
 colormap(parula(256))
hold off
view([39 66])
    drawnow




    if flipkyc

       figure(923)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
hold on
scatter3(KYf,CCf,KXf,'.')
surf(KYGo,CGo,KXGo,'FaceAlpha',0.75)
shading interp
%surf(KYGp,CGp,KXGp,'FaceAlpha',0.75)
%shading interp
%surf(KYGpa,CGpa,KXGpa,'FaceAlpha',0.75)
%shading interp
%surf([KYGp;KYGpa],[CGp;CGpa],[KXGp;KXGpa],'FaceAlpha',0.75)
%shading interp
area([KYpmin(3:end);KYpmax(end:-1:3)],[Cpmax(3:end),Cpmax(end:-1:3)]','FaceAlpha',0.5,'FaceColor',[0 0 0.1])
 xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
 yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    zaxis=zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    set(zaxis,'FontSize',20)
   xlim([0 1.23])
   ylim([0.05 2.1])
 %   zlim([0 max(KX)])
 colormap(parula(256))
hold off
view([125 69])
    drawnow


    end

%%Top down view of interpolated surfaces, with no scatter plot data.
figure(21)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
hold on
surf(KYGo,CGo,KXGo,'FaceAlpha',0.9)
shading interp
%surf([KYGp;KYGpa],[CGp;CGpa],[KXGp;KXGpa],'FaceAlpha',0.75)
%shading interp
area([KYpmin(3:end);KYpmax(end:-1:3)],[Cpmax(3:end),Cpmax(end:-1:3)]','FaceAlpha',0.5,'FaceColor',[0 0 0.1])
plot3(-1,0,0)
plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',30)
plot3(KYI(1:end-50),C1(1:end-50),KX1(1:end-50),'LineWidth',2)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',2)
 xlabel('k_y')
ylabel('c')
zlabel('k_x')
   xlim([0 1.25])
    ylim([0.05 1.92])
 %   zlim([0 max(KX)])
hold off
 xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
 yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    zaxis=zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    colormap(parula(256))
    view([0 90])
    drawnow

if flipkyc

   figure(921)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
hold on
surf(KYGo,CGo,KXGo,'FaceAlpha',0.9)
shading interp
%surf([KYGp;KYGpa],[CGp;CGpa],[KXGp;KXGpa],'FaceAlpha',0.75)
%shading interp
area([KYpmin(3:end);KYpmax(end:-1:3)],[Cpmax(3:end),Cpmax(end:-1:3)]','FaceAlpha',0.5,'FaceColor',[0 0 0.1])
plot3(-1,0,0)
plot3(sqrt(4 + sqrt(3))/2,0.438691,0,'k.','MarkerSize',30)
plot3(KYI(1:end-50),C1(1:end-50),KX1(1:end-50),'LineWidth',2)
    plot3(KYI(iic),C2(iic),KYI(iic)*0,'LineWidth',2)
 xlabel('k_y')
ylabel('c')
zlabel('k_x')
   xlim([0 1.25])
    ylim([0.05 1.92])
 %   zlim([0 max(KX)])
hold off
 xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
 yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    zaxis=zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    colormap(parula(256))
    view([90 90])
    drawnow




end
