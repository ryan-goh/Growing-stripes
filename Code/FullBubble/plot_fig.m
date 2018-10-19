%%%%Plotting figures for zig-zag bubble:
global mu0
%%loads the ky continuation (with c=0) which the bubble emanates from, then
%%clears other data
load('ContData0.mat','PROF')
KX_II = PROF(end-1,:);
KY_II = PROF(end,:);
clear PROF




load('moduli_dat_zz.mat')
%%loads dx,dy,Lx,Ly, mu0
%% CC,KX,KY -> gives the moduli data (combined vectors give triplets of points)
%% C_sn, KX_sn, KY_sn -> the saddle node locations for each point
%% kdat is the moduli data, along with an x-slice of the solution
%% full is stores the "full" solution, core + far-field w + u_p
%% moduli gives the continuation solution: i.e. core solution and k_x, and c:  (w,k_x,c)
%% I is the number of cells in the data storage variables (i.e. I = length(full))
%%KY_I stores the ky values for each continuation run (i.e. KY_I(ii) for
%%moduli{ii})

%%Possible other approach: interpolate a curve for each slice onto a
%%uniform grid, then use surf

ID = 500;      
JD = 9;
CG = [];
KYG = KY_I(1:JD)'*ones(1,ID);
KXG = [];
C_sn1 =[];
KX_sn1 =[];
KY_sn1 = KY_I(1:JD);
%%%%
for iii = 1:JD
    cd = moduli{iii}(end,:);
    [id] = find(cd<0);
    cd(id) = [];
    kxd = moduli{iii}(end-1,:);
    kxd(id) = [];

    %%Find Saddle-node point
    [csn,isn] = max(cd);
    kxsn = kxd(isn);
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
    
    CG = [CG;yy(1,:)];
    KXG = [KXG; yy(2,:)];
end

load('Bubbleplots.mat')
%Preamble
[~,m] = size(eta);
eta = eta(1,:);
C_ch = Cx(1,:);
kzz = 1;
eps_ch = (0:0.0001:0.017)';

%Build Surface
kx_ch = sqrt(2)*kzz*(eps_ch.^(1/2))*eta;
ky_ch = kzz-eps_ch*ones(1,m);
c_ch = 8*(eps_ch.^(3/2))*C_ch;




h=figure(1)
box on
set(gcf, 'Position', [400, 400, 700, 450]);
set(gcf,'PaperPositionMode','auto')
scatter3(KY,CC,KX,'.','LineWidth',2)
hold on
surf(KYG,CG,KXG,'FaceAlpha',0.75)
shading interp
%plot3(KY_sn(1:JD),C_sn(1:JD),KX_sn(1:JD),'g','LineWidth',2)  %%saddle node location
plot3(KY_sn1(1:JD),C_sn1(1:JD),KX_sn1(1:JD),'g','LineWidth',2)  %%saddle node location
plot3(KY_II,KY_II*0,KX_II,'b','LineWidth',2)  %%ky continuation, c=0
surf(ky_ch,c_ch,kx_ch,'FaceAlpha',0.5,'FaceColor',[0.5,0.2,0],'EdgeColor','None')  %%CH-prediction of surface
plot3((kzz-eps_ch)',8*(C_ch(critindex(1))*eps_ch.^(3/2))',sqrt(2)*(eps_ch.^(1/2))'*eta(critindex(1)),'r','LineWidth',2)%%CH-prediction of SN
plot3((kzz-eps_ch)',0*c_ch,sqrt(2)*eps_ch.^(1/2),'k','LineWidth',2)
hold off
xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
zaxis=zlabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
    xlim([min(KY_I) max(KY_I)])
    ylim([0 max(CC)*(1.2)])
    zlim([0 max(KX)])
    view([1 1 0.5])
    drawnow

    

    

    


%%Plot moduli-curves along with the c = 0 curve, and saddle node curve
figure(2)
    scatter3(KY,CC,KX,'.')
    hold on
    plot3(KY_sn1,C_sn1,KX_sn1,'LineWidth',1.5)
    plot3(KY_II,KY_II*0,KX_II,'LineWidth',1.5)
    hold off
    xlabel('ky')
    ylabel('c')
    zlabel('k_x')
    xlim([min(KY_I) max(KY_I)])
    ylim([0 max(CC)])
    zlim([0 max(KX)])
    drawnow
    
    
%%Plot a the saddle node curves projected onto C,KY space  
figure(3)
plot(KY_sn1,C_sn1,'g','LineWidth',2)
hold on
plot3((kzz-eps_ch)',8*(C_ch(critindex(1))*eps_ch.^(3/2))',sqrt(2)*(eps_ch.^(1/2))'*eta(critindex(1)),'r','LineWidth',2)%%CH-prediction of SN
hold off
xlabel('k_y')
ylabel('c')

figure(4)
plot(C_sn1,KX_sn1,'g','LineWidth',2)
hold on
plot(8*(C_ch(critindex(1))*eps_ch.^(3/2))',sqrt(2)*(eps_ch.^(1/2))'*eta(critindex(1)),'r','LineWidth',2)%%CH-prediction of SN
hold off
xlabel('c')
ylabel('k_x')

    
    
    
    
    
%%Plot a series of solutions
%pick ky slice:
jj = 5; 
kky=KY_I(jj)
IND_s = [50, 505, 700 895]; %index of point along continuation curve for plotted solutions 
Nx= floor(Lx/dx)+1;
Ny = floor(Ly/dy)+1;
dxx = Lx/Nx;
dyy = Ly/Ny;
y = dyy*[1:Ny]';
x =[dxx:dxx:Lx]';

h=figure(5)
box on
set(gcf, 'Position', [400, 400, 900, 450]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(moduli{jj}(end,:),moduli{jj}(end-1,:),'LineWidth',2)
plot(moduli{jj}(end,IND_s),moduli{jj}(end-1,IND_s),'o','LineWidth',2)
text(moduli{jj}(end,IND_s(1)),moduli{jj}(end-1,IND_s(1))*1.01,'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}$', 'Interpreter', 'latex','FontSize',18)
text(moduli{jj}(end,IND_s(2)),moduli{jj}(end-1,IND_s(2))*1.015,'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}$', 'Interpreter', 'latex','FontSize',18)
text(moduli{jj}(end,IND_s(3)),moduli{jj}(end-1,IND_s(3))*1.01,'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}$', 'Interpreter', 'latex','FontSize',18)
text(moduli{jj}(end,IND_s(4))-0.0001,moduli{jj}(end-1,IND_s(4)),'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {4}}}$', 'Interpreter', 'latex','FontSize',18)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
hold off

h=figure(6)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii  = 1:length(IND_s)
    %ufull = full{jj}(:,Nx*IND_s(ii)+1:Nx*(IND_s(ii)+1));
    ufull = sprof(x,y,kky,moduli{jj}(:,IND_s(ii)));
    subplot(4,1,ii)
    imagesc(x,y,ufull)
    xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    title([num2str(ii)])
end
hold off
drawnow






% kkx = moduli{jj}(end-1,ii);
% cc = moduli{jj}(end,ii);
% ufull = full{jj}(:,Nx*ii+1:Nx*(ii+1));
% figure(5)
% imagesc(x,y,ufull)
% xlabel('x')
% ylabel('y')
% title([' c = ', num2str(cc),' ky = ', num2str(KY_I(jj)), ', kx = ', num2str(kkx)]) 


% 
% ii = 100;
% jj = 3;
% Nx= floor(Lx/dx)+1;
% Ny = floor(Ly/dy)+1;
% dx = Lx/Nx;
% dy = Ly/Ny;
% y = dy*[1:Ny]';
% x =[dx:dx:Lx]';
% kkx = moduli{jj}(end-1,ii);
% cc = moduli{jj}(end,ii);
% ufull = full{jj}(:,Nx*ii+1:Nx*(ii+1));
% figure(6)
% imagesc(x,y,ufull)
% xlabel('x')
% ylabel('y')
% title([' c = ', num2str(cc),' ky = ', num2str(KY_I(jj)), ', kx = ', num2str(kkx)]) 

% 
% %%%%Run through one continuation run (can put a movie wrapper around this
% %%%%loop to record, also plots the corresponding point on a cross-section
% %%%%of the moduli surface
% jj = 10;
% for ii = 100;%1:3:length(moduli{jj}(1,:))
%     %%Plot a specific solution
% %     Nx= floor(Lx/dx)+1;
% %     Ny = floor(Ly/dy)+1;
% %     y = dy*[1:Ny]';
% %     x =[dx:dx:Lx]';
%     KXdum = moduli{jj}(end-1,:);
%     Cdum = moduli{jj}(end,:);
%     kkx =KXdum(ii);
%     cc = Cdum(ii);
%     kky = KY_I(jj);
%     %%function to output 
%     w = moduli{jj}(:,ii);
%     ufull = sprof(x,y,kky,w);
%     
%    % ufull = full{jj}(:,Nx*ii+1:Nx*(ii+1));
%     figure(3)
%     imagesc(x,y,ufull)
%     xlabel('x')
%     ylabel('y')
%     title([' c = ', num2str(cc),' ky = ', num2str(kky), ', kx = ', num2str(kkx)]) 
%     
% 
%     figure(4)
%     plot(Cdum,KXdum,'-',cc,kkx,'x')
%     xlabel('c')
%     ylabel('k_x')
%     title(['ky = ', num2str(kky)])
%     drawnow
% end

    
    
    
    
    
    
    

