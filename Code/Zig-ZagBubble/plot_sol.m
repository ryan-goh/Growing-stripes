%%%%Plotting figures for zig-zag bubble:
global mu0
%%loads the ky continuation (with c=0) which the bubble emanates from, then
%%clears other data
load('ContData0.mat','PROF')
KX_II = PROF(end-1,:);
KY_II = PROF(end,:);
clear PROF


load('moduli_dat_zz.mat')%,'KY_I')



%%%For making a .avi video of solutions on bifurcaiton curve
makevideo = 1;

%%Plot a series of solutions
%pick ky slice:
jj = 5;
kky=KY_I(jj)
IND_s = [50, 505, 700 900]; %index of point along continuation curve for plotted solutions
Nx= floor(Lx/dx)+1;
Ny = floor(Ly/dy)+1;
dxx = Lx/Nx;
dyy = Ly/Ny;
y = dyy*[1:Ny]';
x =[dxx:dxx:Lx]';

h=figure(5)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(moduli{jj}(end,:),moduli{jj}(end-1,:),'LineWidth',2)
plot(moduli{jj}(end,IND_s),moduli{jj}(end-1,IND_s),'ko','LineWidth',2)
text(moduli{jj}(end,IND_s(1)),moduli{jj}(end-1,IND_s(1))*1.01,'(1)', 'Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}$', 'Interpreter', 'latex','FontSize',18)
text(moduli{jj}(end,IND_s(2)),moduli{jj}(end-1,IND_s(2))*1.02,'(2)', 'Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}$', 'Interpreter', 'latex','FontSize',18)
text(moduli{jj}(end,IND_s(3)),moduli{jj}(end-1,IND_s(3))*1.01,'(3)', 'Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}$', 'Interpreter', 'latex','FontSize',18)
text(moduli{jj}(end,IND_s(4)),moduli{jj}(end-1,IND_s(4))-0.006,'(4)', 'Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {4}}}$', 'Interpreter', 'latex','FontSize',18)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
xlim([0, max(moduli{jj}(end,:))])
title(['$k_y = $',num2str(kky)],'Interpreter', 'latex', 'FontSize', 20)
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
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    set(gca, 'FontSize', 16);
    set(yaxis,'FontSize',20);
    title(['(' ,num2str(ii),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])
end
hold off
xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
 set(xaxis,'FontSize',20);

drawnow















    if makevideo
fname = sprintf('movie_ky%d',kky)
vidfile = VideoWriter(fname);
open(vidfile);
PROF = moduli{jj}(:,:);
 iif = find(PROF(end,:)<PROF(end,1));
PROF =  [PROF(:,iif(end:-1:1)),PROF(:,1:iif(1)-1)];
L2N = vecnorm(PROF(1:end-2,:),2);

for ii = 1:1:length(PROF(end,:))
    znew = PROF(:,ii);
     ufull = sprof(x,y,kky,znew);

    fig = figure(600)
    box on
    set(gcf, 'Position', [400, 400, 1200, 600]);
    set(gcf,'PaperPositionMode','auto')
    subplot(2,2,1)
    plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
    hold on
    plot(PROF(end,ii),PROF(end-1,ii),'ko','LineWidth',2)
    hold off
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
    subplot(2,2,2)
    plot(PROF(end,:),L2N,'LineWidth',2)
    hold on
    plot(PROF(end,ii),L2N(ii),'ko','LineWidth',2)
    hold off
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$\|w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
    subplot(2,2,[3,4])
     imagesc(x,y,ufull)
    xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
        xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
    daspect([1 1 1])
    drawnow;
    Frame = getframe(fig);
    writeVideo(vidfile,Frame);
end






close(vidfile)







end
