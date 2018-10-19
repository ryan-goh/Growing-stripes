%%%Change ky value in filename here for perp stripes and below for oblique
kychoose =1.064602e+00;%%%.78;1.120102+00;9.980035e-01;

fname = sprintf('ContDataPerp%d.mat',kychoose)
load(fname)
PROF1 = PROF;


%%make a movie of solutions along moduli slice
makevideo = 0;



x_higher_order = false;


%Need to formula numerous things to get the jacobian back out
Lx = x(end);
Ly = y(end);
dx = x(2) - x(1);
dy = y(2) - y(1);
Nx = floor(Lx/dx);
Ny = floor(Ly/dy);
N = Nx*Ny

L2N1= vecnorm(PROF(1:end-2,:),2);




%%%%%%Plot some perp solutions before and after the fold curve
IND_s = [20 65 length(PROF(end,:))-10];%[20 65 180];
hb = figure(500)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
plot(PROF(end,IND_s),PROF(end-1,IND_s),'ko','LineWidth',2)
text(PROF(end,IND_s(1)),PROF(end-1,IND_s(1))+0.005,'(1)','Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}$', 'Interpreter', 'latex','FontSize',18)
text(PROF(end,IND_s(2))+0.02,PROF(end-1,IND_s(2))+0.005,'(2)','Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {2}}}$', 'Interpreter', 'latex','FontSize',18)
text(PROF(end,IND_s(3)),PROF(end-1,IND_s(3))+0.005,'(3)','Interpreter', 'latex','FontSize',18);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {3}}}$', 'Interpreter', 'latex','FontSize',18)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
hold off
kky = ky;


hp=figure(501)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii  = 1:length(IND_s)
    %ufull = full{jj}(:,Nx*IND_s(ii)+1:Nx*(IND_s(ii)+1));
    ufull = sprof(x,y,kky,PROF(:,IND_s(ii)));%moduli{jj}(:,IND_s(ii)));
    subplot(4,1,ii)
    imagesc(x,y,ufull)

    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 28);
        xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    set(gca, 'FontSize', 16);

    set(yaxis,'FontSize',28)
    title(['(',num2str(ii),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])
end
hold off

drawnow






%%%%%%%%%%%load the oblique stripesload('ContData9.980035e-01.mat') and
%%%%add to the plot
%load('ContData9.980035e-01.mat')
fname = sprintf('ContData%d.mat',kychoose)
load(fname)

Lx = x(end);
Ly = y(end);
dx = x(2) - x(1);
dy = y(2) - y(1);
Nx = floor(Lx/dx);
Ny = floor(Ly/dy);
N = Nx*Ny
IND_s1 = size(PROF,2)-10;140;

[kxmm,iimx] = min(PROF(end-1,:));
PROF = [PROF(:,iimx:-1:1),PROF(:,iimx+1:end)];

L2N = vecnorm(PROF(1:end-2,:));


figure(500)
hold on
plot(PROF(end,:),PROF(end-1,:),'LineWidth',2,'Color', [0    0.4470    0.7410])
plot(PROF(end,IND_s1),PROF(end-1,IND_s1),'ko','LineWidth',2)
text(PROF(end,IND_s1),PROF(end-1,IND_s1)-0.01,'(4)','Interpreter', 'latex','FontSize',18);%['$\raisebox{.5pt}{\' ...
                   % 'textcircled{\raisebox{-.9pt} {4}}}$'], ...
     %'Interpreter', 'latex','FontSize',18)
xlim([-0.0001 max(PROF(end,:))+0.01])
ylim([-0.015 max(PROF(end-1,:))])
%ylim([-0.015 0.4092])
hold off
title(['$k_y$ = ', num2str(kychoose)],'Interpreter','latex')



kky = ky
figure(501)
hold on
ufull = sprof(x,y,kky,PROF(:,IND_s1));%moduli{jj}(:,IND_s(ii)));
    subplot(4,1,4)
    imagesc(x,y,ufull)

    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    %set(xaxis,'FontSize',28)
    %set(yaxis,'FontSize',28);
    xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    set(gca, 'FontSize', 16);
        xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 28);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 28);


     title(['(',num2str(4),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])
    hold off
%     xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
% set(xaxis,'FontSize',20)
    drawnow



    if makevideo
fname = sprintf('movie_ky%d',ky)
vidfile = VideoWriter(fname);
open(vidfile);
for ii = 1:length(PROF1(end,:))
    znew = PROF1(:,ii);
     ufull = sprof(x,y,kky,znew);

    fig = figure(600)
    box on
    set(gcf, 'Position', [400, 400, 1200, 600]);
    set(gcf,'PaperPositionMode','auto')
    subplot(2,2,1)
    plot(PROF1(end,:),PROF1(end-1,:),'LineWidth',2)
    hold on
    plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
    plot(PROF1(end,ii),PROF1(end-1,ii),'ko','LineWidth',2)
    hold off
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
    subplot(2,2,2)
    plot(PROF1(end,:),L2N1,'LineWidth',2)
    hold on
    plot(PROF(end,:),L2N,'LineWidth',2)
    plot(PROF1(end,ii),L2N1(ii),'ko','LineWidth',2)
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
    daspect([1 1 1]);
    drawnow;
    Frame = getframe(fig);
    writeVideo(vidfile,Frame);
end

for ii = 1:1:length(PROF(end,:))
    znew = PROF(:,ii);
     ufull = sprof(x,y,kky,znew);

    fig = figure(600)
    box on
    set(gcf, 'Position', [400, 400, 1200, 600]);
    set(gcf,'PaperPositionMode','auto')
    subplot(2,2,1)
    plot(PROF1(end,:),PROF1(end-1,:),'LineWidth',2)
    hold on
    plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
    plot(PROF(end,ii),PROF(end-1,ii),'ko','LineWidth',2)
    hold off
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
    subplot(2,2,2)
    plot(PROF1(end,:),L2N1,'LineWidth',2)
    hold on
    plot(PROF(end,:),L2N,'LineWidth',2)
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
     daspect([1 1 1]);
    drawnow;
    Frame = getframe(fig);
    writeVideo(vidfile,Frame);
end






close(vidfile)







end
