%%%Barba Script ---- Figure 8 isola of perpendicular stripes, near fold of
%%%folds at ~.78
clear all
close all

makevideo = 1;

load('ContDataPerp7.800000e-01.mat')
kky = ky;
pPROF = PROF;
pL2 = [];
for ii = 1:length(PROF(end,:))
    pL2 = [pL2,norm(PROF(1:end-2,ii),2)];
end

pI = [ 75 43];
figure(3)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii = 1:2
 ufull = sprof(x,y,kky,pPROF(:,pI(ii)));
subplot(4,1,ii)
imagesc(x,y,ufull)
    %xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    %xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    %set(xaxis,'FontSize',28);
    set(yaxis,'FontSize',28);
    xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    title(['(',num2str(ii),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])
end
hold off
drawnow




load('ContData7.800000e-01.mat')
ipp = find(PROF(end-1,:)<0);
iii = 1:length(PROF(end,:));
iii(ipp) = [];
PROF = PROF(:,iii);  %%%throws out negative speeds
imm = find(abs(PROF(end-1,:))<0.001);
[kdumm,imin] = min(PROF(end-1,1:imm+200));

oPROF = [PROF(:,end:-1:imin+1),PROF(:,1:imin)];




oL2 = [];
for ii = 1:length(oPROF(end,:))
    oL2 = [oL2,norm(oPROF(1:end-2,ii),2)];
    end

oI = [300 525];
figure(3)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii = 1:2
 ufull = sprof(x,y,kky,oPROF(:,oI(ii)));
subplot(4,1,ii+2)
imagesc(x,y,ufull)
%     xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
   % xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
%    set(xaxis,'FontSize',28);
    set(yaxis,'FontSize',28);
    xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    title(['(',num2str(ii+2),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])
end
hold off
  xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 28);
  xt = get(gca, 'XTick');
drawnow

figure(2)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(pPROF(end,:),pPROF(end-1,:),'LineWidth',2)
plot(oPROF(end,:),oPROF(end-1,:),'LineWidth',2)
plot(pPROF(end,pI),pPROF(end-1,pI),'ko','LineWidth',2)
plot(oPROF(end,oI),oPROF(end-1,oI),'ko','LineWidth',2)
for ii = 1:2
    text(pPROF(end,pI(ii)),pPROF(end-1,pI(ii))+0.005,['(',num2str(ii),')'],'Interpreter', 'latex','FontSize',18);
    text(oPROF(end,oI(ii)),oPROF(end-1,oI(ii))+0.005,['(',num2str(ii+2),')'],'Interpreter', 'latex','FontSize',18);
end
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',28)
set(yaxis,'FontSize',28)
ylim([-0.01 0.12])
hold off



figure(4)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(pPROF(end,:),pL2(:),'LineWidth',2)
plot(oPROF(end,:),oL2(:),'LineWidth',2)
plot(pPROF(end,pI),pL2(pI),'ko','LineWidth',2)
plot(oPROF(end,oI),oL2(oI),'ko','LineWidth',2)
for ii = 1:2
    text(pPROF(end,pI(ii)),pL2(pI(ii))+0.015,['(',num2str(ii),')'],'Interpreter', 'latex','FontSize',18);
    text(oPROF(end,oI(ii)),oL2(oI(ii))+0.01,['(',num2str(ii+2),')'],'Interpreter', 'latex','FontSize',18);
end
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',30)
set(yaxis,'FontSize',30)
xlim([0.012 0.045])
hold off



if makevideo
    fname = sprintf('barba_ky%d',ky);
    vidfile = VideoWriter(fname);
    open(vidfile);
    for ii = 1:1:length(pPROF(end,:))
       znew = pPROF(:,ii);
      ufull =  sprof(x,y,kky,znew);

       fig = figure(100)
       box on
       set(gcf,'Position',[400,400,900,600])
        set(gcf,'PaperPositionMode','auto')
        subplot(2,2,[3,4])
        imagesc(x,y,ufull)
        yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
         xt = get(gca, 'XTick');
         set(gca, 'FontSize', 16);
        set(xaxis,'FontSize',28);
      set(yaxis,'FontSize',28);
        xticks(0:20*pi:40*pi)
        xticklabels({'0','20\pi','40\pi'})
        yticks(0:pi:2*pi)
        yticklabels({'0','\pi','2\pi'})
        daspect([ 1 1 1])
        %title(['(',num2str(ii+2),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])

        subplot(2,2,1)

    plot(pPROF(end,:),pPROF(end-1,:),'LineWidth',2)
        hold on
    plot(oPROF(end,:),oPROF(end-1,:),'LineWidth',2)
    plot(pPROF(end,ii),pPROF(end-1,ii),'ko','LineWidth',2)
   % plot(oPROF(end,oI),oPROF(end-1,oI),'ko','LineWidth',2)
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',28)
    set(yaxis,'FontSize',28)
    ylim([-0.01 0.12])
    hold off



        subplot(2,2,2)

    plot(pPROF(end,:),pL2(:),'LineWidth',2)
           hold on
    plot(oPROF(end,:),oL2(:),'LineWidth',2)
    plot(pPROF(end,ii),pL2(ii),'ko','LineWidth',2)
    %plot(oPROF(end,oI),oL2(oI),'ko','LineWidth',2)
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',30)
    set(yaxis,'FontSize',30)
    xlim([0.012 0.045])
    hold off

    drawnow;
    Frame = getframe(fig);
    writeVideo(vidfile,Frame);

    end



 for ii = 1:2:length(oPROF(end,:))
       znew = oPROF(:,ii);
       ufull = sprof(x,y,kky,znew);

       fig = figure(100)

        subplot(2,2,[3,4])
        imagesc(x,y,ufull)
        yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
        xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
         xt = get(gca, 'XTick');
         set(gca, 'FontSize', 16);
        set(xaxis,'FontSize',28);
      set(yaxis,'FontSize',28);
        xticks(0:20*pi:40*pi)
        xticklabels({'0','20\pi','40\pi'})
        yticks(0:pi:2*pi)
        yticklabels({'0','\pi','2\pi'})
        daspect([ 1 1 1])
        %title(['(',num2str(ii+2),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])

        subplot(2,2,1)
    plot(pPROF(end,:),pPROF(end-1,:),'LineWidth',2)
       hold on
    plot(oPROF(end,:),oPROF(end-1,:),'LineWidth',2)
   % plot(pPROF(end,ii),pPROF(end-1,ii),'ko','LineWidth',2)
    plot(oPROF(end,ii),oPROF(end-1,ii),'ko','LineWidth',2)
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',28)
    set(yaxis,'FontSize',28)
    ylim([-0.01 0.12])
    hold off



        subplot(2,2,2)
    plot(pPROF(end,:),pL2(:),'LineWidth',2)
            hold on
    plot(oPROF(end,:),oL2(:),'LineWidth',2)
    %plot(pPROF(end,ii),pL2(ii),'ko','LineWidth',2)
    plot(oPROF(end,ii),oL2(ii),'ko','LineWidth',2)
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',30)
    set(yaxis,'FontSize',30)
    xlim([0.012 0.045])
    hold off
        drawnow;
    Frame = getframe(fig);
    writeVideo(vidfile,Frame);


    end

    close(vidfile)

end







% figure(1)
% subplot(1,2,1)
% plot(pPROF(end,:),pPROF(end-1,:),PROF(end,:),PROF(end-1,:))
% xlabel('c')
% ylabel('k_x')
% subplot(1,2,2)
% plot( pPROF(end,:),pL2,oPROF(end,:),oL2)
% xlabel('c')
% ylabel('| w|_{L^2}')
% title('k_y = ', num2str(ky))


%%%Plot some plain bifurcation diagrams
ky1 = ky;
load('ContDataPerp7.780000e-01 (1).mat')
p2PROF = PROF;
ky2 = ky;
 p2L2 = [];
for ii = 1:length(p2PROF(end,:))
    p2L2 = [p2L2,norm(p2PROF(1:end-2,ii),2)];
end

load('ContDataPerp7.790000e-01 (1).mat')
p3PROF = PROF;
ky3 = ky;
 p3L2 = [];
for ii = 1:length(p3PROF(end,:))
    p3L2 = [p3L2,norm(p3PROF(1:end-2,ii),2)];
end

% %fa = 28;
% figure(7)
% box on
% set(gcf, 'Position', [400, 400, 1300, 600]);
% set(gcf,'PaperPositionMode','auto')
%     subplot(1,3,3)
% plot(pPROF(end,:),pL2,'LineWidth',2)
% xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
% yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 16);
% set(xaxis,'FontSize',28)
% set(yaxis,'FontSize',28)
% title(['$k_y = $ ',num2str(ky1)],'Interpreter','latex','FontSize',24)
% %xlim([0.012 0.045])
%     subplot(1,3,2)
% plot(p3PROF(end,:),p3L2,'LineWidth',2)
% xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
% yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 16);
% set(xaxis,'FontSize',28)
% set(yaxis,'FontSize',28)
% title(['$k_y = $ ',num2str(ky3)],'Interpreter','latex','FontSize',24)
%     subplot(1,3,1)
% plot(p2PROF(end,:),p2L2,'LineWidth',2)
% xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
% yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 16);
% set(xaxis,'FontSize',28)
% set(yaxis,'FontSize',28)
% title(['$k_y = $ ',num2str(ky2)],'Interpreter','latex','FontSize',24)
% %xlim([0.012 0.045])
%
% %xlim([0.012 0.045])



load('ContDataPerp7.779000e-01.mat')
p4PROF = PROF;
ky4 = ky;
 p4L2 = [];
for ii = 1:length(p4PROF(end,:))
    p4L2 = [p4L2,norm(p4PROF(1:end-2,ii),2)];
end

%fa = 28;
figure(7)
box on
set(gcf, 'Position', [400, 400, 1300, 600]);
set(gcf,'PaperPositionMode','auto')
    subplot(1,4,4)
plot(pPROF(end,:),pL2,'LineWidth',2)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',28)
set(yaxis,'FontSize',28)
title(['$k_y = $ ',num2str(ky1)],'Interpreter','latex','FontSize',24)
%xlim([0.012 0.045])
    subplot(1,4,3)
plot(p3PROF(end,:),p3L2,'LineWidth',2)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',28)
set(yaxis,'FontSize',28)
title(['$k_y = $ ',num2str(ky3)],'Interpreter','latex','FontSize',24)
    subplot(1,4,2)
plot(p2PROF(end,:),p2L2,'LineWidth',2)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',28)
set(yaxis,'FontSize',28)
title(['$k_y = $ ',num2str(ky2)],'Interpreter','latex','FontSize',24)
%xlim([0.012 0.045])
   subplot(1,4,1)
plot(p4PROF(end,:),p4L2,'LineWidth',2)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',28)
set(yaxis,'FontSize',28)
title(['$k_y = $ ',num2str(ky4)],'Interpreter','latex','FontSize',24)
%xlim([0.012 0.045])










% load('ContDataPerp7.780000e-01.mat')
% kky = ky;
% iL2 = [];
% for ii = 1:length(PROF(end,:))
%     iL2 = [iL2,norm(PROF(1:end-2,ii),2)];
% end
%
%
%
% figure(6)
% box on
% set(gcf, 'Position', [400, 400, 900, 600]);
% set(gcf,'PaperPositionMode','auto')
% hold on
% %plot(PROF(end,1:100),iL2(1:100),'LineWidth',2)
% scatter(PROF(end,:),iL2(:),'-','LineWidth',2)
% xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
% yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 16);
% set(xaxis,'FontSize',20)
% set(yaxis,'FontSize',20)
% %xlim([0.012 0.045])
% hold off
%
