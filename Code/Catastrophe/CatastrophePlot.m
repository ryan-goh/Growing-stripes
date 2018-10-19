%%%%Short script to make the catastrophe photos:
makevideo = 0;

load('ContDataZZ8.471268e-01.mat') %%For plots of moduli slice before Catastrophe


iii = find(PROF(end,:)>0);
C1 = PROF(end,iii);
KX1 = PROF(end-1,iii);
PROF1 = PROF;

ky1 = .84712;
ip11 = 400;
ip12 = 1200;
znew1 = PROF(:,ip11);
znew2 = PROF(:,ip12);
ufull1 = sprof(x,y,ky1,znew1);
ufull2 = sprof(x,y,ky1,znew2);
figure(1)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
subplot(4,1,1)
imagesc(x,y,ufull1)
%xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
  xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
set(gca, 'FontSize', 16);
%set(xaxis,'FontSize',24)
set(yaxis,'FontSize',24)
title(['(' ,num2str(1),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])
subplot(4,1,2)
imagesc(x,y,ufull2)
%xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
  xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
set(gca, 'FontSize', 16);
%set(xaxis,'FontSize',24)
set(yaxis,'FontSize',28)
title(['(' ,num2str(2),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])


figure(2)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
plot(C1,KX1,'LineWidth',3,'Color',[0    0.4470    0.7410])
hold on
plot(PROF(end,ip11),PROF(end-1,ip11),'ko','MarkerSize',5)
plot(PROF(end,ip12),PROF(end-1,ip12),'ko','MarkerSize',5)
text(PROF(end,ip11)-0.02,PROF(end-1,ip11)-0.02,'(1)', 'Interpreter', 'latex','FontSize',18);
text(PROF(end,ip12)+0.02,PROF(end-1,ip12)+0.02,'(2)', 'Interpreter', 'latex','FontSize',18);
hold off


%%%%Have to import these from c-cont directory
%load('ContData8.454545e-01.mat')
load('ContData8.471268e-01.mat')
iii = find(PROF(end,:)>0);
C1o = PROF(end,iii);
KX1o = PROF(end-1,iii);


iif = find(KX1o<KX1o(1));
C1o = [C1o(iif(end:-1:1)),C1o(1:iif(1)-10)];
KX1o = [KX1o(iif(end:-1:1)),KX1o(1:iif(1)-10)];

ip21 = 1000;
ip22 = 1480;
znew1 = PROF(:,ip21);
znew2 = PROF(:,ip22);
ufull3 = sprof(x,y,ky1,znew1);
ufull4 = sprof(x,y,ky1,znew2);
figure(1)
hold on
subplot(4,1,3)
imagesc(x,y,ufull3)
%xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
  xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
set(gca, 'FontSize', 16);
%set(xaxis,'FontSize',24)
set(yaxis,'FontSize',28)
title(['(' ,num2str(3),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])

subplot(4,1,4)
imagesc(x,y,ufull4)
xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
  xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',28)
set(yaxis,'FontSize',28)
title(['(' ,num2str(4),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])

hold off


figure(2)
hold on
plot(C1o,KX1o,'Linewidth',3,'Color',[0    0.4470    0.7410])
plot(PROF(end,ip21),PROF(end-1,ip21),'ko','MarkerSize',5)
plot(PROF(end,ip22),PROF(end-1,ip22),'ko','MarkerSize',5)
text(PROF(end,ip21)+0.02,PROF(end-1,ip21),'(3)', 'Interpreter', 'latex','FontSize',16);
text(PROF(end,ip22),PROF(end-1,ip22)-0.02,'(4)', 'Interpreter', 'latex','FontSize',16);
hold off
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',24)
set(yaxis,'FontSize',24)
title(['$k_y\approx$ ' ,num2str(ky1),],'Interpreter', 'latex', 'FontSize', 20)





    if makevideo
fname = sprintf('movie_ky%d',ky1)
vidfile = VideoWriter(fname);
open(vidfile);
for ii = 1:3:length(PROF1(end,:))
    znew = PROF1(:,ii);
     ufull = sprof(x,y,ky1,znew);

    fig = figure(600)
    box on
    set(gcf, 'Position', [400, 400, 1200, 600]);
    set(gcf,'PaperPositionMode','auto')
    subplot(2,1,1)
    plot(C1o,KX1o,'Linewidth',3,'Color',[0    0.4470    0.7410])
    hold on
    plot(C1,KX1,'Linewidth',3,'Color',[0    0.4470    0.7410])
    plot(PROF1(end,ii),PROF1(end-1,ii),'ko','LineWidth',2)
    hold off
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
%     subplot(3,1,2)
%     plot(PROF1(end,:),L2N1,'LineWidth',2)
%     hold on
%     plot(PROF(end,:),L2N,'LineWidth',2)
%     plot(PROF1(end,ii),L2N1(ii),'ko','LineWidth',2)
%     hold off
%     xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
%     yaxis=ylabel('$\|w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 16);
%     set(xaxis,'FontSize',20);
%     set(yaxis,'FontSize',20);
    subplot(2,1,2)
     imagesc(x,y,ufull)
    xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
        xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    daspect([1 1 1])
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
    drawnow;
    Frame = getframe(fig);
    writeVideo(vidfile,Frame);
end

for ii = 1:3:length(PROF(end,:))
    znew = PROF(:,ii);
     ufull = sprof(x,y,ky1,znew);

    fig = figure(600)
    box on
    set(gcf, 'Position', [400, 400, 1200, 600]);
    set(gcf,'PaperPositionMode','auto')
    subplot(2,1,1)
   plot(C1o,KX1o,'Linewidth',3,'Color',[0    0.4470    0.7410])
    hold on
    plot(C1,KX1,'Linewidth',3,'Color',[0    0.4470    0.7410])
    plot(PROF(end,ii),PROF(end-1,ii),'ko','LineWidth',2)
    hold off
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
%     subplot(3,1,2)
%     plot(PROF1(end,:),L2N1,'LineWidth',2)
%     hold on
%     plot(PROF(end,:),L2N,'LineWidth',2)
%     plot(PROF(end,ii),L2N(ii),'ko','LineWidth',2)
%     hold off
%     xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
%     yaxis=ylabel('$\|w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 16);
%     set(xaxis,'FontSize',20);
%     set(yaxis,'FontSize',20);
    subplot(2,1,2)
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




















%%%%%%
%%%%Short script to make the catastrophe photos:


%load('ContDataZZ8.355340e-01.mat')
load('ContDataZZ8.419514e-01.mat')

iii = find(PROF(end,:)>1e-4);
C2 = PROF(end,iii);
KX2 = PROF(end-1,iii);
PROF = PROF(:,iii);
PROF1 = PROF;


ky2 = .842;
ip11 = 900;
ip12 = 2200;
znew1 = PROF(:,ip11);
znew2 = PROF(:,ip12)
ufull1 = sprof(x,y,ky2,znew1);
ufull2 = sprof(x,y,ky2,znew2);
figure(3)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
subplot(4,1,1)
imagesc(x,y,ufull1)
%xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
  xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
set(gca, 'FontSize', 16);
%set(xaxis,'FontSize',24)
set(yaxis,'FontSize',28)
title(['(' ,num2str(1),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])

subplot(4,1,2)
imagesc(x,y,ufull2)
%xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
  xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
set(gca, 'FontSize', 16);
%set(xaxis,'FontSize',24)
set(yaxis,'FontSize',28)
title(['(' ,num2str(2),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])


figure(4)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
plot(C2,KX2,'Linewidth',3,'Color',[0    0.4470    0.7410])
hold on
plot(PROF(end,ip11),PROF(end-1,ip11),'ko','MarkerSize',5)
plot(PROF(end,ip12),PROF(end-1,ip12),'ko','MarkerSize',5)
text(PROF(end,ip11),PROF(end-1,ip11)+0.03,'(1)', 'Interpreter', 'latex','FontSize',16);
text(PROF(end,ip12),PROF(end-1,ip12)-0.02,'(2)', 'Interpreter', 'latex','FontSize',16);
hold off


%%%%Have to import these from Perp_Reattachment directory
%load('ContData8.315061e-01.mat')
load('ContData8.410000e-01.mat')
iii = find(PROF(end,:)>0);
C2o = PROF(end,iii);
KX2o = PROF(end-1,iii);

iif = find(KX2o<KX2o(1));
C2o = [C2o(iif(end:-1:1)),C2o(1:iif(1)-10)];
KX2o = [KX2o(iif(end:-1:1)),KX2o(1:iif(1)-10)];

ip21 = 1600;
ip22 = 1050;
znew1 = PROF(:,ip21);
znew2 = PROF(:,ip22);
ufull3 = sprof(x,y,ky2,znew1);
ufull4 = sprof(x,y,ky2,znew2);
figure(3)
hold on
subplot(4,1,3)
imagesc(x,y,ufull3)
%xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
  xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
set(gca, 'FontSize', 16);
%set(xaxis,'FontSize',24)
set(yaxis,'FontSize',28)
title(['(' ,num2str(3),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])

subplot(4,1,4)
imagesc(x,y,ufull4)
xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20)
  xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',28)
set(yaxis,'FontSize',28)
title(['(' ,num2str(4),')'],'Interpreter', 'latex', 'FontSize', 20,'Position',[130,5])



figure(4)
hold on
plot(C2o,KX2o,'Linewidth',3,'Color',[0    0.4470    0.7410])
plot(PROF(end,ip21),PROF(end-1,ip21),'ko','MarkerSize',5)
plot(PROF(end,ip22),PROF(end-1,ip22),'ko','MarkerSize',5)
text(PROF(end,ip21),PROF(end-1,ip21)+0.02,'(3)', 'Interpreter', 'latex','FontSize',16);
text(PROF(end,ip22)+0.02,PROF(end-1,ip22)-0.01,'(4)', 'Interpreter', 'latex','FontSize',16);
hold off
title(['$k_y\approx$ ' ,num2str(ky2)],'Interpreter', 'latex', 'FontSize', 20)
xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20)
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20)
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',24)
set(yaxis,'FontSize',24)
ylim([min(KX2o),max(KX2o)+0.1])


    if makevideo
fname = sprintf('movie_ky%d',ky2)
vidfile = VideoWriter(fname);
open(vidfile);
for ii = 1:3:length(PROF1(end,:))
    znew = PROF1(:,ii);
     ufull = sprof(x,y,ky2,znew);

    fig = figure(600);
    box on
    set(gcf, 'Position', [400, 400, 1200, 600]);
    set(gcf,'PaperPositionMode','auto')
    subplot(2,1,1)
    plot(C2o,KX2o,'Linewidth',3,'Color',[0    0.4470    0.7410])
    hold on
    plot(C2,KX2,'Linewidth',3,'Color',[0    0.4470    0.7410])
    plot(PROF1(end,ii),PROF1(end-1,ii),'ko','LineWidth',2)
    hold off
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
%     subplot(3,1,2)
%     plot(PROF1(end,:),L2N1,'LineWidth',2)
%     hold on
%     plot(PROF(end,:),L2N,'LineWidth',2)
%     plot(PROF1(end,ii),L2N1(ii),'ko','LineWidth',2)
%     hold off
%     xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
%     yaxis=ylabel('$\|w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 16);
%     set(xaxis,'FontSize',20);
%     set(yaxis,'FontSize',20);
    subplot(2,1,2)
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

for ii = 1:3:length(PROF(end,:))
    znew = PROF(:,ii);
     ufull = sprof(x,y,ky2,znew);

    fig = figure(600);
    box on
    set(gcf, 'Position', [400, 400, 1200, 600]);
    set(gcf,'PaperPositionMode','auto')
    subplot(2,1,1)
    plot(C2o,KX2o,'Linewidth',3,'Color',[0    0.4470    0.7410])
    hold on
    plot(C2,KX2,'Linewidth',3,'Color',[0    0.4470    0.7410])
    plot(PROF(end,ii),PROF(end-1,ii),'ko','LineWidth',2)
    hold off
    xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20);
    set(yaxis,'FontSize',20);
%     subplot(3,1,2)
%     plot(PROF1(end,:),L2N1,'LineWidth',2)
%     hold on
%     plot(PROF(end,:),L2N,'LineWidth',2)
%     plot(PROF(end,ii),L2N(ii),'ko','LineWidth',2)
%     hold off
%     xaxis=xlabel('$c_x$', 'Interpreter', 'latex', 'FontSize', 20);
%     yaxis=ylabel('$\|w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 16);
%     set(xaxis,'FontSize',20);
%     set(yaxis,'FontSize',20);
    subplot(2,1,2)
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
