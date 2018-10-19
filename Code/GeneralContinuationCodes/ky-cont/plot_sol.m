%%%%Script to plot a slice of moduli space in ky, with c fixed
close all
%%%Change cchose for different slices of moduli space
%cchoose = 2.948980e-01;
cchoose =9.081633e-02 ;
%%1.064602e+00;1.120102+00;9.980035e-01;
%load('ContData9.980035e-01.mat')
fname = sprintf('ContData%d.mat',cchoose)
load(fname)
PROF1 = PROF;

spp = size(PROF,2);

if cchoose == 9.081633e-02 ;
    iif = find(PROF1(end,:)<0.002);
    PROF1 = PROF1(:,iif(2):end);
end
L2N1 = vecnorm(PROF1(1:end-2,:),2);

c = cchoose;
cc = c;

makevideo = 1;

x_higher_order = false;

%%%plot core-farfield decomposition:
pl = 0;


%Need to formula numerous things to get the jacobian back out
Lx = x(end);
Ly = y(end);
dx = x(2) - x(1);
dy = y(2) - y(1);
Nx = floor(Lx/dx);
Ny = floor(Ly/dy);
N = Nx*Ny




IND_s = [20 spp-450];%<- for cslice.eps
%IND_s = [20 1075 1300];%<-for csmallslice.eps
hb = figure(500)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
plot(PROF(end,IND_s),PROF(end-1,IND_s),'ko','LineWidth',2)
for ii = 1:length(IND_s)
text(PROF(end,IND_s(ii)),PROF(end-1,IND_s(ii))-0.03,['(',num2str(ii),')'],'Interpreter', 'latex','FontSize',16);%'$\raisebox{.5pt}{\textcircled{\raisebox{-.9pt} {1}}}$', 'Interpreter', 'latex','FontSize',18)
end

xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
hold off




hp=figure(501)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii  = 1:length(IND_s)
    %ufull = full{jj}(:,Nx*IND_s(ii)+1:Nx*(IND_s(ii)+1));
    ufull = sprof(x,y,cc,PROF(:,IND_s(ii)),pl);%moduli{jj}(:,IND_s(ii)));
    subplot(4,1,ii)
    imagesc(x,y,ufull)
    %xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
   yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    %xt = get(gca, 'XTick');
    xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    set(gca, 'FontSize', 16);
    %set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',28)
    title(['(',num2str(ii),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])
end
hold off
%xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
%set(xaxis,'FontSize',20)
drawnow


    %%%
    ufull = sprof(x,y,c,PROF(:,IND_s(2)),1);%moduli{jj}(:,IND_s(ii)));
   figure(600)
    imagesc(x,y,ufull)



%%%%%%%%%%%load the perp stripesload('ContData9.980035e-01.mat') and
%%%%add to the plot
%load('ContData9.980035e-01.mat')
fname = sprintf('ContDataPerp%d.mat',cchoose)
load(fname)


Lx = x(end);
Ly = y(end);
dx = x(2) - x(1);
dy = y(2) - y(1);
Nx = floor(Lx/dx);
Ny = floor(Ly/dy);
N = Nx*Ny
sp = size(PROF,2);
 IND_s1 =[450,40]; %<- for cslice.eps   140;
%IND_s1 = 140;% <- for cslicesmall.eps

[kxmm,iimx] = min(PROF(end-1,:));
iif = find(abs(PROF(end,:)-PROF(end,1))<0.001);
if cchoose == 9.081633e-02 ;

    PROF = [PROF(:,iif(5)-1:-1:1),PROF(:,iif(5):end)];

else
    PROF = [PROF(:,iimx:-1:1),PROF(:,iimx+1:end)];

end



L2N = vecnorm(PROF(1:end-2,:),2);

figure(500)
box on
set(gcf, 'Position', [400, 400, 600, 400]);
set(gcf,'PaperPositionMode','auto')
hold on
plot(PROF(end,:),PROF(end-1,:),'LineWidth',2,'Color',[0 .447 .741])
for ii = 1:length(IND_s1)

plot(PROF(end,IND_s1(ii)),PROF(end-1,IND_s1(ii)),'ko','LineWidth',2)
text(PROF(end,IND_s1(ii)),PROF(end-1,IND_s1(ii))-0.03,['(',num2str(ii+length(IND_s)),')'],'Interpreter', 'latex','FontSize',16);%['$\raisebox{.5pt}{\' ...
                   % 'textcircled{\raisebox{-.9pt} {4}}}$'], ...
     %'Interpreter', 'latex','FontSize',18)
%xlim([-0.0001 max(PROF(end,:))+0.01])
%ylim([-0.015 max(PROF(end-1,:))])
%ylim([-0.015 0.4092])
end
hold off
 xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 28);
    yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 28);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',28)
    set(yaxis,'FontSize',28)
title(['$c_x$ = ', num2str(c)],'Interpreter','latex')




figure(501)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii = 1:length(IND_s1)
ufull = sprof(x,y,c,PROF(:,IND_s1(ii)),pl);%moduli{jj}(:,IND_s(ii)));
    subplot(4,1,ii+length(IND_s))
    imagesc(x,y,ufull)
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    xticks(0:20*pi:40*pi)
    xticklabels({'0','20\pi','40\pi'})
    yticks(0:pi:2*pi)
    yticklabels({'0','\pi','2\pi'})
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(yaxis,'FontSize',28)
     title(['(',num2str(ii)+length(IND_s),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])
end
    hold off
xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 28);
set(xaxis,'FontSize',28)

    drawnow



    if makevideo
      pl = 0;
      fname = sprintf('movie_c%d',c)
      vidfile = VideoWriter(fname);
      open(vidfile);


      for ii=1:length(PROF1(end,:))
        znew = PROF1(:,ii);
         ufull = sprof(x,y,c,znew,pl);

        fig = figure(600)
        box on
        set(gcf, 'Position', [400, 400, 1200, 600]);
        set(gcf,'PaperPositionMode','auto')
        subplot(2,2,1)
        plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
                hold on
                plot(PROF1(end,:),PROF1(end-1,:),'LineWidth',2,'Color',[0 .42 .71])
                plot(PROF1(end,ii),PROF1(end-1,ii),'ko','LineWidth',2)
        xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
        yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16);
        set(xaxis,'FontSize',20)
        set(yaxis,'FontSize',20)
        hold off


        subplot(2,2,2)
        plot(PROF(end,:),L2N,'LineWidth',2)
                hold on
                plot(PROF1(end,:),L2N1,'LineWidth',2,'Color',[0 .42 .71])
                plot(PROF1(end,ii),L2N1(ii),'ko','LineWidth',2)
        xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
        yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16);
        set(xaxis,'FontSize',20)
        set(yaxis,'FontSize',20)
        hold off



        subplot(2,2,[3,4])
        imagesc(x,y,ufull);
        yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
        xticks(0:20*pi:40*pi)
        xticklabels({'0','20\pi','40\pi'});
        yticks(0:pi:2*pi)
        yticklabels({'0','\pi','2\pi'});
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16);
        set(yaxis,'FontSize',28);
        daspect([1 1 1])
        drawnow;
        Frame = getframe(fig);
        writeVideo(vidfile,Frame);

      end


      for ii=1:length(PROF(end,:))
        znew = PROF(:,ii);
         ufull = sprof(x,y,c,znew,pl);

        fig = figure(600)
        subplot(2,2,1)
        plot(PROF(end,:),PROF(end-1,:),'LineWidth',2)
                hold on
                plot(PROF1(end,:),PROF1(end-1,:),'LineWidth',2,'Color',[0 .42 .71])
                plot(PROF(end,ii),PROF(end-1,ii),'ko','LineWidth',2)
        xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
        yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16);
        set(xaxis,'FontSize',20)
        set(yaxis,'FontSize',20)
        hold off


        subplot(2,2,2)
        plot(PROF(end,:),L2N,'LineWidth',2)
                hold on
                plot(PROF1(end,:),L2N1,'LineWidth',2,'Color',[0 .42 .71])
                plot(PROF(end,ii),L2N(ii),'ko','LineWidth',2)
        xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
        yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16);
        set(xaxis,'FontSize',20)
        set(yaxis,'FontSize',20)
        hold off



        subplot(2,2,[3,4])
        imagesc(x,y,ufull);
        yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
        xticks(0:20*pi:40*pi)
        xticklabels({'0','20\pi','40\pi'});
        yticks(0:pi:2*pi)
        yticklabels({'0','\pi','2\pi'});
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 16);
        set(yaxis,'FontSize',28);
        daspect([1 1 1])
        drawnow;
        Frame = getframe(fig);
        writeVideo(vidfile,Frame);

      end
      close(vidfile)
    end
