%%%%isola, fixed c above c~.5, continue in ky

clear all
close all

load('ContDataPerp5.400000e-01.mat')
%load('ContDataPerp5.311000e-01.mat')
cc = .54;
pL2 = [];
for ii = 1:length(PROF(end,:))
    pL2 = [pL2,norm(PROF(1:end-2,ii),2)];
end
pPROF = PROF;

figure(1)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on 
plot(pPROF(end,:),pL2(:),'LineWidth',2)
xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
hold off


% 



%%Figures with solution labels
pI = [1 80];
figure(2)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on
for ii = 1:2
 ufull = sprof(x,y,cc,PROF(:,pI(ii)));
subplot(2,1,ii)
imagesc(x,y,ufull)
    xaxis=xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 20);
    yaxis=ylabel('$y$', 'Interpreter', 'latex', 'FontSize', 20);
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16);
    set(xaxis,'FontSize',20)
    set(yaxis,'FontSize',20)
    title(['(',num2str(ii),')'],'Interpreter','latex','FontSize',20,'Position',[130 5])
end
hold off
drawnow

figure(3)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on 
plot(pPROF(end,:),pPROF(end-1,:),'LineWidth',2)
plot(pPROF(end,pI),pPROF(end-1,pI),'ko','LineWidth',2)
for ii = 1:2
    text(pPROF(end,pI(ii)),pPROF(end-1,pI(ii))+0.005,['(',num2str(ii),')'],'Interpreter', 'latex','FontSize',18);
end
xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$k_x$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
ylim([-0.01 0.12])
hold off


figure(4)
box on
set(gcf, 'Position', [400, 400, 900, 600]);
set(gcf,'PaperPositionMode','auto')
hold on 
plot(pPROF(end,:),pL2(:),'LineWidth',2)
plot(pPROF(end,pI),pL2(pI),'ko','LineWidth',2)
for ii = 1:2
    text(pPROF(end,pI(ii)),pL2(pI(ii))+0.015,['(',num2str(ii),')'],'Interpreter', 'latex','FontSize',18);
end
xaxis=xlabel('$k_y$', 'Interpreter', 'latex', 'FontSize', 20);
yaxis=ylabel('$\| w\|_{L^2}$', 'Interpreter', 'latex', 'FontSize', 20);
xt = get(gca, 'XTick');
set(gca, 'FontSize', 16);
set(xaxis,'FontSize',20)
set(yaxis,'FontSize',20)
hold off


% 
% load('ContData5.397959e-01.mat')
% oL2 = [];
% for ii = 1:length(PROF(end,:))
%     oL2 = [oL2,norm(PROF(1:end-2,ii),2)];
% end
% oPROF = PROF;
% 
% 
% figure(1)
% subplot(1,2,1)
% plot(pPROF(end,:),pPROF(end-1,:),oPROF(end,:),oPROF(end-1,:))
% xlabel('k_y')
% ylabel('k_x')
% subplot(1,2,2)
% plot( pPROF(end,:),pL2,oPROF(end,:),oL2)
% xlabel('c')
% ylabel('| w|_{L^2}')
% %title('c = ', num2str(.535))