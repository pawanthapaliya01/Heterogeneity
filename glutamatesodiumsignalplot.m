clc
close all
clear 
%format compact
load Nak1.dat
load Nak2.dat
load Nak3.dat
load Nak4.dat
%load d1.csv
%load d2.csv
%load d3.csv
%load d4.csv
%load Naknmda.dat
%load Naknonmda.dat
load naexp1.csv
load naexp2.csv
load naexp3.csv
load naexp4.csv

X=categorical({'Glu','dap','nbqx+dap','dap+nbqx+tbo'});
X = reordercats(X,{'Glu','dap','nbqx+dap','dap+nbqx+tbo'});

Y = [(max(Nak1(:,2))-14368), 5600;(max(Nak2(:,2))-14368), 2822;(max(Nak3(:,2))-14368), 2000;...
(max(Nak4(:,2))-14368), 2.58];
err=[0,600;0, 100;0, 100;0, 30]*1.0d-3;
ylow = [600 100 100 30]*1.0d-3;

figure
%subplot(1,2,1)
%  plot(Nak1(:,1)-min(Nak1(:,1)),Nak1(:,2)*1.0d-3,'-r','DisplayName','WT ','LineWidth',2)
%  hold on
%  
% plot(Nak2(:,1)-min(Nak2(:,1)),Nak2(:,2)*1.0d-3,'-m','DisplayName','APV','LineWidth',2)
% % 
% 
% plot(Nak3(:,1)-min(Nak3(:,1)),Nak3(:,2)*1.0d-3,'-c','DisplayName','APV+NBQX','LineWidth',2)
%  
% plot(Nak4(:,1)-min(Nak4(:,1)),Nak4(:,2)*1.0d-3,'-k','DisplayName','APV+NBQX+TBOA','LineWidth',2)
%  
% 
% plot(naexp1(:,1)-min(naexp1(:,1))+4.9,naexp1(:,2),'o','MarkerSize',.1,...
%     'MarkerEdgeColor','r',...
%     'LineWidth',2,...
%     'DisplayName','EXP WT')
% % plot(d1(7:end,1),d1(7:end,2),'o','MarkerSize',5,...
% %     'MarkerEdgeColor','r',...
% %     'LineWidth',2,...
% %     'DisplayName','EXP WT')
%  plot(naexp2(:,1)-min(naexp2(:,1))+4.9,naexp2(:,2),'o','MarkerSize',.1,...
%     'MarkerEdgeColor','m',...
%     'LineWidth',2,...
%     'DisplayName','EXP APV')
% 
% % plot(d2(:,1),d2(:,2),'o','MarkerSize',5,...
% %     'MarkerEdgeColor','m',...
% %     'LineWidth',2,...
% %     'DisplayName','EXP APV')
%  plot(naexp3(:,1)-min(naexp3(:,1))+4.9,naexp3(:,2),'o','MarkerSize',.1,...
%     'MarkerEdgeColor','c',...
%    'LineWidth',2,...
%     'DisplayName','EXP APV+NBQX')
% % plot(d3(:,1),d3(:,2),'o','MarkerSize',5,...
% %     'MarkerEdgeColor','c',...
% %    'LineWidth',2,...
% %     'DisplayName','EXP APV+NBQX')
%  plot(naexp4(:,1)-min(naexp4(:,1))+4.9,naexp4(:,2),'o','MarkerSize',.1,...
%     'MarkerEdgeColor','k',...
%     'LineWidth',2,...
%     'DisplayName','EXP APV+NBQX+TBOA')
% % plot(d4(:,1),d4(:,2),'o','MarkerSize',5,...
% %     'MarkerEdgeColor','k',...
% %     'LineWidth',2,...
% %     'DisplayName','EXP APV+NBQX+TBOA')
%  
% ax.YAxis.Exponent = 0;
%  
% xlabel('Time(s)','FontName','Times New Roman','FontSize',12,'FontWeight','bold','Color','k')
% ylabel('[Na]_a(mM)','FontName','Times New Roman','FontSize',12,'FontWeight','bold','Color','k')
% ytickformat('%.0f')
% xtickformat('%.0f')
% legend('FontName','Times New Roman','FontSize',12,'FontWeight','bold','LineWidth',2)
% legend 'boxoff'
% %title('Astrocyte cortex','FontName','Times New Roman','FontSize',12,'FontWeight','bold','Color','black')
% ax=gca;
% %a = get(gca,'XTickLabel');
% %set(gca,'XTickLabel',a,'FontName','Times New Roman','fontsize',18','FontWeight','bold','LineWidth',2.5)
% ax.XAxis.FontSize =12;
% ax.XAxis.FontWeight = 'bold';
% ax.XAxis.Color='black';
% ax.XAxis.LineWidth=2;
% ax.XAxis.FontName='Times new Roman';
% ax.YAxis.LineWidth =2;
% ax.YAxis.FontWeight='bold';
% ax.YAxis.FontSize=12;
% ax.YAxis.FontName='Times new Roman';
% ax.YAxis.Color='black';
% ax.XLim=[0,60];
% ax.YLim=[12,22];
% xticks(0:10:60);
% yticks(12:2:22);
% set(gca, 'box', 'off')
% % hold off

% 
% % subplot(1,2,2)

b=bar( Y*1.0d-3, 'grouped','BarWidth',1,'facecolor','flat','LineWidth',2);
%set(b(1),'FaceColor','r')
b(1).CData(1,:) =[1 1 1];
b(2).CData(1,:) = [17 17 17]; % group 1 2nd bar
b(1).CData(2,:)=[1 1 1];
b(2).CData(2,:)=[17 17 17];
b(1).CData(3,:)=[1 1 1];
b(2).CData(3,:)=[17 17 17];
b(1).CData(4,:)=[1 1 1];
b(2).CData(4,:)=[17 17 17];
b(2).CData(1,:)=[1 0 0];
b(2).CData(2,:)=[1 0 0];
b(2).CData(3,:)=[1 0 0];
b(2).CData(4,:)=[1 0 0];
%set(b(1,1),'FaceColor','r','EdgeColor','None')
%set(b(1,2),'FaceColor','None','EdgeColor','r')

LineArray={ '-' , '-' };
for k=1:2
  set(b(k),'LineStyle',LineArray{k},'LineWidth',2)
end
hold on
[ngroups,nbars] = size(Y);
% Calculating the width for each bar group
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
    %errorbar(x',Y*1.0d-3,err,'k','linestyle','none');
end
%Plot the errorbars
errorbar(x',Y*1.0d-3,err,'k','linestyle','none','MarkerSize',1.5,'LineWidth',1);
hold off

er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
ytickformat('%.0f')

legend('SIM','EXP','FontName','Times New Roman','FontSize',12,'FontWeight','bold');

legend 'boxoff'




ylabel('[\DeltaNa]_a(mM)','FontName','Times New Roman','FontSize',12,'FontWeight','bold','Color','k')
ax=gca;
ax.XAxis.FontSize =12;
ax.XAxis.FontWeight = 'bold';
ax.XAxis.Color='k';
ax.XAxis.FontName='Times new Roman';
ax.XAxis.LineWidth=3;

 a = get(gca,'XTickLabel');
ax.YAxis.LineWidth =2;
ax.YAxis.FontWeight='bold';
ax.YAxis.FontSize=14;
ax.YAxis.FontName='Times new Roman';
ax.YAxis.Color='black';
ax.YLim=[0,8.02];
yticks(0:2:8.02);

set(gca,'XTickLabel',a,'FontName','Times New Roman','FontSize',12,'FontWeight','bold','LineWidth',2)
set(gca, 'XTickLabel', {'Wild Type','APV','NBQX+APV','APV+NBQX+TBOA'});
set(gca, 'box', 'off')
set(gca,'XTickLabelRotation',10)
% saveas(gcf,'Glusodcortex','epsc')
% % 
