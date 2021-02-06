% frames = [25,50,75,85,90,95,100,125,150,175,200,225,250,275,300];
% time1 = ones(length(frames),1);
% for j=1:length(frames)
%     sum = 0;
%     for i=1:frames(j)
%         sum = sum +t(i);
%     end
%     time1(j) = sum;
% end

%% plot times 
function PlotTimes(X1, YMatrix1)
% 创建 figure
figure1 = figure('Color',[1 1 1]);

% 创建 axes
axes1 = axes('Parent',figure1,'Position',[0.1 0.12 0.87 0.85]);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1(:,1:3),'LineWidth',1.5,'Parent',axes1);
plot2 = plot(X1(1:end-1),YMatrix1(1:end-1,4),X1(1:end-3),YMatrix1(1:end-3,5),'LineWidth',1.5,'Parent',axes1);
plot3 = plot(X1,YMatrix1(:,end),'LineWidth',1.5,'Parent',axes1);
plot4 = plot(X1(11),YMatrix1(11,4),'ko','MarkerFaceColor','k');
plot5 = plot(X1(9),YMatrix1(9,5),'ko','MarkerFaceColor','k');
set(plot1(1),'DisplayName','Garg');
set(plot1(2),'DisplayName','Jiang');
set(plot1(3),'DisplayName','Ren');
set(plot2(1),'DisplayName','Wei');
set(plot2(2),'DisplayName','Offline');
set(plot3(1),'DisplayName','Online');
set(plot4,'DisplayName','Out of memory');
set(plot5,'DisplayName','Out of memory');

ylabel('Times');
xlabel('Frames');
box(axes1,'on');
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,'Location','northwest','EdgeColor','none');
end

% figure1 = figure('Color',[1 1 1]);
% 
% % 创建 axes
% axes1 = axes('Parent',figure1,'Position',[0.1 0.12 0.87 0.85]);
% hold(axes1,'on');
% plot(time1(:,1),time1(:,2),'LineWidth',1.5)
% hold on
% plot(time1(:,1),time1(:,3),'LineWidth',1.5)
% hold on
% plot(time1(:,1),time1(:,4),'LineWidth',1.5)
% hold on
% plot(time1(1:end-2,1),time1(1:end-2,5),'LineWidth',1.5)
% hold on
% plot(time1(1:end-3,1),time1(1:end-3,6),'LineWidth',1.5)
% hold on 
% plot(time1(:,1),time1(:,7),'LineWidth',1.5)
% hold on
% plot(time1(10,1),time1(10,5),'ro','MarkerFaceColor','r')
% hold on
% plot(time1(9,1),time1(9,6),'go','MarkerFaceColor','g')
% 
% ylabel('Times');
% xlabel('Frames');
% 
% box(axes1,'on');
% legend1 = legend(axes1,'show');
% set(legend1,'Location','northwest','EdgeColor','none');
