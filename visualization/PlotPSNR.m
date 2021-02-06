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
function PlotPSNR(X1, i1,g1,j1,r1,s1,y1,our1)
% ���� figure
figure1 = figure('Color',[1 1 1]);

% ���� axes
axes1 = axes('Parent',figure1,'Position',[0.1 0.12 0.87 0.85]);
hold(axes1,'on');

w=4;
% ʹ�� plot �ľ������봴������
%figure(1);
set(gca, 'Position', [0.07, 0.08, 0.9, 0.9])
plot(X1,[i1(4:112),s1(2:110),y1,our1(2:110)],'LineWidth',w);
axis tight;

% plot2 = plot(X1,g1,'LineWidth',w);
%plot3 = plot(X1,j1,'LineWidth',w);
%plot4 = plot(X1,r1,'LineWidth',w);
%plot5 = plot(X1,s1,'LineWidth',w);
%plot6 = plot(X1,our1,'LineWidth',w);

xlim([1 120])
ylim([22 36])
set(gca, 'XTick',[1:40:170])
set(gca, 'YTick', [24:4:40])
set(gca, 'fontsize', 22, 'Xcolor', 'black', 'Ycolor','black')
set(gca, 'Color', 'white')
box off
text(55,34.5,'a2', 'fontsize', 30, 'color', 'black')
text(45,23,'Frame No.', 'fontsize', 22, 'color', 'black')
text(-8,35.5,'(dB)', 'fontsize', 22, 'color', 'black')

set(gcf, 'color', 'w');

% ylabel('PSNR(dB)', 'fontsize',28);
% xlabel('Frame No.', 'fontsize', 28);
% box(axes1,'on');
% ���� legend
legend1 = legend(axes1,'show');
% set(legend1,'Location','northwest','EdgeColor','none');
im=frame2im(getframe(gcf));
imshow(im,'border','tight', 'initialmagnification', 'fit');
imwrite(im, 'a2.jpg')
end

% figure1 = figure('Color',[1 1 1]);
% 
% % ���� axes
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
