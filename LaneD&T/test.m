%% 清空现场
clear;
clc;
%% 加载图像-1 & 预处理
RGB = imread('test4.bmp');
frameSize = [size(RGB,2) size(RGB,1)];
refPoint = [floor(frameSize(1)/2),frameSize(2)];
tic;
[~,linesG] = DCHTbandGlobal( RGB );
[ leftOK,rightOK,twoLaneG ] = laneMatchForSearch( linesG,frameSize,refPoint);
t = toc
figure(1);
imshow(RGB);
hold on;
% if (leftOK)
%     xy = [twoLaneG(1).startPoint;twoLaneG(1).endPoint];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end
% if (rightOK)
%     xy = [twoLaneG(2).startPoint;twoLaneG(2).endPoint];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end
deadline = floor(frameSize(2)/2);
[ trackKeyPoint ] = getTrackKeyPoint( twoLaneG,deadline,[size(RGB,2),size(RGB,1)] );
[ drawKeyPoint ] = getDrawKeyPoint( trackKeyPoint,deadline,deadline,[size(RGB,2),size(RGB,1)] );
% 绘制追踪参数
plot([1 frameSize(1)],[deadline deadline],'g-.','LineWidth',2);hold on;
plot([1 frameSize(1)],[frameSize(2) frameSize(2)],'g-.','LineWidth',2);hold on;
plot([trackKeyPoint(1) 1.2*(trackKeyPoint(2)-trackKeyPoint(1))+trackKeyPoint(1)],[frameSize(2) 1.2*(deadline-frameSize(2))+frameSize(2)],'m','LineWidth',2);hold on;
plot([trackKeyPoint(3) 1.2*(trackKeyPoint(4)-trackKeyPoint(3))+trackKeyPoint(3)],[frameSize(2) 1.2*(deadline-frameSize(2))+frameSize(2)],'m','LineWidth',2);hold on;
plot([trackKeyPoint(1) trackKeyPoint(2)],[frameSize(2) deadline],'pb','MarkerSize',16,'MarkerFaceColor','b');hold on;
plot([trackKeyPoint(3) trackKeyPoint(4)],[frameSize(2) deadline],'pb','MarkerSize',16,'MarkerFaceColor','b');hold on;

%
shape2RGB = insertShape(RGB,'FilledPolygon',drawKeyPoint,'Color',[0	255	153],'Opacity',0.2);
figure(2);
imshow(shape2RGB);
