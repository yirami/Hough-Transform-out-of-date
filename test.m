%% 清空现场
clear;
clc;
%% 加载图像-1 & 预处理
RGB = imread('test1.bmp');
I  = rgb2gray(RGB);
img = edge(I,'canny');
img_=uint32(img);

% [H1,~,~]=voteSHT(img_);
% figure(1);
% mesh(H1);

% [H1,~,~]=voteSHTa(img_,20,60);
% figure(1);
% mesh(H1);

% [H_DCHT,theta_DCHT,rho_DCHT]=voteDCHT(img_);
% figure(31);
% mesh(H_DCHT);
% peaks = houghpeaks(double(H_DCHT),8,'Threshold',double(0.3*max(H_DCHT(:))));
% lines = houghlines(img,theta_DCHT,rho_DCHT,peaks);
% figure(32);
% imshow(img);
% hold on;
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end

% [H_DCHT,theta_DCHT,rho_DCHT]=voteDCHTa(img_,20,60);
% figure(31);
% mesh(H_DCHT);
% peaks = houghpeaks(double(H_DCHT),8,'Threshold',double(0.3*max(H_DCHT(:))));
% lines = houghlines(img,theta_DCHT,rho_DCHT,peaks);
% figure(32);
% imshow(img);
% hold on;
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end

% lines = votePPHT(img_,10,50,5);
% figure(4);
% imshow(img);
% hold on;
% for i=1:size(lines,2)
%     xy = [lines(2,i),lines(1,i);lines(4,i),lines(3,i)];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end

% lines = votePPHTa(img_,10,50,5,20,60);
% figure(4);
% imshow(img);
% hold on;
% for i=1:size(lines,2)
%     xy = [lines(2,i),lines(1,i);lines(4,i),lines(3,i)];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end

% lines = voteDCPPHT(img_,10,50,5);
% figure(4);
% imshow(img);
% hold on;
% for i=1:size(lines,2)
%     xy = [lines(2,i),lines(1,i);lines(4,i),lines(3,i)];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end

lines = voteDCPPHTa(img_,10,50,5,20,60);
figure(4);
imshow(img);
hold on;
for i=1:size(lines,2)
    xy = [lines(2,i),lines(1,i);lines(4,i),lines(3,i)];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end

% [H_SHT,theta_SHT,rho_SHT]=hough( img );
% figure(21);
% mesh(H_SHT);
% peaks = houghpeaks(H_SHT,6,'Threshold',0.3*max(H_SHT(:)));
% lines = houghlines(img,theta_SHT,rho_SHT,peaks);
% figure(22);
% imshow(img);
% hold on;
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end

