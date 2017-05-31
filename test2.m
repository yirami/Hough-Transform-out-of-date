%% 清空现场
clear;
clc;
%% 加载图像-1 & 预处理
RGB = imread('test1.bmp');
[~,rect] = imcrop(RGB);

% tic;
[lines] = fixedDCHTband( RGB,rect,-90,90,2,5 );
% t1=toc;
figure(1);
imshow(rgb2gray(RGB));
hold on;
for k = 1:length(lines)
    xy = [lines(k).pointStart;lines(k).pointEnd];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end

% tic;
% I  = rgb2gray(RGB);
% img = edge(I,'canny');
% img_=uint32(img);
%[H_DCHT,theta_DCHT,rho_DCHT]=voteDCHT(img_);
% [H_DCHT,theta_DCHT,rho_DCHT]=voteDCHTa(img_,-90,90);
% [linesInfo] = searchLines(img_,H_DCHT,theta_DCHT,rho_DCHT,10,5);
% t2=toc;
% figure(2);
% imshow(I);
% hold on;
% for k = 1:size(linesInfo,2)
%     rows = size(img,1);
%     xy = [linesInfo(2,k),linesInfo(1,k);linesInfo(4,k),linesInfo(3,k)];
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% end