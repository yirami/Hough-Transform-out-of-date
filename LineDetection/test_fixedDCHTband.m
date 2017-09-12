%% 清空现场
clear;
clc;
%% 加载图像-1 & 预处理
RGB = imread('test4.bmp');
[~,rect] = imcrop(RGB);
angleSpec = [-50,50];
tic;
[lines] = fixedDCHTband( RGB,rect,angleSpec,4,5 );
t = toc
figure(1);
imshow(rgb2gray(RGB));
hold on;
for k = 1:length(lines)
    xy = [lines(k).pointStart;lines(k).pointEnd];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
end