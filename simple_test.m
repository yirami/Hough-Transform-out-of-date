%% 清空现场
clear;
clc;
%% 加载图像-1 & 预处理
RGB = imread('test1.bmp');
I  = rgb2gray(RGB);
img = edge(I,'canny');
figure(1);
%% 标准Hough
[ H,~,~ ] = hough_re( img );
subplot(1,2,1);
mesh(H);
%% 改进Hough
[ H,~,~ ] = hough_en( img );
subplot(1,2,2);
mesh(H);
%% 清空现场
clear;
clc;
%% 加载图像-2 & 预处理
RGB = imread('test2.bmp');
I  = rgb2gray(RGB);
img = edge(I,'canny');
figure(2);
%% 标准Hough
[ H,~,~ ] = hough_re( img );
subplot(1,2,1);
mesh(H);
%% 改进Hough
[ H,~,~ ] = hough_en( img );
subplot(1,2,2);
mesh(H);

