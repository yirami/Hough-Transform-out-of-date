%% ����ֳ�
clear;
clc;
%% ����ͼ��-1 & Ԥ����
RGB = imread('test1.bmp');
I  = rgb2gray(RGB);
img = edge(I,'canny');
figure(1);
%% ��׼Hough
[ H,~,~ ] = hough_re( img );
subplot(1,2,1);
mesh(H);
%% �Ľ�Hough
[ H,~,~ ] = hough_en( img );
subplot(1,2,2);
mesh(H);
%% ����ֳ�
clear;
clc;
%% ����ͼ��-2 & Ԥ����
RGB = imread('test2.bmp');
I  = rgb2gray(RGB);
img = edge(I,'canny');
figure(2);
%% ��׼Hough
[ H,~,~ ] = hough_re( img );
subplot(1,2,1);
mesh(H);
%% �Ľ�Hough
[ H,~,~ ] = hough_en( img );
subplot(1,2,2);
mesh(H);

