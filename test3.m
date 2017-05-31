%% 清空现场
clear;
clc;
%% 加载图像-1 & 预处理
RGB = imread('test1.bmp');
thisGray = rgb2gray(RGB);
thisEdge = edge(thisGray,'canny');
thisEdgeUint32 = uint32(thisEdge);
%%
[thisH,~,~] = fixedDCHTa(thisEdgeUint32,size(RGB,1),size(RGB,2),0,0,-90,90);
[H_DCHT,~,~]=voteDCHT(thisEdgeUint32);
cmp = int32(thisH)-int32(H_DCHT);
% %%
% times = 10;
% angle = floor(180*rand(times,2)-90);
% cmp = 0;
% for i = 1:times
%     
%     [thisH,~,~] = fixedDCHTa(thisEdgeUint32,size(RGB,1),size(RGB,2),0,0,angle(i,1),angle(i,2));
% 
%     %[H_DCHT,~,~]=voteDCHTa(thisEdgeUint32,angle(i,1),angle(i,2));
%     [H_DCHT,~,~]=voteDCHT(thisEdgeUint32);
%     
%     diff = int32(thisH)-int32(H_DCHT);
%     cmp = cmp +diff;
% end
%%
figure(1);
mesh(cmp);