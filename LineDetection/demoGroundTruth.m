%% ����ֳ�
clear;
clc;
%% ����ͼ���б�
testDir = '.\DirectionRate';
picList = dir(testDir);
for i = size(picList,1):-1:1
    if picList(i).isdir
        picList(i) = [];
    end
end
load('benchmark.mat');
for i = 1:size(picList,1)
    thisPic = imread(strcat(picList(i).folder,'\',picList(i).name));
    thisGray = rgb2gray(thisPic);
    thisEdge = edge(thisGray,'canny');
    % ƥ�伯�л�׼�߶ε��ֶ�
    for j=1:size(ref,2)
        if strcmp(ref(j).name,picList(i).name)
            Point = ref(j).Point;
            break;
        end
    end
    % ��ʾԭͼ
    figure(3*(i-1)+1);
    set(gcf,'color','w');
    imshow(thisPic);
    % ��ʾ��Եͼ
    figure(3*(i-1)+2);
    set(gcf,'color','w');
    imshow(thisEdge);
    % ���ƻ�׼
    figure(3*(i-1)+3);
    set(gcf,'color','w');
    imshow(thisEdge);
    hold on;
    for j = 1:length(Point)
       xy = [Point([1 2],j)'; Point([3 4],j)'];
       plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
    end
end