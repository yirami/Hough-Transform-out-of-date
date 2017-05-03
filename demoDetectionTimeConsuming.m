%% 清空现场
clear;
clc;
%% 测试图像列表
testDir = '.\DetectionTimeConsuming';
picList = dir(testDir);
for i = size(picList,1):-1:1
    if picList(i).isdir
        picList(i) = [];
    end
end
%% 统一测试参数
linesMax = 1000;
lineLengthPPHT = 50;
lineGapPPHT = 5;
lineLength = 50;
lineGap = 5;
timeConsuming = zeros(6,size(picList,1));
%% 逐个加载测试图片并采用三种方法检测
for i = 1:size(picList,1)
    thisPic = imread(strcat(picList(i).folder,'\',picList(i).name));
    thisPic = imnoise(thisPic,'salt & pepper',0.0001);
    thisGray = rgb2gray(thisPic);
    tic;
    thisEdge = edge(thisGray,'canny');
    timeConsuming(1,i)=toc;
    uint32thisEdge = uint32(thisEdge);
    %% 累计概率Hough（PPHT）
    tic;
    linesPPHT = votePPHT(uint32thisEdge,linesMax,lineLengthPPHT,lineGapPPHT);
    timeConsuming(2,i)=toc;
    %% 标准Hough（SHT）
    tic;
    [hSHT,thetaSHT,rhoSHT]=voteSHT(uint32thisEdge);
    timeConsuming(3,i)=toc;
    tic;
    peaksSHT = houghpeaks(double(hSHT),linesMax,'Threshold',100);
    linesSHT = houghlines(thisEdge,thetaSHT,rhoSHT,peaksSHT,'MinLength',lineLength,'FillGap',lineGap);
    timeConsuming(4,i)=toc;
    %% 采用方向编码的Hough（DCHT）
    tic;
    [hDCHT,thetaDCHT,rhoDCHT]=voteDCHT(uint32thisEdge);
    timeConsuming(5,i)=toc;
    tic;
    peaksDCHT = houghpeaks(double(hDCHT),linesMax,'Threshold',100);
    linesDCHT = houghlines(thisEdge,thetaDCHT,rhoDCHT,peaksDCHT,'MinLength',lineLength,'FillGap',lineGap);
    timeConsuming(6,i)=toc;
end
Time = mean((timeConsuming'));
figure(1);
set(gcf,'color','w');
subplot(4,1,1);
bar([0 Time(4) Time(6)],'FaceColor',[255 214 197]/255,'EdgeColor',[214 148 148]/255,'LineWidth',1.5);
set(gca,'xlim',[0 4]);set(gca,'ylim',[0 0.105]);
text(4.1,0.05,{'Peaks and Lines';'    Extraction'});
set(gca,'xticklabel',[]);set(gca,'box','off');
set(gca,'position',[0.15 0.71 0.65 0.2]);
subplot(4,1,2);
bar([Time(2) Time(3) Time(5)],'FaceColor',[102 205 170]/255,'EdgeColor',[125 255 180]/255,'LineWidth',1.5);
set(gca,'xlim',[0 4]);set(gca,'ylim',[0 0.035]);
text(4.3,0.02,'Voting');
set(gca,'xticklabel',[]);set(gca,'box','off');
set(gca,'position',[0.15 0.51 0.65 0.2]);
ylabel('Excution Time(s)');
subplot(4,1,3),
bar([Time(1) Time(1) Time(1)],'FaceColor',[169 30 223]/255,'EdgeColor',[1 0 1],'LineWidth',1.5);
set(gca,'xlim',[0 4]);set(gca,'ylim',[0 0.035]);
text(4.1,0.02,'Edge Detection');
set(gca,'xticklabel',[]);set(gca,'box','off');
set(gca,'position',[0.15 0.31 0.65 0.2]);
subplot(4,1,4),
bar([Time(1)+Time(2) Time(1)+Time(3)+Time(4) Time(1)+Time(5)+Time(6)],'FaceColor',[255 116 138]/255,'EdgeColor',[1 0 0],'LineWidth',1.5);
set(gca,'xlim',[0 4]);set(gca,'ylim',[0 0.17]);
text(4.3,0.1,'Total');
set(gca,'position',[0.15 0.11 0.65 0.2]);
set(gca, 'XTick', 0:4);set(gca,'box','off');
set(gca,'XTickLabel',{'','PPHT','SHT','DCHT',''});
xlabel('Algorithm');