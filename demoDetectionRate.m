%% 清空现场
clear;
clc;
%% 测试图像列表
testDir = '.\DirectionRate';
picList = dir(testDir);
for i = size(picList,1):-1:1
    if picList(i).isdir
        picList(i) = [];
    end
end
%% 统一测试参数
load('benchmark.mat');
load('config.mat');
linesMax = [5,10,20,40,80,200,1000,2000];
hitRate=zeros(3,length(linesMax));
hitrateT=zeros(3,size(picList,1));
tolerence=[2 10];
%% 数据收集
subdata = struct('linesMax',{},'linesPPHT',{},'linesSHT',{},'linesDCHT',{},'hitrateT',{});
dataDetectionRate = struct('picName',{},'thisPic',{},'thisGray',{},'thisEdge',{},'benchmark',{},'hSHT',{},'hDCHT',{},'subdata',subdata);
for t=1:length(linesMax)
    %% 逐个加载测试图片并采用三种方法检测
    for i = 1:size(picList,1)
        %% 准备基准线段参数
        for j=1:size(ref,2)
            if strcmp(ref(j).name,picList(i).name)
                Point = ref(j).Point;
                break;
            end
        end
        for j=1:size(Point,2)
                P=[0 0]';
                Q1=Point([1 2],j);
                Q2=Point([3 4],j);
            if (Q2(1)-Q1(1)>0)
                alpha=asind((Q2(2)-Q1(2))/sqrt((Q2(1)-Q1(1))^2+(Q2(2)-Q1(2))^2));
                rho=abs(det([Q2-Q1,P-Q1]))/norm(Q2-Q1);
                Point(5,j)=alpha;
                Point(6,j)=rho;
            else
                alpha=asind((Q1(2)-Q2(2))/sqrt((Q2(1)-Q1(1))^2+(Q2(2)-Q1(2))^2));
                rho=abs(det([Q2-Q1,P-Q1]))/norm(Q2-Q1);
                Point(5,j)=alpha;
                Point(6,j)=rho;
            end
        end
        dataDetectionRate(i).picName=picList(i).name;
        dataDetectionRate(i).benchmark=Point;
        dataDetectionRate(i).subdata(t).linesMax=linesMax(t);
        %% 提取当前图像的配置
        for j=1:size(ref,2)
            if strcmp(config(j).name,picList(i).name)
                lineLengthPPHT=config(j).lineLengthPPHT;
                lineGapPPHT=config(j).lineGapPPHT;
                lineLength=config(j).lineLength;
                lineGap=config(j).lineGap;
                threshold=config(j).threshold;      
                break;
            end
        end
        %% 读取测试图像
        thisPic = imread(strcat(picList(i).folder,'\',picList(i).name));
        thisGray = rgb2gray(thisPic);
        thisEdge = edge(thisGray,'canny');
        dataDetectionRate(i).thisPic=thisPic;
        dataDetectionRate(i).thisGray=thisGray;
        dataDetectionRate(i).thisEdge=thisEdge;
        figure(i);
        set(gcf,'color','w');
        %% 显示原图
        subplot(2,3,1);
        imshow(thisPic);
        title(strcat('当前参数：',num2str(linesMax(t))));
        %% 累计概率Hough（PPHT）
        linesPPHT = votePPHT(uint32(thisEdge),linesMax(t),lineLengthPPHT,lineGapPPHT);
        dataDetectionRate(i).subdata(t).linesPPHT=linesPPHT;   
        subplot(2,3,4); imshow(thisEdge);
        hold on;
        hitListPPHT=plotLines( linesPPHT,tolerence,Point,size(thisEdge) );
        title('PPHT的检测结果');
        %% 标准Hough（SHT）
        [hSHT,thetaSHT,rhoSHT]=voteSHT(uint32(thisEdge));
        peaksSHT = houghpeaks(double(hSHT),linesMax(t),'Threshold',double(0.1*max(hSHT(:))));
        linesSHT = houghlines(thisEdge,thetaSHT,rhoSHT,peaksSHT,'MinLength',lineLength,'FillGap',lineGap);
        dataDetectionRate(i).hSHT=hSHT;
        dataDetectionRate(i).subdata(t).linesSHT=linesSHT;
        subplot(2,3,2); mesh(hSHT);
        title('SHT的Hough空间'); xlabel('\theta'); ylabel('\rho');zlabel('Voting number');
        subplot(2,3,5); imshow(thisEdge);
        hold on;
        hitListSHT=plotLines( linesSHT,tolerence,Point,size(thisEdge) );
        title('SHT的检测结果');
        %% 采用方向编码的Hough（DCHT）
        [hDCHT,thetaDCHT,rhoDCHT]=voteDCHT(uint32(thisEdge));
        peaksDCHT = houghpeaks(double(hDCHT),linesMax(t),'Threshold',double(0.1*max(hDCHT(:))));
        linesDCHT = houghlines(thisEdge,thetaDCHT,rhoDCHT,peaksDCHT,'MinLength',lineLength,'FillGap',lineGap);
        dataDetectionRate(i).hDCHT=hDCHT;
        dataDetectionRate(i).subdata(t).linesDCHT=linesDCHT;
        subplot(2,3,3); mesh(hDCHT);
        title('DCHT的Hough空间'); xlabel('\theta'); ylabel('\rho');zlabel('Voting number');
        subplot(2,3,6); imshow(thisEdge);
        hold on;
        hitListDCHT=plotLines( linesDCHT,tolerence,Point,size(thisEdge) );
        title('DCHT的检测结果');
        %% 计算命中率
        hitrateT(1,i)=length(hitListPPHT)/size(Point,2);
        hitrateT(2,i)=length(hitListSHT)/size(Point,2);
        hitrateT(3,i)=length(hitListDCHT)/size(Point,2);
        dataDetectionRate(i).subdata(t).hitrateT=hitrateT(:,i);
    end
    hitRate(1,t)=mean(hitrateT(1,:));
    hitRate(2,t)=mean(hitrateT(2,:));
    hitRate(3,t)=mean(hitrateT(3,:));
end
figure(10);
set(gcf,'color','w');
Y=0.5:0.5:1.5;
c=[255 0 0;128 0 128;214 148 148;131 18 37;0 191 138;79 144 193;255 116 208;227 159 89]/255;
plot(hitRate(:,1),Y,'-d','color',c(1,:),'linewidth',2,'markersize',10);hold on;
plot(hitRate(:,2),Y,'-s','color',c(2,:),'linewidth',2,'markersize',10);hold on;
plot(hitRate(:,3),Y,'-*','color',c(3,:),'linewidth',2,'markersize',10);hold on;
plot(hitRate(:,4),Y,'-+','color',c(4,:),'linewidth',2,'markersize',10);hold on;
plot(hitRate(:,5),Y,'-p','color',c(5,:),'linewidth',2,'markersize',10);hold on;
plot(hitRate(:,6),Y,'-h','color',c(6,:),'linewidth',2,'markersize',10);hold on;
plot(hitRate(:,7),Y,'-o','color',c(7,:),'linewidth',2,'markersize',10);hold on;
plot(hitRate(:,8),Y,':x','color',c(8,:),'linewidth',2,'markersize',10);
legend('N_{l}=5','N_{l}=10','N_{l}=20','N_{l}=40','N_{l}=80','N_{l}=200','N_{l}=1000','N_{l}\rightarrow\infty','Location',[0 0.6 0.1 0.3]);
axis([0 1 0 2]);
xlabel('Detection Rate');
ylabel('Algorithm');
set(gca, 'XTick', 0:0.1:1);
set(gca, 'YTick', 0:0.5:2);
set(gca,'FontSize',12);
set(gca,'YTickLabel',{'','PPHT','SHT','DCHT',''});
set(gca,'position',[0.15 0.16 0.71 0.8]);
