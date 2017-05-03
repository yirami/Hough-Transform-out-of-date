%% «Âø’œ÷≥°
clear;
clc;
%% Detectionrate
load dataDetectionRate
tolerence=[2 10];
    % Õº1
figure(1);
set(gcf,'color','w');
num = length(dataDetectionRate(1).subdata);
hitRate = zeros(3,num);
for i = 1:num
    temp=[];
    for j = 1:length(dataDetectionRate)
        temp = [temp dataDetectionRate(j).subdata(i).hitrateT];
    end
    hitRate(:,i) = mean(temp,2);
end
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
legend('N_{l}=5','N_{l}=10','N_{l}=20','N_{l}=40','N_{l}=80','N_{l}=200','N_{l}=1000','N_{l}\rightarrow\infty','northeast');
axis([0 1 0 2]);
xlabel('Detection Rate');ylabel('Algorithm');
set(gca, 'XTick', 0:0.1:1);set(gca, 'YTick', 0:0.5:2);
set(gca,'FontSize',12);set(gca,'FontName','Times new roman');
set(gca,'YTickLabel',{'','PPHT','SHT','DCHT',''});
set(gca,'position',[0.15 0.16 0.71 0.8]);
    % Õº2
figure(2);
Position = zeros(12,4);
Position(:,3) = 0.33; Position(:,4) = 0.22;
Position(:,1) = repmat([0 0.335 0.67],1,4);
Position(:,2) = [0.78 0.78 0.78 0.53 0.53 0.53 0.28 0.28 0.28 0.03 0.03 0.03];
txtPosition = zeros(12,2);
txtPosition(:,1) = ones(12,1)*floor(size(dataDetectionRate(1).thisEdge,2)*0.45);
txtPosition(:,2) = ones(12,1)*floor(size(dataDetectionRate(1).thisEdge,1)*1.05);
set(gcf,'color','w');
subplot('position',Position(1,:));
imshow(dataDetectionRate(1).thisPic);
text(txtPosition(1,1),txtPosition(1,2),'( a )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(2,:));
imshow(dataDetectionRate(1).thisEdge);
text(txtPosition(2,1),txtPosition(2,2),'( b )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(3,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
Point = dataDetectionRate(1).benchmark;
for j = 1:size(Point,2)
   xy = [Point([1 2],j)'; Point([3 4],j)'];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end
text(txtPosition(3,1),txtPosition(3,2),'( c )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(4,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(1).linesPPHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(4,1),txtPosition(4,2),'( d )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(5,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(1).linesSHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(5,1),txtPosition(5,2),'( e )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(6,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(1).linesDCHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(6,1),txtPosition(6,2),'( f )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(7,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(3).linesPPHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(7,1),txtPosition(7,2),'( g )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(8,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(3).linesSHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(8,1),txtPosition(8,2),'( h )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(9,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(3).linesDCHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(9,1),txtPosition(9,2),'( i )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(10,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(6).linesPPHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(10,1),txtPosition(10,2),'( j )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(11,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(6).linesSHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(11,1),txtPosition(11,2),'( k )','FontSize',18,'FontName','Times new roman');
subplot('position',Position(12,:));
imshow(dataDetectionRate(1).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(1).subdata(6).linesDCHT,tolerence,Point,size(dataDetectionRate(1).thisEdge) );
text(txtPosition(12,1),txtPosition(12,2),'( l )','FontSize',18,'FontName','Times new roman');
   % Õº3
figure(3);
set(gcf,'color','w');
subplot('position',[0.1 0.1 0.4 0.8]);
mesh(dataDetectionRate(1).hSHT);
subplot('position',[0.6 0.1 0.4 0.8]);
mesh(dataDetectionRate(1).hDCHT);
    % Õº4
figure(4);
set(gcf,'color','w');
Point = dataDetectionRate(2).benchmark;
txtPosition = zeros(5,2);
txtPosition(:,1) = ones(5,1)*floor(size(dataDetectionRate(2).thisEdge,2)*0.45);
txtPosition(:,2) = ones(5,1)*floor(size(dataDetectionRate(2).thisEdge,1)*1.05);
subplot('position',[0 0.4 0.5 0.6]);
imshow(dataDetectionRate(2).thisPic);
text(txtPosition(1,1),txtPosition(1,2),'( a )','FontSize',18,'FontName','Times new roman');
subplot('position',[0.5 0.4 0.5 0.6]);
imshow(dataDetectionRate(2).thisEdge);
text(txtPosition(2,1),txtPosition(2,2),'( b )','FontSize',18,'FontName','Times new roman');
subplot('position',[0 0 0.33 0.4]);
imshow(dataDetectionRate(2).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(2).subdata(4).linesPPHT,tolerence,Point,size(dataDetectionRate(2).thisEdge) );
text(txtPosition(3,1),txtPosition(3,2),'( c )','FontSize',18,'FontName','Times new roman');
subplot('position',[0.335 0 0.33 0.4]);
imshow(dataDetectionRate(2).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(2).subdata(4).linesSHT,tolerence,Point,size(dataDetectionRate(2).thisEdge) );
text(txtPosition(4,1),txtPosition(4,2),'( d )','FontSize',18,'FontName','Times new roman');
subplot('position',[0.67 0 0.33 0.4]);
imshow(dataDetectionRate(2).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(2).subdata(4).linesDCHT,tolerence,Point,size(dataDetectionRate(2).thisEdge) );
text(txtPosition(5,1),txtPosition(5,2),'( e )','FontSize',18,'FontName','Times new roman');
    % Õº5
figure(5);
set(gcf,'color','w');
Point = dataDetectionRate(3).benchmark;
txtPosition = zeros(5,2);
txtPosition(:,1) = ones(5,1)*floor(size(dataDetectionRate(3).thisEdge,2)*0.45);
txtPosition(:,2) = ones(5,1)*floor(size(dataDetectionRate(3).thisEdge,1)*1.05);
subplot('position',[0 0.4 0.5 0.6]);
imshow(dataDetectionRate(3).thisPic);
text(txtPosition(1,1),txtPosition(1,2),'( a )','FontSize',18,'FontName','Times new roman');
subplot('position',[0.5 0.4 0.5 0.6]);
imshow(dataDetectionRate(3).thisEdge);
text(txtPosition(2,1),txtPosition(2,2),'( b )','FontSize',18,'FontName','Times new roman');
subplot('position',[0 0 0.33 0.4]);
imshow(dataDetectionRate(3).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(3).subdata(4).linesPPHT,tolerence,Point,size(dataDetectionRate(3).thisEdge) );
text(txtPosition(3,1),txtPosition(3,2),'( c )','FontSize',18,'FontName','Times new roman');
subplot('position',[0.335 0 0.33 0.4]);
imshow(dataDetectionRate(3).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(3).subdata(4).linesSHT,tolerence,Point,size(dataDetectionRate(3).thisEdge) );
text(txtPosition(4,1),txtPosition(4,2),'( d )','FontSize',18,'FontName','Times new roman');
subplot('position',[0.67 0 0.33 0.4]);
imshow(dataDetectionRate(3).thisEdge);hold on;
[~]=plotLines( dataDetectionRate(3).subdata(4).linesDCHT,tolerence,Point,size(dataDetectionRate(3).thisEdge) );
text(txtPosition(5,1),txtPosition(5,2),'( e )','FontSize',18,'FontName','Times new roman');
%% TimeConsuming
load dataTimeConsuming
figure(6);
set(gcf,'color','w');
subplot(4,1,1);
bar([0 Time(4) Time(6)],'FaceColor',[255 214 197]/255,'EdgeColor',[214 148 148]/255,'LineWidth',1.5);
set(gca,'xlim',[0 4]);set(gca,'ylim',[0 0.105]);
text(4.1,0.05,{'Peaks and Lines';'    Extraction'},'FontSize',12,'FontName','Times new roman');
set(gca,'xticklabel',[]);set(gca,'box','off');
set(gca,'position',[0.15 0.71 0.65 0.2]);
subplot(4,1,2);
bar([Time(2) Time(3) Time(5)],'FaceColor',[102 205 170]/255,'EdgeColor',[125 255 180]/255,'LineWidth',1.5);
set(gca,'xlim',[0 4]);set(gca,'ylim',[0 0.035]);
text(4.3,0.02,'Voting','FontSize',12,'FontName','Times new roman');
set(gca,'xticklabel',[]);set(gca,'box','off');
set(gca,'position',[0.15 0.51 0.65 0.2]);
ylabel('Excution Time (s)','FontSize',12,'FontName','Times new roman');
subplot(4,1,3),
bar([Time(1) Time(1) Time(1)],'FaceColor',[169 30 223]/255,'EdgeColor',[1 0 1],'LineWidth',1.5);
set(gca,'xlim',[0 4]);set(gca,'ylim',[0 0.035]);
text(4.1,0.02,'Edge Detection','FontSize',12,'FontName','Times new roman');
set(gca,'xticklabel',[]);set(gca,'box','off');
set(gca,'position',[0.15 0.31 0.65 0.2]);
subplot(4,1,4),
bar([Time(1)+Time(2) Time(1)+Time(3)+Time(4) Time(1)+Time(5)+Time(6)],'FaceColor',[255 116 138]/255,'EdgeColor',[1 0 0],'LineWidth',1.5);
set(gca,'xlim',[0 4]);set(gca,'ylim',[0 0.17]);
text(4.3,0.1,'Total','FontSize',12,'FontName','Times new roman');
set(gca,'position',[0.15 0.11 0.65 0.2]);
set(gca, 'XTick', 0:4);set(gca,'box','off');
set(gca,'XTickLabel',{'','PPHT','SHT','DCHT',''},'FontName','Times new roman');
xlabel('Algorithm','FontSize',12,'FontName','Times new roman');
