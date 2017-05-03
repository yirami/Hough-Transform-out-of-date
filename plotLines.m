function [ hitList ] = plotLines( lines,tolerence,Point,sizeP )
%plotLines ���Ƽ�����ֱ��
%   ����������е�ֱ�ߵ�����
hitList=[];
for k=1:size(lines,2)
    if isstruct(lines)
        xy = [lines(k).point1; lines(k).point2];
    else
        xy = [lines(2,k),lines(1,k);lines(4,k),lines(3,k)];
    end
    xy=double(xy);
    %   ���з��ء���ɫ�������򡰺�ɫ��
    P=[0 0];
    Q1=double(xy(1,:));
    Q2=double(xy(2,:));
    if (Q2(1)-Q1(1)>0)
        alpha=asind((Q2(2)-Q1(2))/sqrt((Q2(1)-Q1(1))^2+(Q2(2)-Q1(2))^2));
        rho=abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);
    else
        alpha=asind((Q1(2)-Q2(2))/sqrt((Q2(1)-Q1(1))^2+(Q2(2)-Q1(2))^2));
        rho=abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1);
    end
%% �˶λ���ʵ�ʼ�����߶β����϶ϵ�
%     flagC = 'red';
%     for j=1:size(Point,2)
%         obj1 = abs(alpha-Point(5,j));
%         obj2 = abs(rho-Point(6,j));
%         if ((obj1<tolerence(1)||obj1>180-tolerence(1))&&obj2<tolerence(2))
%             flagC='green';
%             hitList=[hitList j];
%             break;
%         end
%     end
%     plot(xy(:,1),xy(:,2),'LineWidth',2,'Color',flagC);
%     plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%     plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','yellow');
%% �˶ζ��ڣ�û�гɹ�ƥ����߶�ֱ�ӻ��Ʋ����϶˵㣨�����ֱ�ߣ�
    flag=false;
    for j=1:size(Point,2)
        obj1 = abs(alpha-Point(5,j));
        obj2 = abs(rho-Point(6,j));
        if ((obj1<tolerence(1)||obj1>180-tolerence(1))&&obj2<tolerence(2))
            flag=true;
            hitList=[hitList j];
            break;
        end
    end
    if ~flag
        % �߶�
        plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','red');
        plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
        plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','yellow');
        % ֱ��
%         kl=(xy(2,2)-xy(1,2))/(xy(2,1)-xy(1,1));
%         bl=xy(1,2)-kl*xy(1,1);
%         klP=sizeP(1)/sizeP(2);
%         if kl<klP
%             px=[0 sizeP(2)];
%             py=kl*px+bl;
%             plot(px,py,'r:','LineWidth',2);
%         else
%             py=[0 sizeP(1)];
%             px=(py-bl)/kl;
%             plot(px,py,'r:','LineWidth',2);
%         end
    end
%%
end
hitList=unique(hitList);
%% �˶ζ��ڣ��ɹ�ƥ����߶λ��Ʊ궨�߶�
for i=1:length(hitList)
    xy = [Point(1,i),Point(2,i);Point(3,i),Point(4,i)];
    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
end
end

