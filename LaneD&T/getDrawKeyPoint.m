function [ drawKeyPoint ] = getDrawKeyPoint( trackKeyPoint,trackDeadLine,drawDeadLine,frameSize )
%getDrawKeyPoint 计算绘图需要的关键点位置信息
%   依据用于追踪的关键点计算帧上相应区域的关键点以便绘图
%   trackKeyPoint   -> [leftBottomColumn,leftTopColumn,rightBottomColumn,rightTopColumn]
%   frameSize       -> [width,height]
%   drawKeyPoint    -> [up,left,right]
    %% [x(c),y(r)]
    refLeftBottom = [trackKeyPoint(1),frameSize(2)];
    refLeftTop = [trackKeyPoint(2),trackDeadLine];
    refRightBottom = [trackKeyPoint(3),frameSize(2)];
    refRightTop = [trackKeyPoint(4),trackDeadLine];
    cLeftTop = floor((refLeftBottom(1)-refLeftTop(1))*(drawDeadLine-refLeftTop(2))/(refLeftBottom(2)-refLeftTop(2))+refLeftTop(1)+0.5);
    cRightTop = floor((refRightBottom(1)-refRightTop(1))*(drawDeadLine-refRightTop(2))/(refRightBottom(2)-refRightTop(2))+refRightTop(1)+0.5);
    %%
    %
    if (cLeftTop<cRightTop)
        upSegment = [cRightTop,drawDeadLine,cLeftTop,drawDeadLine];
    else
        ratio = (refRightBottom(1)-refLeftBottom(1))/(refLeftTop(1)-refRightTop(1)+refRightBottom(1)-refLeftBottom(1));
        crossPoint = ratio*(refLeftTop-refLeftBottom)+refLeftBottom;
        upSegment = crossPoint;
    end
    %
    if (refLeftBottom(1)<1)
        rLeftLeft = floor((refLeftBottom(2)-refLeftTop(2))*(1-refLeftTop(1))/(refLeftBottom(1)-refLeftTop(1))+refLeftTop(2)+0.5);
        leftSegment = [1,rLeftLeft,1,frameSize(2)];
    else
        leftSegment = refLeftBottom;
    end
    %
    if (refRightBottom(1)>frameSize(1))
        rRightRight = floor((refRightBottom(2)-refRightTop(2))*(1-refRightTop(1))/(refRightBottom(1)-refRightTop(1))+refRightTop(2)+0.5);
        rightSegment = [frameSize(1),frameSize(2),frameSize(1),rRightRight];
    else
        rightSegment = refRightBottom;
    end
    drawKeyPoint = [upSegment,leftSegment,rightSegment];
end

