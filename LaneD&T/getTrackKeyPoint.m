function [ trackKeyPoint ] = getTrackKeyPoint( laneInfoStruct,trackDeadLine,frameSize )
%getTrackKeyPoint 计算追踪需要的关键点位置信息
%   对这些关键点进行Kalman滤波
%   laneInfoStruct  -> structure of two lanes
%   frameSize       -> [width,height]
%   trackKeyPoint   -> [leftBottomColumn,leftTopColumn,rightBottomColumn,rightTopColumn]
    %%
    for i=2:-1:1
        thisTheta = laneInfoStruct(i).theta;
        thisStartPoint = laneInfoStruct(i).startPoint;
        thisTan = tand(thisTheta);
        RefB = floor((thisStartPoint(2)-frameSize(2))*thisTan+thisStartPoint(1)+0.5);
        trackRefT = floor((thisStartPoint(2)-trackDeadLine)*thisTan+thisStartPoint(1)+0.5);
        trackKeyPoint(2*i) = trackRefT;
        trackKeyPoint(2*i-1) = RefB;
    end
end

