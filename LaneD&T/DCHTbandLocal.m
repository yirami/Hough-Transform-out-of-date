function [ hasLines,lines ] = DCHTbandLocal( image,DoLeft,DoRight,trackInfo )
%fDCHTband 固定区域DCHT检测直线
%   对于给定区域，给定的角度范围采用DCHT方法进行直线检测并输出最终的直线信息
%   image       -> original image
%   rect        -> specified area
%   angleSpec   -> specified angle start : end
%% Prepare data
frameSize = size(image);
refPoint = [frameSize(1),floor(frameSize(2)/2)];
refLine = floor(frameSize(1)/2);
linesInfo =[];
%% Tracking left lane
if (DoLeft)
    [ leftRegionOK,leftLaneRect ] = calculateRegion( trackInfo(1).startPoint,trackInfo(1).theta,refLine,frameSize );
    if (leftRegionOK)
        targetInfo = [trackInfo(1).theta,trackInfo(1).refRho,trackInfo(1).startPoint(2),trackInfo(1).startPoint(1)];
        [ linesInfo ] = task(image,leftLaneRect,targetInfo,refPoint,trackInfo(1).theta-10,trackInfo(1).theta+10);
        [ linesInfo ] = fixedPoint( linesInfo,leftLaneRect);
    end
end
%% Tracking right lane
if (DoRight)
    [ rightRegionOK,rightLaneRect ] = calculateRegion( trackInfo(2).startPoint,trackInfo(2).theta,refLine,frameSize );
    if (rightRegionOK)
        targetInfo = [trackInfo(2).theta,trackInfo(2).refRho,trackInfo(2).startPoint(2),trackInfo(2).startPoint(1)];
        [ linesInfoTemp ] = task( image,rightLaneRect,targetInfo,refPoint,trackInfo(2).theta-10,trackInfo(2).theta+10);
        [ linesInfo ] = [linesInfo,fixedPoint( linesInfoTemp,rightLaneRect)];
    end
end
%% Format stored
NumOfLines = size(linesInfo,2);
if (NumOfLines)
    hasLines = true;
    lines(NumOfLines) =struct('startPoint',[],'endPoint',[],'theta',[]);
    for i=1:NumOfLines
        lines(i).startPoint = [linesInfo(2,i),linesInfo(1,i)];
        lines(i).endPoint = [linesInfo(4,i),linesInfo(3,i)];
        lines(i).theta = linesInfo(5,i);
    end
else
    hasLines = false;
    lines(1) =struct('startPoint',[],'endPoint',[],'theta',[]);
end

end

function [ isOK,rect ] = calculateRegion( startPoint,theta,refLine,frameSize )
isOK = true;
thisSin = sind(theta);
thisCos = cosd(theta);
rhoG = startPoint(1)*thisCos+startPoint(2)*thisSin;
refCB = floor((rhoG-frameSize(1)*thisSin)/thisCos+0.5);
refCT = floor((rhoG-refLine*thisSin)/thisCos+0.5);
if ((refCB<1 || refCB>frameSize(2)) && (refCT<1 || refCT>frameSize(2)))
    isOK = false;
end
if (refCB<refCT)
    maxC = refCT;
    minC = refCB;
else
    maxC = refCB;
    minC = refCT;
end
if (minC>11)
    leftC = minC-10;
else
    leftC = 1;
end
if (maxC+10<frameSize(2))
    rightC = maxC+10;
else
    rightC = frameSize(2);
end
rect = [leftC refLine rightC-leftC-1 frameSize(1)-refLine-1];
end

function [ linesInfo ] = task( image,rect,targetInfo,refPoint,angleS,angleE)
fixedImage = imcrop(image,rect);
% Note: imcrop -> 'rect' return by imcrop is usually decimal
%              -> rect = [xmin(cLeft),ymin(rTop),width,height];
%              -> start at column 'xmin' with 'width+1'
%              -> start at row 'ymin' with 'height+1'
thisGray = colorFilter(fixedImage);
thisEdge = edge(thisGray,'canny');
thisEdgeUint32 = uint32(thisEdge);
if (angleS<-90)
    angleS = 180+angleS;
end
if (angleE>90)
    angleE=angleE-180;
end
[thisH,thisTheta,thisRho] = voteDCHTa(thisEdgeUint32,angleS,angleE);
[linesInfo] = laneTrack(thisEdgeUint32,thisH,thisTheta,thisRho,1,targetInfo,refPoint);
end

function [ linesInfo ] = fixedPoint( linesInfo,rect)
linesInfo(1,:) = linesInfo(1,:)+rect(2)-1;
linesInfo(2,:) = linesInfo(2,:)+rect(1)-1;
linesInfo(3,:) = linesInfo(3,:)+rect(2)-1;
linesInfo(4,:) = linesInfo(4,:)+rect(1)-1;
end