function [ hasLines,lines ] = DCHTbandGlobal( image )
%DCHTbandGlobal 全局车道线搜索
%   对于给定帧采用DCHT方法进行直线检测并输出最终的直线信息
%   image       -> original image
%% Configuration
defaultArea = [332,151,840,568];
fixedRect = [defaultArea(1)+1,defaultArea(2)+1,defaultArea(3)-1,defaultArea(4)-1];
fixedImage = imcrop(image,fixedRect);
%% Transform the format of image and detect the edge
thisGray = colorFilter(fixedImage);
thisEdge = edge(thisGray,'canny');
thisEdgeUint32 = uint32(thisEdge);
%% Transform
[thisH,thisTheta,thisRho] = voteDCHT(thisEdgeUint32);
%% Search
[linesInfo] = laneSearch(thisEdgeUint32,thisH,thisTheta,thisRho);
%% Format stored
NumOfLines = size(linesInfo,2);
if (NumOfLines)
    hasLines = true;
    lines(NumOfLines) =struct('startPoint',[],'endPoint',[],'theta',[],'rho',[]);
    for i=1:NumOfLines
        lines(i).startPoint = [linesInfo(2,i)+defaultArea(1),linesInfo(1,i)+defaultArea(2)];
        lines(i).endPoint = [linesInfo(4,i)+defaultArea(1),linesInfo(3,i)+defaultArea(2)];
        lines(i).theta = linesInfo(5,i);
        lines(i).rho = linesInfo(6,i);
    end
else
    hasLines = false;
    lines(1) =struct('startPoint',[],'endPoint',[],'theta',[],'rho',[]);
end

end

