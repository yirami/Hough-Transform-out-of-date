function [ lines ] = fixedDCHTband( image,rect,angleSpec,linesMax,linesGap )
%fDCHTband 固定区域DCHT检测直线
%   对于给定区域，给定的角度范围采用DCHT方法进行直线检测并输出最终的直线信息
%   image       -> original image
%   rect        -> specified area
%   angleSpec   -> specified angle start : end
%   linesMax    -> maximum number of lines detected
%   linesGap    -> maximum distance allowed when searching for lines
%% Size of image
sizeRaw = size(image);
%% Fix the rectangle
rect = round(rect);
fixedRect = [rect(1)+1,rect(2)+1,rect(3)-1,rect(4)-1];
% Note: imcrop -> 'rect' return by imcrop is usually decimal
%              -> rect = [xmin(cLeft),ymin(rTop),width,height];
%              -> start at column 'xmin' with 'width+1'
%              -> start at row 'ymin' with 'height+1'
%% Transform the format of image and detect the edge
fixedImage = imcrop(image,fixedRect);
thisGray = colorFilter(fixedImage);
% thisGray = rgb2gray(fixedImage);
thisEdge = edge(thisGray,'canny');
thisEdgeUint32 = uint32(thisEdge);
%% Transform
[thisH,thisTheta,thisRho] = fixedDCHTa(thisEdgeUint32,sizeRaw(1:2),rect(2:-1:1),angleSpec);
%% Search
[linesInfo] = fixedsearchLines(thisEdgeUint32,thisH,thisTheta,thisRho,linesMax,linesGap,rect(2:-1:1),angleSpec);
%% Format stored
SizeOfLines = size(linesInfo,2);
lines(SizeOfLines) =struct('pointStart',[],'pointEnd',[],'theta',[],'rho',[]);
for i=1:SizeOfLines
    lines(i).pointStart = [linesInfo(2,i),linesInfo(1,i)];
    lines(i).pointEnd = [linesInfo(4,i),linesInfo(3,i)];
    lines(i).theta = linesInfo(5,i);
    lines(i).rho = linesInfo(6,i);
end

end

