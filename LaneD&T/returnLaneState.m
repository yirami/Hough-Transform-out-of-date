function [ laneStateInStructure ] = returnLaneState( laneStateInVector,trackDeadLine,frameSize,refPoint )
%returnLaneState �ش���ǰ������״̬
%   ��������ʾ��״̬���ƽṹ���ʾ�ĳ���״̬

leftStartPoint = [laneStateInVector(1),frameSize(2)];
leftEndPoint = [laneStateInVector(2),trackDeadLine];
rightStartPoint = [laneStateInVector(3),frameSize(2)];
rightEndPoint = [laneStateInVector(4),trackDeadLine];
leftTheta = atand((leftStartPoint(1)-leftEndPoint(1))/(leftEndPoint(2)-leftStartPoint(2)));
rightTheta = atand((rightStartPoint(1)-rightEndPoint(1))/(rightEndPoint(2)-rightStartPoint(2)));
laneStateInStructure = struct('startPoint',[],'endPoint',[],'theta',[],'refRho',[]);
laneStateInStructure(1).startPoint = leftStartPoint;
laneStateInStructure(1).endPoint = leftEndPoint;
laneStateInStructure(1).theta = leftTheta;
laneStateInStructure(1).refRho = (leftStartPoint(1)-refPoint(1))*cosd(leftTheta)+(leftStartPoint(2)-refPoint(2))*sind(leftTheta);
laneStateInStructure(2).startPoint = rightStartPoint;
laneStateInStructure(2).endPoint = rightEndPoint;
laneStateInStructure(2).theta = rightTheta;
laneStateInStructure(2).refRho = (rightStartPoint(1)-refPoint(1))*cosd(rightTheta)+(rightStartPoint(2)-refPoint(2))*sind(rightTheta);
end

