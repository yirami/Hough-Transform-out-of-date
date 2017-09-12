function [ leftOK,rightOK,twoLane ] = laneMatchForTrack( priorLines,DoLeft,DoRight,lines,frameSize)
%laneMatchForTrack 从线段结果中匹配左右车道线
%   基于先验车道线参数，匹配当前帧车道线，超过一定阈值则匹配失败
%%
refPoint = [floor(frameSize(1)/2),frameSize(2)];
leftOK = false;
rightOK = false;
failThreshold = 10;
if (DoLeft)
    leftOK = true;
    leftRhoBenchmark = priorLines(1).refRho;
    leftThetaBenchmark = priorLines(1).theta;
    leftGoal = failThreshold;
end
if (DoRight)
    rightOK = true;
    rightRhoBenchmark = priorLines(2).refRho;
    rightThetaBenchmark = priorLines(2).theta;
    rightGoal = failThreshold;
end
%%
twoLane(2) = struct('startPoint',[],'endPoint',[],'theta',[],'refRho',[]);
for i = 1:length(lines)
    thisPs = lines(i).startPoint;
    thisPe = lines(i).endPoint;
    thisTheta = lines(i).theta;
    thisRefRho = (thisPs(1)-refPoint(1))*cosd(thisTheta)+(thisPs(2)-refPoint(2))*sind(thisTheta);
    if (DoLeft)
        goal = 0.1*abs(thisRefRho-leftRhoBenchmark)+abs(thisTheta-leftThetaBenchmark);
        if (goal<leftGoal)
            twoLane(1).startPoint = thisPs;
            twoLane(1).endPoint = thisPe;
            twoLane(1).theta = thisTheta;
            twoLane(1).refRho = thisRefRho;
            leftGoal = goal;
        end
    end
    if (DoRight)
        goal = 0.1*abs(thisRefRho-rightRhoBenchmark)+abs(thisTheta-rightThetaBenchmark);
        if (goal<rightGoal)
            twoLane(2).startPoint = thisPs;
            twoLane(2).endPoint = thisPe;
            twoLane(2).theta = thisTheta;
            twoLane(2).refRho = thisRefRho;
            rightGoal = goal;
        end
    end
end
%%
if (DoLeft && leftGoal==failThreshold)
    leftOK = false;
end
if (DoRight && rightGoal==failThreshold)
    rightOK = false;
end
end
