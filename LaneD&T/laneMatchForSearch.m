function [ leftOK,rightOK,twoLane ] = laneMatchForSearch( lines,frameSize,refPoint)
%laneMatchForSearch ���߶ν����ƥ�����ҳ�����
%   ���ڲο��㣬����ÿ���߶���������rho����ѡ�������С������һ��
%   refPoint    -> [x(c),y(r)]
%%
leftOK = false;
rightOK = false;
hasVague = false;
%%
twoLane(2) = struct('startPoint',[],'endPoint',[],'theta',[],'refRho',[]);
rhoBenchmark = sum(frameSize);
leftBenchmark = -rhoBenchmark;
rightBenchmark = rhoBenchmark;
vagueBenchmark = rhoBenchmark;
for i = 1:length(lines)
    thisPs = lines(i).startPoint;
    thisPe = lines(i).endPoint;
    thisTheta = lines(i).theta;
    thisRho = (thisPs(1)-refPoint(1))*cosd(thisTheta)+(thisPs(2)-refPoint(2))*sind(thisTheta);
    if (~thisTheta)
        if (abs(thisRho)<vagueBenchmark)
            vaguePosition = i;
            vagueBenchmark = abs(thisRho);
        end
        continue;
    end
    if (thisTheta>0 && thisRho<0)
        if (thisRho>leftBenchmark)
            twoLane(1).startPoint = thisPs;
            twoLane(1).endPoint = thisPe;
            twoLane(1).theta = thisTheta;
            twoLane(1).refRho = thisRho;
            leftBenchmark = thisRho;
        end
    elseif (thisTheta<0 && thisRho>0)
        if (thisRho<rightBenchmark)
            twoLane(2).startPoint = thisPs;
            twoLane(2).endPoint = thisPe;
            twoLane(2).theta = thisTheta;
            twoLane(2).refRho = thisRho;
            rightBenchmark = thisRho;
        end
    end
end
%%
if (leftBenchmark~=-rhoBenchmark)
    leftOK = true;
end
if (rightBenchmark~=rhoBenchmark)
    rightOK = true;
end
%%
if ~(leftOK && rightOK) && hasVague
    if xor(leftOK,rightOK)
        if (leftOK)
            twoLane(2).startPoint = lines(vaguePosition).startPoint;
            twoLane(2).endPoint = lines(vaguePosition).endPoint;
            twoLane(2).theta = lines(vaguePosition).theta;
            twoLane(2).refRho = vagueBenchmark;
            rightOK = true;
        else
            twoLane(1).startPoint = lines(vaguePosition).startPoint;
            twoLane(1).endPoint = lines(vaguePosition).endPoint;
            twoLane(1).theta = lines(vaguePosition).theta;
            twoLane(1).refRho = -vagueBenchmark;
            leftOK = true;
        end
    end
end
end

