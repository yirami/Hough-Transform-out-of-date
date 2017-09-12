%% Clear ALL
clear;
clc;
%% Creat Video Object
VideoSrc = vision.VideoFileReader('tracking.mp4','VideoOutputDataType','uint8');
frameSize = VideoSrc.info.VideoSize;
VideoPosition = [100,100,frameSize];
VideoOut = vision.VideoPlayer('Name','Demo of tracking','Position',VideoPosition);
%% Config the visualization parameter
drawDeadLine = floor(frameSize(2)/2);
drawVector = [];
%% Initial the tool of frame rate analysis
thisFigure = figure(1);
thisFigure.Name = '帧率分析';
xlabel('批（50帧）');
ylabel('平均帧速');
title('实时帧速分析');
hold on;
lastPoint = [0,0];
sizeofPool = 50;
timeofPool = 0;
index = sizeofPool;
countofPool = 0;
%% Initial the state variable
% search need
stateOK = false;
refPoint = [floor(frameSize(1)/2),frameSize(2)];
% track need
trackDeadLine = floor(frameSize(2)/2);
DoLeft = true;
DoRight = true;
laneStateInStructure = [];
laneStateInVector = [];
% kalman filter
isKalmanFilterInitial = false;
maxMissDetectionFrame = 20;
restMissDetectionFrame = maxMissDetectionFrame;
thisCase.MotionModel           = 'ConstantAcceleration';
thisCase.InitialEstimateError  = 1E5 * ones(1, 3);
%thisCase.MotionNoise           = [25, 10, 1];
thisCase.MotionNoise           = [0.5, 0.2, 0.02];
thisCase.MeasurementNoise      = 25;
%% Debug
% frameNum = 0;
%% Frame processing loop
while ~isDone(VideoSrc)
    %% Start the stopwatch
    tic;
    %% Load the frame
    rawImage = step(VideoSrc);
    %%
    if (~stateOK)
        [hasLines,lines] = DCHTbandGlobal( rawImage );
        if (hasLines)
            [ leftOK,rightOK,twoLane ] = laneMatchForSearch( lines,frameSize,refPoint);
            if (leftOK && rightOK)
                stateOK = true;
                restMissDetectionFrame = maxMissDetectionFrame;
                laneStateInStructure = twoLane;
                laneStateInVector = getTrackKeyPoint( laneStateInStructure,trackDeadLine,frameSize );
                if (~isKalmanFilterInitial)
                    isKalmanFilterInitial = true;
                    thisCase.InitialLocation = laneStateInVector;
                    thisCase.KalmanFilter = configureKalmanFilter( ...
                        thisCase.MotionModel, ...
                        thisCase.InitialLocation, ...
                        thisCase.InitialEstimateError, ...
                        thisCase.MotionNoise, ...
                        thisCase.MeasurementNoise);
                    drawVector = getDrawKeyPoint( laneStateInVector,trackDeadLine,drawDeadLine,frameSize );
                else
                    predict(thisCase.KalmanFilter);
                    laneStateInVector = correct(thisCase.KalmanFilter,laneStateInVector);
                    drawVector = getDrawKeyPoint( laneStateInVector,trackDeadLine,drawDeadLine,frameSize );
                end
            else
                stateOK = false;
            end
        else
            stateOK = false;
        end
    else
        [ hasLines,lines ] = DCHTbandLocal( rawImage,DoLeft,DoRight,laneStateInStructure );
        if (hasLines)
            [ leftOK,rightOK,twoLane ] = laneMatchForTrack( laneStateInStructure,DoLeft,DoRight,lines,frameSize);
            if (leftOK && rightOK)
                restMissDetectionFrame = maxMissDetectionFrame;
                laneStateInStructure = twoLane;
                laneStateInVector = getTrackKeyPoint( laneStateInStructure,trackDeadLine,frameSize );
                predict(thisCase.KalmanFilter);
                laneStateInVector = correct(thisCase.KalmanFilter,laneStateInVector);
                drawVector = getDrawKeyPoint( laneStateInVector,trackDeadLine,drawDeadLine,frameSize );
            else
                stateOK = false;
            end
        else
            stateOK = false;
        end
        if (~stateOK && restMissDetectionFrame>0)
            stateOK = true;
            restMissDetectionFrame = restMissDetectionFrame-1;
            laneStateInVector = predict(thisCase.KalmanFilter);
            laneStateInStructure = returnLaneState( laneStateInVector,trackDeadLine,frameSize,refPoint );
            drawVector = getDrawKeyPoint( laneStateInVector,trackDeadLine,drawDeadLine,frameSize );
        end
    end

    %% Show this frame
    if (stateOK)
        road2Image = insertShape(rawImage,'FilledPolygon',drawVector,'Color',[0	255	128],'Opacity',0.2);
        step(VideoOut, road2Image);
    else
        step(VideoOut, rawImage);
    end
    %% Stop the stopwatch and Refresh the data 
    timeConsuming = toc;
    if (index>0)
        index =index-1;
        timeofPool = timeofPool + timeConsuming;
    else
        index = sizeofPool-1;
        countofPool = countofPool+1;
        thisPoint = [countofPool,sizeofPool/timeofPool];
        plot([lastPoint(1),thisPoint(1)],[lastPoint(2),thisPoint(2)]);
        hold on;
        lastPoint = thisPoint;
        timeofPool = timeConsuming;
    end
end
%%
release(VideoOut);
release(VideoSrc);