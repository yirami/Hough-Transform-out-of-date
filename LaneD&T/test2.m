clear;
clc;

load Record

thisCase.MotionModel           = 'ConstantAcceleration';
thisCase.InitialLocation       = record(1,:);
thisCase.InitialEstimateError  = 1E5 * ones(1, 3);
thisCase.MotionNoise           = [25, 10, 1];
thisCase.MeasurementNoise      = 25;

thisCase.KalmanFilter = configureKalmanFilter(thisCase.MotionModel,thisCase.InitialLocation,thisCase.InitialEstimateError,thisCase.MotionNoise,thisCase.MeasurementNoise);

for i=1:200
    if (mod(i,10)==0)
        detect(i,:) = record(i,:);
        track(i,:) = predict(thisCase.KalmanFilter);
    else
        detect(i,:) = record(i,:);
        predict(thisCase.KalmanFilter);
        track(i,:) = correct(thisCase.KalmanFilter, detect(i,:));
    end
end

figure(1);
set(gcf,'color','w');
plot(detect(:,1),'r');hold on;
plot(detect(:,2),'r');hold on;
plot(detect(:,3),'r');hold on;
plot(detect(:,4),'r');hold on;
plot(track(:,1),'g');hold on;
plot(track(:,2),'g');hold on;
plot(track(:,3),'g');hold on;
plot(track(:,4),'g');hold on;

thisCase.MotionNoise           = [0.5, 0.2, 0.02];
thisCase.KalmanFilter = configureKalmanFilter(thisCase.MotionModel,thisCase.InitialLocation,thisCase.InitialEstimateError,thisCase.MotionNoise,thisCase.MeasurementNoise);

for i=1:200
    if (mod(i,10)==0)
        detect(i,:) = record(i,:);
        track(i,:) = predict(thisCase.KalmanFilter);
    else
        detect(i,:) = record(i,:);
        predict(thisCase.KalmanFilter);
        track(i,:) = correct(thisCase.KalmanFilter, detect(i,:));
    end
end

figure(2);
set(gcf,'color','w');
plot(detect(:,1),'r');hold on;
plot(detect(:,2),'r');hold on;
plot(detect(:,3),'r');hold on;
plot(detect(:,4),'r');hold on;
plot(track(:,1),'g');hold on;
plot(track(:,2),'g');hold on;
plot(track(:,3),'g');hold on;
plot(track(:,4),'g');hold on;


