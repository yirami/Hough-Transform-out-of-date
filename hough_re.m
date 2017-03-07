function [ H,theta,rho ] = hough_re( bw )
% 标准Hough变换的重写，由于内置Hough已经采用Mex实现了算法，
%                           性能优化充分，不适合用来在M平台比较改进算法的性能；
% 由于仅用于比较，故只保留了精度为1的设置，后续可完善；
% Input[bw] must be a binary image;

% Set the defaults
thetaResolution = 1;
rhoResolution = 1;
theta = [];

% Compute theta and rho
[M,N] = size(bw);

if (isempty(theta))
    theta = linspace(-90, 0, ceil(90/thetaResolution) + 1);
    theta = [theta -fliplr(theta(2:end - 1))];
end

D = sqrt((M - 1)^2 + (N - 1)^2);
q = ceil(D/rhoResolution);
nrho = 2*q + 1;
rho = linspace(-q*rhoResolution, q*rhoResolution, nrho);
% Get the coordinates of nonzero pixels 
[r,c]=find(bw~=0);
% Compute rho
rho_cal=fix(c*cosd(theta)+r*sind(theta));
% Compute the voting position of rho
rho_cal_v=rho_cal+1-rho(1);
% Reconstruct the voting position of rho and theta
theta_v=theta+1-theta(1);
theta_v=repmat(theta_v,size(rho_cal_v,1),1);
rho_cal_v=reshape(rho_cal_v,[],1);
theta_v=reshape(theta_v,[],1);
% Prepare data
subs=[rho_cal_v theta_v];
val=ones(length(subs),1);
% Vote
H=accumarray(subs,val,[length(rho) length(theta)]);
end
