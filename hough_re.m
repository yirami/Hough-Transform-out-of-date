function [ H,theta,rho ] = hough_re( bw )
% ��׼Hough�任����д����������Hough�Ѿ�����Mexʵ�����㷨��
%                           �����Ż���֣����ʺ�������Mƽ̨�ȽϸĽ��㷨�����ܣ�
% ���ڽ����ڱȽϣ���ֻ�����˾���Ϊ1�����ã����������ƣ�
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
