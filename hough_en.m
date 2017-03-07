function [ H,theta,rho ] = hough_en( img )
%hough_en Implement the improved Hough transform
%   Using a convolution kernel containing directional encoding information
%   to determine the straight line texture at the pixel and to conduct 
%   a targeted vote

%% Set the defaults
thetaResolution = 1;
rhoResolution = 1;
kernal=[1000 100 10;1 0 1;10 100 1000;];
%% Compute theta and rho
[M,N] = size(img);
D = sqrt((M - 1)^2 + (N - 1)^2);
q = ceil(D/rhoResolution);
nrho = 2*q + 1;
rho = linspace(-q*rhoResolution, q*rhoResolution, nrho);                    % rho is constructed;
theta=linspace(-90,90,4*ceil(45/thetaResolution)+1);
theta(end)=[];                                                              % theta is constructed;
%% Compute matrix of judge
judge=zeros(size(img));
img_=double(img);
judge(2:end-1,2:end-1)=conv2(img_,kernal,'valid');
judge=img.*judge;
judge=sparse(judge);
%% Compute votes
[votes] = compute_votes (judge,theta,rho(1));
val=ones(size(votes,1),1);
%% Vote & Outout matrix of Hough
H=accumarray(votes,val,[length(rho) length(theta)]);
end
%% Module of Sub-function
function [votes] = compute_votes (jMat,theta,rho_min)
% Prepare data
vote_n90n45=[]; vote_n45n0=[];  vote_n0p45=[];  vote_p45p90=[];
theta_=reshape(theta,[],4);                                                 % theta is reconstructed;
theta_v=theta_-theta_(1)+1;                                                 % generate voting sequence of theta;
% Judge & Compute the rho of each point
[r_n90n45,c_n90n45]=find(jMat==1002|jMat==2001|jMat==1001|jMat==2002|jMat==1101|jMat==2000|jMat==1000|jMat==1011|jMat==2|jMat==1|jMat==2101|jMat==2011|jMat==1201|jMat==1021|jMat==1102|jMat==1012);
if (~isempty(r_n90n45))
    rho_v=fix(c_n90n45*cosd(theta_(:,1)')+r_n90n45*sind(theta_(:,1)'))-rho_min+1;
    vote_n90n45=[reshape(rho_v,[],1) reshape(repmat(theta_v(:,1)',length(r_n90n45),1),[],1)];
end
[r_n45n0,c_n45n0]=find(jMat==2100|jMat==1200|jMat==1100|jMat==2200|jMat==1101|jMat==2000|jMat==1000|jMat==1110|jMat==200|jMat==100|jMat==2110|jMat==2101|jMat==1210|jMat==1201|jMat==1120|jMat==1102);
if (~isempty(r_n45n0))
    rho_v=fix(c_n45n0*cosd(theta_(:,2)')+r_n45n0*sind(theta_(:,2)'))-rho_min+1;
    vote_n45n0=[reshape(rho_v,[],1) reshape(repmat(theta_v(:,2)',length(r_n45n0),1),[],1)];
end
[r_n0p45,c_n0p45]=find(jMat==210|jMat==120|jMat==110|jMat==220|jMat==1110|jMat==200|jMat==100|jMat==111|jMat==20|jMat==10|jMat==2110|jMat==1210|jMat==211|jMat==1120|jMat==121|jMat==112);
if (~isempty(r_n0p45))
    rho_v=fix(c_n0p45*cosd(theta_(:,3)')+r_n0p45*sind(theta_(:,3)'))-rho_min+1;
    vote_n0p45=[reshape(rho_v,[],1) reshape(repmat(theta_v(:,3)',length(r_n0p45),1),[],1)];
end
[r_p45p90,c_p45p90]=find(jMat==21|jMat==12|jMat==11|jMat==22|jMat==111|jMat==20|jMat==10|jMat==1011|jMat==2|jMat==1|jMat==2011|jMat==211|jMat==1021|jMat==121|jMat==1012|jMat==112);
if (~isempty(r_p45p90))
    rho_v=fix(c_p45p90*cosd(theta_(:,4)')+r_p45p90*sind(theta_(:,4)'))-rho_min+1;
    vote_p45p90=[reshape(rho_v,[],1) reshape(repmat(theta_v(:,4)',length(r_p45p90),1),[],1)];
end
votes=[vote_n90n45;vote_n45n0;vote_n0p45;vote_p45p90];
end