function EnKFLorenz95
clc
clear all
close all

%model parameters
N_s=40; %No. of states
Fo=8;

dt=.01; %integration sample time
tvec=0:dt:160*dt; %time vector
ini_model_per=.008; %initial model perturbation
dev_ini_per=.05; %standard deviation for model perturbations
obs_percent=50; %percent for No. of locations of sensors
N_obs=round(N_s*(obs_percent/100)); %No. of observations
obs_error_dev=.5; %observations error standard deviation


R=obs_error_devˆ2*eye(N_obs); %observations noise covriance

H=zeros(N_obs,N_s); %initial value for measurement matrix
temp_loc=randperm(N_s);
sen_loc=sort(temp_loc(1:N_obs)); %sensor locations randomly generated

for j=1:N_obs %generating measurement matrix
    H(j,sen_loc(j))=1;
end

%truth and process noise
truth=ones(N_s,1)*Fo;
truth(round(end/2))=truth(round(end/2))+ini_model_per;
randn('seed',0)
ymeas(:,1)=truth(sen_loc(:),1)+sqrtm(R)*randn(N_obs,1);

%ensembles
NE=200; %No. of ensembles
%generate initial ensemble
