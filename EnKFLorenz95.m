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

P_a=dev_ini_perˆ2*eye(N_s);

xinit=truth(:,1)+dev_ini_perˆ2*randn(N_s,1); %analysis at t=0
Af(:,:,1)=repmat(xinit,1,NE)+sqrtm(P_a)*randn(N_s,NE);

Q=P_a; %Q: process noise covariance matrix

for indt=1:length(tvec)
    
    A=Af(:,:,indt);
    Abar=repmat(mean(A,2),1,NE);
    Aprime=A-Abar;
    Pf=Aprime*Aprime'/(NE-1);
    E=sqrtm(R)*randn(length(R),NE);
    D=repmat(ymeas(:,indt),1,NE)+E;
    Re=E*E'/(NE-1);
    K_gain=Pf*H'*pinv(H*Pf*H'+Re);
    Aa=A+K_gain*(D-H*A);
    Pa=(eye(N_s)-K_gain*H)*Pf*(eye(N_s)-K_gain*H)'+K_gain*Re*K_gain';
    
    truth(:)=truth(:)+rk4(truth(:),dt,Fo,N_s)+sqrtm(Q)*randn(N_s,1);
    ymeas(:,indt+1)=truth(sen_loc(:),1)+sqrtm(R)*randn(N_obs,1);
    
    for indE=1:NE
        A(:,indE)=Aa(:,indE)+rk4(Aa(:,indE),dt,Fo,N_s)+sqrtm(Q)*randn(N_s,1);
    end
    
    Af(:,:,indt+1)=A;
    Afmean(:,indt+1)=mean(A,2);
    Aamean(:,indt+1)=mean(Aa,2);
    tt(:,indt)=truth(:);
end

figure
plot(truth(:),'r','Linewidth',2);
hold on
plot(sen_loc,ymeas(:,end),'g*');
hold on
plot(mean(Aa,2),'b--','Linewidth',2);
legend('truth','measurement','analysis');
set(gca,'fontsize',20);

rms1=rms(truth(:)-mean(Aa,2)) %root mean squared error

keyboard


function k=L95(X,F,N_s)

    k=zeros(N_s,1);
    k(1)=(X(2)-X(N_s-1))*X(N_s)-X(1);
    k(2)=(X(3)-X(N_s))*X(1)-X(2);
    k(N_s)=(X(1)-X(N_s-2))*X(N_s-1)-X(N_s);
    for j=3:N_s-1
        k(j)=(X(j+1)-X(j-2))*X(j-1)-X(j);
    end

    k=k+F;


%runge kutta 4th order to solve the ode
function dx=rk4(Xprev,dt,F,N_s)

    k1=L95(Xprev,F,N_s);
    k2=L95(Xprev+.5*dt*k1,F,N_s);
    k3=L95(Xprev+.5*dt*k2,F,N_s);
    k4=L95(Xprev+dt*k3,F,N_s);
    dx=1/6*dt*(k1+2*k2+2*k3+k4);
    dx=1/6*dt*(k1+2*k2+2*k3+k4);
