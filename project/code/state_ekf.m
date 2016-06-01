%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% state_ekf.m
% Author: Andrew Bylard
% Stanford University
% May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

ns = 6;     % Dimension of state
np = 4;     % Dimension of parameters
m = 6;      % Dimension of measurements

% Free-flyer parameters
mf = 10;            % [kg] Mass
rf = 0.15;          % [m] Radius
Jf = 2/5*mf*rf^2;   % [kg*m^2] Moment of inertia
L = 0.05;           % [m] Gripper length

% Grasped object parameters
mo = 5;             % [kg] Mass
ro = 0.30;          % [m] Radius
Jo = 2/5*mo*ro^2;   % [kg*m^2] Moment of inertia
Rcm = rf+L+ro;      % [m] Distance between centers of mass
alpha = 0;          % [rad] Object center of mass offset

mtot = mf+mo;       % [kg] Total mass
Jtot = Jf+Jo;       % [kg*m^2] Summed moments of inertia
J1 = Jtot*mtot+mf*mo*Rcm^2;     % [kg^2*m^2] Common expression

T = 70;
dt = 0.005;
Nsteps = T/dt;
t = 0:dt:T;

%% State EKF Setup
x0 = zeros(ns, 1);
Qs = diag([0.001*ones(3,1); 0.01*ones(3,1)]);
R = diag([0.001*ones(3,1); 0.01*ones(3,1)]);
mu0 = x0;
Sigma0 = 0.1*eye(ns);
    
mu_up = zeros(ns, Nsteps+1);
mu_pr = zeros(ns, Nsteps+1);
Sigma_up = zeros(ns,ns,Nsteps+1);
Sigma_pr = zeros(ns,ns,Nsteps+1);
mu_up(:,1) = mu0;
mu_pr(:,1) = mu0;
Sigma_up(:,:,1) = Sigma0;
Sigma_pr(:,:,1) = Sigma0;
y = zeros(m,Nsteps+1);
x = zeros(ns, Nsteps+1);
x(:,1) = x0;
y_up = zeros(m,Nsteps+1);

At = zeros(ns,ns);
Ct = zeros(m,ns);

for i = 1:Nsteps
    % Input
    a = 0.5;
    Fx = a*sin(t(i)/5);
    Fy = a*cos(t(i)/10);
    tau = a*cos(t(i)/7);
    
    % Predict next state
    mu_pxd = mu_up(1,i);
    mu_pyd = mu_up(2,i);
    mu_thd = mu_up(3,i);
    mu_px = mu_up(4,i);
    mu_py = mu_up(5,i);
    mu_th = mu_up(6,i);
    cta = cos(mu_th+alpha);
    sta = sin(mu_th+alpha);
    mu_pr(:,i+1) = [mu_pxd+dt*(Fx/mtot+mo*Rcm*mu_thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1))
                    mu_pyd+dt*(Fy/mtot+mo*Rcm*mu_thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1))
                    mu_thd+dt*(mtot*tau+mo*Rcm*(Fy*cta+Fx*sta))/J1
                    mu_pxd+dt*mu_pxd
                    mu_pyd+dt*mu_pyd
                    mu_thd+dt*mu_thd];
    At = eye(ns);
    At(1,3) = 2*dt/mtot*mo*Rcm*mu_thd*cta;
    At(2,3) = 2*dt/mtot*mo*Rcm*mu_thd*sta;
    At(1,6) = -dt*mo*Rcm*(-Fy*mo*Rcm*cos(2*(mu_th+alpha))+J1*mu_thd^2*sta-cta*(mtot*tau+2*Fx*mo*Rcm*sta))/(mtot*J1);
    At(2,6) =  dt*mo*Rcm*(-Fx*mo*Rcm*cos(2*(mu_th+alpha))+J1*mu_thd^2*sta+sta*(mtot*tau+2*Fy*mo*Rcm*cta))/(mtot*J1);    
    At(3,6) =  dt*mo*Rcm/J1*(Fx*cta-Fy*sta);
    Sigma_pr(:,:,i+1) = At*Sigma_up(:,:,i)*At'+Qs;
     
    % Simulate the satellite
    pxd = x(1,i); pyd = x(2,i); thd = x(3,i); px = x(4,i); py = x(5,i); th = x(6,i);
    cta = cos(th+alpha); sta = sin(th+alpha);
    x(:,i+1) = [pxd+dt*(Fx/mtot+mo*Rcm*thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1))
                pyd+dt*(Fy/mtot+mo*Rcm*thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1))
                thd+dt*(mtot*tau+mo*Rcm*(Fy*cta+Fx*sta))/J1
                px+dt*pxd
                py+dt*pyd
                th+dt*thd] + mvnrnd(zeros(ns,1),Qs*dt^2)';
            
    % Take measurement
    pxd = x(1,i+1); pyd = x(2,i+1); thd = x(3,i+1); px = x(4,i+1); py = x(5,i+1); th = x(6,i+1);
    cta = cos(th+alpha); sta = sin(th+alpha);
    y(:,i+1) = [px
                py
                th
                Fx/mtot+mo*Rcm*thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1)
                Fy/mtot+mo*Rcm*thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1)
                thd] + mvnrnd(zeros(m,1),R)';
    
    % Update state estimate
    mu_pxd = mu_pr(1,i+1); mu_pyd = mu_pr(2,i+1); mu_thd = mu_pr(3,i+1); mu_px = mu_pr(4,i+1); mu_py = mu_pr(5,i+1); mu_th = mu_pr(6,i+1);
    cta = cos(th+alpha); sta = sin(th+alpha);
    Ct(1,4) = 1; Ct(2,5) = 1; Ct(3,6) = 1; Ct(6,3) = 1;
    Ct(4,6) = -mo*Rcm*(-Fy*mo*Rcm*cos(2*(mu_th+alpha))+J1*mu_thd^2*sta-cta*(mtot*tau+2*Fx*mo*Rcm*sta))/(mtot*J1);
    Ct(5,6) =  mo*Rcm*(-Fx*mo*Rcm*cos(2*(mu_th+alpha))+J1*mu_thd^2*sta+sta*(mtot*tau+2*Fy*mo*Rcm*cta))/(mtot*J1);    
    mu_up(:,i+1) = mu_pr(:,i+1)+Sigma_pr(:,:,i+1)*Ct'*inv(Ct*Sigma_pr(:,:,i+1)*Ct'+R)*(y(:,i+1)-Ct*mu_pr(:,i+1));
    Sigma_up(:,:,i+1) = Sigma_pr(:,:,i+1)-Sigma_pr(:,:,i+1)*Ct'*inv(Ct*Sigma_pr(:,:,i+1)*Ct'+R)*Ct*Sigma_pr(:,:,i+1);
    
    % Predict measurement
    mu_pxd = mu_up(1,i+1); mu_pyd = mu_up(2,i+1); mu_thd = mu_up(3,i+1); mu_px = mu_up(4,i+1); mu_py = mu_up(5,i+1); mu_th = mu_up(6,i+1);
    cta = cos(mu_th+alpha); sta = sin(mu_th+alpha);
    y_up(:,i+1) = [mu_px
                   mu_py
                   mu_th
                   Fx/mtot+mo*Rcm*mu_thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1)
                   Fy/mtot+mo*Rcm*mu_thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1)
                   mu_thd];
end

for i = 1:ns
    figure
    plot(t,x(i,:)); hold on
    plot(t,mu_up(i,:),'r')
    title(['Estimating state ' num2str(i)],'fontsize',14)
end

for i = 1:ns
    figure
    plot(t,y(i,:)); hold on
    plot(t,y_up(i,:),'r')
    title(['Estimating measurement ' num2str(i)],'fontsize',14)
end