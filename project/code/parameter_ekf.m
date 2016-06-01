%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameter_ekf.m
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

T = 50;
dt = 0.005;
Nsteps = T/dt;
t = 0:dt:T;

%% Parameter EKF Setup
xp0 = zeros(np, 1);
Qp = diag([0.001*ones(4,1)]);
R = diag([0.001*ones(3,1); 0.01*ones(3,1)]);
mup0 = xp0;
Sigmap0 = 0.1*eye(np);
    
mup_up = zeros(np, Nsteps+1);
mup_pr = zeros(np, Nsteps+1);
Sigmap_up = zeros(np,np,Nsteps+1);
Sigmap_pr = zeros(np,np,Nsteps+1);
mup_up(:,1) = mup0;
mup_pr(:,1) = mup0;
Sigmap_up(:,:,1) = Sigmap0;
Sigmap_pr(:,:,1) = Sigmap0;
y = zeros(m,Nsteps+1);
xp = zeros(np, Nsteps+1);
xp(:,1) = xp0;
y_up = zeros(m,Nsteps+1);

Apt = zeros(np,np);
Cpt = zeros(m,np);

for i = 1:Nsteps
    % Input
    Fx = 0.2*sin(t(i)/5);
    Fy = 0.2*cos(t(i)/10);
    tau = 0.2*cos(t(i)/7);
    
    % Predict next state
    mu_mo = mup_up(1,i);
    mu_jo = mup_up(2,i);
    mu_rcm = mup_up(3,i);
    mu_alp = mup_up(4,i);
    th = 
    cta = cos(mu_th+alpha);
    sta = sin(mu_th+alpha);
    mup_pr(:,i+1) = [mu_pxd+dt*(Fx/mtot+mo*Rcm*mu_thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1))
                    mu_pyd+dt*(Fy/mtot+mo*Rcm*mu_thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1))
                    mu_thd+dt*(mtot*tau+mo*Rcm*(Fy*cta+Fx*sta))/J1
                    mu_pxd+dt*mu_pxd
                    mu_pyd+dt*mu_pyd
                    mu_thd+dt*mu_thd];
    Apt = eye(np);
    Apt(1,3) = 2*dt/mtot*mo*Rcm*mu_thd*cta;
    Apt(2,3) = 2*dt/mtot*mo*Rcm*mu_thd*sta;
    Apt(1,6) = -dt*mo*Rcm*(-Fy*mo*Rcm*cos(2*(mu_th+alpha))+J1*mu_thd^2*sta-cta*(mtot*tau+2*Fx*mo*Rcm*sta))/(mtot*J1);
    Apt(2,6) =  dt*mo*Rcm*(-Fx*mo*Rcm*cos(2*(mu_th+alpha))+J1*mu_thd^2*sta+sta*(mtot*tau+2*Fy*mo*Rcm*cta))/(mtot*J1);    
    Apt(3,6) =  dt*mo*Rcm/J1*(Fx*cta-Fy*sta);
    Sigmap_pr(:,:,i+1) = Apt*Sigmap_up(:,:,i)*Apt'+Qp;
     
    % Simulate the satellite
    pxd = xp(1,i); pyd = xp(2,i); thd = xp(3,i); px = xp(4,i); py = xp(5,i); th = xp(6,i);
    cta = cos(th+alpha); sta = sin(th+alpha);
    xp(:,i+1) = [pxd+dt*(Fx/mtot+mo*Rcm*thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1))
                pyd+dt*(Fy/mtot+mo*Rcm*thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1))
                thd+dt*(mtot*tau+mo*Rcm*(Fy*cta+Fx*sta))/J1
                pxd+dt*pxd
                pyd+dt*pyd
                thd+dt*thd] + mvnrnd(zeros(np,1),Qp*dt^2)';
            
    % Take measurement
    pxd = xp(1,i+1); pyd = xp(2,i+1); thd = xp(3,i+1); px = xp(4,i+1); py = xp(5,i+1); th = xp(6,i+1);
    cta = cos(th+alpha); sta = sin(th+alpha);
    y(:,i+1) = [px
                py
                th
                Fx/mtot+mo*Rcm*thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1)
                Fy/mtot+mo*Rcm*thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1)
                thd] + mvnrnd(zeros(m,1),R)';
    
    % Update state estimate
    mu_pxd = mup_pr(1,i+1); mu_pyd = mup_pr(2,i+1); mu_thd = mup_pr(3,i+1); mu_px = mup_pr(4,i+1); mu_py = mup_pr(5,i+1); mu_th = mup_pr(6,i+1);
    cta = cos(th+alpha); sta = sin(th+alpha);
    Cpt(1,4) = 1; Cpt(2,5) = 1; Cpt(3,6) = 1; Cpt(6,3) = 1;
    Cpt(4,6) = -mo*Rcm*(-Fy*mo*Rcm*cos(2*(mu_th+alpha))+J1*mu_thd^2*sta-cta*(mtot*tau+2*Fx*mo*Rcm*sta))/(mtot*J1);
    Cpt(5,6) =  mo*Rcm*(-Fx*mo*Rcm*cos(2*(mu_th+alpha))+J1*mu_thd^2*sta+sta*(mtot*tau+2*Fy*mo*Rcm*cta))/(mtot*J1);    
    mup_up(:,i+1) = mup_pr(:,i+1)+Sigmap_pr(:,:,i+1)*Cpt'*inv(Cpt*Sigmap_pr(:,:,i+1)*Cpt'+R)*(y(:,i+1)-Cpt*mup_pr(:,i+1));
    Sigmap_up(:,:,i+1) = Sigmap_pr(:,:,i+1)-Sigmap_pr(:,:,i+1)*Cpt'*inv(Cpt*Sigmap_pr(:,:,i+1)*Cpt'+R)*Cpt*Sigmap_pr(:,:,i+1);
    
    % Predict measurement
    mu_pxd = mup_up(1,i+1); mu_pyd = mup_up(2,i+1); mu_thd = mup_up(3,i+1); mu_px = mup_up(4,i+1); mu_py = mup_up(5,i+1); mu_th = mup_up(6,i+1);
    cta = cos(mu_th+alpha); sta = sin(mu_th+alpha);
    y_up(:,i+1) = [mu_px
                   mu_py
                   mu_th
                   Fx/mtot+mo*Rcm*mu_thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1)
                   Fy/mtot+mo*Rcm*mu_thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1)
                   mu_thd];
end

for i = 1:np
    figure
    plot(t,xp(i,:)); hold on
    plot(t,mup_up(i,:),'r')
    title(['Estimating state ' num2str(i)],'fontsize',14)
end

for i = 1:np
    figure
    plot(t,y(i,:)); hold on
    plot(t,y_up(i,:),'r')
    title(['Estimating measurement ' num2str(i)],'fontsize',14)
end