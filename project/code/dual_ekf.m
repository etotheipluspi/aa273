%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project_workspace.m
% Author: Andrew Bylard
% Stanford University
% May 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

ns = 6;     % Dimension of state
np = 4;     % Dimension of parameters
m = 6;      % Dimension of measurements
% m = 7;      % Dimension of measurements

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

mo_min = 0.1;
Jo_min = 0.05;
Rcm_min = rf+L;

mtot = mf+mo;       % [kg] Total mass
Jtot = Jf+Jo;       % [kg*m^2] Summed moments of inertia
J1 = Jtot*mtot+mf*mo*Rcm^2;     % [kg^2*m^2] Common expression

T = 70;
dt = 0.005;
Nsteps = T/dt;
t = 0:dt:T;

%% State EKF Setup
xs0 = zeros(ns, 1);
Qs = diag([0.001*ones(3,1); 0.01*ones(3,1)]);
R = diag([0.1*ones(3,1); 0.1*ones(3,1)]);
mus0 = xs0;
Sigmas0 = 0.1*eye(ns);
    
mus_up = zeros(ns, Nsteps+1);
mus_pr = zeros(ns, Nsteps+1);
Sigmas_up = zeros(ns,ns,Nsteps+1);
Sigmas_pr = zeros(ns,ns,Nsteps+1);
mus_up(:,1) = mus0;
mus_pr(:,1) = mus0;
Sigmas_up(:,:,1) = Sigmas0;
Sigmas_pr(:,:,1) = Sigmas0;
xs = zeros(ns, Nsteps+1);
xs(:,1) = xs0;
ys_up = zeros(m,Nsteps+1);

Ast = zeros(ns,ns);
Cst = zeros(m,ns);

%% Parameter EKF Setup
xp0 = [mo
       Jo
       Rcm
       0];
Qp = diag([0.05*ones(4,1)]);
mup0 = [mf
       Jf
       rf+L
       0];
mup0 = xp0;
Sigmap0 = 0.9*eye(np);
    
mup_up = zeros(np, Nsteps+1);
mup_pr = zeros(np, Nsteps+1);
Sigmap_up = zeros(np,np,Nsteps+1);
Sigmap_pr = zeros(np,np,Nsteps+1);
mup_up(:,1) = mup0;
mup_pr(:,1) = mup0;
Sigmap_up(:,:,1) = Sigmap0;
Sigmap_pr(:,:,1) = Sigmap0;
xp = zeros(np, Nsteps+1);
xp(:,1) = xp0;
yp_up = zeros(m,Nsteps+1);

Apt = zeros(np,np);
Cpt = zeros(m,np);

y = zeros(m,Nsteps+1);

%% Dual EKF
for i = 1:Nsteps
    % Parameter estimation flag
    dekf_on = 1;
    
    % Input
    a = 0.5;
    Fx = a*sin(t(i)/5);
    Fy = a*cos(t(i)/10);
    tau = a*cos(t(i)/7);
    
    % Predict parameters
    mu_pxd = mus_up(1,i); mu_pyd = mus_up(2,i); mu_thd = mus_up(3,i); mu_px = mus_up(4,i); mu_py = mus_up(5,i); mu_th = mus_up(6,i);
    mu_mo = mup_up(1,i); mu_jo = mup_up(2,i); mu_rcm = mup_up(3,i); mu_alp = mup_up(4,i);
    cta = cos(mu_th+mu_alp); sta = sin(mu_th+mu_alp); c2ta = cos(mu_th+mu_alp);
    mu_mtot = mf+mu_mo;
    mu_Jtot = Jf+mu_jo;
    mu_J1 = mu_Jtot*mu_mtot+mf*mu_mo*mu_rcm^2;
    mu_J2 = mu_Jtot*mu_mtot-mf*mu_mo*mu_rcm^2;
    if dekf_on
        mup_pr(:,i+1) = mup_up(:,i);
        Apt(1,1) = dt*((-Fx+mf*mu_rcm*mu_thd^2*cta)/mu_mtot^2+mu_Jtot*mf*mu_rcm*tau*sta/mu_J1^2+mf*mu_mo*mu_rcm^2*(mu_J1+mu_Jtot*mu_mtot)*sta*(Fy*cta+Fx*sta)/(mu_mtot^2*mu_J1^2));
        Apt(2,1) = dt*((-Fy+mf*mu_rcm*mu_thd^2*sta)/mu_mtot^2-mu_Jtot*mf*mu_rcm*tau*cta/mu_J1^2-mf*mu_mo*mu_rcm^2*(mu_J1+mu_Jtot*mu_mtot)*cta*(Fy*cta+Fx*sta)/(mu_mtot^2*mu_J1^2));
        Apt(1,2) = -dt*mu_mo*mu_rcm*sta/mu_J1^2*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta));
        Apt(2,2) = -dt*mu_mo*mu_rcm*cta/mu_J1^2*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta));
        Apt(1,3) = dt*mu_mo*( (mu_J2*tau+2*Fx*mu_Jtot*mu_mo*mu_rcm*sta)*sta/mu_J1^2+(mu_thd^2+2*Fy*mu_Jtot*mu_mtot*mu_mo*mu_rcm*sta/mu_J1^2)*cta/mu_mtot);
        Apt(2,3) = dt*mu_mo*(-(mu_J2*tau+2*Fx*mu_Jtot*mu_mo*mu_rcm*sta)*cta/mu_J1^2+(mu_thd^2*sta+2*Fy*mu_Jtot*mu_mtot*mu_mo*mu_rcm*cta^2/mu_J1^2)/mu_mtot);
        Apt(1,4) = -dt*mu_mo*mu_rcm/(mu_mtot*mu_J1)*(-Fy*mu_mo*mu_rcm*c2ta+mu_J1*mu_thd^2*sta-(mu_mtot*tau+2*Fx*mu_mo*mu_rcm*sta)*cta);
        Apt(2,4) =  dt*mu_mo*mu_rcm/(mu_mtot*mu_J1)*(-Fx*mu_mo*mu_rcm*c2ta+mu_J1*mu_thd^2*cta+(mu_mtot*tau+2*Fx*mu_mo*mu_rcm*cta)*sta);
        Apt(3,1) = -dt*mf*mu_rcm/mu_J1^2*(mf*mu_rcm*tau-Fy*mu_Jtot*cta-Fx*mu_Jtot*sta);
        Apt(3,2) = -dt*mu_mtot/mu_J1^2*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta));
        Apt(3,3) = dt*mu_mo/mu_J1^2*(mu_J1*(Fy*cta+Fx*sta)-2*mf*mu_rcm*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta)));
        Apt(3,4) = dt*mu_mo*mu_rcm/mu_J1*(Fx*cta-Fy*sta);
        Sigmap_pr(:,:,i+1) = Apt*Sigmap_up(:,:,i)*Apt'+Qp;

        mup_pr(1,i+1) = max(mo_min, mup_pr(1,i+1));
        mup_pr(2,i+1) = max(Jo_min, mup_pr(2,i+1));
        mup_pr(3,i+1) = max(Rcm_min, mup_pr(3,i+1));
        mup_pr(4,i+1) = 0;
    else
        mup_pr(1,i+1) = mo;
        mup_pr(2,i+1) = Jo;
        mup_pr(3,i+1) = Rcm;
        mup_pr(4,i+1) = 0;
    end
    
    % Predict next state
    mus_pr(:,i+1) = [mu_pxd+dt*(Fx/mu_mtot+mu_mo*mu_rcm*mu_thd^2*cta/mu_mtot+mu_mo*mu_rcm*tau*sta/mu_J1+mu_mo^2*mu_rcm^2*(Fy*cta*sta+Fx*sta^2)/(mu_mtot*mu_J1))
                     mu_pyd+dt*(Fy/mu_mtot+mu_mo*mu_rcm*mu_thd^2*sta/mu_mtot-mu_mo*mu_rcm*tau*cta/mu_J1-mu_mo^2*mu_rcm^2*(Fx*cta*sta+Fy*cta^2)/(mu_mtot*mu_J1))
                     mu_thd+dt*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta))/mu_J1
                     mu_px+dt*mu_pxd
                     mu_py+dt*mu_pyd
                     mu_th+dt*mu_thd];
    Ast = eye(ns);
    Ast(1,3) = 2*dt/mu_mtot*mu_mo*mu_rcm*mu_thd*cta;
    Ast(2,3) = 2*dt/mu_mtot*mu_mo*mu_rcm*mu_thd*sta;
    Ast(1,6) = -dt*mu_mo*mu_rcm*(-Fy*mu_mo*mu_rcm*cos(2*(mu_th+mu_alp))+mu_J1*mu_thd^2*sta-cta*(mu_mtot*tau+2*Fx*mu_mo*mu_rcm*sta))/(mu_mtot*mu_J1);
    Ast(2,6) =  dt*mu_mo*mu_rcm*(-Fx*mu_mo*mu_rcm*cos(2*(mu_th+mu_alp))+mu_J1*mu_thd^2*sta+sta*(mu_mtot*tau+2*Fy*mu_mo*mu_rcm*cta))/(mu_mtot*mu_J1);    
    Ast(3,6) =  dt*mu_mo*mu_rcm/mu_J1*(Fx*cta-Fy*sta);
    Sigmas_pr(:,:,i+1) = Ast*Sigmas_up(:,:,i)*Ast'+Qs;
    
    % Simulate the free-flyer
    pxd = xs(1,i); pyd = xs(2,i); thd = xs(3,i); px = xs(4,i); py = xs(5,i); th = xs(6,i);
    cta = cos(th+alpha); sta = sin(th+alpha);
    xs(:,i+1) = [pxd+dt*(Fx/mtot+mo*Rcm*thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1))
                 pyd+dt*(Fy/mtot+mo*Rcm*thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1))
                 thd+dt*(mtot*tau+mo*Rcm*(Fy*cta+Fx*sta))/J1
                 px+dt*pxd
                 py+dt*pyd
                 th+dt*thd] + mvnrnd(zeros(ns,1),Qs*dt^2)';
    xp(:,i+1) = xp0;
            
    % Take measurement
    pxd = xs(1,i+1); pyd = xs(2,i+1); thd = xs(3,i+1); px = xs(4,i+1); py = xs(5,i+1); th = xs(6,i+1);
    cta = cos(th+alpha); sta = sin(th+alpha);
    y(:,i+1) = [px
                py
                th
                Fx/mtot+mo*Rcm*thd^2*cta/mtot+mo*Rcm*tau*sta/J1+mo^2*Rcm^2*(Fy*cta*sta+Fx*sta^2)/(mtot*J1)
                Fy/mtot+mo*Rcm*thd^2*sta/mtot-mo*Rcm*tau*cta/J1-mo^2*Rcm^2*(Fx*cta*sta+Fy*cta^2)/(mtot*J1)
                thd] + mvnrnd(zeros(m,1),R)';
%                 (mtot*tau+mo*Rcm*(Fy*cta+Fx*sta))/J1] + mvnrnd(zeros(m,1),R)';
    
    % Update parameter estimates
    mu_pxd = mus_pr(1,i+1); mu_pyd = mus_pr(2,i+1); mu_thd = mus_pr(3,i+1); mu_px = mus_pr(4,i+1); mu_py = mus_pr(5,i+1); mu_th = mus_pr(6,i+1);
    mu_mo = mup_pr(1,i+1); mu_jo = mup_pr(2,i+1); mu_rcm = mup_pr(3,i+1); mu_alp = mup_pr(4,i+1);
    cta = cos(mu_th+mu_alp); sta = sin(mu_th+mu_alp); c2ta = cos(mu_th+mu_alp);
    mu_mtot = mf+mu_mo;
    mu_Jtot = Jf+mu_jo;
    mu_J1 = mu_Jtot*mu_mtot+mf*mu_mo*mu_rcm^2;
    
    if dekf_on
        Cpt(4,1) = ((-Fx+mf*mu_rcm*mu_thd^2*cta)/mu_mtot^2+mu_Jtot*mf*mu_rcm*tau*sta/mu_J1^2+mf*mu_mo*mu_rcm^2*(mu_J1+mu_Jtot*mu_mtot)*sta*(Fy*cta+Fx*sta)/(mu_mtot^2*mu_J1^2));
        Cpt(5,1) = ((-Fy+mf*mu_rcm*mu_thd^2*sta)/mu_mtot^2-mu_Jtot*mf*mu_rcm*tau*cta/mu_J1^2-mf*mu_mo*mu_rcm^2*(mu_J1+mu_Jtot*mu_mtot)*cta*(Fy*cta+Fx*sta)/(mu_mtot^2*mu_J1^2));
        Cpt(4,2) = -mu_mo*mu_rcm*sta/mu_J1^2*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta));
        Cpt(5,2) = -mu_mo*mu_rcm*cta/mu_J1^2*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta));
        Cpt(4,3) = mu_mo*( (mu_J2*tau+2*Fx*mu_Jtot*mu_mo*mu_rcm*sta)*sta/mu_J1^2+(mu_thd^2+2*Fy*mu_Jtot*mu_mtot*mu_mo*mu_rcm*sta/mu_J1^2)*cta/mu_mtot);
        Cpt(5,3) = mu_mo*(-(mu_J2*tau+2*Fx*mu_Jtot*mu_mo*mu_rcm*sta)*cta/mu_J1^2+(mu_thd^2*sta+2*Fy*mu_Jtot*mu_mtot*mu_mo*mu_rcm*cta^2/mu_J1^2)/mu_mtot);
        Cpt(4,4) = -mu_mo*mu_rcm/(mu_mtot*mu_J1)*(-Fy*mu_mo*mu_rcm*c2ta+mu_J1*mu_thd^2*sta-(mu_mtot*tau+2*Fx*mu_mo*mu_rcm*sta)*cta);
        Cpt(5,4) =  mu_mo*mu_rcm/(mu_mtot*mu_J1)*(-Fx*mu_mo*mu_rcm*c2ta+mu_J1*mu_thd^2*cta+(mu_mtot*tau+2*Fx*mu_mo*mu_rcm*cta)*sta);
%         Cpt(7,1) = -mf*mu_rcm/mu_J1^2*(mf*mu_rcm*tau-Fy*mu_Jtot*cta-Fx*mu_Jtot*sta);
%         Cpt(7,2) = -mu_mtot/mu_J1^2*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta));
%         Cpt(7,3) = mu_mo/mu_J1^2*(mu_J1*(Fy*cta+Fx*sta)-2*mf*mu_rcm*(mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta)));
%         Cpt(7,4) = mu_mo*mu_rcm/mu_J1*(Fx*cta-Fy*sta);
        Cst(1,4) = 1; Cst(2,5) = 1; Cst(3,6) = 1; Cst(6,3) = 1;
        Cst(4,6) = -mu_mo*mu_rcm*(-Fy*mu_mo*mu_rcm*cos(2*(mu_th+mu_alp))+mu_J1*mu_thd^2*sta-cta*(mu_mtot*tau+2*Fx*mu_mo*mu_rcm*sta))/(mu_mtot*mu_J1);
        Cst(5,6) =  mu_mo*mu_rcm*(-Fx*mu_mo*mu_rcm*cos(2*(mu_th+mu_alp))+mu_J1*mu_thd^2*sta+sta*(mu_mtot*tau+2*Fy*mu_mo*mu_rcm*cta))/(mu_mtot*mu_J1);
%         Cst(7,6) =  mu_mo*mu_rcm/mu_J1*(Fx*cta-Fy*sta);
        mup_up(:,i+1) = mup_pr(:,i+1)+Sigmap_pr(:,:,i+1)*Cpt'*inv(Cpt*Sigmap_pr(:,:,i+1)*Cpt'+R)*(y(:,i+1)-[Cst Cpt]*[mus_pr(:,i+1); mup_pr(:,i+1)]);
        if i > 463
            i;
        end
        Sigmap_up(:,:,i+1) = Sigmap_pr(:,:,i+1)-Sigmap_pr(:,:,i+1)*Cpt'*inv(Cpt*Sigmap_pr(:,:,i+1)*Cpt'+R)*Cpt*Sigmap_pr(:,:,i+1);

        mup_up(1,i+1) = max(mo_min, mup_up(1,i+1));
        mup_up(2,i+1) = max(Jo_min, mup_up(2,i+1));
        mup_up(3,i+1) = max(Rcm_min, mup_up(3,i+1));
        mup_up(4,i+1) = 0;
    else
        mup_up(1,i+1) = mo;
        mup_up(2,i+1) = Jo;
        mup_up(3,i+1) = Rcm;
        mup_up(4,i+1) = 0;
    end
    
    % Update state estimates
    mu_mo = mup_up(1,i+1); mu_jo = mup_up(2,i+1); mu_rcm = mup_up(3,i+1); mu_alp = mup_up(4,i+1);
    cta = cos(mu_th+mu_alp); sta = sin(mu_th+mu_alp); c2ta = cos(mu_th+mu_alp);
    mu_mtot = mf+mu_mo;
    mu_Jtot = Jf+mu_jo;
    mu_J1 = mu_Jtot*mu_mtot+mf*mu_mo*mu_rcm^2;
    
    Cst(1,4) = 1; Cst(2,5) = 1; Cst(3,6) = 1; Cst(6,3) = 1;
    Cst(4,6) = -mu_mo*mu_rcm*(-Fy*mu_mo*mu_rcm*cos(2*(mu_th+mu_alp))+mu_J1*mu_thd^2*sta-cta*(mu_mtot*tau+2*Fx*mu_mo*mu_rcm*sta))/(mu_mtot*mu_J1);
    Cst(5,6) =  mu_mo*mu_rcm*(-Fx*mu_mo*mu_rcm*cos(2*(mu_th+mu_alp))+mu_J1*mu_thd^2*sta+sta*(mu_mtot*tau+2*Fy*mu_mo*mu_rcm*cta))/(mu_mtot*mu_J1);
%     Cst(7,6) =  mu_mo*mu_rcm/mu_J1*(Fx*cta-Fy*sta);
    mus_up(:,i+1) = mus_pr(:,i+1)+Sigmas_pr(:,:,i+1)*Cst'*inv(Cst*Sigmas_pr(:,:,i+1)*Cst'+R)*(y(:,i+1)-Cst*mus_pr(:,i+1));
    Sigmas_up(:,:,i+1) = Sigmas_pr(:,:,i+1)-Sigmas_pr(:,:,i+1)*Cst'*inv(Cst*Sigmas_pr(:,:,i+1)*Cst'+R)*Cst*Sigmas_pr(:,:,i+1);
    
    % Predict measurement
    mu_pxd = mus_up(1,i+1); mu_pyd = mus_up(2,i+1); mu_thd = mus_up(3,i+1); mu_px = mus_up(4,i+1); mu_py = mus_up(5,i+1); mu_th = mus_up(6,i+1);
    cta = cos(mu_th+mu_alp); sta = sin(mu_th+mu_alp); c2ta = cos(mu_th+mu_alp);
    ys_up(:,i+1) = [mu_px
                    mu_py
                    mu_th
                    Fx/mu_mtot+mu_mo*mu_rcm*mu_thd^2*cta/mu_mtot+mu_mo*mu_rcm*tau*sta/mu_J1+mu_mo^2*mu_rcm^2*(Fy*cta*sta+Fx*sta^2)/(mu_mtot*mu_J1)
                    Fy/mu_mtot+mu_mo*mu_rcm*mu_thd^2*sta/mu_mtot-mu_mo*mu_rcm*tau*cta/mu_J1-mu_mo^2*mu_rcm^2*(Fx*cta*sta+Fy*cta^2)/(mu_mtot*mu_J1)
                    mu_thd];
%                     (mu_mtot*tau+mu_mo*mu_rcm*(Fy*cta+Fx*sta))/mu_J1];
end

for i = 1:ns
    figure
    plot(t,xs(i,:)); hold on
    plot(t,mus_up(i,:),'r')
    title(['Estimating state ' num2str(i)],'fontsize',14)
end

for i = 1:ns
    figure
    plot(t,y(i,:)); hold on
    plot(t,ys_up(i,:),'r')
    title(['Estimating measurement ' num2str(i)],'fontsize',14)
end

if dekf_on
    for i = 1:np
        figure
        plot(t,xp(i,:)); hold on
        plot(t,mup_up(i,:),'r')
        title(['Estimating parameter ' num2str(i)],'fontsize',14)
    end
end