% MAIN - Main script to solve the Optimal Control Problem
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK
% ICLOCS (Imperial College London Optimal Control) Version 2.0
% 1 May 2018
% iclocs@imperial.ac.uk


%--------------------------------------------------------

clear all; close all; clc; format compact;


[problem,guess]=UAVChain;  		% Fetch the problem definition

options= settings_h(100);                  % Get options and solver settings
% options= settings_hp(5,4);                  % Get options and solver settings
%options= settings_hp([4 5 3],[-1 0.3 0.4 1]);                  % Get options and solver settings

[infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
[solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
[solution]=output(problem,solution,options,data,4);          % Output solutions

%% Sanity check
q1 = solution.X(:,4);
q2 = solution.X(:,5);
q3 = solution.X(:,6);

p2 = solution.U(:,2);
p1 = solution.U(:,1);
p3 = solution.U(:,3);

AP = data.data.auxdata(12);
Disp = data.data.auxdata(13);
alt = data.data.auxdata(4);
B       = data.data.auxdata(1);
K       = data.data.auxdata(2);
alpha   = data.data.auxdata(3);
ss      = data.data.auxdata(5);
%distance 1
dist1 = sqrt((q1-q2).^2 + Disp^2);
dist2 = sqrt((q2-q3).^2 + Disp^2);
dist3 = sqrt((q3-AP).^2 + alt^2 + Disp^2);

figure;
hold on;
plot(dist1,'LineWidth',2);
plot(dist2,'LineWidth',2);
plot(dist3,'LineWidth',2);


ff1 = 100-(B)*log2(1 + (K.*p1)./(ss.*((q1-q2).^2 + Disp^2).^alpha));
ff2 = 100-(B)*log2(1 + (K.*p2)./(ss.*((q2-q3).^2 + Disp^2).^alpha));
ff3 = 100-(B)*log2(1 + (K.*p3)./(ss.*((q3-AP).^2 + alt^2 + Disp^2).^alpha));

figure; plot(ff1); hold on; plot(ff2); plot(ff3);
%% Plot States and Inputs
timePlot=1;
if(~timePlot)
    figure('name','Transmission Powers','position', [200, 200, 450, 300]);
    h1 = plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Up,1,solution.T),'b','Linewidth',2);  hold on; grid on;
    h2 = plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Up,2,solution.T),'r','Linewidth',2);
    h3 = plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Up,3,solution.T),'g','Linewidth',2);
    plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Up,1,solution.T),'bo','Linewidth',1);
    plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Up,2,solution.T),'ro','Linewidth',1);
    plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Up,3,solution.T),'go','Linewidth',1);
    ylabel('Transmission Powers ($W$)','Interpreter','latex', 'FontSize',18);
    xlabel('Position ($km$)','Interpreter','latex', 'FontSize',18);
    legend([h1,h2,h3],{'$p_1$','$p_2$','$p_3$'},'Interpreter','latex', 'FontSize',18);
    
    figure('name','Acceleration','position', [200, 200, 450, 300]);
    h1 = plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Up,4,solution.T),'b','Linewidth',2);  hold on; grid on;
    h2 = plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Up,5,solution.T),'r','Linewidth',2);
    h3 = plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Up,6,solution.T),'g','Linewidth',2);
    plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Up,4,solution.T),'bo','Linewidth',1);
    plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Up,5,solution.T),'ro','Linewidth',1);
    plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Up,6,solution.T),'go','Linewidth',1);
    ylabel('Acceleration $ms^{-2}$','Interpreter','latex', 'FontSize',18);
    xlabel('Position ($km$)','Interpreter','latex', 'FontSize',18);
    legend([h1,h2,h3],{'$a_1$','$a_2$','$a_3$'},'Interpreter','latex', 'FontSize',18);
    
    figure('name','Velocities','position', [200, 200, 450, 300]);
    h1 = plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Xp,1,solution.T),'b','Linewidth',2);  hold on; grid on;
    h2 = plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Xp,2,solution.T),'r','Linewidth',2);
    h3 = plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Xp,3,solution.T),'g','Linewidth',2);
    plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Xp,1,solution.T),'bo','Linewidth',1);
    plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Xp,2,solution.T),'ro','Linewidth',1);
    plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Xp,3,solution.T),'go','Linewidth',1);
    ylabel('UAV Velocities $ms^{-1}$','Interpreter','latex', 'FontSize',18);
    xlabel('Position ($km$)','Interpreter','latex', 'FontSize',18);
    legend([h1,h2,h3],{'$v_1$','$v_2$','$v_3$'},'Interpreter','latex', 'FontSize',18);
else
    figure('name','Transmission Powers','position', [200, 200, 450, 300]);
    h1 = plot(solution.T(:,1),speval(solution.Up,1,solution.T),'Linewidth',2);  hold on; grid on;
    h2 = plot(solution.T(:,1),speval(solution.Up,2,solution.T),'Linewidth',2);
    h3 = plot(solution.T(:,1),speval(solution.Up,3,solution.T),'Linewidth',2);
%     plot(solution.T(:,1),speval(solution.Up,1,solution.T),'bo','Linewidth',1);
%     plot(solution.T(:,1),speval(solution.Up,2,solution.T),'ro','Linewidth',1);
%     plot(solution.T(:,1),speval(solution.Up,3,solution.T),'go','Linewidth',1);
    ylabel('Transmission Powers ($W$)','Interpreter','latex', 'FontSize',18);
    xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
    legend([h1,h2,h3],{'$p_1$','$p_2$','$p_3$'},'Interpreter','latex', 'FontSize',18);
    
    figure('name','Acceleration','position', [200, 200, 450, 300]);
    h1 = plot(solution.T(:,1),speval(solution.Up,4,solution.T),'Linewidth',2);  hold on; grid on;
    h2 = plot(solution.T(:,1),speval(solution.Up,5,solution.T),'Linewidth',2);
    h3 = plot(solution.T(:,1),speval(solution.Up,6,solution.T),'Linewidth',2);
%     plot(solution.T(:,1),speval(solution.Up,4,solution.T),'bo','Linewidth',1);
%     plot(solution.T(:,1),speval(solution.Up,5,solution.T),'ro','Linewidth',1);
%     plot(solution.T(:,1),speval(solution.Up,6,solution.T),'go','Linewidth',1);
    ylabel('Acceleration ($ms^{-2}$)','Interpreter','latex', 'FontSize',18);
    xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
    legend([h1,h2,h3],{'$a_1$','$a_2$','$a_3$'},'Interpreter','latex', 'FontSize',18);
    
    figure('name','Velocities','position', [200, 200, 450, 300]);
    h1 = plot(solution.T(:,1),speval(solution.Xp,1,solution.T),'Linewidth',2);  hold on; grid on;
    h2 = plot(solution.T(:,1),speval(solution.Xp,2,solution.T),'Linewidth',2);
    h3 = plot(solution.T(:,1),speval(solution.Xp,3,solution.T),'Linewidth',2);
%     plot(solution.T(:,1),speval(solution.Xp,1,solution.T),'bo','Linewidth',1);
%     plot(solution.T(:,1),speval(solution.Xp,2,solution.T),'ro','Linewidth',1);
%     plot(solution.T(:,1),speval(solution.Xp,3,solution.T),'go','Linewidth',1);
    ylabel('UAV Velocities ($ms^{-1}$)','Interpreter','latex', 'FontSize',18);
    xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
    legend([h1,h2,h3],{'$v_1$','$v_2$','$v_3$'},'Interpreter','latex', 'FontSize',18);
   
   figure('name','Positions','position', [200, 200, 450, 300]);
    h1 = plot(solution.T(:,1),speval(solution.Xp,4,solution.T)./1000,'Linewidth',2);  hold on; grid on;
    h2 = plot(solution.T(:,1),speval(solution.Xp,5,solution.T)./1000,'Linewidth',2);
    h3 = plot(solution.T(:,1),speval(solution.Xp,6,solution.T)./1000,'Linewidth',2);
%     plot(solution.T(:,1),speval(solution.Xp,1,solution.T),'bo','Linewidth',1);
%     plot(solution.T(:,1),speval(solution.Xp,2,solution.T),'ro','Linewidth',1);
%     plot(solution.T(:,1),speval(solution.Xp,3,solution.T),'go','Linewidth',1);
    ylabel('UAV Positions ($km$)','Interpreter','latex', 'FontSize',18);
    xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
    legend([h1,h2,h3],{'$v_1$','$v_2$','$v_3$'},'Interpreter','latex', 'FontSize',18); 
    
end




