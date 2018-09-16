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


[problem,guess]=UAVSim;  		% Fetch the problem definition

options= settings_h(30);                  % Get options and solver settings
% options= settings_hp(5,4);                  % Get options and solver settings
%options= settings_hp([4 5 3],[-1 0.3 0.4 1]);                  % Get options and solver settings

[infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
[solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
[solution]=output(problem,solution,options,data,4);          % Output solutions


figure('name','Transmission Powers','position', [200, 200, 700, 300]);
plot(solution.T(:,1),speval(solution.Up,1,solution.T),'-o','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Up,2,solution.T),'-o','Linewidth',2);
plot(solution.T(:,1),speval(solution.Up,3,solution.T),'-o','Linewidth',2);
plot(solution.T(:,1),speval(solution.Up,4,solution.T),'-o','Linewidth',2);
ylabel('Transmission Powers ($W$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'$p_1$','$p_2$','$p_3$','$p_4$'},'Interpreter','latex', 'FontSize',18);

figure('name','Transmission Rates','position', [200, 200, 700, 300]);
plot(solution.T(:,1),speval(solution.Up,5,solution.T),'-o','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Up,6,solution.T),'-o','Linewidth',2);
plot(solution.T(:,1),speval(solution.Up,7,solution.T),'-o','Linewidth',2);
plot(solution.T(:,1),speval(solution.Up,8,solution.T),'-o','Linewidth',2);
ylabel('Transmission Rates ($bps$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'$r_1$','$r_2$','$r_3$','$r_4$'},'Interpreter','latex', 'FontSize',18);

figure('name','Acceleration','position', [200, 200, 700, 300]);
plot(solution.T(:,1),speval(solution.Up,9,solution.T),'-o','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Up,10,solution.T),'-o','Linewidth',2);
ylabel('Acceleration ($ms^{-2}$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'UAV$_1$','UAV$_2$'},'Interpreter','latex', 'FontSize',18);

figure('name','Storage Buffers','position', [200, 200, 700, 300]);
plot(solution.T(:,1),speval(solution.Xp,1,solution.T),'-o','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Xp,2,solution.T),'-o','Linewidth',2);
plot(solution.T(:,1),speval(solution.Xp,3,solution.T),'-o','Linewidth',2);
ylabel('Storage Buffers ($B$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'Source','UAV$_1$','UAV$_2$'},'Interpreter','latex', 'FontSize',18);

figure('name','UAV Velocities','position', [200, 200, 700, 300]);
plot(solution.T(:,1),speval(solution.Xp,4,solution.T),'-o','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Xp,5,solution.T),'-o','Linewidth',2);
ylabel('Velocities ($ms^{-1}$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'UAV$_1$','UAV$_2$'},'Interpreter','latex', 'FontSize',18);

figure('name','UAV Positions','position', [200, 200, 700, 300]);
plot(solution.T(:,1),speval(solution.Xp,6,solution.T),'-o','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Xp,7,solution.T),'-o','Linewidth',2);
ylabel('Position ($km$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'UAV$_1$','UAV$_2$'},'Interpreter','latex', 'FontSize',18);

