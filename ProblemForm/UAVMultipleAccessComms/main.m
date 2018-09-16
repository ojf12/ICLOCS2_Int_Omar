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


[problem,guess]=UAVMultipleAccessComms;  		% Fetch the problem definition

 options= settings_h(20);                  % Get options and solver settings 
% options= settings_hp(5,4);                  % Get options and solver settings 
%options= settings_hp([4 5 3],[-1 0.3 0.4 1]);                  % Get options and solver settings 

[infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
[solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
[solution]=output(problem,solution,options,data,4);          % Output solutions


%% Plot States
xx=linspace(solution.T(1,1),solution.T(end,1),1000);
figure; 
plot(solution.T(:,1),speval(solution.Xp,1,solution.T(:,1)),'ro' ); hold on
plot(solution.T(:,1),speval(solution.Xp,2,solution.T(:,1)),'bo' )
plot(xx,speval(solution.Xp,1,xx),'r-' )
plot(xx,speval(solution.Xp,2,xx),'b-' )
xlabel('Time [s]')
ylabel('Storage Buffer States')
legend('s1','s2')
grid on

%% Plot Inputs
figure
plot(solution.T(:,1),speval(solution.Up,1,solution.T),'-bo' ); hold on
plot(solution.T(:,1),speval(solution.Up,3,solution.T),'-ko' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)],'k-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(3), problem.inputs.uu(3)],'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
grid on
ylabel('Control Input Transmission Powers')
legend('p1','p2')

figure
plot(solution.T(:,1),speval(solution.Up,2,solution.T),'-bo' ); hold on
plot(solution.T(:,1),speval(solution.Up,4,solution.T),'-ko' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(2), problem.inputs.ul(2)],'k-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(4), problem.inputs.uu(4)],'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
grid on
ylabel('Control Input Transmission Rates')
legend('r1','r2')