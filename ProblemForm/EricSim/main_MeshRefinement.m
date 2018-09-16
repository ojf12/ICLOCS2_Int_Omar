% MAIN_MESHREFINEMENT - Main script to solve the Optimal Control Problem
% with mesh refinement
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

clear all;close all; clc; format compact;

global sol;  
sol=[];                             % Initialize solution structure

[problem,guess]=UAVSim;  		% Fetch the problem definition


%%%% Choice among the ones belwo %%%%
% options= settings_auto(40);                  
options= settings_h(30);                
% options= settings_hp(5,4);                 

% Declare variables for storage
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9; % the initial (before start) maximum absolute local eeror
i=1; % starting iteration
imax=3; % maximum number of mesh refinement iterations

while (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax    
    [infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
    [solution,infoNLP1,data] = solveNLP(infoNLP,data);      % Solve the NLP
    [solution]=output(problem,solution,options,data,4);         % Output solutions
    
    maxAbsError=max(abs(solution.Error));
    maxAbsConstraintError=max(solution.ConstraintError);
    errorHistory(i,:)=maxAbsError;
    ConstraintErrorHistory(i,:)=maxAbsConstraintError;
    timeHistory(i)=solution.computation_time;
    solutionHistory{i}=solution;
    
    if (any(maxAbsError>problem.states.xErrorTol) || any(maxAbsConstraintError>problem.constraintErrorTol)) && i<=imax
        [ options, guess ] = doMeshRefinement( options, problem, guess, data, solution, i );
    end
    i=i+1;
end

MeshRefinementHistory.errorHistory=errorHistory;
MeshRefinementHistory.timeHistory=timeHistory;
MeshRefinementHistory.ConstraintErrorHistory=ConstraintErrorHistory;



figure('name','Transmission Powers','position', [200, 200, 700, 300]);set(gcf,'color','w');
plot(solution.T(:,1),speval(solution.Up,1,solution.T),'-','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Up,2,solution.T),'-','Linewidth',2);
plot(solution.T(:,1),speval(solution.Up,3,solution.T),'-','Linewidth',2);
plot(solution.T(:,1),speval(solution.Up,4,solution.T),'-','Linewidth',2);
ylabel('Transmission Powers ($W$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'$p_1$','$p_2$','$p_3$','$p_4$'},'Interpreter','latex', 'FontSize',18);

figure('name','Transmission Rates','position', [200, 200, 700, 300]);set(gcf,'color','w');
plot(solution.T(:,1),speval(solution.Up,5,solution.T),'-','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Up,6,solution.T),'-','Linewidth',2);
plot(solution.T(:,1),speval(solution.Up,7,solution.T),'-','Linewidth',2);
plot(solution.T(:,1),speval(solution.Up,8,solution.T),'-','Linewidth',2);
ylabel('Transmission Rates ($bps$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'$r_1$','$r_2$','$r_3$','$r_4$'},'Interpreter','latex', 'FontSize',18);

figure('name','Acceleration','position', [200, 200, 700, 300]);set(gcf,'color','w');
plot(solution.T(:,1),speval(solution.Up,9,solution.T),'-','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Up,10,solution.T),'-','Linewidth',2);
ylabel('Acceleration ($ms^{-2}$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'UAV$_1$','UAV$_2$'},'Interpreter','latex', 'FontSize',18);

figure('name','Storage Buffers','position', [200, 200, 700, 300]);set(gcf,'color','w');
plot(solution.T(:,1),speval(solution.Xp,1,solution.T),'-','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Xp,2,solution.T),'-','Linewidth',2);
plot(solution.T(:,1),speval(solution.Xp,3,solution.T),'-','Linewidth',2);
ylabel('Storage Buffers ($B$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'Source','UAV$_1$','UAV$_2$'},'Interpreter','latex', 'FontSize',18);

figure('name','UAV Velocities','position', [200, 200, 700, 300]);set(gcf,'color','w');
plot(solution.T(:,1),speval(solution.Xp,4,solution.T),'-','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Xp,5,solution.T),'-','Linewidth',2);
ylabel('Velocities ($ms^{-1}$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'UAV$_1$','UAV$_2$'},'Interpreter','latex', 'FontSize',18);

figure('name','UAV Positions','position', [200, 200, 700, 300]);set(gcf,'color','w');
plot(solution.T(:,1),speval(solution.Xp,6,solution.T),'-','Linewidth',2);  hold on; grid on;
plot(solution.T(:,1),speval(solution.Xp,7,solution.T),'-','Linewidth',2);
ylabel('Position ($km$)','Interpreter','latex', 'FontSize',18);
xlabel('Time ($s$)','Interpreter','latex', 'FontSize',18);
legend({'UAV$_1$','UAV$_2$'},'Interpreter','latex', 'FontSize',18);

