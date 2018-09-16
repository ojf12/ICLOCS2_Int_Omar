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

[problem,guess]=UAVChain;  		% Fetch the problem definition


%%%% Choice among the ones belwo %%%%
% options= settings_auto(40);                  
options= settings_h(20);                
% options= settings_hp(5,4);                 

% Declare variables for storage
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9; % the initial (before start) maximum absolute local eeror
i=1; % starting iteration
imax=50; % maximum number of mesh refinement iterations

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

%% Plot States and Inputs
figure('name','Transmission Powers','position', [200, 200, 700, 300]);
h1 = plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Up,1,solution.T),'b','Linewidth',2);  hold on; grid on;
h2 = plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Up,2,solution.T),'r','Linewidth',2);
h3 = plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Up,3,solution.T),'g','Linewidth',2);
plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Up,1,solution.T),'bo','Linewidth',1);
plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Up,2,solution.T),'ro','Linewidth',1);
plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Up,3,solution.T),'go','Linewidth',1);
ylabel('Transmission Powers (W)','Interpreter','latex', 'FontSize',18);
xlabel('Position (km)','Interpreter','latex', 'FontSize',18);
legend([h1,h2,h3],{'$p_1$','$p_2$','$p_3$'},'Interpreter','latex', 'FontSize',18);

figure('name','Acceleration','position', [200, 200, 700, 300]);
h1 = plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Up,4,solution.T),'b','Linewidth',2);  hold on; grid on;
h2 = plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Up,5,solution.T),'r','Linewidth',2);
h3 = plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Up,6,solution.T),'g','Linewidth',2);
plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Up,4,solution.T),'bo','Linewidth',1);
plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Up,5,solution.T),'ro','Linewidth',1);
plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Up,6,solution.T),'go','Linewidth',1);
ylabel('Acceleration $ms^{-2}$','Interpreter','latex', 'FontSize',18);
xlabel('Position (km)','Interpreter','latex', 'FontSize',18);
legend([h1,h2,h3],{'$a_1$','$a_2$','$a_3$'},'Interpreter','latex', 'FontSize',18);

figure('name','Velocities','position', [200, 200, 700, 300]);
h1 = plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Xp,1,solution.T),'b','Linewidth',2);  hold on; grid on;
h2 = plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Xp,2,solution.T),'r','Linewidth',2);
h3 = plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Xp,3,solution.T),'g','Linewidth',2);
plot(10^-3.*speval(solution.Xp,4,solution.T),speval(solution.Xp,1,solution.T),'bo','Linewidth',1);
plot(10^-3.*speval(solution.Xp,5,solution.T),speval(solution.Xp,2,solution.T),'ro','Linewidth',1);
plot(10^-3.*speval(solution.Xp,6,solution.T),speval(solution.Xp,3,solution.T),'go','Linewidth',1);
ylabel('UAV Velocities $ms^{-1}$','Interpreter','latex', 'FontSize',18);
xlabel('Position (km)','Interpreter','latex', 'FontSize',18);
legend([h1,h2,h3],{'$v_1$','$v_2$','$v_3$'},'Interpreter','latex', 'FontSize',18);



figure('name','Velocities2','position', [200, 200, 700, 300]);
h1 = plot(solution.T(:,1)./60,speval(solution.Xp,1,solution.T),'Linewidth',2);  hold on; grid on;
h2 = plot(solution.T(:,1)./60,speval(solution.Xp,2,solution.T),'Linewidth',2);
h3 = plot(solution.T(:,1)./60,speval(solution.Xp,3,solution.T),'Linewidth',2);
ylabel('UAV Velocities $ms^{-1}$','Interpreter','latex', 'FontSize',17);
xlabel('Time (min)','Interpreter','latex', 'FontSize',17);
legend([h1,h2,h3],{'$v_1$','$v_2$','$v_3$'},'Interpreter','latex', 'FontSize',14);

figure('name','Transmission Powers2','position', [200, 200, 700, 300]);
h1 = plot(solution.T(:,1)./60,speval(solution.Up,1,solution.T),'Linewidth',2);  hold on; grid on;
h2 = plot(solution.T(:,1)./60,speval(solution.Up,2,solution.T),'Linewidth',2);
h3 = plot(solution.T(:,1)./60,speval(solution.Up,3,solution.T),'Linewidth',2);
ylabel('Transmission Powers (W)','Interpreter','latex', 'FontSize',17);
xlabel('Time (min)','Interpreter','latex', 'FontSize',17);
legend([h1,h2,h3],{'$p_1$','$p_2$','$p_3$'},'Interpreter','latex', 'FontSize',14);


