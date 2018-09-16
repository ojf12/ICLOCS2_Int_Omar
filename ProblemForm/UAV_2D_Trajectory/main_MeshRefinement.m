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

simNum=4;
if(simNum==1)
    [problem,guess]=UAVSim1;  		% Fetch the problem definition
elseif(simNum==2)
    [problem,guess]=UAVSim2;  		% Fetch the problem definition
elseif(simNum==3)
    [problem,guess]=UAVSim3; 
elseif(simNum==4)
    [problem,guess]=UAVSim4; 
end

%%%% Choice among the ones belwo %%%%
% options= settings_auto(40);                  
options= settings_h(100);                
% options= settings_hp(5,4);                 

% Declare variables for storage
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9; % the initial (before start) maximum absolute local eeror
i=1; % starting iteration
imax=20; % maximum number of mesh refinement iterations

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


if(simNum==1 || simNum==4)
    figure('name','Transmission Powers','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Up,1,solution.T),'-o','Linewidth',2);  hold on; grid on;
    
    figure('name','Transmission Rates','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Up,2,solution.T),'-o','Linewidth',2);  hold on; grid on;
    
    figure('name','Acceleration','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Up,3,solution.T),'-o','Linewidth',2);  hold on; grid on;
    
    figure('name','Storage Buffers','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Xp,1,solution.T),'-o','Linewidth',2);  hold on; grid on;
    
    figure('name','UAV position','position', [200, 200, 700, 300]);
    plot(speval(solution.Xp,2,solution.T),speval(solution.Xp,3,solution.T),'-o','Linewidth',2); hold on; grid on;
    
elseif(simNum==2 || simNum==3)
    figure('name','Transmission Powers','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Up,1,solution.T),'-o','Linewidth',2);  hold on; grid on;
    plot(solution.T(:,1),speval(solution.Up,3,solution.T),'-o','Linewidth',2);  hold on; grid on;

    figure('name','Transmission Rates','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Up,2,solution.T),'-o','Linewidth',2);  hold on; grid on;
    plot(solution.T(:,1),speval(solution.Up,4,solution.T),'-o','Linewidth',2);  hold on; grid on;

    figure('name','Acceleration','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Up,5,solution.T),'-o','Linewidth',2);  hold on; grid on;
    
    
    figure('name','Storage Buffers','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Xp,1,solution.T),'-o','Linewidth',2);  hold on; grid on;
    plot(solution.T(:,1),speval(solution.Xp,2,solution.T),'-o','Linewidth',2);  hold on; grid on;

    figure('name','UAV position','position', [200, 200, 700, 300]);
    plot(speval(solution.Xp,3,solution.T),speval(solution.Xp,4,solution.T),'-o','Linewidth',2); hold on; grid on;
end
if(simNum==4)
    figure('name','UAV Velocity','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Xp,6,solution.T),'-o','Linewidth',2);  hold on; grid on;
end
