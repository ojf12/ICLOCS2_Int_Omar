% main_h_MeshRefinement - Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% Aircraft go around in the present of wind-shear problem
%
% The problem was originally presented by: 
% R. Bulirsch, F. Montrone, and H. Pesch, "Abort landing in the presence of windshear as a minimax optimal control problem, part 1: Necessary conditions", Journal of Optimization Theory and Applications, 70(1), pp 1-23, 1991.
% R. Bulirsch, F. Montrone, and H. Pesch, "Abort landing in the presence of windshear as a minimax optimal control problem, part 2: Multiple shooting and homotopy", Journal of Optimization Theory and Applications, 70(2), pp 223-254, 1991.
% This example was formulated originally by:
% J. Betts, "Practical Methods for Optimal Control and Estimation Using Nonlinear Programming: Second Edition," Advances in Design and Control, Society for Industrial and Applied Mathematics, 2010.
% This implementation contains modifications
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

clear all;close all;format compact;
global sol;  
sol=[];                             % Initialize solution structure

options= settings_h(20);                  % Get options and solver settings 
[problem,guess]=WindshearGoAround;          % Fetch the problem definition
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9;
i=1; imax=1;

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


%%
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,3,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,3,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Velocity [ft/s]')
grid on


figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,4,[solution.T(:,1); solution.tf])*180/pi,'bo' )
hold on
plot(xx,speval(solution.Xp,4,xx)*180/pi,'b-' )
hold on
plot([solution.T(1,1); solution.tf],[problem.states.xl(4), problem.states.xl(4)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(4), problem.states.xu(4)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,1,xx),'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Position [ft]')
grid on

%%
figure
plot(solution.X(:,1),solution.X(:,2),'bo-')
hold on
h1=plot([solution.X(1,1); solution.X(end,1)],[solution.p, solution.p],'k--' );
xlabel('Position [ft]','Interpreter','latex', 'FontSize',13);
ylabel('Altitude [ft]','Interpreter','latex', 'FontSize',13);
xlim([solution.X(1,1); solution.X(end,1)])
legend(h1,'minimum altitude')
grid on

%%
figure
plot(solution.T,speval(solution.Up,1,solution.T)*180/pi,'bo' )
hold on
plot(xx,speval(solution.Up,1,xx)*180/pi,'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu, problem.inputs.uu]*180/pi,'r--' )
xlim([0 solution.tf])
xlabel('Time [s]','Interpreter','latex', 'FontSize',13);
ylabel('Angle of attack [deg]','Interpreter','latex', 'FontSize',13);
grid on

figure
plot(solution.T(1:end-1,1),solution.dU*180/pi,'-bo')
hold on
plot([solution.T(1,1); solution.tf],[problem.inputs.url, problem.inputs.url]*180/pi,'r--' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uru, problem.inputs.uru]*180/pi,'r--' )
xlim([0 solution.tf])
xlabel('Time [s]','Interpreter','latex', 'FontSize',13);
ylabel('Angle of attack rate [deg/s]','Interpreter','latex', 'FontSize',13);
grid on

