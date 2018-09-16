% main_h_MeshRefinement - Main script to solve the Optimal Control Problem with h-typed mesh and refinement
%
% F50MinTimeFlight - Minimum Time Flight Profile for Commercial Aircraft
%
% The aerodynamic and propulsion data are obtained from
% "Performance model Fokker 50", Delft University of Technology, 2010
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

options= settings_h(40);                  % Get options and solver settings 
[problem,guess]=F50MinTimeFlight;          % Fetch the problem definition
errorHistory=zeros(2,length(problem.states.x0));
npsegmentHistory=zeros(2,1);
ConstraintErrorHistory=zeros(2,length(problem.constraintErrorTol));
timeHistory=zeros(1,2);
solutionHistory=cell(1,2);

maxAbsError=1e9;
i=1; imax=50;

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


%%
MeshRefinementHistory.errorHistory=errorHistory;
MeshRefinementHistory.timeHistory=timeHistory;
MeshRefinementHistory.ConstraintErrorHistory=ConstraintErrorHistory;


%%
xx=linspace(solution.T(1,1),solution.tf,1000);

figure
center=[data.data.auxdata.obs_epos data.data.auxdata.obs_npos];
obspos = [center-data.data.auxdata.obs_r 2*data.data.auxdata.obs_r 2*data.data.auxdata.obs_r];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
hold on
plot(solution.X(:,3),solution.X(:,2),'b-','LineWidth',2)
xlabel('East Position [m]')
ylabel('North Position [m]')
grid on
plot(problem.states.x0(3)-1000, problem.states.x0(2)-1000, '.k', 'MarkerSize',20)
plot(problem.states.xfl(3)-1000, problem.states.xfl(2)-1000, '.k', 'MarkerSize',20)
text(problem.states.x0(3)+20000,problem.states.x0(2),'ORG')
text(problem.states.xfl(3)+10000,problem.states.xfl(2)+10000,'DES')
xlim([-1 10]*10^5)
ylim([-1 9]*10^5)

figure
center=[data.data.auxdata.obs_epos data.data.auxdata.obs_npos];
obspos = [center-data.data.auxdata.obs_r 2*data.data.auxdata.obs_r 2*data.data.auxdata.obs_r];
rectangle('Position',obspos,'Curvature',[1 1], 'FaceColor', 'red', 'Edgecolor','none');
hold on
plot(solution.X(:,3),solution.X(:,2),'b-','LineWidth',2)
xlabel('East Position [m]')
ylabel('North Position [m]')
grid on
plot(problem.states.x0(3)-1000, problem.states.x0(2)-1000, '.k', 'MarkerSize',20)
plot(problem.states.xfl(3)-1000, problem.states.xfl(2)-1000, '.k', 'MarkerSize',20)
text(problem.states.x0(3)+10000,problem.states.x0(2),'ORG')
text(problem.states.xfl(3)-10000,problem.states.xfl(2)-5000,'DES')
text(8.25e05, 7.45e05,'NO FLIGHT ZONE','Color','white','FontSize',14)
xlim([8.1 9.1]*10^5)
ylim([7.1 8.1]*10^5)
%%
figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,1,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,1,xx),'b-' )
% plot([solution.T(1,1); solution.tf],[problem.states.xl(1), problem.states.xl(1)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(1), problem.states.xu(1)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Altitude [m]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,4,[solution.T(:,1); solution.tf]),'bo' )
hold on
plot(xx,speval(solution.Xp,4,xx),'b-' )
plot([solution.T(1,1); solution.tf],[problem.states.xl(4), problem.states.xl(4)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.states.xu(4), problem.states.xu(4)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('True Airspeed [m/s]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,5,[solution.T(:,1); solution.tf])*180/pi,'bo' )
hold on
plot(xx,speval(solution.Xp,5,xx)*180/pi,'b-' )
hold on
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Flight Path Angle [deg]')
grid on


figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,6,[solution.T(:,1); solution.tf])*180/pi,'bo' )
hold on
plot(xx,speval(solution.Xp,6,xx)*180/pi,'b-' )
hold on
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Tracking Angle [deg]')
grid on

figure
plot([solution.T(:,1); solution.tf],speval(solution.Xp,7,[solution.T(:,1); solution.tf])/9.81,'bo' )
hold on
plot(xx,speval(solution.Xp,7,xx)/9.81,'b-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Aircraft Mass [kg]')
grid on

figure
plot(solution.T(:,1),speval(solution.Up,1,solution.T)*180/pi,'bo')
hold on
plot(xx,speval(solution.Up,1,xx)*180/pi,'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(1), problem.inputs.ul(1)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(1), problem.inputs.uu(1)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Angle of attack [deg]')
grid on


figure
plot(solution.T(:,1),speval(solution.Up,2,solution.T)*180/pi,'bo')
hold on
plot(xx,speval(solution.Up,2,xx)*180/pi,'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(2), problem.inputs.ul(2)]*180/pi,'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(2), problem.inputs.uu(2)]*180/pi,'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Roll angle [deg]')
grid on

figure
plot(solution.T(:,1),speval(solution.Up,3,solution.T),'bo')
hold on
plot(xx,speval(solution.Up,3,xx),'b-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.ul(3), problem.inputs.ul(3)],'r-' )
plot([solution.T(1,1); solution.tf],[problem.inputs.uu(3), problem.inputs.uu(3)],'r-' )
xlim([0 solution.tf])
xlabel('Time [s]')
ylabel('Throttle Setting [-]')
grid on