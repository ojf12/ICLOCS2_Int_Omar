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

simNum=5;
if(simNum==1)
    [problem,guess]=UAVSim1;  		% Fetch the problem definition
elseif(simNum==2)
    [problem,guess]=UAVSim2;  		% Fetch the problem definition
elseif(simNum==3)
    [problem,guess]=UAVSim3; 
elseif(simNum==4)
    [problem,guess]=UAVSim4; 
elseif(simNum==5)
    [problem,guess]=UAVSim5; 
end
options= settings_h(100);                  % Get options and solver settings
% options= settings_hp(5,4);                  % Get options and solver settings
%options= settings_hp([4 5 3],[-1 0.3 0.4 1]);                  % Get options and solver settings

[infoNLP,data,options]=transcribeOCP(problem,guess,options); % Format for NLP solver
[solution,status,data] = solveNLP(infoNLP,data);      % Solve the NLP
[solution]=output(problem,solution,options,data,4);          % Output solutions

close all;
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
    th = 0:pi/50:2*pi;
    r = 500;
    x1 = r * cos(th) + 1000;     x2 = r * cos(th) - 1000;
    yunit = r * sin(th);
    plot(x1, yunit,'LineWidth',2);
    plot(x2, yunit,'LineWidth',2);
    
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
elseif(simNum==5)
    figure('name','UAV position','position', [200, 200, 700, 300]);
    plot(speval(solution.Xp,1,solution.T),speval(solution.Xp,2,solution.T),'-o','Linewidth',2); hold on; grid on;
    plot([0,-1000],[-1000,1000],'rx','Linewidth',3);
    figure('name','Lambda','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Up,3,solution.T),'-o','Linewidth',2);  hold on; grid on;
    
    figure('name','Integrator lambda','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Xp,5,solution.T),'-o','Linewidth',2);  hold on; grid on;
end
if(simNum==4)
    figure('name','UAV Velocity','position', [200, 200, 700, 300]);
    plot(solution.T(:,1),speval(solution.Xp,6,solution.T),'-o','Linewidth',2);  hold on; grid on;
end

function circle(x,y,r)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp);
end
