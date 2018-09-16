function [problem,guess] = myProblem
%myProblem - Template file for optimal control problem definition
%
%Syntax:  [problem,guess] = myProblem
%
% Outputs:
%    problem - Structure with information on the optimal control problem
%    guess   - Guess for state, control and multipliers.
%
% Other m-files required: none
% MAT-files required: none
%
% Copyright (C) 2018 Yuanbo Nie, Omar Faqir, and Eric Kerrigan. All Rights Reserved.
% The contribution of Paola Falugi, Eric Kerrigan and Eugene van Wyk for the work on ICLOCS Version 1 (2010) is kindly acknowledged.
% This code is published under the BSD License.
% Department of Aeronautics and Department of Electrical and Electronic Engineering,
% Imperial College London London  England, UK 
% ICLOCS (Imperial College London Optimal Control) Version 2.0 
% 1 May 2018
% iclocs@imperial.ac.uk

%------------- BEGIN CODE --------------
%
%%  Data Parameter Definition
% Conversion from bits to MB
MB_b = 8e6;
% Bandwidth allocation B_r=B_a=B Hz
B = 10^5;
% Antenna and channel gain parameter
K=1;
% Max transmit power of all nodes
P_max = 100;
% noise power at receivers
ss=10^(-10);
% UAV altitude
alt = 2000;
% Path loss exponent
alpha = 1.5;
% Parameters for UAV drag model
c_1 = 9.26*10^-4;
c_2 = 2250;
% Starting data load, scaled to 100
D_init = 5*MB_b;
D_scaling = 100/D_init;
D_init = 100;
r_Des = 10000;
% Relay Memory Constraint
M = 1.5*D_init;
% UAV mass (kg)
mass = 1.5;
% UAV minimum speed v_min=30km/h maximum speed v_max=100km/h
v_min = 12; 
v_max = 27.7778; 
v_avg = (v_max+v_min)/2;
% Minimum and maximum UAV thrust constraints
F_min = -50;
F_max = 50;
% UAV starting point, and source positions
numUAV=3;
Dist = v_avg*(20*60);

Q_init = zeros(numUAV,2);
Q_init(1,1) = -Dist;
Q_init(1,2) = 0;
Q_init(2,:) = Q_init(1,:) + Dist;
Q_init(3,:) = Q_init(1,:) + 2*Dist;
Disp = 2000;
%% Problem Formulation

%Initial Time. t0<tf
problem.time.t0=0;

% Final time. Let tf_min=tf_max if tf is fixed.
fixed_final_time = 20*60;
problem.time.tf_min=fixed_final_time;     
problem.time.tf_max=fixed_final_time; 
guess.tf=fixed_final_time;

% Parameters bounds. pl=< p <=pu
problem.parameters.pl=[];
problem.parameters.pu=[];
guess.parameters=[];

% Initial conditions for system.
problem.states.x0=[v_avg,v_avg,v_avg, Q_init(1,1), Q_init(2,1),Q_init(3,1)];

% Initial conditions for system. Bounds if x0 is free s.t. x0l=< x0 <=x0u
problem.states.x0l=[v_avg,v_avg,v_avg, Q_init(1,1), Q_init(2,1), Q_init(3,1)];
problem.states.x0u=[v_avg,v_avg,v_avg, Q_init(1,1), Q_init(2,1), Q_init(3,1)];

% State bounds. xl=< x <=xu
problem.states.xl=[v_min,v_min,v_min, Q_init(1,1), Q_init(2,1), Q_init(3,1) ];
problem.states.xu=[v_max,v_max,v_max, Q_init(1,2), Q_init(2,2), Q_init(3,2)];

% State rate bounds. xrl=< x <=xru
problem.states.xrl=[-inf -inf -inf -inf -inf -inf]; 
problem.states.xru=[ inf  inf inf inf inf inf]; 

% State error bounds
problem.states.xErrorTol=[1 1 1 10 10 10]; 

% State constraint error bounds
problem.states.xConstraintTol = [1 1 1 5 5 5];
problem.states.xrConstraintTol= [1 1 1 5 5 5];

% Terminal state bounds. xfl=< xf <=xfu
problem.states.xfl=[v_min,v_min,v_min, Q_init(1,2), Q_init(2,2), Q_init(3,2)]; 
problem.states.xfu=[v_max,v_max,v_max, Q_init(1,2), Q_init(2,2), Q_init(3,2)];

% Guess the state trajectories with [x0 ... xf]
% guess.time=[];
guess.states=[];

% Number of control actions N 
% Set problem.inputs.N=0 if N is equal to the number of integration steps.  
% Note that the number of integration steps defined in settings.m has to be divisible 
% by the  number of control actions N whenever it is not zero.
problem.inputs.N=0;       
      
% Input bounds
problem.inputs.ul=[0, 0,0,             -inf,   -inf, -inf];
problem.inputs.uu=[P_max, P_max, P_max, inf, inf, inf];

% Bounds on the first control action
problem.inputs.u0l=[0, 0, 0,            -inf, -inf, -inf];
problem.inputs.u0u=[P_max, P_max, P_max, inf, inf, inf];

% Input rate bounds
problem.inputs.url=[-inf -inf -inf -5 -5 -5];
problem.inputs.uru=[ inf inf inf 5 5 5];

% Input constraint error bounds
problem.inputs.uConstraintTol  = [1 1 1 1 1 1];
problem.inputs.urConstraintTol = [1 1 1 1 1 1];

% Guess the input sequences with [u0 ... uf]
guess.inputs = [];

% Bounds for path constraint function gl =< g(x,u,p,t) =< gu
problem.constraints.gl=[-inf,  -inf,   -inf, F_min, F_min, F_min];
problem.constraints.gu=[0,0,0,F_max, F_max, F_max];

% Path constraint error bounds
problem.constraints.gTol=[ 1 1 1 6 6 6];

% Bounds for boundary constraints bl =< b(x0,u0,p,t0,tf) =< bu
problem.constraints.bl=[];
problem.constraints.bu=[];

problem.setpoints.states = [];
problem.setpoints.inputs = [];

% store the necessary problem parameters used in the functions
problem.data.auxdata=[B, K, alpha, alt, ss, P_max, M, mass, D_scaling, c_1, c_2, Q_init(end), Disp];

% Plant model name, for use with Adigator
problem.data.plantmodel = '...';

% For algebraic variable rate constraint
problem.data.xrl=problem.states.xrl;
problem.data.xru=problem.states.xru;
problem.data.xrConstraintTol=problem.states.xrConstraintTol;
problem.data.url=problem.inputs.url;
problem.data.uru=problem.inputs.uru;
problem.data.urConstraintTol=problem.inputs.urConstraintTol;

% Get function handles and return to Main.m
problem.functions={@L,@E,@f,@g,@avrc,@b};
problem.functions_unscaled={@L_unscaled,@E_unscaled,@f_unscaled,@g_unscaled,@avrc_unscaled,@b_unscaled};
problem.constraintErrorTol=[problem.constraints.gTol,problem.constraints.gTol,problem.states.xConstraintTol,problem.states.xConstraintTol,problem.inputs.uConstraintTol,problem.inputs.uConstraintTol];

%------------- END OF CODE --------------

function stageCost=L_unscaled(x,u,p,t,data)


% L - Returns the stage cost.
% The function must be vectorized and
% xi, ui are column vectors taken as x(:,i) and u(:,i) (i denotes the i-th
% variable)
% 
% Syntax:  stageCost = L(x,xr,u,ur,p,t,data)
%
% Inputs:
%    x  - state vector
%    xr - state reference
%    u  - input
%    ur - input reference
%    p  - parameter
%    t  - time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    stageCost - Scalar or vectorized stage cost
%
%  Remark: If the stagecost does not depend on variables it is necessary to multiply
%          the assigned value by t in order to have right vector dimesion when called for the optimization. 
%          Example: stageCost = 0*t;
%
%------------- BEGIN CODE --------------
% Get relevant problem data
c_1 = data.auxdata(10);
c_2 = data.auxdata(11);

% Define relevant input variables
p1 = u(:,1);
p2 = u(:,2);
p3 = u(:,3);

% Define relevant state variables
v1  = x(:,1);
v2  = x(:,2);
v3  = x(:,3);

stageCost = (p1 + p2 + p3) + Drag(v1,c_1,c_2).*v1 + Drag(v2,c_1,c_2).*v2 + + Drag(v3,c_1,c_2).*v3;

%------------- END OF CODE --------------


function boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,data) 

% E - Returns the boundary value cost
%
% Syntax:  boundaryCost=E(x0,u0,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    boundaryCost - Scalar boundary cost
%
%------------- BEGIN CODE --------------

% Get relevant problem data
mass = data.auxdata(8);

% Define relevant state variables
v10 = x0(1);  v1f = xf(1);
v20 = x0(2);  v2f = xf(2);
v30 = x0(3);  v3f = xf(3);

% boundary cost is change in kinetic energy
boundaryCost = (mass/2)*(v1f^2 - v10^2) + ...
               (mass/2)*(v2f^2 - v20^2) + ...
               (mass/2)*(v3f^2 - v30^2);

%------------- END OF CODE --------------


function dx = f_unscaled(x,u,p,t,data)
% f - Returns the ODE right hand side where x'= f(x,u,p,t)
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% state corresponds to one column of dx.
% 
% 
% Syntax:  dx = f(x,u,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%    data-structured variable containing the values of additional data used inside
%          the function 
%
% Output:
%    dx - time derivative of x
%
%  Remark: If the i-th ODE right hand side does not depend on variables it is necessary to multiply
%          the assigned value by a vector of ones with the same length  of t  in order 
%          to have  a vector with the right dimesion  when called for the optimization. 
%          Example: dx(:,i)= 0*ones(size(t,1)); 
%
%------------- BEGIN CODE --------------

% Define relevant input variables
a1 = u(:,4);
a2 = u(:,5);
a3 = u(:,6);

% Define relevant state variables
v1  = x(:,1);
v2  = x(:,2);
v3  = x(:,3);

q1  = x(:,4);
q2  = x(:,5);
q3  = x(:,6);

%Define ODE right-hand side
dx(:,1) = a1;
dx(:,2) = a2;
dx(:,3) = a3;

dx(:,4) = v1;
dx(:,5) = v2;
dx(:,6) = v3;

%------------- END OF CODE --------------


function c=g_unscaled(x,u,p,t,data)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% constraint corresponds to one column of c
% 
% Syntax:  c=g(x,u,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%   data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    c - constraint function
%
%
%------------- BEGIN CODE --------------

% Get relevant problem data
B       = data.auxdata(1);
K       = data.auxdata(2);
alpha   = data.auxdata(3);
alt     = data.auxdata(4);
ss      = data.auxdata(5);
mass    = data.auxdata(8);
scaling = data.auxdata(9);
c_1     = data.auxdata(10);
c_2     = data.auxdata(11);
AP      = data.auxdata(12);
Disp    = data.auxdata(13);
% Define relevant input variables
p1 = u(:,1);
p2 = u(:,2);
p3 = u(:,3);

a1 = u(:,4);
a2 = u(:,5);
a3 = u(:,6);

% Define relevant state variables
v1  = x(:,1);
v2  = x(:,2);
v3  = x(:,3);

q1  = x(:,4);
q2  = x(:,5);
q3  = x(:,6);


% Second Transmission Constraint
c(:,1)= 10000-(B)*log2(1 + (K.*p1)./(ss.*((q1-q2).^2 + Disp^2).^alpha));

% Third Transmission Constraint
c(:,2)= 10000-(B)*log2(1 + (K.*p2)./(ss.*((q2-q3).^2 + Disp^2).^alpha));

% Fourth Transmission Constraint
c(:,3)= 10000-(B)*log2(1 + (K.*p3)./(ss.*((q3-AP).^2 + Disp^2).^alpha));

% First Newtonian force constraint
c(:,4)=  mass*a1 + Drag(v1,c_1,c_2);

% Second Newtonian force constraint
c(:,5)=  mass*a2 + Drag(v2,c_1,c_2);

% Third Newtonian force constraint
c(:,6)=  mass*a3 + Drag(v3,c_1,c_2);

%------------- END OF CODE --------------

function cr=avrc_unscaled(x,u,p,t,data)
% avrc_unscaled - Returns the rate constraint algebraic function where [xrl url] =<
% avrc(x,u,p,t) =< [xru uru]
% The function must be vectorized and
% xi, ui, pi are column vectors taken as x(:,i), u(:,i) and p(:,i). Each
% constraint corresponds to one column of c
% 
% Syntax:  cr=avrc_unscaled(x,u,p,t,data)
%
% Inputs:
%    x  - state vector
%    u  - input
%    p  - parameter
%    t  - time
%   data- structured variable containing the values of additional data used inside
%          the function
%
% Output:
%    cr - constraint function
%
%
%------------- BEGIN CODE --------------

% when having rate constraints (i.e. not all equal to +/-inf)
[ cr ] = addRateConstraint( x,u,p,t,data );

% when not having rate constraints
% [ cr ] = [];

%------------- END OF CODE --------------

function bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin)

% b_unscaled - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
%
% Syntax:  bc=b(x0,xf,u0,uf,p,tf,data)
%
% Inputs:
%    x0  - state at t=0
%    xf  - state at t=tf
%    u0  - input at t=0
%    uf  - input at t=tf
%    p   - parameter
%    tf  - final time
%    data- structured variable containing the values of additional data used inside
%          the function
%
%          
% Output:
%    bc - column vector containing the evaluation of the boundary function 
%
% Leave it here
varargin=varargin{1};
%------------- BEGIN CODE --------------
bc = [];
%------------- END OF CODE --------------
% When adpative time interval add constraint on time
if length(varargin)==2
    options=varargin{1};
    t_segment=varargin{2};
    if strcmp(options.transcription,'hpLGR') && options.adaptseg==1 
        if size(t_segment,1)>size(t_segment,2)
            bc=[bc;diff(t_segment)];
        else
            bc=[bc,diff(t_segment)];
        end
    end
end

%------------- END OF CODE --------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Leave the following unchanged! %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stageCost=L(x,xr,u,ur,p,t,vdat)

% L - Returns the stage cost.
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    if ~isempty(xr)
        xr=scale_variables_back( xr, vdat.Xscale, vdat.Xshift );
    end
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if ~isempty(ur)
        ur=scale_variables_back( ur, vdat.Uscale, vdat.Ushift );
    end
    stageCost=L_unscaled(x,u,p,t,vdat);
else
    stageCost=L_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------


function boundaryCost=E(x0,xf,u0,uf,p,t0,tf,vdat) 
% E - Returns the boundary value cost
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x0=scale_variables_back( x0', vdat.Xscale, vdat.Xshift )';
    xf=scale_variables_back( xf', vdat.Xscale, vdat.Xshift )';
    u0=scale_variables_back( u0', vdat.Uscale, vdat.Ushift )';
    uf=scale_variables_back( uf', vdat.Uscale, vdat.Ushift )';
    boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat);
else
    boundaryCost=E_unscaled(x0,xf,u0,uf,p,t0,tf,vdat);
end


%------------- END OF CODE --------------


function dx = f(x,u,p,t,vdat)
% f - Returns the ODE right hand side where x'= f(x,u,p,t)
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    dx = f_unscaled(x,u,p,t,vdat);
    dx= scale_variables( dx, vdat.Xscale, 0 );
else
    dx = f_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------


function c=g(x,u,p,t,vdat)

% g - Returns the path constraint function where gl =< g(x,u,p,t) =< gu
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    c = g_unscaled(x,u,p,t,vdat);
else
    c = g_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------


function cr=avrc(x,u,p,t,vdat)
% avrc - Returns the rate constraint algebraic function where [xrl url] =< avrc(x,u,p,t) =< [xru uru]
% Warp function
%------------- BEGIN CODE --------------

if isfield(vdat,'Xscale')
    x=scale_variables_back( x, vdat.Xscale, vdat.Xshift );
    u=scale_variables_back( u, vdat.Uscale, vdat.Ushift );
    if isfield(vdat,'Pscale')
        p=scale_variables_back( p, vdat.Pscale, vdat.Pshift );
    end
    cr = avrc_unscaled(x,u,p,t,vdat);
else
    cr = avrc_unscaled(x,u,p,t,vdat);
end

%------------- END OF CODE --------------

function bc=b(x0,xf,u0,uf,p,t0,tf,vdat,varargin)
% b - Returns a column vector containing the evaluation of the boundary constraints: bl =< bf(x0,xf,u0,uf,p,t0,tf) =< bu
% Warp function
%------------- BEGIN CODE --------------
bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin);
if isfield(vdat,'Xscale')
    if ~isempty(bc)
        x0=scale_variables_back( x0', vdat.Xscale, vdat.Xshift );
        xf=scale_variables_back( xf', vdat.Xscale, vdat.Xshift );
        u0=scale_variables_back( u0', vdat.Uscale, vdat.Ushift );
        uf=scale_variables_back( uf', vdat.Uscale, vdat.Ushift );
        if isfield(vdat,'Pscale')
            p=scale_variables_back( p', vdat.Pscale, vdat.Pshift );
        end
        bc=b_unscaled(x0,xf,u0,uf,p,t0,tf,vdat,varargin);
    end
end
%------------- END OF CODE --------------