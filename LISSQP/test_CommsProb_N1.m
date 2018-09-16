clear all;
close all;
clc;

%------------------------------------------------------------------------%
%% Setupt Simulation Environment
global P nt  Vars M Int K M T B D Dt;
nt = 100;
Vars = zeros(4,1);

% Define variable parameters and set global parameters
Vt = 10;
DMrat = 0.75;
sigSq_t = 10e-11*ones(nt+1,1);
setGlobalParametersN1(sigSq_t, DMrat, Vt,[],[]);
Dt = D*P/M;
%------------------------------------------------------------------------%

%% Function definitions
funcs.objective         = @(x)objective(x);
funcs.constraints       = @(x)constraints(x);
funcs.gradient          = @(x)gradient(x);
funcs.jacobian          = @(x)jacobian(x);
funcs.hessian           = @(x,v)hessian(x,v);

%% LISSQP Solver Options
options.max_iter=500;
options.vars.eta = 0.1; % between 0 and 0.5
options.vars.tau = 0.75; % between 0 and 1
options.vars.delta = 1e4; % big enough
options.vars.rho = 0.5; %
options.lp_options = optimoptions(@linprog,'disp','off');
options.realHess = 1;
options.linElasticStart = 0;

% QP SOLVER OPTIONS
% Uncomment to use Quadprog
qp_options = optimoptions(@quadprog,'disp','off');
qpSolver = @(H,f,A,b,Aeq,beq,lb,ub,x0,Ac)qps_ip(H,f,A,b,Aeq,beq,lb,ub,x0,qp_options,Ac);

% Uncomment to use qpOASES
% qpSolver = @(H,f,A,b,Aeq,beq,lb,ub,x0,Ac)qps_qpOASES(H,f,A,b,Aeq,beq,lb,ub,x0,Ac);

%% Lower and Upper bounds on decision variables and constraints
bounds.xl = zeros(3*(nt+1),1);
bounds.xu = ones(3*(nt+1),1).*P;
bounds.xu(2:3:end-1) = inf;

bounds.cl = zeros(nt+1,1);
bounds.cu = inf*ones(nt+1,1);

bounds.bl = [Dt;zeros(nt+1,1)];
bounds.bu = [inf;zeros(nt+1,1)];


x0 = ones(3*(nt+1),1);

tmp1 = kron(eye(nt+1), [ 0 0 1]);
tmp2 = [zeros(1,3*(nt+1)) ; ...
    kron(eye(nt), [0 1 -1]) , zeros(nt,3)];
A   = [ tmp1 + tmp2 ; ...
    zeros(1,3*nt),0,-1,1];


%% Call SQP solver
tic;
x = LISSQP(funcs,A,bounds, qpSolver, options, x0);
time=toc;

%% Plot solution
figure;
subplot(3,1,1); plot(x(1:3:end-2)); hold on; plot(x(1:3:end-2),'b*');
subplot(3,1,2); plot(x(2:3:end-1)); hold on; plot(x(2:3:end-1),'b*');
subplot(3,1,3); plot(x(3:3:end-0)); hold on; plot(x(3:3:end-0),'b*');
f_sol = objective(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% QP solver definition

function [x, lambdaQP, flag,Ac] = qps_ip(H,f,A,b,Aeq,beq,lb,ub,x0,opt,Ac)
[x,t1,flag,t2,l_QuadProg] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,opt);
lambdaQP = [ l_QuadProg.eqlin; l_QuadProg.ineqlin; l_QuadProg.lower; l_QuadProg.upper];
Ac = [];
end

function [x, lambdaQP, flag,Ac] = qps_qpOASES(H,f,A,b,Aeq,beq,lb,ub,x0,Ac)

At  = [Aeq; A];
ubA = [beq ; b];
lbA = [beq ; -1*inf*ones(size(b))];

m1 = size(Aeq,1);
m2 = size(A,1);
m3 = length(lb);

WorkingSetGuess = qpOASES_auxInput('guessedWorkingSetB',Ac(1:length(x0)),'guessedWorkingSetC',Ac(length(x0)+1:end));
[x,fval,flag,iter,lambdaQP,auxOut] = qpOASES(H,f,At,lb,ub,lbA,ubA,qpOASES_options,WorkingSetGuess);

% -1 indicates lower bound is active
% 0 indicates that no bound is active
% 1 indicates that the upper bound is active 
W1 = auxOut.workingSetB;
W2 = auxOut.workingSetC;

% Generate masks for lagrange multipliers from working set
W1_l = -1*W1;
W1_l(W1_l<0) = 0;

W1_u= W1;
W1_u(W1_u<0) = 0;

% Get multipliers and apply active set mask for inequalities
lambdaQP = [-1*lambdaQP(m3+1:m3+m1); ...
    -1*lambdaQP(m3+m1+1 : end) .* (W2(m1+1:m1+m2)); ...
    [lambdaQP(1:m3) ; lambdaQP(1:m3)].*[W1_l ; -1*W1_u]];
Ac = [W1;W2];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function definitions
function [obj] = objective(x)
global Int ;
p = x(1:3:end-2);
obj = Int'*p;
end

function [grad]   = gradient(x)
global Int nt;

grad = zeros(1,3*(nt+1));
grad(1:3:end-2) =  Int;
end

function [H] = hessian(x,lambda_C)

global nt K B M Int P;

p = x(1:3:end-2);
coeff = Int*(B*P/(M*log(2)));

% Get appropriate lambdas
dcp = (lambda_C.*(coeff.*(K./(1+K.*p)).^2))';
dcr = zeros(1,nt+1);
dcS = zeros(1,nt+1);

H = diag(reshape( [dcp ; dcr ; dcS],[3*(nt+1),1]));

end

function [c] = constraints(x)
% Sum data rate onstraint
global  Int K B M P;
p = x(1:3:end-2);
r = x(2:3:end-1);

%data rate constraints
c = -(r - Int.*((B*P/M)*log2(1+K.*p)));
end

function [jac]   = jacobian(x)
global B  K nt Int M P;
p = x(1:3:end-2);

dcp = diag((B*P/(M*log(2)))*Int.*(K./(1+K.*p)));
dcr = -1*eye(nt+1);
dcS = zeros(nt+1);

jac = reshape( [dcp , dcr , dcS]',[nt+1,3*(nt+1)]);
jac = sparse(jac);
end