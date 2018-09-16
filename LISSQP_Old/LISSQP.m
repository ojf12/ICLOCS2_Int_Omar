function [ x, output ] = LISSQP( funcs, A, bounds,qpSolver, options, data, x_init)
% LISSQP  Line search sequential quadratic programming.
%   X = LISSQP(funcs,A,bounds,qpSolver, options, x_init) attempts to solve
%   the constrained nonlinear programming problem:
%       min funcs.objective(x)
%   s.t.
%       bounds.bl <= A*x                  <= bounds.bu,
%       bounds.cl <= funcs.constraints(x) <= bounds.cu,
%       bounds.xl <= x                    <= bounds.xu.
%
%   FUNCS
%   The first input is a structure containing function handles for the
%   following MATLAB routines.
%
%   * funcs.objective(required)
%       Calculates the value of the objective function at the current
%       iterate x. Returns a scalar value. 
%
%   * funcs.gradient(required)
%       Calculates the gradient of the ojective function at the current
%       iterate x and returns the result in a row vector of length n =
%       length(x).
%
%   * funcs.Hessian(not required)
%
%   * funcs.constraints(required)
%       Calculates the value of the nonlinear constraints at the current
%       iterate x and stores the result in a column vector of length m,
%       where m is the number of constraints.
%
%   * funcs.jacobian(required)
%       Calculates the value of the jacobian of the nonlinear constraints
%       at the current iterate x and returns a matrix of dimensions m x n.
%
%   A
%   The second input is a matrix A. This matrix corresponds to the jacobian
% of the linear constraints.
%
%   BOUNDS
%   The third input is a structure containing upper and lower bounds for
%   the decision variables x, linear constraints Ax and nonlinear
%   constraints c(x) such that: 
%       bounds.bl <= A*x                  <= bounds.bu,
%       bounds.cl <= funcs.constraints(x) <= bounds.cu,
%       bounds.xl <= x                    <= bounds.bu.
%   If the upper and lower bounds for a particular constraint are
%   equal then this constraint will be interpreted as an equality
%   constraint. If there does not exist a lower or upper bound on a 
%   particular constraint then set the corresponding entry to -+ Inf
%   respectively.
%
%   QPSOLVER
%   This is a function handle for a QP solver that generates a solution to
%   the convex programming problem: (subprob 1)
%       min x'Hx + f'x
%   s.t.
%       Aeq*x = beq     (l1),
%       A*x >= b        (l2),
%       lb<=x<=ub       (l3),
%   where the variables in parentheses correspond to the dual variables of
%   each constraint.
%   This function is called at each iteration of the SQP to solve the
%   arising QP subproblem. To interface with the SQP framework this
%   function handle must be of the following form:
%
%   [x, lambdaQP, flag,Ac] = QPSOLVER(H,f,A,b,Aeq,beq,lb,ub,x0,Ac)
%
%   The first 8 inputs correspond to vectors and matrices in (subprob 1).
%   x0 is an initial guess of the solution. Ac is an estimate of the active
%   set at the solution, in a format appropriate to the solver of choice. 
%   If QPSOLVER is a function handle to an active set QP solver, then 
%   supplying an initial quess of the active set (warm starting) may reduce 
%   solution time of the inner SQP. If the QP solver of choice does not 
%   accept warm starting, then set AC=[].
%
%   QPSOLVER returns the primal and dual solutions x and lambdaQP. To
%   interface with the outer SQP loop the variable lambdaQP must be ordered
%   such that 
%       lambdaQP = [l1 ; l2 ; l3]
%
%   FLAG is an integer value, -2 if the QP subproblem does not converge to
%   a feasible point, 0 otherwise.
%
%   AC is the true value of the active set at the solution of the QP
%   subproblem. 
%
%   Here is an example of a function handle to Matlab's interior point
%   QUADPROG solver.
%
%       function [x, lambdaQP, flag,Ac] = qpsolver(H,f,A,b,Aeq,beq,lb,ub,x0,opt,Ac)
%           [x,~,flag,~,l_QuadProg] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,opt);
%           lambdaQP = [ l_QuadProg.eqlin; l_QuadProg.ineqlin; l_QuadProg.lower; l_QuadProg.upper];
%           Ac = [];
%       end
%
%   OPTIONS
%
%     *options.max_iter=500;
%       The maximum number of major iterations of the SQP algorithm.
%       Options of the QP solver used should be set outside of the SQP call
%       and will be assumed constant over the major iterations.
%
%     *options.vars.eta = 0.1;
%       Eta takes a value between 0<eta<0.5. It is the weighting parameter
%       used in the iterative line search procedure to determine a suitable
%       step size that result in a descent direciton. At iteration k, with
%       point (x_k), the new candidate point x_(k+1) = x_(k)+alpha_(k)p_(k)
%       must satisfy:
%           L1(x_(k+1)) > L1(x_(k)) + eta * D_L1(x_(k),p)
%       where L1(x) is the l1-merit function evaluated at point x.
%       D_L1(x,p) is the directional derivative of the l1-merit function
%       at point x in direction p. Qualitatively, a smaller eta results 
%       in a smaller step length. 
%
%     *options.vars.tau = 0.75; 
%       Tau takes a value between 0<tau<1. This is a scaling parameter used
%       in the iterative line search procedure. Candidate step lengths
%       alpha are scaled down by this value iteratively until a descent
%       direction is found. Smaller values of tau result in fewer
%       iterations of this procedure but also result in sparser sampling of
%       the merit function along direciton p and therefore may generate
%       more conservative step lengths.
%
%     *options.vars.delta = 1e4; 
%       Delta is a weighting parameter used in the l1-merit function.
%       Loosely speaking, delta relates to the dynamics of the penalty
%       weighting for constraint violation. A larger delta results in high
%       er constraint penalties as the algorithm iterates. Delta must be
%       sufficiently large.
%
%     *options.vars.rho = 0.5; 
%       Rho is a weighting parameter used in the l1-merit function, taking
%       values 0<Rho<1. Rho is used for generating a lower bound for the
%       penalty weighting associated with constraint violation. 
%
%     *options.lp_options = optimoptions(@linprog,'disp','off');
%       If options.linElasticStart is selected, then lissqp solves a linear
%       elastic problem to check feasibility of the nonlinear constraints, 
%       using user defined x0 as a starting point.
%       If a feasible point is found the SQP is started with this point,
%       otherwise if the linear constraints are found infeasible the
%       algorithm terminates without further function evaluations. Note
%       that the feasible point generated by this procedure may not satisfy
%       the nonlinear constraints and may in general not be a good starting
%       point. User defined starting points may result in better
%       performance. 
%
%     *options.realHess = 0;
%       Set to 1 if the real hessian is supplied, defualt is 0. If the real
%       Hessian is not supplied then a Damped BFGS update is used to
%       generate an approximation of the hessian of the Lagrangian. 
%
%     *options.linElasticStart = 0;
%       Set to 1 if the initial point x0 is to be generated by a linear
%       feaisbility problem, described above. Set to 0 if the user defined
%       initial guess of x0 is to be used. 
%
%     X0 is the user supplied starting point

%% Initialization and Error checking
% Save dimensions of x
n = length(x_init);
nSave = n;

% Initialise constant solver parameters
eta     = options.vars.eta;
tau     = options.vars.tau;
delta   = options.vars.delta;
rho     = options.vars.rho;

x = x_init;

%% Parse equality and inequality linear constraints.
eq_ind = find(bounds.bl==bounds.bu);
in_ind = find(bounds.bl~=bounds.bu);

Aeq   = A(eq_ind,:);
beq   = bounds.bl(eq_ind);

Aineq = [A(in_ind,:); -1*A(in_ind,:)];
bineq = [bounds.bl(in_ind); -1*bounds.bu(in_ind)];
Aineq = Aineq(~isinf(bineq),:);
bineq = bineq(~isinf(bineq));

% Get dimensions of linear constraint matrices
mAEq   = size(Aeq,1);
mAIneq = size(Aineq,1);
mA     = mAEq + mAIneq;

%Format linear elastic program to check linear constraint feasibility
if(options.linElasticStart)
    [xt,fval] = linElasticProg(Aeq,beq,Aineq,bineq,mAEq,mAIneq,mA,bounds,n,x,options.lp_options);
    if(fval > 0)
        disp('Infeasible linear constraints, abort')
        return;
    else
        disp('All is swell')
        x = xt(1:n);
    end
end


%% Parse nonlinear constraints, initialise matrix buffs for Hessian approx
usefulMat = [ eye(n) ; -1*eye(n) ];

funcsLoc = funcs;
buff = 1;
buff_new = 2;

AeqX=[];
if(mAEq>0)
   Aeq  = sparse(Aeq);
   AeqX = Aeq*x - beq; 
end

AineqX=[];
if(mAIneq>0)
   Aineq  = sparse(Aineq);
   AineqX = Aineq*x-bineq;
end

% Parse nonlinear equality and inequality constraints
constInd.ceq = find(bounds.cl==bounds.cu);
cbEq = bounds.cl(constInd.ceq);

c_in_index = find(bounds.cl~=bounds.cu);
cl = bounds.cl(c_in_index);
cu = bounds.cu(c_in_index);
cbIn = [cl ; -1*cu];
cbIn = cbIn(~isinf(cbIn));
constInd.cl = c_in_index(~isinf(cl));
constInd.cu = c_in_index(~isinf(cu));

% Constraint buffer
tmp1 = funcsLoc.constraints(x,data);
const(:,1) = [AeqX ; (tmp1(constInd.ceq) - cbEq) ; AineqX ; ([tmp1(constInd.cl);-1*tmp1(constInd.cu)] - cbIn)];
m = length(const);
const(:,2) = zeros(m,1);

% Get dimensions of nonlinear constraint matrices
mNLEq = length(constInd.ceq);
mNLIn = length(constInd.cu)+length(constInd.cl);
mNL = m - mA;
mEq = mAEq + mNLEq;

% Jacobian buffer
tmp2 = funcsLoc.jacobian(x,data);
jac = zeros(mNL,n,2);
jac(:,:,1) = [tmp2(constInd.ceq,:);tmp2(constInd.cl,:);-1*tmp2(constInd.cu,:)];

% Objective buffer
obj = zeros(1,2);
obj(buff) = funcsLoc.objective(x,data);

% Gradiant buffer
grad = zeros(n,2);
grad(:,buff) = funcsLoc.gradient(x,data);

% Initialize remaining Variables
lambda = rand(m+2*n,1);
p      = zeros(n,1);
mu     = 1;
numInf = 0;
Ac     = zeros(n+m,1);

% Compute initial positive definite Hessian Approx or exact Hessian
if(options.realHess)
%     B = funcs.hessian(x_init,getLambda_C(lambda,constInd,mNL,mAEq,mEq,mAIneq,mNLIn,mNLEq));
    B = funcs.hessian(x_init,1,getLambda_C(lambda,constInd,mNL,mAEq,mEq,mAIneq,mNLIn,mNLEq),data);
else
    B = eye(length(x_init));
end

%% Start major iteration
disp('Entering major iterations');
iter = 1;
while(iter < options.max_iter)

    % Check if termination test is satisfied, SNOPT
    if((iter > 1) && (exit_flag~=-2) && terminationTest(const(:,buff),lambda,grad(:,buff),allJac,bounds,x,mEq))
        disp('Termination criteria satisfied');
        if( (n>nSave) && (norm(x(nSave+1:end))>1e-5) )
           disp('Nonzero slack variables detected, nonlinear constraints infeasible'); 
        end
        break;
    end
    
    % Solve quadratic subproblem for p - minor iterations
    
    jac=full(jac);
    [p,l_qp,exit_flag,Ac] = ...
        qpSolver(B, grad(:,buff), ... % Hessian and gradient info
        -[Aineq ; jac(mNLEq+1:end,:,buff)] , const(mEq+1:end,buff), ... % Inequality constraints
        -[Aeq   ; jac(1:mNLEq,:,buff)]     , const(1:mEq,buff),... % Equality constraints
        bounds.xl-x,bounds.xu-x,p,Ac); % Variable bounds

    if((exit_flag==-2))
        numInf=numInf+1;
        if(numInf>1)
            disp('Problem found to have infeasible linear constraints, exiting major iterations');
            break;
        end
        % Add slack variables if nonlinear constraints are infeasible
        gamma = 100*norm(grad(:,buff));
        
        disp('Problem found to be inconsistent, introducing slack variables');
        
        disp('->Formatting Elastic problem functions')
        funcsLoc.objective         = @(x)(funcs.objective(x(1:nSave)) + gamma*sum(x(nSave+1:end)));
        funcsLoc.constraints       = @(x)(funcs.constraints(x(1:nSave)) - x(nSave+1:nSave+mNL) + x(nSave+mNL+1:end));
        funcsLoc.gradient          = @(x)([funcs.gradient(x(1:nSave)) , gamma*ones(1,2*mNL)]);
        funcsLoc.jacobian          = @(x)([funcs.jacobian(x(1:nSave)), -1*eye(mNL), eye(mNL)]);
        
        disp('->Formatting Elastic problem matrices and bounds')
        Aeq   = [Aeq ,   zeros(mAEq,2*mNL)];
        Aineq = [Aineq , zeros(mAIneq,2*mNL)];
        usefulMat =  [ eye(n+2*mNL) ; -1*eye(n+2*mNL) ];
        
        if(options.realHess)   
            grad     = [grad ; gamma*ones(2*mNL,1)];
            slackJac = [-1*eye(mNL), eye(mNL)];
            jac      = [jac, slackJac];
            
            funcsLoc.hessian = @(x,lambda_C)([funcs.hessian(x(1:nSave),lambda_C), zeros(n,2*mNL) ; ...
                                    zeros(2*mNL,n) , eye(2*mNL)]);
        else
            grad     = [grad ; gamma*ones(2*mNL,2)];
            slackJac = [-1*eye(mNL), eye(mNL)];
            jac      = cat(2,jac,cat(3,slackJac,slackJac));
        
            B = [B , zeros(n,2*mNL) ; zeros(2*mNL,n) , eye(2*mNL)];
        end
        
        bounds.xl = [bounds.xl; zeros(2*mNL,1)];
        bounds.xu = [bounds.xu; inf*ones(2*mNL,1)];
        
        n = n + 2*mNL;
        p = [p;zeros(2*mNL,1)];
        x = [x;zeros(2*mNL,1)];
        lambda = [lambda ; zeros(4*mNL,1)];
        
        disp('->Solving Elastic problem');
    else
        
        % Get lambda step directions
        p_lambda = l_qp-lambda;
        
        %% Perform Line search based on merit function
        %Choose good mu
        mu = choose_mu(mu,p,B,grad(:,buff),const(:,buff),rho,delta);
        
        alpha = 1;
        tau_alpha = 0.5;
        
        D_phi   = D_merit(x,const(:,buff),grad(:,buff),mEq,bounds,mu,p);
        merit_x = merit(x,funcsLoc,Aeq,beq,Aineq,bineq,bounds,constInd,mEq, cbIn, cbEq, mu,data);
        while(merit((x+alpha*p),funcsLoc,Aeq,beq,Aineq,bineq,bounds,constInd,mEq, cbIn, cbEq, mu,data) > ...
                (merit_x + eta*alpha*D_phi))
            alpha = alpha*tau_alpha;
            tau_alpha = tau_alpha*0.5;
        end
        
        %% Update function evaluation and buffers/Hessian approximation if necessary
        x = x + alpha*p;
        lambda = lambda + alpha*p_lambda;
        
        % Update buffers
        if(mAEq>0)
            AeqX=Aeq*x-beq;
        end
        if(mAIneq>0)
            AineqX=Aineq*x-bineq;
        end
     
        if(options.realHess)
            % No need for circular buffer in this case 
            
%             B = funcs.hessian(x,getLambda_C(lambda,constInd,mNL,mAEq,mEq,mAIneq,mNLIn,mNLEq));
            B = funcs.hessian(x,1,getLambda_C(lambda,constInd,mNL,mAEq,mEq,mAIneq,mNLIn,mNLEq),data);
            
            tmp2 = funcsLoc.jacobian(x,data);
            jac = [tmp2(constInd.ceq,:);tmp2(constInd.cl,:);-1*tmp2(constInd.cu,:)];
            allJac = [Aeq ; jac(1:mNLEq,:) ; Aineq ; jac(mNLEq+1:end,:) ; usefulMat];

            tmp1 = funcsLoc.constraints(x,data);
            const = [AeqX ; (tmp1(constInd.ceq) - cbEq) ; AineqX ; ([tmp1(constInd.cl);-1*tmp1(constInd.cu)] - cbIn)];

            obj   = funcsLoc.objective(x,data);
            grad  = funcsLoc.gradient(x,data)';
        else
            % Put updates in circular buffer
            
            buff_new = mod(buff,2)+1;
            
            tmp2 = funcsLoc.jacobian(x,data);
            jac(:,:,buff_new) = [tmp2(constInd.ceq,:);tmp2(constInd.cl,:);-1*tmp2(constInd.cu,:)];
            
            tmp1 = funcsLoc.constraints(x,data);
            const(:,buff_new) = [AeqX ; (tmp1(constInd.ceq) - cbEq) ; AineqX ; ([tmp1(constInd.cl);-1*tmp1(constInd.cu)] - cbIn)];
            
            obj(buff_new)     = funcsLoc.objective(x,data);
            grad(:,buff_new)  = funcsLoc.gradient(x,data);
            
            allJac = [Aeq ; jac(1:mNLEq,:,buff_new) ; Aineq ; jac(mNLEq+1:end,:,buff_new) ; usefulMat];
            y = (grad(:,buff_new) - (lambda'*allJac)') - ...
                (grad(:,buff) - (lambda'*[Aeq ; jac(1:mNLEq,:,buff) ; Aineq ; jac(mNLEq+1:end,:,buff) ; usefulMat])');
            B = Damped_BFGS(alpha*p,y,B);
            
            buff = buff_new;
        end
        
        
        disp(['Iteration ',num2str(iter),', primal step size norm ',num2str(norm(p)*alpha)]);
        iter = iter+1;
     end

end

output.iter = iter;
output.lambda = lambda;
output.obj = obj(buff); 
x = x(1:nSave);

% Function end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lambda_C = getLambda_C(lambda,constInd,mNL,mAEq,mEq,mAIneq,mNLIn,mNLEq)
% Get lambdas corresponding to equality and inequality constraints. This is
% for use if the user wishes to supply the real hessian of the Lagrangian.

    tmp = length([constInd.cl;setdiff(constInd.cu,constInd.cl)]);
if(mNL>0)
    lambda_C = zeros(mNLEq + tmp,1);

    % Lambda corresponding to equality constraints is directly available
    lambda_C(constInd.ceq) = -1*lambda(mAEq+1 : mEq);
    
    % Parse lambda corresponding to inequality constriants
    lambda_lu = lambda(mEq + mAIneq +1 : mEq + mAIneq + mNLIn);
    lambda_l  = lambda_lu(1:length(constInd.cl));
    lambda_u  = lambda_lu(length(constInd.cl)+1:end);
    
    % Lambda corresponding to inequality constraints   
    lambda_C(constInd.cl) = lambda_C(constInd.cl) + lambda_l;
    lambda_C(constInd.cu) = lambda_C(constInd.cu) - lambda_u;
    else
       lambda_C=[]; 
    end
end

function [xt,fval] = linElasticProg(Aeq,beq,Aineq,bineq,mAEq,mAIneq,mA,bounds,n,x,op)
% Format and solve a linear elastic program to check feasibility of the
% linear constraints.
%           min_(x,v,w) sum(v + w)
%           s.t.
%           xl <= x <= xu
%           bl <= Ax -v + w <= bu
% If constraints are feasible, then a point is selected that satisfies
% these constraints, and therefore each iteration of the SQP algorithm will
% also satisfy the linear constraints.
% If linear constraints are infeasible, SQP terminates without entering
% major iterations.

% Format matrices
f = [zeros(1,n),ones(1,2*mA)];
Ateq   = -1*[Aeq , -1*eye(mAEq), eye(mAEq), zeros(mAEq,2*mAIneq)];
Atineq = -1*[Aineq , zeros(mAIneq,2*mAEq), -1*eye(mAIneq), eye(mAIneq)];

% Call linprog
[xt,fval] = linprog(f,Atineq,-bineq,Ateq,-beq, ...
    [bounds.xl;zeros(2*mA,1)],[bounds.xu;inf*ones(2*mA,1)],...
    [x;zeros(2*mA,1)],op);

end
    
function t = terminationTest(const,lambda,grad,allJac,bounds,x,mEq)
% Perform an SNOPT-like termination test to evaluate satisfaction of the
% first order conditions of optimality
% SNOPT lilke termination test

t=0;

tx = 1e-5*(1+max(abs(x)));
tl = 1e-5*(1+max(abs(lambda)));

allConst = [const ; x-bounds.xl ; -x+bounds.xu];

ni = ~isinf(allConst);

cond1 = all(allConst(ni) >= -tx);
cond2 = all(lambda(mEq+1:end) >= -tl);
cond3 = all(allConst(ni).*lambda(ni) <= tl);
cond4 = all(abs(grad - allJac'*lambda) <= tl);

if(cond1 && cond2 && cond3 && cond4)
    t=1;
end

end

function phi = merit(x,funcs,Aeq,beq,Aineq,bineq,bounds,constInd,mEq, cbIn, cbEq,mu,data)
% Compute the li merit function of the problem at the point x

% If linear constraint matrices are non-empty then compute constraintmp1t value
if(size(Aeq,1)>0)
   AeqX=Aeq*x-beq; 
else
   AeqX=[]; 
end
if(size(Aineq,1)>0)
   AineqX=Aineq*x-bineq;
else
   AineqX=[];
end

% Evaluate and form full constraint matrix
tmp1 = funcs.constraints(x,data);
const = -1*[AeqX ; (tmp1(constInd.ceq) - cbEq) ; AineqX ; ([tmp1(constInd.cl);-1*tmp1(constInd.cu)] - cbIn) ;...
            x-bounds.xl ; bounds.xu-x];
tmp = const(mEq+1:end);

% Calculate merit function
phi = funcs.objective(x,data) + mu.*sum(abs(const(1:mEq))) + ...
    mu.*( sum(tmp(tmp>0)));

end

function d_phi = D_merit(x,const,grad,mEq,bounds,mu,p)
% Compute the directional derivitive of the l1 merit function at point x in
% direction p

tmp = -1*[const(mEq+1:end) ; x-bounds.xl ; bounds.xu-x];

d_phi = grad'*p + mu*sum(abs(const(1:mEq))) + mu.*sum(tmp(tmp>0));


end

function mu = choose_mu(mu,p,B,grad, const,rho,delta)
% Compute mu, to be used in the merit function, to ensure convergence of
% the SQP. Practical lower bound on mu taken from e (18.36), Numerical
% Optimization, Nocedal and Wright

% Can also do this with a max function
tmp = (grad'*p + (1/2)*p'*B*p)/((1-rho)*sum(abs(const)));
if(mu < tmp)
    mu = mu + delta;
end

end

function B_k = Damped_BFGS(s,y,B_k)
% Perform BFGS update of the Hessian approximation if the exact Hessian of
% the Lagrangian is not supplied. Inclusion of damping terms ensures
% positive definiteness of the Hessian. See procedure (18.2), Numerical
% Optimization, Nocedal and Wright

tmp1 = B_k*s;
tmp2 = s'*tmp1;
% Determine interpolation amount
theta = 1;
if( (s'*y) < 0.2*(tmp2))
    theta = (0.8*(tmp2)) / ( (tmp2) - (s'*y) );
end

% Get interpolated value of approximation of gradient of lagrangian
r = theta*y + (1-theta)*tmp1;

% Damped BFGS update
B_k = B_k - (tmp1*tmp1')/(tmp2) + (r*r')/(s'*r);

end