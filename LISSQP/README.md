# LISSQP
# help lissqp

  LISSQP  Line search sequential quadratic programming.
    X = LISSQP(funcs,A,bounds,qpSolver, options, x_init) attempts to solve
    the constrained nonlinear programming problem:
        
        min funcs.objective(x)
    
    s.t.
    
      bounds.bl <= A*x                  <= bounds.bu,
      
      bounds.cl <= funcs.constraints(x) <= bounds.cu,
      
      bounds.xl <= x                    <= bounds.bu.
 
    FUNCS
    The first input is a structure containing function handles for the
    following MATLAB routines.
 
    * funcs.objective(required)
        Calculates the value of the objective function at the current
        iterate x. Returns a scalar value. 
 
    * funcs.gradient(required)
        Calculates the gradient of the ojective function at the current
        iterate x and returns the result in a row vector of length n =
        length(x).
 
    * funcs.Hessian(not required)
 
    * funcs.constraints(required)
        Calculates the value of the nonlinear constraints at the current
        iterate x and stores the result in a column vector of length m,
        where m is the number of constraints.
 
    * funcs.jacobian(required)
        Calculates the value of the jacobian of the nonlinear constraints
        at the current iterate x and returns a matrix of dimensions m x n.
 
    A
    The second input is a matrix A. This matrix corresponds to the jacobian
  of the linear constraints.
 
    BOUNDS
    The third input is a structure containing upper and lower bounds for
    the decision variables x, linear constraints Ax and nonlinear
    constraints c(x) such that: 
        
        bounds.bl <= A*x                  <= bounds.bu,
        
        bounds.cl <= funcs.constraints(x) <= bounds.cu,
        
        bounds.xl <= x                    <= bounds.bu.
        
    If the upper and lower bounds for a particular constraint are
    equal then this constraint will be interpreted as an equality
    constraint. If there does not exist a lower or upper bound on a 
    particular constraint then set the corresponding entry to -+ Inf
    respectively.
 
    QPSOLVER
    This is a function handle for a QP solver that generates a solution to
    the convex programming problem: (subprob 1)
        
        min x'Hx + f'x
    
    s.t.
    
      Aeq*x = beq     (l1),
    
      A*x >= b        (l2),
    
      lb<=x<=ub       (l3),
      
    where the variables in parentheses correspond to the dual variables of
    each constraint.
    This function is called at each iteration of the SQP to solve the
    arising QP subproblem. To interface with the SQP framework this
    function handle must be of the following form:
 
    [x, lambdaQP, flag,Ac] = QPSOLVER(H,f,A,b,Aeq,beq,lb,ub,x0,Ac)
 
    The first 8 inputs correspond to vectors and matrices in (subprob 1).
    x0 is an initial guess of the solution. Ac is an estimate of the active
    set at the solution, in a format appropriate to the solver of choice. 
    If QPSOLVER is a function handle to an active set QP solver, then 
    supplying an initial quess of the active set (warm starting) may reduce 
    solution time of the inner SQP. If the QP solver of choice does not 
    accept warm starting, then set AC=[].
 
    QPSOLVER returns the primal and dual solutions x and lambdaQP. To
    interface with the outer SQP loop the variable lambdaQP must be ordered
    such that 
        lambdaQP = [l1 ; l2 ; l3]
 
    FLAG is an integer value, -2 if the QP subproblem does not converge to
    a feasible point, 0 otherwise.
 
    AC is the true value of the active set at the solution of the QP
    subproblem. 
 
    Here is an example of a function handle to Matlab's interior point
    QUADPROG solver.
 
        function [x, lambdaQP, flag,Ac] = qpsolver(H,f,A,b,Aeq,beq,lb,ub,x0,opt,Ac)
            [x,~,flag,~,l_QuadProg] = quadprog(H,f,A,b,Aeq,beq,lb,ub,x0,opt);
            lambdaQP = [ l_QuadProg.eqlin; l_QuadProg.ineqlin; l_QuadProg.lower; l_QuadProg.upper];
            Ac = [];
        end
 
    OPTIONS
 
      *options.max_iter=500;
        The maximum number of major iterations of the SQP algorithm.
        Options of the QP solver used should be set outside of the SQP call
        and will be assumed constant over the major iterations.
 
      *options.vars.eta = 0.1;
        Eta takes a value between 0<eta<0.5. It is the weighting parameter
        used in the iterative line search procedure to determine a suitable
        step size that result in a descent direciton. At iteration k, with
        point (x_k), the new candidate point x_(k+1) = x_(k)+alpha_(k)p_(k)
        must satisfy:
            L1(x_(k+1)) > L1(x_(k)) + eta * D_L1(x_(k),p)
        where L1(x) is the l1-merit function evaluated at point x.
        D_L1(x,p) is the directional derivative of the l1-merit function
        at point x in direction p. Qualitatively, a smaller eta results 
        in a smaller step length. 
 
      *options.vars.tau = 0.75; 
        Tau takes a value between 0<tau<1. This is a scaling parameter used
        in the iterative line search procedure. Candidate step lengths
        alpha are scaled down by this value iteratively until a descent
        direction is found. Smaller values of tau result in fewer
        iterations of this procedure but also result in sparser sampling of
        the merit function along direciton p and therefore may generate
        more conservative step lengths.
 
      *options.vars.delta = 1e4; 
        Delta is a weighting parameter used in the l1-merit function.
        Loosely speaking, delta relates to the dynamics of the penalty
        weighting for constraint violation. A larger delta results in high
        er constraint penalties as the algorithm iterates. Delta must be
        sufficiently large.
 
      *options.vars.rho = 0.5; 
        Rho is a weighting parameter used in the l1-merit function, taking
        values 0<Rho<1. Rho is used for generating a lower bound for the
        penalty weighting associated with constraint violation. 
 
      *options.lp_options = optimoptions(@linprog,'disp','off');
        If options.linElasticStart is selected, then lissqp solves a linear
        elastic problem to check feasibility of the nonlinear constraints, 
        using user defined x0 as a starting point.
        If a feasible point is found the SQP is started with this point,
        otherwise if the linear constraints are found infeasible the
        algorithm terminates without further function evaluations. Note
        that the feasible point generated by this procedure may not satisfy
        the nonlinear constraints and may in general not be a good starting
        point. User defined starting points may result in better
        performance. 
 
      *options.realHess = 0;
        Set to 1 if the real hessian is supplied, defualt is 0. If the real
        Hessian is not supplied then a Damped BFGS update is used to
        generate an approximation of the hessian of the Lagrangian. 
 
      *options.linElasticStart = 0;
        Set to 1 if the initial point x0 is to be generated by a linear
        feaisbility problem, described above. Set to 0 if the user defined
        initial guess of x0 is to be used. 
 
      X0 is the user supplied starting point
