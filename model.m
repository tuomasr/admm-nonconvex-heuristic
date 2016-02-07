%% About this script
% This script implements an ADMM-heuristic for solving a mixed-integer
% quadratic problem for embedded devices
% The heuristic is based on Takapoui et al. (2015)
% http://stanford.edu/~boyd/papers/pdf/miqp_admm.pdf
% Implementation of the heuristic for this sample problem
% Tuomas Rintamäki 2015
% tuomas.rintamaki@aalto.fi
%% Setup data for a sample problem
% The model represents a n-node chain. Each of the nodes has a load and a generator.
% Depending on the cost parameters of the generators, the nodes can import
% from or export to neighbouring nodes.

% Mathematically, the model solves the following problem
% min (1/2)*x'Qx + c'x + alpha
% s.t. Ax = b
% the first n variables in x contain generation decisions (continuous).
% the following n-1 variables contain flows in the chain (continuous)
% the following n variables are slack variables (continuous) for the generation, i.e.,
% generation <= (max generation)*(generator on/off), which is reformulated as
% generation + slack = max generation
% the last n variables are generator on/off decisions (binary)

clear
n = 50;             % the number of nodes
m = 2*n;            % the number of equality constraints (nodal balance & generation capacity)
v = [randi([1 10],n,1);zeros(3*n-1,1)]; % cost parameters are sampled
Q = diag(v);
c = zeros(4*n-1,1);
A = zeros(m,4*n-1);
alpha = 0;

no_total = 4*n-1;
no_conts = 3*n-1;
no_ints = n;

p_max = 2;

% build the constraint matrix A
for i=1:n
    A(i,i) = 1;      % nodal balance (x)
    A(n+i,i) = 1;    % maximum generation (x)

    if i==1
        A(i,n+i) = 1;
    elseif i>1 && i<n
        A(i,n+i-1) = -1;
        A(i,n+i) = 1;
    elseif i==n
        A(i,n+i-1) = -1;
    end

    A(n+i,2*n-1+i) = 1;         % maximum generation (s)
    A(n+i,3*n-1+i) = -p_max;    % maximum generation (z)
end

d = 1;

if p_max < d
    disp('error: p_max is less than d');
end

b = [ones(n,1)*d;zeros(n,1)];
%% Solve the problem with Gurobi (exactly)
clear model;

model.modelsense = 'min';
model.Q = sparse(Q);
model.A = sparse(A);
model.obj = c;
model.objcon = alpha;
model.sense = '=';
model.rhs = b;
cont = 'C'; int = 'B';
model.vtype = [repmat(cont,no_conts,1);repmat(int,no_ints,1)];
f_min = -2;
f_max = 2;
lb = [zeros(n,1);ones(n-1,1)*f_min;zeros(n,1);zeros(n,1)];
ub = [ones(n,1)*Inf;ones(n-1,1)*f_max;ones(n,1)*Inf;ones(n,1)];
model.lb = lb;
model.ub = ub;

gurobi_write(model, 'qp.lp');
results = gurobi(model);

results.status
xs = results.x
objval = results.objval
%% Solve the problem with non-convex ADMM (approximate)
% set parameters
rho = 5*n;
e_tol = 10^(-7);
%% Preconditioning. This boosts the performance of the algorithm
row_norm = sqrt(sum(A.^2,2));
E = diag(1./row_norm);
F = eye(no_total);
%% Run the algorithm with random restarts
% set the initial best
admm_x_best = zeros(no_total,1);
admm_y_best = Inf;

tic
for i=1:n
    % sample a random starting point. It may not be feasible so its objective value is infinity
    admm_x0 = [p_max.*rand(n,1);(f_max-f_min).*rand(n-1,1)+f_min;p_max.*rand(n,1);ones(n,1)];
    admm_x0((2*n):(3*n-1))=p_max-admm_x0(1:n);
    admm_y0 = Inf;
    
    % optimize
    [admm_x, admm_y] = admm_sparse(admm_x0,admm_y0,m,no_total,Q,c,alpha,A,E,F,b,no_conts,lb,ub,rho,e_tol);
    
    % update the best value
    if admm_y < admm_y_best
        admm_x_best = admm_x;
        admm_y_best = admm_y;
    end
end
toc

% print objective values
admm_y_best
objval
%% Compare Gurobi and ADMM results by plotting them side-by-side
figure
plot(admm_x_best);
hold on
plot(xs)
legend('admm','optimal')