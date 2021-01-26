%% Student infos %%
% Name              - Matriculation number
% Daniel Döhring    - 366448
% Huu Duc Nguyen    - 405242

function [f_lb, f_ub] = convex_bound(n, m, c, H, Q, A, b, lb, ub)
%% Step 1: Check input properties: Dimensions, symmetry...
% Dimensions:
assert(isequal(size(c), [n 1]),   'Dimensions of c are not as expected.');
assert(isequal(size(H), [n n]),   'Dimensions of H are not as expected.');
assert(isequal(size(Q), [m*n n]), 'Dimensions of Q are not as expected.');
assert(isequal(size(A), [m n]),   'Dimensions of A are not as expected.');
assert(isequal(size(b), [m 1]),   'Dimensions of b are not as expected.');
assert(isequal(size(lb), [n 1]),  'Dimensions of lb are not as expected.');
assert(isequal(size(ub), [n 1]),  'Dimensions of ub are not as expected.');

% Matrix properties
% Symmetry of Qi
for i = 0:(m-1)
    Qi = Q(i * n + 1: (i + 1) * n, :);
    assert(issymmetric(Qi), 'Q%d is not symmetric.', i + 1);
end

% Positive semi-definiteness and symmetry of H
assert(all(eig(H)>=0), 'H is not positive semi-definite');
assert(issymmetric(H), 'H is not symmetric.');

%% Step 2: Generate auxiliary variables
% Since Qi is symmetric, we will only add new variables that arise
% considering only upper/lower half of all Qi
Q_Summed = zeros(n);
for i = 0:(m-1)
    Qi = Q(i * n + 1 : (i + 1) * n, :);
    Q_Summed = Q_Summed + abs(Qi); % Use abs to avoid Qi cancels out each others
end

[ind_row, ind_col] = find(triu(Q_Summed) ~= 0);
no_aux_vars = size(ind_row,1); % No. of auxiliary variables

%% Inputs for quadprog
%{
[x,fval] = quadprog(H,f,A,b,Aeq,beq,lb,ub)
min 1/2.x'.H.x + f'.x
s.t.:
A.x <= b
Aeq.x = beq
lb <= x <= ub
%}

% Initialize inputs
N = n + no_aux_vars; % Shape of new optimization variables, x = [x1, .. xn, .. wij]
% Inputs for Objective Functions
H_qp = zeros(N); f_qp = zeros(N,1);
% Inputs for Linear Inequality Constraints
A_qp = zeros(4*no_aux_vars,N); b_qp = zeros(4*no_aux_vars,1);
% Inputs for Linear Equality Constraints
Aeq_qp = zeros(m,N); %beq_qp = zeros(m,1);
% Inputs for Upper/Lower Bounds
lb_qp = zeros(N,1); ub_qp = zeros(N,1);

%% Step 3: Generate matrix B implementing the inequality constraints of the auxiliary variables
% The required matrix B is the matrix A_qp here.
% Inputs for inequality constraints: A_qp, b_qp
for i=1:no_aux_vars
    i_row = ind_row(i);
    i_col = ind_col(i);
    A_qp(4*i-3,[i_row, i_col, n+i]) = [ lb(i_col)  lb(i_row) -1];%  xjL.xi + xiL.xj - wij <=  xiL.xjL
    A_qp(4*i-2,[i_row, i_col, n+i]) = [ ub(i_col)  ub(i_row) -1];%  xjU.xi + xiU.xj - wij <=  xiU.xjU
    A_qp(4*i-1,[i_row, i_col, n+i]) = [-lb(i_col) -ub(i_row)  1];% -xjL.xi - xiU.xj + wij <= -xiU.xjL
    A_qp(4*i,  [i_row, i_col, n+i]) = [-ub(i_col) -lb(i_row)  1];% -xjU.xi - xiL.xj + wij <= -xiL.xjU
    b_qp(4*i-3:4*i) = [lb(i_row)*lb(i_col);
                       ub(i_row)*ub(i_col);
                      -ub(i_row)*lb(i_col);
                      -lb(i_row)*ub(i_col)];
end

%% Step 4: Compose quadratic program with linear constraints
% Inputs for objective function: H_qp, f_qp
H_qp(1:n, 1:n) = 2 * diag(diag(H)); % Matlab expects QP's of the form 0.5*x'*H*x
f_qp(1:n, 1) = c;
for i = 1:no_aux_vars
    f_qp(n+i) = 2*H(ind_row(i),ind_col(i));
end

% Inputs for equality constraints: Aeq_qp, beq_qp
Aeq_qp(:, 1:n) = A;
for i = 1:m
    for k=1:no_aux_vars
        Qi = Q(n*(i-1)+1:n*i,:);
        if ind_row(k) == ind_col(k)
            Aeq_qp(i,n+k) = Qi(ind_row(k),ind_col(k));
        else
            % Every auxiliary variable is visited only once 
            % Due to symmetry of Qi: Off-diagonal components need to be doubled
            Aeq_qp(i,n+k) = 2*Qi(ind_row(k),ind_col(k));
        end
    end
end
beq_qp = b;

% Extended bounds including auxiliary variables
lb_qp(1:n) = lb;
ub_qp(1:n) = ub;
% Set bounds for auxiliary variables as default (-Inf, Inf)
% Sufficient, since all ready have 4 inequalities for each w_ij
for i = 1:no_aux_vars
   lb_qp(n + i) = -Inf; % -Inf for lower bound
   ub_qp(n + i) = Inf;  % Inf for upper bound
end

% Run quadratic programming solver
[~,f_lb] = quadprog(H_qp,f_qp,A_qp,b_qp,Aeq_qp,beq_qp,lb_qp,ub_qp);
% We get f_lb as the underestimation of f

%% Step 5: Apply 'fmincon' to original problem
% Input for objective function
objective = @(x) x'*H*x + c'*x;
% Inputs for linear constraints (both inequality and equality)
A_nlp = []; b_nlp = [];
Aeq_nlp = []; beq_nlp = [];
% Input for initial point
x_0 = 0.5 * (lb + ub); % Use average as initial value
% Input for Non-Linear Constraints
function [c, c_eq] = NonLinearConstraints(x)
    c = []; % No nonlinear inequalities
    c_eq = zeros(m,1);
    for i = 1:m % Loop over the matrices defining the constraints
        Qi = Q((i-1)*n+1:i*n,:);
        c_eq(i)= x'*Qi*x + A(i,:)*x - b(i);
    end
end
nonlcon = @NonLinearConstraints;
%options = optimoptions('fmincon','Display','iter');

% Run fmincon solver for constrained nonlinear optimization problem
[~, f_ub] = fmincon(objective, x_0, A_nlp, b_nlp, Aeq_nlp, beq_nlp, lb, ub, nonlcon);
% We get f_ub as the upper bound of f

%% Display the results
disp('Result for the problem:');
fprintf('The lower bound is: %g.\n', f_lb);
fprintf('The upper bound is: %g.\n', f_ub);

end