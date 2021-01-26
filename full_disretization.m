%% Student infos %%

% Name                               - Matriculation number
% Daniel Döhring                     - 366448
% Huu Duc Nguyen                     - 405242
% Philipus Benizi Angga Aditya Putra - 402726

function optVal = full_disretization(phi_x, der_x, xmin, xmax, umin, umax, x0, tf, M, u0)

%{
Reformulate given dynamic optimization problem to nonlinear optimization problem

fmincon:
min f(x)
s.t.
c(x) <= 0; ceq(x) <= 0;
A.x <= b;  Aeq.x = beq;
lb <= x <= ub;

New optimization variable of size 4*M
x = [u1, x1, u2, x2, ..., uM, xM]
%}


%% Input for objective function
obj = @(x) phi_x( x(4*M-2:4*M) );

%% Input for initial point
x_0 = repmat([u0; x0], M, 1);

%% Inputs for linear inequality & equality constraints
A = []; b = [];
Aeq = []; beq = [];

%% Inputs for lower / upper bounds
lb = repmat([umin; xmin], M, 1);
ub = repmat([umax; xmax], M, 1);

%% Input for Non-Linear Constraints
dt = tf / M;
function [c, ceq] = NonLinearConstraints(x)
    c = []; % No nonlinear inequalities
    ceq = zeros(3*M,1);
    for i = 1:M
        xcur = x(4*i-2:4*i); % x(k)
        if i==1
            xprev = x0; % use value of x0 for constraints on x(t=1) 
        else
            xprev = x(4*i-6:4*i-4); % x(k-1)
        end
        uk = x(4*i-3); % u(k)
        ceq(3*i-2:3*i) = xcur - xprev - dt * der_x(xcur, uk);
        % c = x(k) - x(k-1) - dt . f(x(k+1),u(k+1)) = 0
    end
end
nonlcon = @NonLinearConstraints;

%% Apply 'fmincon' to discretized problem
% Run fmincon solver
options = optimoptions('fmincon','MaxFunctionEvaluations',27000);
[~, optVal] = fmincon(obj, x_0, A, b, Aeq, beq, lb, ub, nonlcon,options);

%options = optimoptions('fmincon','Display','iter','MaxFunctionEvaluations',27000);
%[~, optVal] = fmincon(obj, x_0, A, b, Aeq, beq, lb, ub, nonlcon, options);

%clc;
fprintf('The optimal value is: %g.\n', optVal);

end