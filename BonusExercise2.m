%% Student infos %%
% Name              - Matriculation number
% Daniel Döhring    - 366448
% Huu Duc Nguyen    - 405242

%% Generic optimization problem formulation
%{
min x'*H*x + c'*x
s.t. x'*Qi*x + ai*x = bi with i = 1, 2, .., m (m < n)
xL <= x <= xU
%}

%% Problem 1 %%

clear; clc;
% Inputs for Objective Function
H = zeros(3, 3); H(3, 3) = 1;
c = [1 1 0]';

% Inputs for Non-linear Constraints
Q1 = zeros(3, 3); Q1(1, 2) = 0.5; Q1(2, 1) = 0.5;
Q2 = zeros(3, 3); Q2(3, 2) = 0.5; Q2(2, 3) = 0.5;
Q = [Q1; Q2];
A = zeros(2, 3); A(1, 3) = 1;
b = [8 15]';

% Inputs for Upper/Lower Bounds
lb = [0 0 0]'; ub = [10 10 10]';

[m,n] = size(A);
% Find underestimation and upper bounds
[f_lb, f_ub] = convex_bound(n, m, c, H, Q, A, b, lb, ub);

%% Problem 2 %%

clear; clc;
% Inputs for Objective Function
H = zeros(4); H(3,3) = 1; H(4,4) = 1;
c = [1 1 0 0]';

% Inputs for Non-linear Constraints
Q1 = zeros(4);
Q1(1, 2) = 0.5; Q1(2, 1) = 0.5;
Q1(2, 3) = 0.5; Q1(3, 2) = 0.5;
Q2 = zeros(4); Q2(1, 2) = 0.5; Q2(2, 1) = 0.5;
Q3 = zeros(4); Q3(2, 3) = 0.5; Q3(3, 2) = 0.5;
Q = [Q1; Q2; Q3];
A = zeros(3,4); A(2,4) = 1; A(3,1) = 1;
b = [2 3 5]';

% Inputs for Upper/Lower Bounds
lb = [0 0 0 0]'; ub = [10 4 10 10]';

[m,n] = size(A);
% Find underestimation and upper bounds
[f_lb, f_ub] = convex_bound(n, m, c, H, Q, A, b, lb, ub);