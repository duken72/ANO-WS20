%% Student infos %%

% Name                               - Matriculation number
% Daniel Döhring                     - 366448
% Huu Duc Nguyen                     - 405242
% Philipus Benizi Angga Aditya Putra - 402726

clear;
clc;

%% Generic optimization problem formulation
%{
min PHI(x) = phi(x(tf))
s.t.
der(x) = f(x, u)
x(0) = x0
xmin <= x <= xmax
umin <= u <= umax
%}

%% Problem %%
phi_x = @(x) x(3);
der_x = @(x, u) [(1-x(2)^2)*x(1) - x(2) + u;
                 x(1);
                 x(1)^2 + x(2)^2 + u^2];
x0 = [0, 1, 0]';
xmin = [-0.4 -Inf -Inf]';
xmax = [Inf Inf Inf]';
umin = -0.3; umax = 1;
tf = 5;

%% Solve Optimization problem with Full Discretization
M = 50;
u0 = 0;

optVal = full_disretization(phi_x, der_x, xmin, xmax, umin, umax, x0, tf, M, u0);

%% Single Shooting Approach
%{
Objective value in general would decrease??
%}