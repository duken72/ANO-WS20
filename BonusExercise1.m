%% Student infos %%
%%%%%%%%%%%%%%%%%%%

% Name           - Matriculation number
% Daniel Döhring                     - 366448
% Huu Duc Nguyen                     - 405242
% Philipus Benizi Angga Aditya Putra - 402726

clear;
clc;

%% System Constants and Parameters %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turbine data
% Power limitations
Min_P1 = 2500;     % [kW]
Max_P1 = 6250;     % [kW]
Min_P2 = 3000;     % [kW]
Max_P2 = 9000;     % [kW]
% Flow limitations
Max_I1   = 87000;  % [kg/h]
Max_I2   = 110000; % [kg/h]
Max_LPS  = 64000;  % [kg/h]
Max_Cond = 28000;  % [kg/h] Condensate flow
Max_Int  = 60000;  % [kg/h] Internal flow

% Vapor data (specific enthalpies)
h_HPS = 3163; % [kJ/kg]
h_MPS = 2949; % [kJ/kg]
h_LPS = 2911; % [kJ/kg]
h_CW  = 449;  % [kJ/kg] CW ^= Condensed Water

% System requirement data
Min_MPS   = 123000; % [kg/h]
Min_LPS   = 45000;  % [kg/h]
Min_ElecP = 24500;  % [kW] Electrical Power

% Energy data
C_Fuel   = 1.5;   % [€/TJ = € / (10^6 kJ)]
eta_Evap = 0.75;  % [1]      Evaporator
L_Cond   = 0.008; % [€/kg*h] Loss (of) Condensate 
C_Prod   = 0.02;  % [€/kW*h] Cost of Produced  Power
C_Purch  = 0.05;  % [€/kW*h] Cost of PURCHASED Power 
C_Penal  = 0.001; % [€/kW*h] Demand Penalty
P_Basic  = 12000; % [kW]     Basic power

%% Variable Conventions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% x_1  = P1, Power of Turbine 1           [kW]
% x_2  = P2, Power of Turbine 2           [kW]
% x_3  = PP, Produced Power               [kW]
% x_4  = EP, Electrical Power             [kW]
% x_5  = Power                            [kW]     
% x_6  = Fuel (Energy)                    [kJ]
% x_7  = C, Condensate Flow Rate          [kg/h]
% x_8  = I1, Inflow Turbine 1             [kg/h]
% x_9  = I2, Inflow Turbine 2             [kg/h]
% x_10 = HE1, Outflow rate Turbine 1 high [kg/h]
% x_11 = HE2, Outflow rate Turbine 2 high [kg/h]
% x_12 = LE1, Outflow rate Turbine 1 low  [kg/h]
% x_13 = LE2, Outflow rate Turbine 2 low  [kg/h]
% x_14 = BF1, By-pass flow rate valve 1   [kg/h]
% x_15 = BF2, By-pass flow rate valve 2   [kg/h]
% x_16 = HPS, High   pressure flow rate   [kg/h]
% x_17 = MPS, Medium pressure flow rate   [kg/h]
% x_18 = LPS, low    pressure flow rate   [kg/h]

%% Equality Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%

A_Eq = zeros(10, 18); % Dimensions: #Eq. Constraints * #Variables
b_Eq = zeros(10, 1);

%% Mass Balances 
% MB1 - transformed such that I1 + I2 + BF1 - HPS = 0
A_Eq(1, [8, 9, 14, 16]) =    [1,   1,   1,   -1];

% MB2 - transformed such that HE1 + LE1 + C - I1 = 0
A_Eq(2, [10, 12, 7, 8]) =    [1,    1,    1, -1];

% MB3 - transformed such that HE2 + LE2 - I2 = 0
A_Eq(3, [11, 13, 9]) =       [1,    1,   -1];

% MB4 - transformed such that    HE1 + HE2 + BF1 - BF2 - MPS = 0
A_Eq(4, [10, 11, 14, 15, 17]) = [1,   1,    1,   -1,   -1];

% MB5 - transformed such that LE1 + LE2 + BF2 - LPS = 0
A_Eq(5, [12, 13, 15, 18]) =  [1,    1,    1,   -1];

%% Energy Balances
% EB1 - transformed such that h_MPS*HE1 + h_LPS*LE1 + h_CW*C + 3600*P1 - h_MPS*I1= 0
A_Eq(6, [10, 12, 7, 1, 8]) = [h_MPS,      h_LPS,      h_CW,    3600,   - h_HPS];

% EB2 - transformed such that h_MPS*HE2 + h_LPS*LE2 + 3600*P2 - h_MPS*I2 = 0
A_Eq(7, [11, 13, 2, 9]) =    [h_MPS,      h_LPS,      3600,   - h_HPS];

% EB3 - transformed such that P1 + P2 - PP = 0
A_Eq(8, [1, 2, 3]) =         [1,  1,  -1];

% EB4 - transformed such that PP + EP - Power = 0
A_Eq(9, [3, 4, 5]) =         [1,   1,  -1];

% EB5 - transformed such that h_HPS * HPS - eta_Evap * Fuel = 0
A_Eq(10, [16, 6]) =          [h_HPS,      - eta_Evap]; 

%% Inequality Constraints %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% "True" inequality constraints (not simply bounds)
A_InEq = zeros(8, 18); % Dimensions: #InEq. Constraints * #Variables
b_InEq = zeros(8, 1);

% IC1, upper bound: P1 <= Max_P1
A_InEq(1, 1) = 1;        % +P1
b_InEq(1, 1) = Max_P1;   % +6250

% IC2: I1 <= 87000
A_InEq(2, 8) = 1;        % +I1
b_InEq(2, 1) = Max_I1;   % +87000

% IC3: C <= 28000
A_InEq(3, 7) = 1;        % +C
b_InEq(3, 1) = Max_Cond; % +28000

% IC4: I1 - HE1 <= Max_Int
A_InEq(4, [8, 10]) = [1, -1]; % +I1 - HE1
b_InEq(4, 1)  = Max_Int;      % +60000[kg/h]

% IC5, upper bound: P2 <= Max_P2
A_InEq(5, 2) = 1;        % +P2
b_InEq(5, 1) = Max_P2;   % +9000

% IC6: I2 <= Max_I2
A_InEq(6, 9) = 1;        % +I1
b_InEq(6, 1) = Max_I2;   % +110000

% IC7: LE2 <= Max_LPS
A_InEq(7, 13) = 1;       % +LE2
b_InEq(7, 1)  = Max_LPS; % +64000

% IC8, IC9: Handled by lower bounds

% IC10: -P1 -P2 -EP <= -24500
A_InEq(8, [1, 2, 4]) = [-1, -1, -1];
b_InEq(8, 1) = -Min_ElecP; % - 24500

% Lower bound (0) for most optimization variables (IC11) 
LowerBound = zeros(1, 18);

% IC1, lower bound: 2500 <= P1
LowerBound(1) = Min_P1;

% IC5, lower bound: 3000 <= P2
LowerBound(2) = Min_P2;

% IC8: MPS >= 123000
LowerBound(17) = Min_MPS;

% IC9: LPS >= 45000
LowerBound(18) = Min_LPS;

%% The Optimization Cases %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case a: EP <= 12000 %
Case_a_constraint    = zeros(1, 18);
Case_a_constraint(4) = 1;
A_InEq_a = [A_InEq; Case_a_constraint]; % Append 
b_InEq_a = [b_InEq; P_Basic];

f_a               = zeros(1, 18);
f_a([3, 4, 6, 7]) = [C_Prod, C_Purch - C_Penal, C_Fuel * 10^(-6), L_Cond];

% Case b: 12000 <= EP %
LowerBound_b    = LowerBound;
LowerBound_b(4) = P_Basic;

f_b               = zeros(1, 18);
f_b([3, 4, 6, 7]) = [C_Prod, C_Purch, C_Fuel * 10^(-6), L_Cond];

%% Run the optimization solver %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_a    = linprog(f_a, A_InEq_a, b_InEq_a, A_Eq, b_Eq, LowerBound, []);
Cost_a = f_a * x_a + C_Penal * P_Basic;

x_b    = linprog(f_b, A_InEq, b_InEq, A_Eq, b_Eq, LowerBound_b, []);
Cost_b = f_b * x_b;

%% Display results %%
%%%%%%%%%%%%%%%%%%%%%

clc;
disp('Optimal cost for each case:');
fprintf('Case a: %d €/hour.\n', round(Cost_a));
fprintf('Case b: %d €/hour.\n', round(Cost_b));
fprintf('Thus, final optimal cost value is: %d €/hour.\n', min(round(Cost_a), round(Cost_b)));
if Cost_a < Cost_b
    disp('The Demand Penalty is ACTIVE');
else
    disp('The Demand Penalty is NOT ACTIVE');
end