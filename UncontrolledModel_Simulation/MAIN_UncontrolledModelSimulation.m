
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% MATLAB Code for simulating the un-controlled COVID-19 University epidemic evolution based on the designed University model
% in the work
%
% Ranking the Effectiveness of Non-pharmaceutical Interventions to Counter COVID-19 in Universities with Vaccinated Population
% by Zirui Niu and Giordano Scarciotti
% 
%
% Zirui Niu, October 21, 2021
% Contact: z.niu@outlook.com
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%


%% Initialization
clc; clear;

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%       Model parameter initialization with vaccinated population         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%    
param.beta_yy = 0.0678;
param.p = 1.48;
param.beta_ss = 0.1006;
param.k = 1.5;
param.beta_cy = 0.0901;
param.beta_cs = 0.1338;
param.epsilon_y = 0.2;
param.epsilon_s = 0.1;
param.xi_y = 0.0882;
param.xi_s = 0.0176;
param.eta_y = 0.012;
param.eta_s = 0.0212;
param.rho_y = 1/8;
param.rho_s = 1/12;
param.phi_y = 1/14;
param.phi_s = 1/14;
param.mu_y = 0.3;
param.mu_s = 0.3;
param.delta = 0.7;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%              Set initial conditions of the optimal problem              %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Set the total number of studnets and staffs in the University department
param.Pop_y = 1200;
param.Pop_s = 150;
param.Pop = param.Pop_y + param.Pop_s;

% Number of exposed students at day 1
E_yinit = 5;
S_yinit = param.Pop_y - E_yinit;

% Number of exposed staffs at day 1
E_sinit = 2;
S_sinit = param.Pop_s - E_sinit;

% Get the fraction of susceptible and exposed compartments
param.S_yinit = S_yinit/param.Pop;
param.E_yinit = E_yinit/param.Pop;
param.S_sinit = S_sinit/param.Pop;
param.E_sinit = E_sinit/param.Pop;

% Set time duration (250 days) and number of steps (1000)
numstep = 1000;
tRange = [0 250];


%% Solve the ODE and run the simulation
% Define the initial values of states
x0 = [param.S_yinit; param.S_sinit; param.E_yinit; param.E_sinit; 0; 0; 0; 0; 0; 0; 0];

% No non-pharmaceutical control measures are conducted
u0 = [0; 0; 0; 1; 1];

% Compute initial value of R0
R0 = getR0(param, u0)

% Sove the ODE and get the epidemic evolution trajectories for 250 days
[tSol,YSol] = ode4(@(t, x)SEIQR(t, x, param), tRange, x0, numstep);

% Plot the trajectories of student compartments and staff compartments
plotTraj(1, tSol, YSol', param);

% Compute the trajectory of R0 across the 250 days
R0_time = getR0(param, [0; 0; 0; 1; 1], YSol');

% Plot the R0 trajectory
figure(4)
plot(tSol, R0_time, 'LineWidth', 1.5);
grid on;
xlabel('Time (days)')
ylabel('R0')
title('Variation of R0')
set(gca,'FontSize',26);
set(get(gca,'XLabel'),'FontSize',28);
set(get(gca,'YLabel'),'FontSize',28);
set(get(gca,'ZLabel'),'FontSize',28);
set(get(gca,'title'),'FontSize',30);

