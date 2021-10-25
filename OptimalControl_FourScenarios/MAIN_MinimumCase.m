
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% MATLAB Code for deriving and plotting optimal trajectories under minimum
% case scenario based on the designed University COVID-19 pandemic model 
% in the work
%
% Ranking the Effectiveness of Non-pharmaceutical Interventions to counter COVID-19 in Universities with Vaccinated Population
% by Zirui Niu and Giordano Scarciotti
% 
%
% Zirui Niu, October 21, 2021
% Contact: z.niu@outlook.com
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

%% Initialization
clc; clear;
addpath('OptimTraj-master')

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
%%% First Phase: incubation, last 14 days
tRangeP1 = [0 14];
numstepP1 = 50;

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

%%% Get un-controlled 14-day trajectories
% Define the initial values of states
x0_p1 = [param.S_yinit; param.S_sinit; param.E_yinit; param.E_sinit; 0; 0; 0; 0; 0; 0; 0];

% No non-pharmaceutical control measures are conducted
u0 = [0; 0; 0; 1; 1];

% Sove the ODE and get the epidemic evolution trajectories for 14 days
[tSolP1,YSolP1] = ode4(@(t, x)SEIQR(t, x, param), tRangeP1, x0_p1, numstepP1);
x0_p2 = YSolP1(end,:)';

% Compute initial value of R0
R0 = getR0(param, u0)

%%% Second Phase: control, last 120 days
tRangeP2 = [0 120];
numstepP2 = 600;

% Get un-controlled 120-day trajectories: from day 14 to day 134
[~,Yguess] = ode4(@(t, x)SEIQR(t, x, param), tRangeP2, x0_p2, numstepP2);

% Set initial guess for the optimal control problem
guessTime = 360;
guess = Yguess(1:guessTime,:)';


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                       Set the cost function                             %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%%% Set coefficient matrices of the quadratic cost function
% Maximize the number of susceptible
H = diag([-100, -100, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
Q = diag([-100; -100; 0; 0; 0; 0; 0; 0; 0; 0; 0]);

% To guarantee normal convergence
R = diag([5e-4; 5e-4; 5e-4; 5e-4; 5e-4]);

% Define the system dynamics and the cost function
problem1.func.dynamics = @(t,x,u)( StaffStudentDynamics(x,u,param) );
problem1.func.pathObj = @(t,x,u)( pathObj(x, u, Q, R) );
problem1.func.bndObj = @(t0,x0,tF,xF)( xF'*H*xF );


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                      Set the problem constraints                        %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Total control period: 120 days
problem1.bounds.initialTime.low = 0;
problem1.bounds.initialTime.upp = 0;
problem1.bounds.finalTime.low = tRangeP2(2);
problem1.bounds.finalTime.upp = tRangeP2(2);

% Start to control at the 14-th day
problem1.bounds.initialState.low = x0_p2;
problem1.bounds.initialState.upp = x0_p2;

% Set terminal constraints for states
problem1.bounds.finalState.low = zeros(11,1);
problem1.bounds.finalState.upp = ones(11,1);

% Set path constraints for states
problem1.bounds.state.low = zeros(11,1);
problem1.bounds.state.upp = ones(11,1);

% Set path constraints for input control variables. The upper bounds are 
% set to limit the maximum effectiveness because the actual enforcement is 
% not ideal
problem1.bounds.control.low = [0; 0; 0; 1; 1];
problem1.bounds.control.upp = [0.7; 0.7; 0.7; 0.5/param.eta_y; 0.5/param.eta_s];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                      Set feasible initial guess                         %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
% Use simulation results of un-controlled model as initial guesses
problem1.guess.time = 0:1:guessTime;
problem1.guess.state = [x0_p2, guess];
problem1.guess.control = [0.5*ones(3,guessTime+1); 5*ones(2,guessTime+1)];


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                           Solver options                                %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
problem1.options.nlpOpt = optimset(...
            'Display','iter',...
            'TolFun',1e-6,...
            'MaxFunEvals',1e6);
problem1.options.verbose = 3;
problem1.options.method = 'trapezoid';
problem1.options.trapezoid.nGrid = 30;


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%            Solve the problem and obtain optimal trajectories            %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
soln1 = optimTraj(problem1);

% Unpack the simulation
tSolP2 = linspace(soln1.grid.time(1), soln1.grid.time(end), guessTime);
xSolP2 = soln1.interp.state(tSolP2);
uSolP2 = soln1.interp.control(tSolP2);


%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%
%                         Plot optimal trajectories                       %
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%

% Plot the optimal trajectories including the first 14-day un-controlled 
% epidemic evolution
tSol = [tSolP1(1:numstepP1)', (tSolP2 + 14)];
plotOptOut(1, tSol, [YSolP1(1:numstepP1, :)', xSolP2], uSolP2, numstepP1, param)

% Compute variation of R0 in 134 days
uSol2 = [kron(ones(1, numstepP1), [0;0;0;1;1]), uSolP2];
R0 = getR0(param, uSol2, [YSolP1(1:numstepP1, :)', xSolP2]);

% Plot the R0 trajectory during the 134 days
figure(3)
plot(tSol, R0, 'LineWidth', 1.5);
grid on;
xlabel('Time (days)')
ylabel('R0')
title('Variation of R0')
set(gca,'FontSize',26);
set(get(gca,'XLabel'),'FontSize',28);
set(get(gca,'YLabel'),'FontSize',28);
set(get(gca,'ZLabel'),'FontSize',28);
set(get(gca,'title'),'FontSize',30);