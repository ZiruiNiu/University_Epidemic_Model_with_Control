function [t,x] = ode4(f, trange, xinitial, numstep)
% The function is defined to solve the different equation beased on the
% fixed-step fourth order Runge-Kutta method
%
% Inputs: 
%       f: the ODE function handle
%       trange: vector containing the initial and final value of time
%       xinital: initial value of the independent variable x
%       numstep: number of steps to solve the ODE by RK4
%
% Outputs:
%       t: resulting time vector
%       x: a matrix containing values of independent variables with respect
%       to each time in the vector "t"
       
% RK4
    A = [0 0 0 0; 0.5  0 0 0; 0 0.5 0 0; 0 0 1 0];
    b = [1/6 1/3 1/3 1/6];
    c = [0 0.5 0.5 1];

% Initialize values of x
    x = [];
    x(:,1) = xinitial;

% Calculate the step size
    h = trange(end)/numstep;

% Define the time discretisation
    t = 0:h:trange(end);

% Figure out the value of s from the size of matrices A, b, or c
    s = length(b);
    F = zeros(s,length(xinitial));

% Update t and x according to the formulas for RK method using matrices A, b, c
    for ii = 1:(length(t)-1)
        for kk = 1:s
            F(kk,:) = (f((t(ii) + c(kk) * h), (x(:,ii) + h * (A(kk,:) * F)')))';
        end
        x(:,ii+1) = x(:,ii) + h * (b * F)';
    end
    t = t';
    x = x';
end
