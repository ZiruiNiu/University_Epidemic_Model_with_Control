function dYdt = SEIQR2(t,Y, param)  
% The function represents the original un-controlled epidemic model.
% Having 11 compartments, the system is composed of 11 differential 
% equations
%
% Inputs:
%       t: time vectior
%       Y: state vector
%       param: the structure variable containing model parameter values.
%
% Outputs:
%       dYdt: the differential of state vector

    % Expand the state variable
    S_y = Y(1);     % susceptible students
    S_s = Y(2);     % susceptible staffs
    E_y = Y(3);     % exposed students
    E_s = Y(4);     % exposed staffs
    I_y = Y(5);     % infected students
    I_s = Y(6);     % infected staffs
    Q_y = Y(7);     % quarantined students
    Q_s = Y(8);     % quarantined staffs
    R_y = Y(9);     % recovered students
    R_s = Y(10);    % recovered staffs
    C = Y(11);      % viral environment concentration
    
    % Get the value of model parameters
    beta_yy = param.beta_yy;
    beta_ss = param.beta_ss;
    beta_cy = param.beta_cy;
    beta_cs = param.beta_cs;
    k = param.k;
    epsilon_y = param.epsilon_y;
    epsilon_s = param.epsilon_s;
    xi_y = param.xi_y;
    xi_s = param.xi_s;
    eta_y = param.eta_y;
    eta_s = param.eta_s;
    rho_y = param.rho_y;
    rho_s = param.rho_s;
    phi_y = param.phi_y;
    phi_s = param.phi_s;
    mu_y = param.mu_y;
    mu_s = param.mu_s;
    delta = param.delta;
    
    % Define the differential equations
    dSydt = -(beta_yy*(E_y + k*I_y + E_s + k*I_s) + beta_cy*C)*S_y;
    dSsdt = -(beta_ss*(E_s + k*I_s + E_y + k*I_y) + beta_cs*C)*S_s;
    dEydt = (beta_yy*(E_y + k*I_y + E_s + k*I_s) + beta_cy*C)*S_y - (epsilon_y + xi_y)*E_y;
    dEsdt = (beta_ss*(E_s + k*I_s + E_y + k*I_y) + beta_cs*C)*S_s - (epsilon_s + xi_s)*E_s;
    dIydt = epsilon_y*E_y - (eta_y + rho_y)*I_y;
    dISdt = epsilon_s*E_s - (eta_s + rho_s)*I_s;
    dQydt = eta_y*I_y - phi_y*Q_y;
    dQsdt = eta_s*I_s - phi_s*Q_s;
    dRydt = xi_y*E_y + rho_y*I_y + phi_y*Q_y;
    dRsdt = xi_s*E_s + rho_s*I_s + phi_s*Q_s;
    dCdt  = mu_y*(E_y + k*I_y)+ mu_s*(E_s + k*I_s) - delta*C;

    % Create output column vector representing differential of state variable
    dYdt = [dSydt; dSsdt; dEydt; dEsdt; dIydt; dISdt; dQydt; dQsdt; dRydt; dRsdt; dCdt];
    
end