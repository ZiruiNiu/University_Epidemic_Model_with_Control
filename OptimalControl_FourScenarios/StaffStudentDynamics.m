function Y = StaffStudentDynamics(x, u, param)
% The function describe the dynamics of the overall University epidemic 
% model under control of different non-pharmaceutical intervention measures.
%
% Inputs:
%       x: state variable matrix
%       u: input variable matrix
%       param: the structure variable containing model parameter values.
%
% Outputs:
%       Y: the differential of state vector

    % Expand the state variable
    S_y = x(1,:);       % susceptible students
    S_s = x(2,:);       % susceptible staffs
    E_y = x(3,:);       % exposed students
    E_s = x(4,:);       % exposed staffs
    I_y = x(5,:);       % infected students
    I_s = x(6,:);       % infected staffs
    Q_y = x(7,:);       % quarantined students
    Q_s = x(8,:);       % quarantined staffs
    R_y = x(9,:);       % recovered students
    R_s = x(10,:);      % recovered staffs
    C = x(11,:);        % viral environment concentration
    
    k_m = u(1,:);       % effectiveness of mask wearing
    k_d = u(2,:);       % effectiveness of social-distancing
    k_e = u(3,:);       % effectiveness of environmental disinfection
    u_qy = u(4,:);       % effectiveness of quarantine on infected students
    u_qs = u(5,:);       % effectiveness of quarantine on infected staffs
    
    % convert effectiveness of each measure to the degree of promotion or reduction on
    % the affected model parameter
    u_p = 1 - (1 - 0.4*k_m).*(1 - 0.65*k_d);
    u_m = 0.6*k_m;
    u_e = 1 + 0.3*k_e;
    
    % Express the system dynamics by differential equations
    Y(1,:) = -((1 - u_p).*param.beta_yy.*(E_y + param.k.*I_y + E_s + param.k.*I_s) + (1 - u_m).*param.beta_cy.*C).*S_y;
    Y(2,:) = -((1 - u_p).*param.beta_ss.*(E_s + param.k.*I_s + E_y + param.k.*I_y) + (1 - u_m).*param.beta_cs.*C).*S_s;
    Y(3,:) = ((1 - u_p).*param.beta_yy.*(E_y + param.k.*I_y + E_s + param.k.*I_s) + (1 - u_m).*param.beta_cy.*C).*S_y - (param.epsilon_y + param.xi_y).*E_y;
    Y(4,:) = ((1 - u_p).*param.beta_ss.*(E_s + param.k.*I_s + E_y + param.k.*I_y) + (1 - u_m).*param.beta_cs.*C).*S_s - (param.epsilon_s + param.xi_s).*E_s;
    Y(5,:) = param.epsilon_y.*E_y - (u_qy.*param.eta_y + param.rho_y).*I_y;
    Y(6,:) = param.epsilon_s.*E_s - (u_qs.*param.eta_s + param.rho_s).*I_s;
    Y(7,:) = u_qy.*param.eta_y.*I_y - param.phi_y.*Q_y;
    Y(8,:) = u_qs.*param.eta_s.*I_s - param.phi_s.*Q_s;
    Y(9,:) = param.xi_y.*E_y + param.phi_y.*Q_y + param.rho_y.*I_y;
    Y(10,:) = param.xi_s.*E_s + param.phi_s.*Q_s + param.rho_s.*I_s;
    Y(11,:) = (1 - u_m)*param.mu_y.*(E_y + param.k.*I_y)+ (1 - u_m)*param.mu_s.*(E_s + param.k.*I_s) - u_e * param.delta.*C;
    
end
        