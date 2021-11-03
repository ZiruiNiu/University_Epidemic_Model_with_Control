function R0 = getR0(param, u, x)
% This function is defined to compute the value of R0 by designated model 
% parameter values and the effectiveness of non-pharmaceutical control 
% measures
%
% Inputs:
%       param: the structure variable containing model parameter values.
%       u: input variable matrix
%       x: state variable matrix
%       
% Outputs:
%       R0: a value or a vector of R0

    % Get value of each input control variable
    k_m = u(1,:);
    k_d = u(2,:);
    k_e = u(3,:);
    k_qy = u(4,:);
    k_qs = u(5,:);
    
    % Compute the degree of pomotion or reduction on the affected model 
    % parameters
    u_p = 1 - (1 - 0.4*k_m).*(1 - 0.65*k_d);
    u_m = 0.6*k_m;
    u_e = 1 + 0.3*k_e;
    
    % Compute the resulting model parameter valaues under control of 
    % non-pharmaceutical intervention measures
    beta_yy = param.beta_yy * (1-u_p);
    p = param.p;
    k = param.k;
    epsilon_y = param.epsilon_y;
    epsilon_s = param.epsilon_s;
    eta_y = (0.5 - param.eta_y)*k_qy + param.eta_y;
    eta_s = (0.5 - param.eta_s)*k_qs + param.eta_s;
    sigma_y = eta_y + param.rho_y;
    sigma_s = eta_s + param.rho_s;
    xi_y = param.xi_y;
    xi_s = param.xi_s;
    mu_y = param.mu_y * (1-u_m);
    beta_cy = param.beta_cy * (1-u_m);
    delta = param.delta * u_e;
    
    % Two input arguments: compute R0 at a point 
    if nargin == 2
        S_yinit = param.S_yinit;
        S_sinit = param.S_sinit;
    % Three input arguments: compute R0 over a timeline
    elseif nargin == 3
        S_yinit = x(1,:);
        S_sinit = x(2,:);
    else
        error('Wrong number of inputs');
    end
    
    % Compute R0 based on the derived equation
    a = beta_yy.*(S_yinit.*(sigma_y + epsilon_y.*k)./(sigma_y.*(epsilon_y+xi_y)) + S_sinit.*p.*(sigma_s + epsilon_s.*k)./(sigma_s.*(epsilon_s+xi_s)));
    c = mu_y.*beta_cy./(delta.*beta_yy);
    R0 = (a/2 .* (1+sqrt(1+4.*c./a)))';
end