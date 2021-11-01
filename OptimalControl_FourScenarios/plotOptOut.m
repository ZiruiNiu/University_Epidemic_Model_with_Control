function plotOptOut(plotNum, tSol, YSol, uSol, numstepP1, param)
% The function is used to plot the derived optimal trajectories from the 
% "OptimTraj" library
%
% Inputs:
%       plotNum: figure number
%       tSol: the row vector that containing time values
%       YSol: a matrix containing values of the state variables
%       uSol: a matrix containing values of the input variables
%       numstepP1: the number of steps when sovling the epidemic evolution
%       in the first 14 days
%       param: the structure that contains all model parameter values

    % Get trajectory of each compartment
    S_y = YSol(1,:);     % susceptible students
    S_s = YSol(2,:);     % susceptible staffs
    E_y = YSol(3,:);     % exposed students
    E_s = YSol(4,:);     % exposed staffs
    I_y = YSol(5,:);     % infected students
    I_s = YSol(6,:);     % infected staffs
    Q_y = YSol(7,:);     % quarantined students
    Q_s = YSol(8,:);     % quarantined staffs
    R_y = YSol(9,:);     % recovered students
    R_s = YSol(10,:);    % recovered staffs
    C = YSol(11,:);      % viral environment concentration
    
    % For student compartments, convert the ratio of total population to the 
    % ratio of student population
    S_yp = param.Pop * S_y / param.Pop_y;
    E_yp = param.Pop * E_y / param.Pop_y;
    I_yp = param.Pop * I_y / param.Pop_y;
    Q_yp = param.Pop * Q_y / param.Pop_y;
    R_yp = param.Pop * R_y / param.Pop_y;
    SEIQR_yp = [S_yp; E_yp; I_yp; Q_yp; R_yp];

    % For staff compartments, convert the ratio of total population to the
    % ratio of staff population
    S_sp = param.Pop * S_s / param.Pop_s;
    E_sp = param.Pop * E_s / param.Pop_s;
    I_sp = param.Pop * I_s / param.Pop_s;
    Q_sp = param.Pop * Q_s / param.Pop_s;
    R_sp = param.Pop * R_s / param.Pop_s;
    SEIQR_sp = [S_sp; E_sp; I_sp; Q_sp; R_sp];
    
    % Plot the trajectories of student compartments
    figure(plotNum)
    subplot(2,2,1)
    plot(tSol, SEIQR_yp, 'LineWidth', 1.2);
    axis([0, tSol(end), 0, 1])
    grid on;
    legend('S_y(t)','E_y(t)','I_y(t)', 'Q_y(t)', 'R_y(t)')
    xlabel('Time (days)')
    ylabel('Student Case Proportion')
    title('Trajectories of Student Group')
    set(legend,'FontSize',18);
    set(gca,'FontSize',18);
    set(get(gca,'XLabel'),'FontSize',20);
    set(get(gca,'YLabel'),'FontSize',20);
    set(get(gca,'ZLabel'),'FontSize',20);
    set(get(gca,'title'),'FontSize',22);
    
    % Plot the trajectories of staff compartments
    subplot(2,2,2)
    plot(tSol, SEIQR_sp, 'LineWidth', 1.2);
    axis([0, tSol(end), 0, 1])
    grid on;
    legend('S_s(t)','E_s(t)','I_s(t)', 'Q_s(t)', 'R_s(t)')
    xlabel('Time (days)')
    ylabel('Staff Case Proportion')
    title('Trajectories of Staff Group')
    set(legend,'FontSize',18);
    set(gca,'FontSize',18);
    set(get(gca,'XLabel'),'FontSize',20);
    set(get(gca,'YLabel'),'FontSize',20);
    set(get(gca,'ZLabel'),'FontSize',20);
    set(get(gca,'title'),'FontSize',22);

    % Plot trajectories of the first three non-quarantine control variables
    subplot(2,2,3)
    plot(tSol, [zeros(3, numstepP1), uSol(1:3,:)], 'LineWidth', 1.2);
    grid on;
    legend('\kappa_m(t)','\kappa_d(t)','\kappa_e(t)');
    axis([0, tSol(end), 0, 1])
    xlabel('Time (days)')
    ylabel('Control Variable Values')
    title('\kappa_m, \kappa_d and \kappa_e')
    set(legend,'FontSize',18);
    set(gca,'FontSize',18);
    set(get(gca,'XLabel'),'FontSize',20);
    set(get(gca,'YLabel'),'FontSize',20);
    set(get(gca,'ZLabel'),'FontSize',20);
    set(get(gca,'title'),'FontSize',22);

    % Plot trajectories of the two quarantine control variables
    subplot(2,2,4)
    plot(tSol, [zeros(2, numstepP1), uSol(4:5,:)], 'LineWidth', 1.2);
    grid on;
    legend('\kappa_{qy}(t)', '\kappa_{qs}(t)');
    axis([0, tSol(end), 0, 1])
    xlabel('Time (days)')
    ylabel('Control Variable Values')
    title('\kappa_{qy} and \kappa_{qs}')
    set(legend,'FontSize',18);
    set(gca,'FontSize',18);
    set(get(gca,'XLabel'),'FontSize',20);
    set(get(gca,'YLabel'),'FontSize',20);
    set(get(gca,'ZLabel'),'FontSize',20);
    set(get(gca,'title'),'FontSize',22);
    
    % Plot trajectories of all five control variables in one figure
    figure(plotNum+1)
    plot(tSol, [zeros(5, numstepP1), uSol(1:5,:)], 'LineWidth', 1.2);
    grid on;
    legend('\kappa_m(t)','\kappa_d(t)','\kappa_e(t)','\kappa_{qy}(t)', '\kappa_{qs}(t)');
    axis([0, tSol(end), 0, 1])
    xlabel('Time (days)')
    ylabel('Control Variable Values')
    title('Effectivenesses of Different Measures')
    set(legend,'FontSize',20);
    set(gca,'FontSize',18);
    set(get(gca,'XLabel'),'FontSize',20);
    set(get(gca,'YLabel'),'FontSize',20);
    set(get(gca,'ZLabel'),'FontSize',20);
    set(get(gca,'title'),'FontSize',22);
    
end