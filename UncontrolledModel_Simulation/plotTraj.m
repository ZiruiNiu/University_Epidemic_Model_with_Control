function plotTraj(plotNum, tSol, YSol, param)
% The function is used to plot the simulation results of the original
% un-controlled epidemic model
%
% Inputs:
%       plotNum: figure number
%       tSol: the row vector that containing time values
%       YSol: a matrix containing values of the independent variables
%       param: the structure that contains all model parameter values

    % Expand the state variable
    S_y = YSol(1,:)';     % susceptible students
    S_s = YSol(2,:)';     % susceptible staffs
    E_y = YSol(3,:)';     % exposed students
    E_s = YSol(4,:)';     % exposed staffs
    I_y = YSol(5,:)';     % infected students
    I_s = YSol(6,:)';     % infected staffs
    Q_y = YSol(7,:)';     % quarantined students
    Q_s = YSol(8,:)';     % quarantined staffs
    R_y = YSol(9,:)';     % recovered students
    R_s = YSol(10,:)';    % recovered staffs
    C = YSol(11,:)';      % viral environment concentration
    
    % For student compartments, convert the ratio of total population to the 
    % ratio of student population
    S_yp = param.Pop * S_y / param.Pop_y;
    E_yp = param.Pop * E_y / param.Pop_y;
    I_yp = param.Pop * I_y / param.Pop_y;
    Q_yp = param.Pop * Q_y / param.Pop_y;
    R_yp = param.Pop * R_y / param.Pop_y;

    % For staff compartments, convert the ratio of total population to the
    % ratio of staff population
    S_sp = param.Pop * S_s / param.Pop_s;
    E_sp = param.Pop * E_s / param.Pop_s;
    I_sp = param.Pop * I_s / param.Pop_s;
    Q_sp = param.Pop * Q_s / param.Pop_s;
    R_sp = param.Pop * R_s / param.Pop_s;
    
    % Plot the trajectories of student compartments
    figure(plotNum)
    subplot(2,1,1)
    plot(tSol, S_yp, 'LineWidth', 1.5);
    hold on;
    plot(tSol, E_yp, 'LineWidth', 1.5);
    hold on;
    plot(tSol, I_yp, 'LineWidth', 1.5);
    hold on;
    plot(tSol, Q_yp, 'LineWidth', 1.5);
    hold on;
    plot(tSol, R_yp, 'LineWidth', 1.5);
    hold off;
    axis([0, tSol(end), 0, 1])
    grid on;
    legend('S_y(t)','E_y(t)','I_y(t)', 'Q_y(t)', 'R_y(t)')
    xlabel('Time (days)')
    ylabel('Student Case Proportion')
    title('Trajectories of Student Group')
    set(legend,'FontSize',20);
    set(gca,'FontSize',18);
    set(get(gca,'XLabel'),'FontSize',20);
    set(get(gca,'YLabel'),'FontSize',20);
    set(get(gca,'ZLabel'),'FontSize',20);
    set(get(gca,'title'),'FontSize',22);
    
    % Plot the trajectories of staff compartments
    subplot(2,1,2)
    plot(tSol, S_sp, 'LineWidth', 1.5);
    hold on;
    plot(tSol, E_sp, 'LineWidth', 1.5);
    hold on;
    plot(tSol, I_sp, 'LineWidth', 1.5);
    hold on;
    plot(tSol, Q_sp, 'LineWidth', 1.5);
    hold on;
    plot(tSol, R_sp, 'LineWidth', 1.5);
    hold off;
    axis([0, tSol(end), 0, 1])
    grid on;
    legend('S_s(t)','E_s(t)','I_s(t)', 'Q_s(t)', 'R_s(t)')
    xlabel('Time (days)')
    ylabel('Staff Case Proportion')
    title('Trajectories of Staff Group')
    set(legend,'FontSize',20);
    set(gca,'FontSize',18);
    set(get(gca,'XLabel'),'FontSize',20);
    set(get(gca,'YLabel'),'FontSize',20);
    set(get(gca,'ZLabel'),'FontSize',20);
    set(get(gca,'title'),'FontSize',22);
    
    % Plot the trajectory of environmental concentration compartment
    figure(plotNum + 2)
    plot(tSol, C, 'LineWidth', 1.5);
    grid on;
    xlabel('Time (days)')
    ylabel('Concentration C(t)')
    title('Concentration of Virus in the Environment')
    set(gca,'FontSize',26);
    set(get(gca,'XLabel'),'FontSize',28);
    set(get(gca,'YLabel'),'FontSize',28);
    set(get(gca,'ZLabel'),'FontSize',28);
    set(get(gca,'title'),'FontSize',30);
end