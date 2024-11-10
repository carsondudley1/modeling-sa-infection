clear
close all

% parameter values
tspan = [0 1]; % I just want to look at a short section of time


sa0 = 1e9; % initial pneumococcal cells 1e3-1e5
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
init1 = [sa0, 0, 0, 0, 8e5, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
[t1, state1] = ode15s(@(t,state)sa_perit_rhs_hill(t, state), tspan, init1);

init2 = [sa0, 0, 0, 0, 8e5, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
init3 = [sa0, 0, 0, 0, 8e5, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

[t2, state2] = ode15s(@(t,state)sa_saf_rhs_hill(t, state), tspan, init2);

[t3, state3] = ode15s(@(t,state)sa_leukfg_rhs_hill(t, state), tspan, init3);

init4 = [sa0, 8e5, 0, 0];
[t4, state4] = ode15s(@(t,state)sa_nofg_rhs(t, state), tspan, init4);


%% Plotting the results

figure(1)
semilogy(t4, state4(:, 1), t3, state3(:, 1), t2, state2(:, 1), t1, state1(:, 1), 'LineWidth', 4)
ylim([10^(6) 10^10])



% hold on
title('SA Count in CFUs vs. Time in Hours')
xlabel('Time (Hours)')
ylabel('Bacteria (CFUs)')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf, 'Color', 'w');
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex')
set(findall(gcf,'-property','FontName'),'FontName','CMU Serif')
legend('Reduced Model', 'Model 2', 'Model 3', 'Complete Model')
