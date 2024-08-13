clear
close all

% parameter values
tspan = [0 4]; % I just want to look at a short section of time


sa0 = 1e9; % initial pneumococcal cells 1e3-1e5
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9);
init1 = [sa0, 0, 0, 0, 8e5, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
[t1, state1] = ode15s(@(t,state)sa_perit_rhs(t, state), tspan, init1);

init2 = [sa0, 0, 0, 0, 8e5, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
[t2, state2] = ode15s(@(t,state)sa_saf_rhs(t, state), tspan, init2);

init3 = [sa0, 0, 0, 0, 8e5, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
[t3, state3] = ode15s(@(t,state)sa_leukfg_rhs(t, state), tspan, init3);

init4 = [sa0, 8e5, 0, 0];
[t4, state4] = ode15s(@(t,state)sa_nofg_rhs(t, state), tspan, init4);


%% Plotting the results

figure(1)
% semilogy(t4, state4(:,1), t3, state3(:, 1), t2, state2(:, 1), t1, state1(:, 1), 'LineWidth', 4)
% ylim([10^(6) 10^10])
% plot(t, state(:,10)+state(:,11)+2*state(:,12)+3*state(:,13)+4*state(:,14)+5*state(:,15)+6*state(:,16)+7*state(:,17)+8*state(:,18)+9*state(:,19)+10*state(:,20)+state(:,23)+state(:,24))
% plot(t, state(:, 11) + 2*state(:,12)+3*state(:,13)+4*state(:,14)+5*state(:,15)+6*state(:,16)+7*state(:,17)+8*state(:,18)+9*state(:,19)+10*state(:,20)+state(:,23)+state(:,24))
plot(t3, state3(:,6) + state3(:, 9), t3, state3(:,5) + state3(:, 8) + state3(:,6) + state3(:, 9), 'LineWidth', 4)
% plot(t, state(:, 26))
% plot(t, state(:, 8) + state(:,7) + state(:, 10) + state(:, 11))
% plot(t, state(:, 23) + state(:, 24))
%legend('Fibrin(ogen) concentration')

%addpath 'iCloud Drive/Desktop/Research'


% hold on
title('Leukocyte Population vs Time in Hours')
xlabel('Time (Hours)')
ylabel('Leukocytes')
set(findall(gcf,'-property','FontSize'),'FontSize',20)
set(gcf, 'Color', 'w');
set(findall(gcf,'-property','Interpreter'),'Interpreter','Latex')
set(findall(gcf,'-property','FontName'),'FontName','CMU Serif')
legend('Bound leukocytes', 'Total leukocytes') %, 'Model 3', 'Complete Model')

export_fig('leuknoFG_result.pdf', '-a1', '-pdf')