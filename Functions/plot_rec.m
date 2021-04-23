function plot_rec(T_vec, n_B, R_nah, steps, error_times, color)
% Calculate mean ans std of recovered state
std_R_nah = std(R_nah);
avg_R_nah = mean(R_nah);

% Population dynamics when no rewiring
col = color; % 127,255,0 0 100/255 0 Old green
plot(T_vec(steps), avg_R_nah(steps), 'color', col, 'linewidth', 2.5)
hold on
errorbar(T_vec(error_times), avg_R_nah(error_times), std_R_nah(error_times), 'LineStyle','none','LineWidth',2.5, 'color', col,'HandleVisibility','off')
xlabel('Days')
ylabel('Recovered people')
xlim([0 T_vec(end)])
ylim([0 n_B])
set(gca, 'fontsize', 14)
%legend(['R ' label_rewiring ], 'location', 'northeast')
%legend box off
end