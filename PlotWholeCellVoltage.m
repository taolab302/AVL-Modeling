function PlotWholeCellVoltage(strain, time)

% Initialization
I_span = -4:2:12;
if strcmp(strain, 'wt') || strcmp(strain, 'exp-2;shl-1')
    I_star = 10;
elseif strcmp(strain, 'shl-1') || strcmp(strain, 'exp-2(lf)')
    I_star = 8;
elseif strcmp(strain, 'exp-2(gf)')
    I_star = 12;
end


[params, initPos] = AVLParameters(strain, 'Voltage');
now = [nowtime() '_' strain];
RecordNowParams(params, now);
Folder = ['Whole Cell Voltage\' now '\'];
mkdir(Folder);

step = 0.01;
t_span = 0:step:time;
colors = [081 081 081 ;
          241 064 064 ;
          026 111 223;
          129 196 154;
          177 119 222;
          204 153 000;
          000 203 204;
          125 078 078;
          251 101 001;
          ]/255;
      
% Simulation
tic;
h = figure; 
set(h, 'position',[0,0,750,2160]);
for i = 1:length(I_span)
    disp(['Processing: I = ' num2str(I_span(i))]);

    I_EXT = set_constant_current_sequence(length(t_span), I_span(i));
    [t, track] = ode15s(@AVLModel, t_span, initPos, [], I_EXT, step, params);
    V_m = track(:, 1); 
    
    subplot(12,1,1); box off; hold on;
    plot(t/1000, I_EXT, 'color', colors(i,:), 'LineWidth', 1.3);
    axis off; box off;
    
    if I_span(i) == I_star
        i_star = i;
        t_star = t;
        track_star = track;
        save([Folder strain 'I_' num2str(I_span(i)) '.mat'], 'I_EXT', 't', 'V_m');
        continue;
    end
    
    % plot V over t
    subplot(12,1,3:6); box off; hold on;
    plot(t/1000, V_m, 'color', [0.7451 0.7451 0.7451], 'LineWidth', 1.3);
    set(gca,'tickdir', 'out', 'Fontsize', 11, 'Fontname', 'Arial');
    axis([0, time/1000, -80, 60]);
    ylabel('Vm (mV)', 'Fontsize', 12, 'FontName', 'Arial');
    if strcmp(strain, 'wt')
        title('WT', 'Fontsize', 15, 'FontName', 'Arial');
    else
        title(['\it' strain], 'Fontsize', 15, 'FontName', 'Arial');
    end
    box off;
    
    save([Folder strain '_I_' num2str(I_span(i)) '.mat'], 'I_EXT', 't', 'V_m');
end

subplot(12,1,1);
text(0.37, -5, '\it{I_{inj}}', 'Fontsize', 13, 'Fontname', 'Arial');
text(-0.033, 5, '2 pA/step', 'Fontsize', 13, 'Fontname', 'Arial');

subplot(12,1,3:6);
V_m = track_star(:, 1);
plot(t/1000, V_m, 'color', colors(i_star,:), 'LineWidth', 1.3); 

subplot(12,1,8:9); box off; hold on;
% plot Ca2+ conductance over t
gg_UNC2 = track_star(:, 2).^2.* track_star(:, 3);
if params.g_UNC2 == 0
    gg_UNC2 = gg_UNC2 * 0;
end

gg_EGL19 = track_star(:, 4).* track_star(:, 5);
if params.g_EGL19 == 0
    gg_EGL19 = gg_EGL19 * 0;
end

gg_CCA1 = track_star(:, 6).^2.* track_star(:, 7);
if params.g_CCA1 == 0
    gg_CCA1 = gg_CCA1 * 0;
end

plot(t_star/1000, gg_UNC2 , 'k', 'LineWidth', 1.3); 
plot(t_star/1000, gg_EGL19, 'r', 'LineWidth', 1.3); 
plot(t_star/1000, gg_CCA1, 'color', colors(3,:), 'LineWidth', 1.3); 
set(gca,'tickdir', 'out', 'Fontsize', 11, 'Fontname', 'Arial');
axis([0, time/1000, 0, 1]);
legend({'UNC-2', 'EGL-19', 'CCA-1'}, 'Fontsize', 10, 'Box', 'off');
ylabel('G_{Ca} / g_{Ca}', 'Fontsize', 12, 'FontName', 'Arial');
box off; hold off;

subplot(12,1,11:12); box off; hold on;
% plot K+ conductance over t
gg_SHL1 = track_star(:, 8).^3 .* (0.7*track_star(:, 9) + 0.3*track_star(:, 10));
if params.g_SHL1== 0
    gg_SHL1 = gg_SHL1 * 0; 
end

gg_EGL36 = 0.31*track_star(:, 11) + 0.36*track_star(:, 12) + 0.39*track_star(:,13);
if params.g_EGL36 == 0
    gg_EGL36 = gg_EGL36 * 0;  
end

gg_EXP2 = track_star(:, 17);
if params.g_EXP2 == 0
    gg_EXP2 = gg_EXP2 * 0;
end

plot(t_star/1000, gg_SHL1, 'k', 'LineWidth', 1.3); 
plot(t_star/1000, gg_EGL36, 'r', 'LineWidth', 1.3); 
plot(t_star/1000, gg_EXP2, 'color', colors(3,:), 'LineWidth', 1.3);
set(gca,'tickdir', 'out', 'Fontsize', 11, 'Fontname', 'Arial');
axis([0, time/1000, 0, 1]);
legend({'SHL-1', 'EGL-36', 'EXP-2'}, 'Fontsize', 10, 'Box', 'off');
xlabel('Time (s)', 'Fontsize', 13, 'FontName', 'Arial');
ylabel('G_K / g_k', 'Fontsize', 12, 'FontName', 'Arial');
box off; 

save([Folder strain '_I_' num2str(I_star) '_channel.mat'], 'gg_UNC2', 'gg_EGL19', 'gg_CCA1',...
                                                          'gg_SHL1', 'gg_EGL36', 'gg_EXP2');

print(h, [Folder strain '.jpg'], '-djpeg', '-r600');
savefig(h, [Folder strain '.fig']);

timespend = toc;
disp(['Total time cost: ' num2str(timespend) ' s']);
