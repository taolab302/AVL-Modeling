function PlotSpikePeakValue(time, param_name_1, param_span_1, param_name_2, param_span_2, positive_spike_thres)

% Initialization
param_name_strip_1 = strrep(param_name_1, '{', ''); 
param_name_strip_1 = strrep(param_name_strip_1, '}', '');
param_name_latex_1 = ParamNameToLaTex(param_name_1);
param_name_strip_2 = strrep(param_name_2, '{', ''); 
param_name_strip_2 = strrep(param_name_strip_2, '}', '');
param_name_latex_2 = ParamNameToLaTex(param_name_2);

strain = 'exp-2(lf)_rescued';
[params, initPos] = AVLParameters(strain, 'Voltage');
now = [nowtime() '_' strain];
RecordNowParams(params, now);
Folder = ['Spike Peak Value\' now '_' param_name_strip_1 '_' param_name_strip_2 '\']; 
if ~exist(Folder, 'dir')
    mkdir(Folder);
end

step = 0.01;
t_span = 0:step:time;
pos_spike_peak_values = zeros(length(param_span_1), length(param_span_2));
neg_spike_peak_values = zeros(length(param_span_1), length(param_span_2));
pos_neg_spike_diff = zeros(length(param_span_1), length(param_span_2));


% Simulation
tic;
for i = 1:length(param_span_1)
    params = setfield(params, param_name_strip_1, param_span_1(i));
    % apply current
    for j = 1:length(param_span_2)
        disp(['Processing: ' param_name_strip_1 ' = ' num2str(param_span_1(i))...
                        ', ' param_name_strip_2 ' = ' num2str(param_span_2(j))]);
                    
        params = setfield(params, param_name_strip_2, param_span_2(j));
        % apply current
        [I, I_start, I_end] = set_constant_current_sequence(length(t_span), 8.5, 30000);
        [~, track] = ode15s(@AVLModel, t_span, initPos, [], I, step, params);
        
        [pos_spike_peak_values(i, j), pos_spike_peak_idx]  = max(track(I_start:I_end, 1));
        [neg_spike_peak_values(i, j), neg_spike_peak_idx]  = min(track(I_start + pos_spike_peak_idx:I_end, 1));
        if pos_spike_peak_values(i, j) < positive_spike_thres
            pos_neg_spike_diff(i, j) = nan;
            neg_spike_peak_values(i, j) = nan;
            pos_spike_peak_values(i, j) = nan;
        else            
            pos_neg_spike_diff(i, j) = neg_spike_peak_idx / 100;
        end
    end
end
save([Folder param_name_1 '_' param_name_2 '.mat'], 'pos_spike_peak_values',...
                                                    'neg_spike_peak_values',...
                                                    'pos_neg_spike_diff',...
                                                    'param_name_1',...
                                                    'param_span_1',...
                                                    'param_name_2',...
                                                    'param_span_2');
                                                
tick_num = 5;
tickskip = [floor(length(param_span_1)/tick_num) floor(length(param_span_2)/tick_num)];
tickskip(~tickskip) = 1;
ticklabel_1 = cellstr(num2str(param_span_1(:)));
ticklabel_2 = cellstr(num2str(param_span_2(:)));

% positive spike peak
h = figure(100);
imagesc(pos_spike_peak_values); colormap(jet); cb = colorbar();
set(cb, 'ticklabelinterpreter','latex');
set(get(cb,'Title'),'string','$\mathrm{mV}$', 'Interpreter', 'LaTex');
set(gca, 'XDir', 'normal', 'YDir', 'normal');
set(gca,'ticklabelinterpreter','latex','tickdir','out');
set(gca, 'xtick', 1:tickskip(2):length(param_span_2));
set(gca, 'ytick', 1:tickskip(1):length(param_span_1));
set(gca, 'Xticklabel', ticklabel_2(1:tickskip(2):length(param_span_2)));
set(gca, 'Yticklabel', ticklabel_1(1:tickskip(1):length(param_span_1)));
xlabel(['$' param_name_latex_2 '(\mathrm{nS})$'], 'Fontsize', 14, 'Interpreter', 'LaTex');
ylabel(['$' param_name_latex_1 '(\mathrm{nS})$'], 'Fontsize', 14, 'Interpreter', 'LaTex');
title('Amplitude of positive spike', 'Interpreter', 'Latex')
print(h, [Folder 'positive_spike_' param_name_1 '_' param_name_2 '_imagesc.jpg'], '-djpeg', '-r600');

% negative spike peak
h = figure(101);
imagesc(neg_spike_peak_values); colormap(jet); cb = colorbar;
set(cb, 'ticklabelinterpreter','latex');
set(get(cb,'Title'),'string','$\mathrm{mV}$', 'Interpreter', 'LaTex');
set(gca, 'XDir', 'normal', 'YDir', 'normal');
set(gca,'ticklabelinterpreter','latex','tickdir','out');
set(gca, 'xtick', 1:tickskip(2):length(param_span_2));
set(gca, 'ytick', 1:tickskip(1):length(param_span_1));
set(gca, 'Xticklabel', ticklabel_2(1:tickskip(2):length(param_span_2)));
set(gca, 'Yticklabel', ticklabel_1(1:tickskip(1):length(param_span_1)));
xlabel(['$' param_name_latex_2 '(\mathrm{nS})$'], 'Fontsize', 14, 'Interpreter', 'LaTex');
ylabel(['$' param_name_latex_1 '(\mathrm{nS})$'], 'Fontsize', 14, 'Interpreter', 'LaTex');
title('Amplitude of negative spike', 'Interpreter', 'Latex')
print(h, [Folder 'negative_spike_' param_name_1 '_' param_name_2 '_imagesc.jpg'], '-djpeg', '-r600');

% positive negative spike peak time diffrence
h = figure(102);
imagesc(pos_neg_spike_diff); colormap(jet); cb = colorbar;
set(cb, 'ticklabelinterpreter','latex');
set(get(cb,'Title'),'string','$\mathrm{ms}$', 'Interpreter', 'LaTex');
set(gca, 'XDir', 'normal', 'YDir', 'normal');
set(gca,'ticklabelinterpreter','latex','tickdir','out');
set(gca, 'xtick', 1:tickskip(2):length(param_span_2));
set(gca, 'ytick', 1:tickskip(1):length(param_span_1));
set(gca, 'Xticklabel', ticklabel_2(1:tickskip(2):length(param_span_2)));
set(gca, 'Yticklabel', ticklabel_1(1:tickskip(1):length(param_span_1)));
xlabel(['$' param_name_latex_2 '(\mathrm{nS})$'], 'Fontsize', 14, 'Interpreter', 'LaTex');
ylabel(['$' param_name_latex_1 '(\mathrm{nS})$'], 'Fontsize', 14, 'Interpreter', 'LaTex');
title('Time diffrence between positive and negative spike', 'Interpreter', 'Latex')
print(h, [Folder 'pos_neg_spike_diff_' param_name_1 '_' param_name_2 '_imagesc.jpg'], '-djpeg', '-r600');

% isoline
anchor_point = [0.31, 0.5];
idx1 = find(param_span_1 == anchor_point(1));
idx2 = find(param_span_2 == anchor_point(2));
target_neg_spike = neg_spike_peak_values(idx1, idx2);
target_time_diff = pos_neg_spike_diff(idx1, idx2);
delta_neg_spike = abs(neg_spike_peak_values - target_neg_spike);
delta_time_diff = abs(pos_neg_spike_diff - target_time_diff);
[neg_spike_y, neg_spike_x] = find(delta_neg_spike < 0.09);
[time_diff_y, time_diff_x] = find(delta_time_diff < 0.25);
disp(['target_neg_spike = ' num2str(target_neg_spike) ', target_time_diff = ' num2str(target_time_diff)]);
% upper bound
anchor_point_u = [0.31, 0.55];
idx1 = find(param_span_1 == anchor_point_u(1));
idx2 = find(param_span_2 == anchor_point_u(2));
target_neg_spike_u = neg_spike_peak_values(idx1, idx2);
target_time_diff_u = pos_neg_spike_diff(idx1, idx2);
delta_neg_spike_u = abs(neg_spike_peak_values - target_neg_spike_u);
delta_time_diff_u = abs(pos_neg_spike_diff - target_time_diff_u);
[neg_spike_y_u, neg_spike_x_u] = find(delta_neg_spike_u < 0.03);
[time_diff_y_u, time_diff_x_u] = find(delta_time_diff_u < 0.25);
disp(['target_neg_spike_u = ' num2str(target_neg_spike_u) ', target_time_diff_u = ' num2str(target_time_diff_u)]);
% lower bound
anchor_point_l = [0.31, 0.45];
idx1 = find(param_span_1 == anchor_point_l(1));
idx2 = find(param_span_2 == anchor_point_l(2));
target_neg_spike_l = neg_spike_peak_values(idx1, idx2);
target_time_diff_l = pos_neg_spike_diff(idx1, idx2);
delta_neg_spike_l = abs(neg_spike_peak_values - target_neg_spike_l);
delta_time_diff_l = abs(pos_neg_spike_diff - target_time_diff_l);
[neg_spike_y_l, neg_spike_x_l] = find(delta_neg_spike_l < 0.09);
[time_diff_y_l, time_diff_x_l] = find(delta_time_diff_l < 0.25);
disp(['target_neg_spike_l = ' num2str(target_neg_spike_l) ', target_time_diff_l = ' num2str(target_time_diff_l)]);

h = figure(100); hold on
l1 = plot(neg_spike_x, neg_spike_y, 'Color', [026 111 223] / 255, 'LineWidth', 2);
fill([neg_spike_x_u; flipud(neg_spike_x_l)], [neg_spike_y_u; flipud(neg_spike_y_l)],...
     [026 111 223] / 255, 'LineStyle', 'None', 'FaceAlpha', 0.5);
l2 = plot(time_diff_x, time_diff_y, 'Color', [000 203 204] / 255, 'LineWidth', 2);
fill([time_diff_x_u; flipud(time_diff_x_l)], [time_diff_y_u; flipud(time_diff_y_l)],...
     [000 203 204] / 255, 'LineStyle', 'None', 'FaceAlpha', 0.5);
legend([l1, l2], {'isoline of negative spike', 'isoline of time diffrence'}, 'Interpreter', 'LaTex');

timespend = toc;
disp(['Total time cost: ' num2str(timespend) ' s']);
end

function param_name_latex = ParamNameToLaTex(param_name)
prefix_loc = find(param_name == '_');
if ~isempty(prefix_loc)
    if isstrprop(param_name(prefix_loc+2), 'upper')
        param_name = [param_name(1:prefix_loc) '\mathrm' param_name(prefix_loc+1:end)];
    elseif strcmp(param_name(1:3), 'phi')
        param_name = ['\' param_name];
    end
end
param_name_latex = param_name;
end 