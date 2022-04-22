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
    if isfield(params, param_name_strip_1)
        params = setfield(params, param_name_strip_1, param_span_1(i));
    end
    % apply current
    if strcmp(param_name_1, 'I')
        [I, I_start, I_end] = set_constant_current_sequence(length(t_span), param_span_1(i));
    end
    for j = 1:length(param_span_2)
        disp(['Processing: ' param_name_strip_1 ' = ' num2str(param_span_1(i))...
                        ', ' param_name_strip_2 ' = ' num2str(param_span_2(j))]);
                    
        if isfield(params, param_name_strip_2)
            params = setfield(params, param_name_strip_2, param_span_2(j));
        end
        % apply current
        if strcmp(param_name_2, 'I')
            [I, I_start, I_end] = set_constant_current_sequence(length(t_span), param_span(j));
        elseif ~strcmp(param_name_1, 'I')
            [I, I_start, I_end] = set_constant_current_sequence(length(t_span));
        end
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
h = figure;
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
h = figure;
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
h = figure;
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

timespend = toc;
disp(['Total time cost: ' num2str(timespend) ' s']);