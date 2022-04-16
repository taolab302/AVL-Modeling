function PlotTailPeakIV(time)

now = nowtime();
Folder = ['Tail Peak IV\' now '\'];
if ~exist(Folder, 'dir')
    mkdir(Folder);
end
step = 0.01;
t_span = 0:step:time;
V_span = -120:10:70;
peak_radius = 8;
[params_wt, initPos] = AVLParameters('wt', 'Current', now);
params_gf = AVLParameters('exp-2(gf)', 'Current', now);
params_lf = AVLParameters('exp-2(lf)', 'Current', now);

tail_peak_wt     = zeros(1, length(V_span));
tail_peak_exp2lf = zeros(1, length(V_span));
tail_peak_exp2gf = zeros(1, length(V_span));

colors = [000 000 000;
          255 000 000;
          026 111 223;
          ]/255;

tic;
for i = 1:length(V_span)
    disp(['Processing: V = ' num2str(V_span(i))]);
    V = set_constant_voltage_sequence(V_span(i), length(t_span), -30);
    [~, track_wt] = ode15s(@WholeCellCurrent, t_span, initPos, [],  V, step, params_wt);
    [~, track_gf] = ode15s(@WholeCellCurrent, t_span, initPos, [],  V, step, params_gf);
    [~, track_lf] = ode15s(@WholeCellCurrent, t_span, initPos, [],  V, step, params_lf);
    G_UNC2_wt  = params_wt.g_UNC2 .* track_wt(:,1).^2.* track_wt(:,2);
    G_EGL19_wt = params_wt.g_EGL19.* track_wt(:,3)   .* track_wt(:,4);
    G_CCA1_wt  = params_wt.g_CCA1 .* track_wt(:,5).^2.* track_wt(:,6);
    G_SHL1_wt  = params_wt.g_SHL1 .* track_wt(:,7).^3.* (0.7*track_wt(:,8)+0.3*track_wt(:,9));
    G_SHL1_lf  = params_lf.g_SHL1 .* track_lf(:,7).^3.* (0.7*track_lf(:,8)+0.3*track_lf(:,9));
    G_EGL36_wt = params_wt.g_EGL36.* (0.31*track_wt(:,10) + 0.36*track_wt(:,11) + 0.39*track_wt(:,12));
    G_EGL36_lf = params_lf.g_EGL36.* (0.31*track_lf(:,10) + 0.36*track_lf(:,11) + 0.39*track_lf(:,12));
    G_EXP2_wt  = params_wt.g_EXP2 .* track_wt(:,16);
    G_EXP2_gf  = params_gf.g_EXP2 .* track_gf(:,16);
    G_NCA_wt   = params_wt.g_NCA;
    G_L_wt     = params_wt.g_L;
    
    I_UNC2_wt  = G_UNC2_wt   .* (V' - params_wt.v_Ca);
    I_EGL19_wt = G_EGL19_wt  .* (V' - params_wt.v_Ca);
    I_CCA1_wt  = G_CCA1_wt   .* (V' - params_wt.v_Ca);
    I_SHL1_wt  = G_SHL1_wt   .* (V' - params_wt.v_K);
    I_SHL1_lf  = G_SHL1_lf   .* (V' - params_lf.v_K);
    I_EGL36_wt = G_EGL36_wt.* (V' - params_wt.v_K);
    I_EGL36_lf = G_EGL36_lf.* (V' - params_lf.v_K);
    I_EXP2_wt  = G_EXP2_wt .* (V' - params_wt.v_K); 
    I_EXP2_gf  = G_EXP2_gf .* (V' - params_gf.v_K); 
    I_NCA_wt   = G_NCA_wt   * (V' - params_wt.v_Na); 
    I_L_wt     = G_L_wt     * (V' - params_wt.v_L); 
    
    I_wt       = (I_UNC2_wt + I_EGL19_wt + I_CCA1_wt + I_SHL1_wt + I_EGL36_wt + I_NCA_wt + I_L_wt + I_EXP2_wt);
    I_exp2gf   = (I_UNC2_wt + I_EGL19_wt + I_CCA1_wt + I_SHL1_wt + I_EGL36_wt + I_NCA_wt + I_L_wt + I_EXP2_gf);
    I_exp2lf   = (I_UNC2_wt + I_EGL19_wt + I_CCA1_wt + I_SHL1_lf + I_EGL36_lf + I_NCA_wt + I_L_wt);
    
    offset = length(t_span)-10000;
    [~, peak_ind] = max(I_wt(offset+peak_radius+1:end-peak_radius));
    peak_interval = (-peak_radius:peak_radius) + peak_ind + offset;
    tail_peak_wt(i)     = mean(I_wt(peak_interval));
    tail_peak_exp2lf(i) = mean(I_exp2lf(peak_interval));
    tail_peak_exp2gf(i) = mean(I_exp2gf(peak_interval));
end

tail_peak_max = max([tail_peak_wt tail_peak_exp2gf]);
tail_peak_ratio_wt = tail_peak_wt/tail_peak_max;
tail_peak_ratio_exp2lf = tail_peak_exp2lf/max(tail_peak_max);
tail_peak_ratio_exp2gf = tail_peak_exp2gf/tail_peak_max;

wt = tail_peak_ratio_wt(8:end)';
gf = tail_peak_ratio_exp2gf(8:end)';
v  = V_span(8:end)';
startPoint = [0 1 8 -15;
              0 1 6 -10];

Boltzmann_fit = fittype('(A1-A2)/(1+exp((x-x0)/dx))+A2','independent','x','coefficients',{'A1','A2','x0','dx'});
fit_wt     = fit(v, wt, Boltzmann_fit, 'Start', startPoint(1,:));
fit_exp2gf = fit(v, gf, Boltzmann_fit, 'Start', startPoint(2,:));

f = figure; hold on;
h1 = scatter(V_span, tail_peak_ratio_wt,     100, colors(1,:), 's', 'filled');
h2 = scatter(V_span, tail_peak_ratio_exp2lf, 100, colors(2,:), '^', 'filled');
h3 = scatter(V_span, tail_peak_ratio_exp2gf, 70, colors(3,:), 'o', 'filled');

vv = v(1):0.01:v(end);
plot(vv, fit_wt(vv),     'color', colors(1,:), 'LineWidth', 1.3);
plot(vv, fit_exp2gf(vv), 'color', colors(3,:), 'LineWidth', 1.3);

legend([h1, h2, h3], {'WT','\itexp-2 lf','\itexp-2 gf'},...
       'Fontsize', 12, 'FontName', 'Arial', 'Location', 'northwest', 'Box', 'off');
set(gca,'YAxisLocation','origin','XAxisLocation','origin','tickdir','out','Fontsize', 11,'Fontname','Arial');
ylabel('I/I_{max}', 'Fontsize', 13, 'Fontname', 'Arial', 'Rotation', 90,'position', [-15 1.05]);
xlabel('Vm (mV)', 'Fontsize', 13, 'Fontname', 'Arial');
axis([-125 80 -0.1 1.1]);
print(f, [Folder 'IV_simulation.jpg'], '-djpeg', '-r600');
savefig(f, [Folder 'IV_simulation.fig']);
save([Folder 'IV_simulation.mat'], 'V_span', 'tail_peak_wt', 'tail_peak_exp2lf', 'tail_peak_exp2gf',...
                                   'tail_peak_ratio_wt', 'tail_peak_ratio_exp2lf', 'tail_peak_ratio_exp2gf',...
                                   'fit_wt', 'fit_exp2gf');

timespend = toc;
disp(['Total time cost: ' num2str(timespend) ' s']);