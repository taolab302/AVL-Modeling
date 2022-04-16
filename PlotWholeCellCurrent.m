function PlotWholeCellCurrent(time)

V_span = -120:10:70;

% make folder
now = nowtime();
Folder = ['Whole Cell Current\' now '\'];
if ~exist(Folder, 'dir')
    mkdir(Folder);
end

[params_wt, initPos] = AVLParameters('wt', 'Current', now);
params_exp2lf = AVLParameters('exp-2(lf)', 'Current', now);
params_exp2gf = AVLParameters('exp-2(gf)', 'Current', now);
step = 0.01;
t_span = 0:step:time;

colors = [000 000 000;
          127 112 135;
          176 147 197;
          116 074 158;
          184 163 153;
          241 227 134;
          193 123 062;
          172 138 118;
          158 201 224;
          073 147 195;
          077 152 166
          170 212 139;
          104 186 083;
          093 158 067;
          219 154 135;
          240 094 094;
          228 045 038;
          247 148 090;
          250 167 072;
          245 127 032;
          ]/255;
      
V_data = zeros(length(t_span), length(V_span));
I_data_wt = zeros(length(t_span), length(V_span));
I_UNC2_data_wt = zeros(length(t_span), length(V_span));
I_EGL19_data_wt = zeros(length(t_span), length(V_span));
I_CCA1_data_wt = zeros(length(t_span), length(V_span));
I_SHL1_data_wt = zeros(length(t_span), length(V_span));
I_EGL36_data_wt = zeros(length(t_span), length(V_span));
I_EXP2_data_wt = zeros(length(t_span), length(V_span));
I_data_exp2lf = zeros(length(t_span), length(V_span));
I_UNC2_data_exp2lf = zeros(length(t_span), length(V_span));
I_EGL19_data_exp2lf = zeros(length(t_span), length(V_span));
I_CCA1_data_exp2lf = zeros(length(t_span), length(V_span));
I_SHL1_data_exp2lf = zeros(length(t_span), length(V_span));
I_EGL36_data_exp2lf = zeros(length(t_span), length(V_span));
I_EXP2_data_exp2lf = zeros(length(t_span), length(V_span));
I_data_exp2gf = zeros(length(t_span), length(V_span));
I_UNC2_data_exp2gf = zeros(length(t_span), length(V_span));
I_EGL19_data_exp2gf = zeros(length(t_span), length(V_span));
I_CCA1_data_exp2gf = zeros(length(t_span), length(V_span));
I_SHL1_data_exp2gf = zeros(length(t_span), length(V_span));
I_EGL36_data_exp2gf = zeros(length(t_span), length(V_span));
I_EXP2_data_exp2gf = zeros(length(t_span), length(V_span));

tic;
h = figure;
set(h, 'position',[0,0,1500,2000]);
for i = 1:length(V_span)
    disp(['Processing: V = ' num2str(V_span(i))]);
    V = set_constant_voltage_sequence(V_span(i), length(t_span), -60);
    V_data(:, i) = V;
    
    % WT Voltage Clamp
    [t, track] = ode15s(@WholeCellCurrent, t_span, initPos, [],  V, step, params_wt);
    G_UNC2  = params_wt.g_UNC2 .* track(:, 1).^2.* track(:, 2);
    G_EGL19 = params_wt.g_EGL19.* track(:, 3)   .* track(:, 4);
    G_CCA1  = params_wt.g_CCA1 .* track(:, 5).^2.* track(:, 6);
    G_SHL1  = params_wt.g_SHL1 .* track(:, 7).^3.* (0.7.*track(:, 8)+0.3.*track(:, 9));
    G_EGL36 = params_wt.g_EGL36.* (0.31*track(:, 10) + 0.36*track(:, 11) + 0.39*track(:, 12));
    G_EXP2  = params_wt.g_EXP2 .* track(:, 16);
    G_NCA   = params_wt.g_NCA;
    G_L     = params_wt.g_L;
    
    I_UNC2  = G_UNC2 .* (V' - params_wt.v_Ca);
    I_EGL19 = G_EGL19.* (V' - params_wt.v_Ca);
    I_CCA1  = G_CCA1 .* (V' - params_wt.v_Ca);
    I_SHL1  = G_SHL1 .* (V' - params_wt.v_K);
    I_EGL36 = G_EGL36.* (V' - params_wt.v_K);
    I_EXP2  = G_EXP2 .* (V' - params_wt.v_K); 
    I_NCA   = G_NCA   * (V' - params_wt.v_Na); 
    I_L     = G_L     * (V' - params_wt.v_L); 
    I       = (I_UNC2 + I_EGL19 + I_CCA1 + I_EXP2 + I_SHL1 + I_EGL36 + I_NCA + I_L); 
    
    I_data_wt(:,i) = I;
    I_UNC2_data_wt(:,i) = I_UNC2;
    I_EGL19_data_wt(:,i) = I_EGL19;
    I_CCA1_data_wt(:,i) = I_CCA1;
    I_SHL1_data_wt(:,i) = I_SHL1;
    I_EGL36_data_wt(:,i) = I_EGL36;
    I_EXP2_data_wt(:,i) = I_EXP2;
    
    subplot(131);hold on;
    plot(t, I, 'color', colors(i,:), 'LineWidth', 1.05);
    axis([0 650 -200 1500]);
    title('WT', 'Fontsize', 15, 'FontName', 'Arial');
    axis off;

    % exp-2(lf) Voltage Clamp
    [t, track] = ode15s(@WholeCellCurrent, t_span, initPos, [],  V, step, params_exp2lf);
    G_UNC2  = params_exp2lf.g_UNC2 .* track(:, 1).^2.* track(:, 2);
    G_EGL19 = params_exp2lf.g_EGL19.* track(:, 3)   .* track(:, 4);
    G_CCA1  = params_exp2lf.g_CCA1 .* track(:, 5).^2.* track(:, 6);
    G_SHL1  = params_exp2lf.g_SHL1 .* track(:, 7).^3.* (0.7*track(:, 8)+0.3*track(:, 9));
    G_EGL36 = params_exp2lf.g_EGL36.* (0.31*track(:, 10) + 0.36*track(:, 11) + 0.39*track(:, 12));
    G_EXP2  = params_exp2lf.g_EXP2 .* track(:, 16);
    G_NCA   = params_exp2lf.g_NCA;
    G_L     = params_exp2lf.g_L;
    
    I_UNC2  = G_UNC2 .* (V' - params_exp2lf.v_Ca);
    I_EGL19 = G_EGL19.* (V' - params_exp2lf.v_Ca);
    I_CCA1  = G_CCA1 .* (V' - params_exp2lf.v_Ca);
    I_SHL1  = G_SHL1 .* (V' - params_exp2lf.v_K);
    I_EGL36 = G_EGL36.* (V' - params_exp2lf.v_K);
    I_EXP2  = G_EXP2 .* (V' - params_exp2lf.v_K); 
    I_NCA   = G_NCA   * (V' - params_exp2lf.v_Na); 
    I_L     = G_L     * (V' - params_exp2lf.v_L); 
    I       = (I_UNC2 + I_EGL19 + I_CCA1 + I_EXP2 + I_SHL1 + I_EGL36 + I_NCA + I_L); 
    
    I_data_exp2lf(:,i) = I;
    I_UNC2_data_exp2lf(:,i) = I_UNC2;
    I_EGL19_data_exp2lf(:,i) = I_EGL19;
    I_CCA1_data_exp2lf(:,i) = I_CCA1;
    I_SHL1_data_exp2lf(:,i) = I_SHL1;
    I_EGL36_data_exp2lf(:,i) = I_EGL36;
    I_EXP2_data_exp2lf(:,i) = I_EXP2;
    
    subplot(132);hold on;
    plot(t, I, 'color', colors(i,:), 'LineWidth', 1.05);
    axis([0 650 -200 1500]);
    title('\itexp-2 lf', 'Fontsize', 15, 'FontName', 'Arial');
    axis off;
    
     % exp-2(gf) Voltage Clamp
    [t, track] = ode15s(@WholeCellCurrent, t_span, initPos, [],  V, step, params_exp2gf);
    G_UNC2  = params_exp2gf.g_UNC2 .* track(:, 1).^2.* track(:, 2);
    G_EGL19 = params_exp2gf.g_EGL19.* track(:, 3)   .* track(:, 4);
    G_CCA1  = params_exp2gf.g_CCA1 .* track(:, 5).^2.* track(:, 6);
    G_SHL1  = params_exp2gf.g_SHL1 .* track(:, 7).^3.* (0.7*track(:, 8)+0.3*track(:, 9));
    G_EGL36 = params_exp2gf.g_EGL36.* (0.31*track(:, 10) + 0.36*track(:, 11) + 0.39*track(:, 12));
    G_EXP2  = params_exp2gf.g_EXP2 .* track(:, 16);
    G_NCA   = params_exp2gf.g_NCA;
    G_L     = params_exp2gf.g_L;
    
    I_UNC2  = G_UNC2 .* (V' - params_exp2gf.v_Ca);
    I_EGL19 = G_EGL19.* (V' - params_exp2gf.v_Ca);
    I_CCA1  = G_CCA1 .* (V' - params_exp2gf.v_Ca);
    I_SHL1  = G_SHL1 .* (V' - params_exp2gf.v_K);
    I_EGL36 = G_EGL36.* (V' - params_exp2gf.v_K);
    I_EXP2  = G_EXP2 .* (V' - params_exp2gf.v_K); 
    I_NCA   = G_NCA   * (V' - params_exp2gf.v_Na); 
    I_L     = G_L     * (V' - params_exp2gf.v_L); 
    I       = (I_UNC2 + I_EGL19 + I_CCA1 + I_EXP2 + I_SHL1 + I_EGL36 + I_NCA + I_L); 
    
    I_data_exp2gf(:,i) = I;
    I_UNC2_data_exp2gf(:,i) = I_UNC2;
    I_EGL19_data_exp2gf(:,i) = I_EGL19;
    I_CCA1_data_exp2gf(:,i) = I_CCA1;
    I_SHL1_data_exp2gf(:,i) = I_SHL1;
    I_EGL36_data_exp2gf(:,i) = I_EGL36;
    I_EXP2_data_exp2gf(:,i) = I_EXP2;
    
    subplot(133);hold on;
    plot(t, I, 'color', colors(i,:), 'LineWidth', 1.05);
    axis([0 650 -200 1500]);
    title('\itexp-2 gf', 'Fontsize', 15, 'FontName', 'Arial');
    axis off;
end

% scale
scale_xy = [500, 600];
subplot(131);
plot([scale_xy(1) scale_xy(1)+100], [scale_xy(2) scale_xy(2)], 'k', 'LineWidth', 1.05);
plot([scale_xy(1) scale_xy(1)], [scale_xy(2) scale_xy(2)+100], 'k', 'LineWidth', 1.05);
text(scale_xy(1)-10, scale_xy(2)-40, '100 ms', 'Fontsize', 10, 'FontName', 'Arial');
text(scale_xy(1)-130, scale_xy(2)+50, '100 pA', 'Fontsize', 10, 'FontName', 'Arial');
subplot(132);
plot([scale_xy(1) scale_xy(1)+100], [scale_xy(2) scale_xy(2)], 'k', 'LineWidth', 1.05);
plot([scale_xy(1) scale_xy(1)], [scale_xy(2) scale_xy(2)+100], 'k', 'LineWidth', 1.05);
text(scale_xy(1)-10, scale_xy(2)-40, '100 ms', 'Fontsize', 10, 'FontName', 'Arial');
text(scale_xy(1)-130, scale_xy(2)+50, '100 pA', 'Fontsize', 10, 'FontName', 'Arial');
subplot(133);
plot([scale_xy(1) scale_xy(1)+100], [scale_xy(2) scale_xy(2)], 'k', 'LineWidth', 1.05);
plot([scale_xy(1) scale_xy(1)], [scale_xy(2) scale_xy(2)+100], 'k', 'LineWidth', 1.05);
text(scale_xy(1)-10, scale_xy(2)-40, '100 ms', 'Fontsize', 10, 'FontName', 'Arial');
text(scale_xy(1)-130, scale_xy(2)+50, '100 pA', 'Fontsize', 10, 'FontName', 'Arial');



V          = V_data;
I_wt       = I_data_wt;
I_UNC2_wt  = I_UNC2_data_wt;
I_EGL19_wt = I_EGL19_data_wt;
I_CCA1_wt  = I_CCA1_data_wt;
I_SHL1_wt  = I_SHL1_data_wt;
I_EGL36_wt = I_EGL36_data_wt;
I_EXP2_wt  = I_EXP2_data_wt;

I_exp2lf       = I_data_exp2lf;
I_UNC2_exp2lf  = I_UNC2_data_exp2lf;
I_EGL19_exp2lf = I_EGL19_data_exp2lf;
I_CCA1_exp2lf  = I_CCA1_data_exp2lf;
I_SHL1_exp2lf  = I_SHL1_data_exp2lf;
I_EGL36_exp2lf = I_EGL36_data_exp2lf;
I_EXP2_exp2lf  = I_EXP2_data_exp2lf;

I_exp2gf       = I_data_exp2gf;
I_UNC2_exp2gf  = I_UNC2_data_exp2gf;
I_EGL19_exp2gf = I_EGL19_data_exp2gf;
I_CCA1_exp2gf  = I_CCA1_data_exp2gf;
I_SHL1_exp2gf  = I_SHL1_data_exp2gf;
I_EGL36_exp2gf = I_EGL36_data_exp2gf;
I_EXP2_exp2gf  = I_EXP2_data_exp2gf;

wt_Folder = [Folder 'wt\'];
mkdir(wt_Folder);
save([wt_Folder 'I_wt.mat'], 't_span', 'V', 'I_wt');
save([wt_Folder 'I_UNC2_wt.mat'], 't_span', 'V', 'I_UNC2_wt');
save([wt_Folder 'I_EGL19_wt.mat'], 't_span', 'V', 'I_EGL19_wt');
save([wt_Folder 'I_CCA1_wt.mat'], 't_span', 'V', 'I_CCA1_wt');
save([wt_Folder 'I_SHL1_wt.mat'], 't_span', 'V', 'I_SHL1_wt');
save([wt_Folder 'I_EGL36_wt.mat'], 't_span', 'V', 'I_EGL36_wt');
save([wt_Folder 'I_EXP2_wt.mat'], 't_span', 'V', 'I_EXP2_wt');

exp2lf_Folder = [Folder 'exp2lf\'];
mkdir(exp2lf_Folder);
save([exp2lf_Folder 'I_exp2lf.mat'], 't_span', 'V', 'I_exp2lf');
save([exp2lf_Folder 'I_UNC2_exp2lf.mat'], 't_span', 'V', 'I_UNC2_exp2lf');
save([exp2lf_Folder 'I_EGL19_exp2lf.mat'], 't_span', 'V', 'I_EGL19_exp2lf');
save([exp2lf_Folder 'I_CCA1_exp2lf.mat'], 't_span', 'V', 'I_CCA1_exp2lf');
save([exp2lf_Folder 'I_SHL1_exp2lf.mat'], 't_span', 'V', 'I_SHL1_exp2lf');
save([exp2lf_Folder 'I_EGL36_exp2lf.mat'], 't_span', 'V', 'I_EGL36_exp2lf');
save([exp2lf_Folder 'I_EXP2_exp2lf.mat'], 't_span', 'V', 'I_EXP2_exp2lf');

exp2gf_Folder = [Folder 'exp2gf\'];
mkdir(exp2gf_Folder);
save([exp2gf_Folder 'I_exp2gf.mat'], 't_span', 'V', 'I_exp2gf');
save([exp2gf_Folder 'I_UNC2_exp2gf.mat'], 't_span', 'V', 'I_UNC2_exp2gf');
save([exp2gf_Folder 'I_EGL19_exp2gf.mat'], 't_span', 'V', 'I_EGL19_exp2gf');
save([exp2gf_Folder 'I_CCA1_exp2gf.mat'], 't_span', 'V', 'I_CCA1_exp2gf');
save([exp2gf_Folder 'I_SHL1_exp2gf.mat'], 't_span', 'V', 'I_SHL1_exp2gf');
save([exp2gf_Folder 'I_EGL36_exp2gf.mat'], 't_span', 'V', 'I_EGL36_exp2gf');
save([exp2gf_Folder 'I_EXP2_exp2gf.mat'], 't_span', 'V', 'I_EXP2_exp2gf');

print(h, [Folder 'V.jpg'], '-djpeg', '-r600');
savefig(h, [Folder 'V.fig']);

timespend = toc;
disp(['Total time cost: ' num2str(timespend) ' s']);
end
