function [params, initPos] = AVLParameters(strain, PlotOption, now)

if strcmp(strain, 'wt')
    params = struct('C', 5, 'g_L', 0.27, 'v_Ca', 100, 'v_K', -84, 'v_L', -60,...
                    ...
                    'g_NCA', 0.02, 'v_Na', 30,...
                    ...
                    'g_UNC2', 1,...
                     'm_a', 0.1, 'm_b', 25, 'm_c', 10, 'm_d', 0.4, 'm_e', -25, 'm_f', 18,...
                     'h_a', 0.01, 'h_b', -50, 'h_c', 10, 'h_d', 0.03, 'h_e', -17, 'h_f', 17,...
                    ...
                    'g_EGL19', 0.75,...
                    's_1', 2.9, 's_2', -4.8, 's_3', 6.0,...
                    's_4', 1.9, 's_5', -8.6, 's_6', 30.0, 's_7', 2.3,...
                    's_8', 0.4, 's_9', 44.6, 's_10', -33.0, 's_11', 5.0,...
                    's_12', 36.4, 's_13', 18.7, 's_14', 3.7, 's_15', 43.1,...
                    'q_1', -4.4, 'q_2', 7.5,...
                    'q_3', 1.43, 'q_4', 14.9, 'q_5', 12, 'q_6', 0.14,...
                    'q_7', 5.96, 'q_8',  -20.5, 'q_9', 8.1, 'q_10', 0.6,...
                    ...
                    'g_CCA1',0.25,...
                    'c_1', -43.32, 'c_2', 7.6,...
                    'c_3', 40.0, 'c_4', -62.5, 'c_5', -12.6, 'c_6', 0.7,...
                    'd_1', -58, 'd_2', 7.0,...
                    'd_3', 280, 'd_4', -60.7, 'd_5', 8.5, 'd_6', 19.8,...
                    ...
                    'g_SHL1', 4,...
                    'v_1', -6.8, 'v_2', 14.1,...
                    'a_m', 1.4, 'b_m', -17.5, 'c_m', 12.9, 'd_m', -3.7, 'e_m', 6.5, 'f_m', 0.2,...
                    'v_3', -33.1, 'v_4', 8.3,...
                    'a_hf', 5.9, 'b_hf', -8.2, 'c_hf', 2.9, 'd_hf', 2.73,...   
                    'a_hs', 84.2, 'b_hs', -7.7, 'c_hs', 2.4, 'd_hs', 11.9,...
                    ...
                    'g_EGL36', 1.35,...
                    'e_1', 63, 'e_2', 28.5,...
                    't_f', 13, 't_m', 63, 't_s', 355,...
                    ...
                    'g_EXP2', 3.1,...
                    'p_1', 0.0241, 'p_2', 0.0408, 'p_3', 0.0091, 'p_4', 0.030,...
                    'p_5', 0.0372, 'p_6', 0.31, 'p_7', 0.0376, 'p_8', 0.0472,...
                    'p_9', 0.0015, 'p_10', 0.0703, 'p_11', 0.2177, 'p_12', 0.03,...
                    'p_13', 0.0313, 'p_14', 0.1418, 'p_15', 8.7204e-06, 'p_16', 1.4011e-06);
elseif strcmp(strain, 'exp-2(gf)')
    params = struct('C', 5, 'g_L', 0.27, 'v_Ca', 100, 'v_K', -84, 'v_L', -60,...
                    ...
                    'g_NCA', 0.02, 'v_Na', 30,...
                    ...
                    'g_UNC2', 1,...
                     'm_a', 0.1, 'm_b', 25, 'm_c', 10, 'm_d', 0.4, 'm_e', -25, 'm_f', 18,...
                     'h_a', 0.01, 'h_b', -50, 'h_c', 10, 'h_d', 0.03, 'h_e', -17, 'h_f', 17,...
                    ...
                    'g_EGL19', 0.75,...
                    's_1', 2.9, 's_2', -4.8, 's_3', 6.0,...
                    's_4', 1.9, 's_5', -8.6, 's_6', 30.0, 's_7', 2.3,...
                    's_8', 0.4, 's_9', 44.6, 's_10', -33.0, 's_11', 5.0,...
                    's_12', 36.4, 's_13', 18.7, 's_14', 3.7, 's_15', 43.1,...
                    'q_1', -4.4, 'q_2', 7.5,...
                    'q_3', 1.43, 'q_4', 14.9, 'q_5', 12, 'q_6', 0.14,...
                    'q_7', 5.96, 'q_8',  -20.5, 'q_9', 8.1, 'q_10', 0.6,...
                    ...
                    'g_CCA1',0.25,...
                    'c_1', -43.32, 'c_2', 7.6,...
                    'c_3', 40.0, 'c_4', -62.5, 'c_5', -12.6, 'c_6', 0.7,...
                    'd_1', -58, 'd_2', 7.0,...
                    'd_3', 280, 'd_4', -60.7, 'd_5', 8.5, 'd_6', 19.8,...
                    ...
                    'g_SHL1', 4,...
                    'v_1', -6.8, 'v_2', 14.1,...
                    'a_m', 1.4, 'b_m', -17.5, 'c_m', 12.9, 'd_m', -3.7, 'e_m', 6.5, 'f_m', 0.2,...
                    'v_3', -33.1, 'v_4', 8.3,...
                    'a_hf', 5.9, 'b_hf', -8.2, 'c_hf', 2.9, 'd_hf', 2.73,...   
                    'a_hs', 84.2, 'b_hs', -7.7, 'c_hs', 2.4, 'd_hs', 11.9,...
                    ...
                    'g_EGL36', 1.35,...
                    'e_1', 63, 'e_2', 28.5,...
                    't_f', 13, 't_m', 63, 't_s', 355,...
                    ...
                    'g_EXP2', 3.1,...
                    'p_1', 0.0097, 'p_2', 0.0367, 'p_3', 0.002, 'p_4', 8.200e-04,...
                    'p_5', 0.0479, 'p_6', 0.31, 'p_7', 0.0811, 'p_8',  0.0367,...
                    'p_9', 0.0012, 'p_10', 0.0511, 'p_11', 0.1182, 'p_12', 0.006,...
                    'p_13', 0.0313, 'p_14', 0.1418, 'p_15', 8.7204e-06, 'p_16', 1.4011e-06);
elseif strcmp(strain, 'exp-2(lf)')
    params = struct('C', 5, 'g_L', 0.27, 'v_Ca', 100, 'v_K', -84, 'v_L', -60,...
                    ...
                    'g_NCA', 0.02, 'v_Na', 30,...
                    ...
                    'g_UNC2', 1,...
                     'm_a', 0.1, 'm_b', 25, 'm_c', 10, 'm_d', 0.4, 'm_e', -25, 'm_f', 18,...
                     'h_a', 0.01, 'h_b', -50, 'h_c', 10, 'h_d', 0.03, 'h_e', -17, 'h_f', 17,...
                    ...
                    'g_EGL19', 0.75,...
                    's_1', 2.9, 's_2', -4.8, 's_3', 6.0,...
                    's_4', 1.9, 's_5', -8.6, 's_6', 30.0, 's_7', 2.3,...
                    's_8', 0.4, 's_9', 44.6, 's_10', -33.0, 's_11', 5.0,...
                    's_12', 36.4, 's_13', 18.7, 's_14', 3.7, 's_15', 43.1,...
                    'q_1', -4.4, 'q_2', 7.5,...
                    'q_3', 1.43, 'q_4', 14.9, 'q_5', 12, 'q_6', 0.14,...
                    'q_7', 5.96, 'q_8',  -20.5, 'q_9', 8.1, 'q_10', 0.6,...
                    ...
                    'g_CCA1',0.25,...
                    'c_1', -43.32, 'c_2', 7.6,...
                    'c_3', 40.0, 'c_4', -62.5, 'c_5', -12.6, 'c_6', 0.7,...
                    'd_1', -58, 'd_2', 7.0,...
                    'd_3', 280, 'd_4', -60.7, 'd_5', 8.5, 'd_6', 19.8,...
                    ...
                    'g_SHL1', 4*1.35,...
                    'v_1', -6.8, 'v_2', 14.1,...
                    'a_m', 1.4, 'b_m', -17.5, 'c_m', 12.9, 'd_m', -3.7, 'e_m', 6.5, 'f_m', 0.2,...
                    'v_3', -33.1, 'v_4', 8.3,...
                    'a_hf', 5.9, 'b_hf', -8.2, 'c_hf', 2.9, 'd_hf', 2.73,...   
                    'a_hs', 84.2, 'b_hs', -7.7, 'c_hs', 2.4, 'd_hs', 11.9,...
                    ...
                    'g_EGL36', 1.35*1.35,...
                    'e_1', 63, 'e_2', 28.5,...
                    't_f', 13, 't_m', 63, 't_s', 355,...
                    ...
                    'g_EXP2', 0,...
                    'p_1', 0.0241, 'p_2', 0.0408, 'p_3', 0.0091, 'p_4', 0.030,...
                    'p_5', 0.0372, 'p_6', 0.31, 'p_7', 0.0376, 'p_8', 0.0472,...
                    'p_9', 0.0015, 'p_10', 0.0703, 'p_11', 0.2177, 'p_12', 0.03,...
                    'p_13', 0.0313, 'p_14', 0.1418, 'p_15', 8.7204e-06, 'p_16', 1.4011e-06);
elseif strcmp(strain, 'shl-1')
    params = struct('C', 5, 'g_L', 0.27, 'v_Ca', 100, 'v_K', -84, 'v_L', -60,...
                    ...
                    'g_NCA', 0.02, 'v_Na', 30,...
                    ...
                    'g_UNC2', 1,...
                     'm_a', 0.1, 'm_b', 25, 'm_c', 10, 'm_d', 0.4, 'm_e', -25, 'm_f', 18,...
                     'h_a', 0.01, 'h_b', -50, 'h_c', 10, 'h_d', 0.03, 'h_e', -17, 'h_f', 17,...
                    ...
                    'g_EGL19', 0.75,...
                    's_1', 2.9, 's_2', -4.8, 's_3', 6.0,...
                    's_4', 1.9, 's_5', -8.6, 's_6', 30.0, 's_7', 2.3,...
                    's_8', 0.4, 's_9', 44.6, 's_10', -33.0, 's_11', 5.0,...
                    's_12', 36.4, 's_13', 18.7, 's_14', 3.7, 's_15', 43.1,...
                    'q_1', -4.4, 'q_2', 7.5,...
                    'q_3', 1.43, 'q_4', 14.9, 'q_5', 12, 'q_6', 0.14,...
                    'q_7', 5.96, 'q_8',  -20.5, 'q_9', 8.1, 'q_10', 0.6,...
                    ...
                    'g_CCA1',0.25,...
                    'c_1', -43.32, 'c_2', 7.6,...
                    'c_3', 40.0, 'c_4', -62.5, 'c_5', -12.6, 'c_6', 0.7,...
                    'd_1', -58, 'd_2', 7.0,...
                    'd_3', 280, 'd_4', -60.7, 'd_5', 8.5, 'd_6', 19.8,...
                    ...
                    'g_SHL1', 0,...
                    'v_1', -6.8, 'v_2', 14.1,...
                    'a_m', 1.4, 'b_m', -17.5, 'c_m', 12.9, 'd_m', -3.7, 'e_m', 6.5, 'f_m', 0.2,...
                    'v_3', -33.1, 'v_4', 8.3,...
                    'a_hf', 5.9, 'b_hf', -8.2, 'c_hf', 2.9, 'd_hf', 2.73,...   
                    'a_hs', 84.2, 'b_hs', -7.7, 'c_hs', 2.4, 'd_hs', 11.9,...
                    ...
                    'g_EGL36', 1.35*1.5,...
                    'e_1', 63, 'e_2', 28.5,...
                    't_f', 13, 't_m', 63, 't_s', 355,...
                    ...
                    'g_EXP2', 3.1*1.5,...
                    'p_1', 0.0241, 'p_2', 0.0408, 'p_3', 0.0091, 'p_4', 0.030,...
                    'p_5', 0.0372, 'p_6', 0.31, 'p_7', 0.0376, 'p_8', 0.0472,...
                    'p_9', 0.0015, 'p_10', 0.0703, 'p_11', 0.2177, 'p_12', 0.03,...
                    'p_13', 0.0313, 'p_14', 0.1418, 'p_15', 8.7204e-06, 'p_16', 1.4011e-06);
elseif strcmp(strain, 'exp-2;shl-1')
    params = struct('C', 5, 'g_L', 0.27, 'v_Ca', 100, 'v_K', -84, 'v_L', -60,...
                    ...
                    'g_NCA', 0.02, 'v_Na', 30,...
                    ...
                    'g_UNC2', 1,...
                     'm_a', 0.1, 'm_b', 25, 'm_c', 10, 'm_d', 0.4, 'm_e', -25, 'm_f', 18,...
                     'h_a', 0.01, 'h_b', -50, 'h_c', 10, 'h_d', 0.03, 'h_e', -17, 'h_f', 17,...
                    ...
                    'g_EGL19', 0.75,...
                    's_1', 2.9, 's_2', -4.8, 's_3', 6.0,...
                    's_4', 1.9, 's_5', -8.6, 's_6', 30.0, 's_7', 2.3,...
                    's_8', 0.4, 's_9', 44.6, 's_10', -33.0, 's_11', 5.0,...
                    's_12', 36.4, 's_13', 18.7, 's_14', 3.7, 's_15', 43.1,...
                    'q_1', -4.4, 'q_2', 7.5,...
                    'q_3', 1.43, 'q_4', 14.9, 'q_5', 12, 'q_6', 0.14,...
                    'q_7', 5.96, 'q_8',  -20.5, 'q_9', 8.1, 'q_10', 0.6,...
                    ...
                    'g_CCA1',0.25,...
                    'c_1', -43.32, 'c_2', 7.6,...
                    'c_3', 40.0, 'c_4', -62.5, 'c_5', -12.6, 'c_6', 0.7,...
                    'd_1', -58, 'd_2', 7.0,...
                    'd_3', 280, 'd_4', -60.7, 'd_5', 8.5, 'd_6', 19.8,...
                    ...
                    'g_SHL1', 0,...
                    'v_1', -6.8, 'v_2', 14.1,...
                    'a_m', 1.4, 'b_m', -17.5, 'c_m', 12.9, 'd_m', -3.7, 'e_m', 6.5, 'f_m', 0.2,...
                    'v_3', -33.1, 'v_4', 8.3,...
                    'a_hf', 5.9, 'b_hf', -8.2, 'c_hf', 2.9, 'd_hf', 2.73,...   
                    'a_hs', 84.2, 'b_hs', -7.7, 'c_hs', 2.4, 'd_hs', 11.9,...
                    ...
                    'g_EGL36', 1.35*2,...
                    'e_1', 63, 'e_2', 28.5,...
                    't_f', 13, 't_m', 63, 't_s', 355,...
                    ...
                    'g_EXP2', 0,...
                    'p_1', 0.0241, 'p_2', 0.0408, 'p_3', 0.0091, 'p_4', 0.0345,...
                    'p_5', 0.0372, 'p_6', 0.31, 'p_7', 0.0376, 'p_8', 0.0472,...
                    'p_9', 0.0011, 'p_10', 0.0703, 'p_11', 0.2177, 'p_12', 0.03,...
                    'p_13', 0.0313, 'p_14', 0.1418, 'p_15', 8.7204e-06, 'p_16', 1.4011e-06);
end
            
if strcmp(PlotOption, 'Voltage')
    initPos = [-60,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0];
elseif strcmp(PlotOption, 'Current')
    initPos = [0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0];
else
    error('Command Error!');
end

if nargin == 3
    RecordNowParams(params, [now '_' strain]);
end
