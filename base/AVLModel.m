function dy = AVLModel(t, y, I, step, params)

% Basic parameters
C     = params.C;
g_L   = params.g_L;
v_Ca  = params.v_Ca;
v_K   = params.v_K;
v_L   = params.v_L;

% NCA-1/2 channel
g_NCA = params.g_NCA;
v_Na  = params.v_Na;

% UNC-2 channel
g_UNC2 = params.g_UNC2;
m_a    = params.m_a;
m_b    = params.m_b;
m_c    = params.m_c;
m_d    = params.m_d;
m_e    = params.m_e;
m_f    = params.m_f;
h_a    = params.h_a;
h_b    = params.h_b;
h_c    = params.h_c;
h_d    = params.h_d;
h_e    = params.h_e;
h_f    = params.h_f;

% EGL-19 channel
g_EGL19 = params.g_EGL19;
s_1     = params.s_1;
s_2     = params.s_2;
s_3     = params.s_3;
s_4     = params.s_4;
s_5     = params.s_5;
s_6     = params.s_6;
s_7     = params.s_7;
s_8     = params.s_8;
s_9     = params.s_9;
s_10    = params.s_10;
s_11    = params.s_11;
s_12    = params.s_12;
s_13    = params.s_13;
s_14    = params.s_14;
s_15    = params.s_15;
q_1     = params.q_1;
q_2     = params.q_2;
q_3     = params.q_3;
q_4     = params.q_4;
q_5     = params.q_5;
q_6     = params.q_6;
q_7     = params.q_7;
q_8     = params.q_8;
q_9     = params.q_9;
q_10    = params.q_10;

% CCA-1 channel
g_CCA1 = params.g_CCA1;
c_1    = params.c_1;
c_2    = params.c_2;
c_3    = params.c_3;
c_4    = params.c_4;
c_5    = params.c_5;
c_6    = params.c_6;
d_1    = params.d_1;
d_2    = params.d_2;
d_3    = params.d_3;
d_4    = params.d_4;
d_5    = params.d_5;
d_6    = params.d_6;

% EGL36 channel
g_EGL36 = params.g_EGL36;
e_1    = params.e_1;
e_2    = params.e_2;
t_f    = params.t_f;
t_m    = params.t_m;
t_s    = params.t_s;

% SHL-1 channel
g_SHL1    = params.g_SHL1;
a_m    = params.a_m;
b_m    = params.b_m;
c_m    = params.c_m;
d_m    = params.d_m;
e_m    = params.e_m;
f_m    = params.f_m;
a_hf   = params.a_hf;
b_hf   = params.b_hf;
c_hf   = params.c_hf;
d_hf   = params.d_hf;
a_hs   = params.a_hs;
b_hs   = params.b_hs;
c_hs   = params.c_hs;
d_hs   = params.d_hs;
v_1    = params.v_1;
v_2    = params.v_2;
v_3    = params.v_3;
v_4    = params.v_4;

% EXP-2 channel
g_EXP2 = params.g_EXP2;
p_1    = params.p_1;
p_2    = params.p_2;
p_3    = params.p_3;
p_4    = params.p_4;
p_5    = params.p_5;
p_6    = params.p_6;
p_7    = params.p_7;
p_8    = params.p_8;
p_9    = params.p_9;
p_10   = params.p_10;
p_11   = params.p_11;
p_12   = params.p_12;
p_13   = params.p_13;
p_14   = params.p_14;
p_15   = params.p_15;
p_16   = params.p_16;

% Variables
V        = y(1);
m_UNC2   = y(2);
h_UNC2   = y(3);
m_EGL19  = y(4);
h_EGL19  = y(5);
m_CCA1   = y(6);
h_CCA1   = y(7);
m_SHL1   = y(8);
hf_SHL1  = y(9);
hs_SHL1  = y(10);
mf_EGL36 = y(11);
mm_EGL36 = y(12);
ms_EGL36 = y(13);
C1_EXP2  = y(14);
C2_EXP2  = y(15);
C3_EXP2  = y(16);
O_EXP2   = y(17);
I_EXP2   = 1 - C1_EXP2 - C2_EXP2 - C3_EXP2 - O_EXP2;

% Whole cell voltage
dvdt = (I(int64(t/step+1))...
        - g_UNC2  * m_UNC2^2 * h_UNC2  * (V-v_Ca)...
        - g_EGL19 * m_EGL19  * h_EGL19 * (V-v_Ca)...
        - g_CCA1  * m_CCA1^2 * h_CCA1  * (V-v_Ca)...
        - g_SHL1  * m_SHL1^3 * (0.7*hf_SHL1 + 0.3*hs_SHL1) * (V-v_K)...    
        - g_EGL36 * (0.31*mf_EGL36 + 0.36*mm_EGL36 + 0.39*ms_EGL36) * (V-v_K)...
        - g_EXP2  * O_EXP2   * (V-v_K)...
        - g_NCA * (V-v_Na)...
        - g_L   * (V-v_L)...
       ) / C;

% UNC-2 channel
m_alpha = m_a * (V-m_b) / (1 - exp(-(V-m_b)/m_c));
m_beta  = m_d * exp(-(V-m_e)/m_f);
h_alpha = h_a * exp(-(V-h_b)/h_c);
h_beta  = h_d / (1 + exp(-(V-h_e)/h_f));
dmUNC2dt = m_alpha * (1 - m_UNC2) - m_beta * m_UNC2;
dhUNC2dt = h_alpha * (1 - h_UNC2) - h_beta * h_UNC2;

% EGL-19 channel
tau_m_EGL19 = s_1 * exp(-((V-s_2)/s_3)^2) + s_4 * exp(-((V-s_5)/s_6)^2) + s_7;
tau_h_EGL19 = s_8 * (s_9 / (1 + exp((V-s_10)/s_11)) + s_12 / (1 + exp((V-s_13)/s_14)) + s_15);
m_EGL19_inf = 1 / (1 + exp(-(V-q_1)/q_2));
h_EGL19_inf = (q_3 / (1 + exp(-(V-q_4)/q_5)) + q_6) * (q_7 / (1 + exp((V-q_8)/q_9)) + q_10);
dmEGL19dt = (m_EGL19_inf - m_EGL19) / tau_m_EGL19;
dhEGL19dt = (h_EGL19_inf - h_EGL19) / tau_h_EGL19;

% CCA-1 channel
m_CCA1_inf = 1 / (1 + exp(-(V-c_1)/c_2));
h_CCA1_inf = 1 / (1 + exp( (V-d_1)/d_2));
tau_m_CCA1 = c_3 / (1 + exp(-(V-c_4)/c_5)) + c_6;
tau_h_CCA1 = d_3 / (1 + exp( (V-d_4)/d_5)) + d_6;
dmCCA1dt   = (m_CCA1_inf - m_CCA1) / tau_m_CCA1;
dhCCA1dt   = (h_CCA1_inf - h_CCA1) / tau_h_CCA1;

% SHL-1 channel
tau_m_SHL1  = a_m / (exp(-(V-b_m)/c_m) + exp((V-d_m)/e_m)) + f_m;
tau_hf_SHL1 = a_hf / (1 + exp((V-b_hf)/c_hf)) + d_hf;
tau_hs_SHL1 = a_hs / (1 + exp((V-b_hs)/c_hs)) + d_hs;
m_SHL1_inf  = 1 / (1 + exp(-(V-v_1)/v_2));
h_SHL1_inf = 1 / (1 + exp( (V-v_3)/v_4));
dmSHL1dt    = (m_SHL1_inf - m_SHL1) / tau_m_SHL1;
dhfSHL1dt   = (h_SHL1_inf - hf_SHL1) / tau_hf_SHL1;
dhsSHL1dt   = (h_SHL1_inf - hs_SHL1) / tau_hs_SHL1;

% EGL-36 channel
m_EGL36_inf = 1 / (1 + exp(-(V-e_1)/e_2));
tau_mf_EGL36 = t_f;
tau_mm_EGL36 = t_m;
tau_ms_EGL36 = t_s;
dmfEGL36dt   = (m_EGL36_inf - mf_EGL36) / tau_mf_EGL36;
dmmEGL36dt   = (m_EGL36_inf - mm_EGL36) / tau_mm_EGL36;
dmsEGL36dt   = (m_EGL36_inf - ms_EGL36) / tau_ms_EGL36;

% EXP-2 channel
alpha_1   = p_1   * exp( p_2  * V);
beta_1    = p_3   * exp(-p_4  * V);
K_f       = p_5;
K_b       = p_6;  
alpha_2   = p_7   * exp( p_8  * V);
beta_2    = p_9   * exp(-p_10 * V);
alpha_i   = p_11  * exp( p_12 * V);
beta_i    = p_13  * exp(-p_14 * V);
alpha_i2  = p_15  * exp( p_16 * V);
psi       = beta_2 * beta_i * alpha_i2 / (alpha_2 * alpha_i);
dC1EXP2dt = beta_1 * C2_EXP2 - alpha_1 * C1_EXP2;
dC2EXP2dt = alpha_1 * C1_EXP2 + K_b * C3_EXP2 - (beta_1 + K_f) * C2_EXP2;
dC3EXP2dt = K_f * C2_EXP2 + psi * I_EXP2 + beta_2 * O_EXP2 - (K_b + alpha_i2 + alpha_2) * C3_EXP2;
dOEXP2dt  = beta_i * I_EXP2 + alpha_2 * C3_EXP2 - (beta_2 + alpha_i) * O_EXP2;

dy = [dvdt;
      dmUNC2dt; 
      dhUNC2dt; 
      dmEGL19dt;
      dhEGL19dt;
      dmCCA1dt;
      dhCCA1dt;
      dmSHL1dt; 
      dhfSHL1dt; 
      dhsSHL1dt; 
      dmfEGL36dt;
      dmmEGL36dt;
      dmsEGL36dt;
      dC1EXP2dt; 
      dC2EXP2dt; 
      dC3EXP2dt;
      dOEXP2dt];
