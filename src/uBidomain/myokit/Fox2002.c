/*
fox_2002
Generated on 2022-10-19 13:22:05
*/


#define N_STATE 13

/* Declare intermediary, temporary and system variables */
static realtype AC_C_sc;
static realtype AC_P_Ca;
static realtype AC_P_CaK;
static realtype AV_i_Ca;
static realtype AV_i_CaK;
static realtype AC_i_Ca_half;
static realtype AV_i_Ca_max;
static realtype AV_L_type_Ca_current_d_gate_E0_m;
static realtype AV_d_infinity;
static realtype AV_tau_d;
static realtype AC_K_mfCa;
static realtype AV_f_Ca_infinity;
static realtype AC_tau_f_Ca;
static realtype AV_f_infinity;
static realtype AV_tau_f;
static realtype AC_K_NaCa;
static realtype AC_K_mCa;
static realtype AC_K_mNa;
static realtype AC_K_sat;
static realtype AC_eta;
static realtype AV_i_NaCa;
static realtype AV_E_Ca;
static realtype AC_g_Cab;
static realtype AV_i_Ca_b;
static realtype AC_A_Cap;
static realtype AC_CMDN_tot;
static realtype AC_CSQN_tot;
static realtype AV_J_leak;
static realtype AV_J_rel;
static realtype AV_J_up;
static realtype AC_K_mCMDN;
static realtype AC_K_mCSQN;
static realtype AC_K_mup;
static realtype AC_P_leak;
static realtype AC_P_rel;
static realtype AC_V_SR;
static realtype AC_V_myo;
static realtype AC_V_up;
static realtype AV_beta_SR;
static realtype AV_beta_i;
static realtype AV_calcium_dynamics_gamma;
static realtype AV_time;
static realtype AC_E_Na;
static realtype AC_g_Na;
static realtype AV_i_Na;
static realtype AV_alpha_h;
static realtype AV_beta_h;
static realtype AC_shift_h;
static realtype AV_alpha_j;
static realtype AV_beta_j;
static realtype AC_shift_j;
static realtype AV_fast_sodium_current_m_gate_E0_m;
static realtype AV_alpha_m;
static realtype AV_beta_m;
static realtype AC_F;
static realtype AC_R;
static realtype AC_T;
static realtype AV_i_Stim;
static realtype AV_pace;
static realtype AC_stim_amplitude;
static realtype AC_stim_end;
static realtype AC_g_Kp;
static realtype AV_i_Kp;
static realtype AV_Kp_V;
static realtype AC_E_K;
static realtype AV_R_V;
static realtype AC_g_Kr;
static realtype AV_i_Kr;
static realtype AV_X_kr_inf;
static realtype AV_tau_X_kr;
static realtype AC_K_mpCa;
static realtype AC_i_pCa_max;
static realtype AV_i_p_Ca;
static realtype AC_E_Ks;
static realtype AC_g_Ks;
static realtype AV_i_Ks;
static realtype AV_X_ks_infinity;
static realtype AV_tau_X_ks;
static realtype AC_g_Nab;
static realtype AV_i_Na_b;
static realtype AC_K_mKo;
static realtype AC_K_mNai;
static realtype AV_f_NaK;
static realtype AV_i_NaK;
static realtype AC_i_NaK_max;
static realtype AC_sigma;
static realtype AC_Ca_o;
static realtype AC_K_i;
static realtype AC_K_o;
static realtype AC_Na_i;
static realtype AC_Na_o;
static realtype AC_K_mK1;
static realtype AC_g_K1;
static realtype AV_i_K1;
static realtype AV_K1_infinity;
static realtype AC_g_to;
static realtype AV_i_to;
static realtype AV_alpha_X_to;
static realtype AV_beta_X_to;
static realtype AV_alpha_Y_to;
static realtype AV_beta_Y_to;

/* Set values of constants */
static void
updateConstants(void)
{
    /* L_type_Ca_current_f_Ca_gate */
    AC_K_mfCa = 0.18;
    AC_tau_f_Ca = 30.0;
    
    /* fast_sodium_current_h_gate */
    AC_shift_h = 0.0;
    
    /* fast_sodium_current_j_gate */
    AC_shift_j = 0.0;
    
    /* sarcolemmal_calcium_pump */
    AC_K_mpCa = 0.05;
    AC_i_pCa_max = 0.05;
    
    /* standard_ionic_concentrations */
    AC_Ca_o = 2000.0;
    AC_K_i = 149.4;
    AC_K_o = 4.0;
    AC_Na_i = 10.0;
    AC_Na_o = 138.0;
    
    /* L_type_Ca_current */
    AC_C_sc = 1.0;
    AC_P_Ca = 2.26e-05;
    AC_P_CaK = 5.79e-07;
    AC_i_Ca_half = (-0.265);
    
    /* Na_Ca_exchanger */
    AC_K_NaCa = 1500.0;
    AC_K_mCa = 1380.0;
    AC_K_mNa = 87.5;
    AC_K_sat = 0.2;
    AC_eta = 0.35;
    
    /* calcium_background_current */
    AC_g_Cab = 0.0003842;
    
    /* calcium_dynamics */
    AC_A_Cap = 0.0001534;
    AC_CMDN_tot = 10.0;
    AC_CSQN_tot = 10000.0;
    AC_K_mCMDN = 2.0;
    AC_K_mCSQN = 600.0;
    AC_K_mup = 0.32;
    AC_P_leak = 1e-06;
    AC_P_rel = 6.0;
    AC_V_SR = 2e-06;
    AC_V_myo = 2.584e-05;
    AC_V_up = 0.1;
    
    /* fast_sodium_current */
    AC_g_Na = 12.8;
    
    /* membrane */
    AC_F = 96.5;
    AC_R = 8.314;
    AC_T = 310.0;
    AC_stim_amplitude = (-80.0);
    AC_stim_end = 9000.0;
    
    /* plateau_potassium_current */
    AC_g_Kp = 0.002216;
    
    /* rapid_activating_delayed_rectifiyer_K_current */
    AC_E_K = AC_R * AC_T / AC_F * log(AC_K_o / AC_K_i);
    AC_g_Kr = 0.0136;
    
    /* slow_activating_delayed_rectifiyer_K_current */
    AC_E_Ks = AC_R * AC_T / AC_F * log((AC_K_o + 0.01833 * AC_Na_o) / (AC_K_i + 0.01833 * AC_Na_i));
    AC_g_Ks = 0.0245;
    
    /* sodium_background_current */
    AC_g_Nab = 0.0031;
    
    /* sodium_potassium_pump */
    AC_K_mKo = 1.5;
    AC_K_mNai = 10.0;
    AC_i_NaK_max = 0.693;
    AC_sigma = 1.0 / 7.0 * (exp(AC_Na_o / 67.3) - 1.0);
    
    /* time_independent_potassium_current */
    AC_K_mK1 = 13.0;
    AC_g_K1 = 2.8;
    
    /* transient_outward_potassium_current */
    AC_g_to = 0.23815;
    
    /* *remaining* */
    AC_E_Na = AC_R * AC_T / AC_F * log(AC_Na_o / AC_Na_i);
    
}

/* Right-hand-side function of the model ODE */
static int rhs(realtype t, const N_Vector y, N_Vector ydot, void *f_data)
{
    /* L_type_Ca_current_d_gate */
    auto AV_L_type_Ca_current_d_gate_E0_m = NV_Ith_S(y, 0) + 40.0;
    auto AV_d_infinity = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 10.0) / (-6.24)));
    auto AV_tau_d = 1.0 / (0.25 * exp((-0.01) * NV_Ith_S(y, 0)) / (1.0 + exp((-0.07) * NV_Ith_S(y, 0))) + 0.07 * exp((-0.05) * AV_L_type_Ca_current_d_gate_E0_m) / (1.0 + exp(0.05 * AV_L_type_Ca_current_d_gate_E0_m)));
    NV_Ith_S(ydot, 9) = (AV_d_infinity - NV_Ith_S(y, 9)) / AV_tau_d;
    
    /* L_type_Ca_current_f_Ca_gate */
    auto AV_f_Ca_infinity = 1.0 / (1.0 + pow(NV_Ith_S(y, 11) / AC_K_mfCa, 3.0));
    NV_Ith_S(ydot, 10) = (AV_f_Ca_infinity - NV_Ith_S(y, 10)) / AC_tau_f_Ca;
    
    /* L_type_Ca_current_f_gate */
    auto AV_f_infinity = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 12.5) / 5.0));
    auto AV_tau_f = 30.0 + 200.0 / (1.0 + exp((NV_Ith_S(y, 0) + 20.0) / 9.5));
    NV_Ith_S(ydot, 8) = (AV_f_infinity - NV_Ith_S(y, 8)) / AV_tau_f;
   
    
    /* fast_sodium_current_h_gate */
    auto AV_alpha_h = 0.135 * exp((NV_Ith_S(y, 0) + 80.0 - AC_shift_h) / (-6.8));
    auto AV_beta_h = 7.5 / (1.0 + exp((-0.1) * (NV_Ith_S(y, 0) + 11.0 - AC_shift_h)));
    NV_Ith_S(ydot, 2) = AV_alpha_h * (1.0 - NV_Ith_S(y, 2)) - AV_beta_h * NV_Ith_S(y, 2);
    
    /* fast_sodium_current_j_gate */
    auto AV_alpha_j = 0.175 * exp((NV_Ith_S(y, 0) + 100.0 - AC_shift_j) / (-23.0)) / (1.0 + exp(0.15 * (NV_Ith_S(y, 0) + 79.0 - AC_shift_j)));
    auto AV_beta_j = 0.3 / (1.0 + exp((-0.1) * (NV_Ith_S(y, 0) + 32.0 - AC_shift_j)));
    NV_Ith_S(ydot, 3) = AV_alpha_j * (1.0 - NV_Ith_S(y, 3)) - AV_beta_j * NV_Ith_S(y, 3);
    
    /* fast_sodium_current_m_gate */
    auto AV_fast_sodium_current_m_gate_E0_m = NV_Ith_S(y, 0) + 47.13;
    auto AV_beta_m = 0.08 * exp((-NV_Ith_S(y, 0)) / 11.0);
    auto AV_alpha_m = 0.32 * AV_fast_sodium_current_m_gate_E0_m / (1.0 - exp((-0.1) * AV_fast_sodium_current_m_gate_E0_m));
    NV_Ith_S(ydot, 1) = AV_alpha_m * (1.0 - NV_Ith_S(y, 1)) - AV_beta_m * NV_Ith_S(y, 1);
    
    /* plateau_potassium_current_Kp_gate */
    auto AV_Kp_V = 1.0 / (1.0 + exp((7.488 - NV_Ith_S(y, 0)) / 5.98));
    
    /* rapid_activating_delayed_rectifiyer_K_current_X_kr_gate */
    auto AV_X_kr_inf = 1.0 / (1.0 + exp((-2.182) - 0.1819 * NV_Ith_S(y, 0)));
    auto AV_tau_X_kr = 43.0 + 1.0 / (exp((-5.495) + 0.1691 * NV_Ith_S(y, 0)) + exp((-7.677) - 0.0128 * NV_Ith_S(y, 0)));
    NV_Ith_S(ydot, 4) = (AV_X_kr_inf - NV_Ith_S(y, 4)) / AV_tau_X_kr;
    
    /* sarcolemmal_calcium_pump */
    auto AV_i_p_Ca = AC_i_pCa_max * NV_Ith_S(y, 11) / (AC_K_mpCa + NV_Ith_S(y, 11));
    
    /* slow_activating_delayed_rectifiyer_K_current_X_ks_gate */
    auto AV_X_ks_infinity = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) - 16.0) / (-13.6)));
    auto AV_tau_X_ks = 1.0 / (7.19e-05 * (NV_Ith_S(y, 0) - 10.0) / (1.0 - exp((-0.148) * (NV_Ith_S(y, 0) - 10.0))) + 0.000131 * (NV_Ith_S(y, 0) - 10.0) / (exp(0.0687 * (NV_Ith_S(y, 0) - 10.0)) - 1.0));
    NV_Ith_S(ydot, 5) = (AV_X_ks_infinity - NV_Ith_S(y, 5)) / AV_tau_X_ks;
    
    /* transient_outward_potassium_current_X_to_gate */
    auto AV_alpha_X_to = 0.04516 * exp(0.03577 * NV_Ith_S(y, 0));
    auto AV_beta_X_to = 0.0989 * exp((-0.06237) * NV_Ith_S(y, 0));
    NV_Ith_S(ydot, 6) = AV_alpha_X_to * (1.0 - NV_Ith_S(y, 6)) - AV_beta_X_to * NV_Ith_S(y, 6);
    
    /* transient_outward_potassium_current_Y_to_gate */
    auto AV_alpha_Y_to = 0.005415 * exp((NV_Ith_S(y, 0) + 33.5) / (-5.0)) / (1.0 + 0.051335 * exp((NV_Ith_S(y, 0) + 33.5) / (-5.0)));
    auto AV_beta_Y_to = 0.005415 * exp((NV_Ith_S(y, 0) + 33.5) / 5.0) / (1.0 + 0.051335 * exp((NV_Ith_S(y, 0) + 33.5) / 5.0));
    NV_Ith_S(ydot, 7) = AV_alpha_Y_to * (1.0 - NV_Ith_S(y, 7)) - AV_beta_Y_to * NV_Ith_S(y, 7);
    
    /* calcium_dynamics */
    auto AV_calcium_dynamics_gamma = 1.0 / (1.0 + pow(2000.0 / NV_Ith_S(y, 12), 3.0));
    auto AV_J_leak = AC_P_leak * (NV_Ith_S(y, 12) - NV_Ith_S(y, 11));
    auto AV_J_rel = AC_P_rel * NV_Ith_S(y, 8) * NV_Ith_S(y, 9) * NV_Ith_S(y, 10) * (AV_calcium_dynamics_gamma * NV_Ith_S(y, 12) - NV_Ith_S(y, 11)) / (1.0 + 1.65 * exp(NV_Ith_S(y, 0) / 20.0));
    auto AV_J_up = AC_V_up / (1.0 + pow(AC_K_mup / NV_Ith_S(y, 11), 2.0));
    auto AV_beta_SR = 1.0 / (1.0 + AC_CSQN_tot * AC_K_mCSQN / pow(AC_K_mCSQN + NV_Ith_S(y, 12), 2.0));
    auto AV_beta_i = 1.0 / (1.0 + AC_CMDN_tot * AC_K_mCMDN / pow(AC_K_mCMDN + NV_Ith_S(y, 11), 2.0));
    NV_Ith_S(ydot, 12) = AV_beta_SR * (AV_J_up - AV_J_leak - AV_J_rel) * AC_V_myo / AC_V_SR;
   
    
    /* rapid_activating_delayed_rectifiyer_K_current */
    auto AV_R_V = 1.0 / (1.0 + 2.5 * exp(0.1 * (NV_Ith_S(y, 0) + 28.0)));
    auto AV_i_Kr = AC_g_Kr * AV_R_V * NV_Ith_S(y, 4) * sqrt(AC_K_o / 4.0) * (NV_Ith_S(y, 0) - AC_E_K);
    
    /* slow_activating_delayed_rectifiyer_K_current */
    auto AV_i_Ks = AC_g_Ks * pow(NV_Ith_S(y, 5), 2.0) * (NV_Ith_S(y, 0) - AC_E_Ks);
    
    /* sodium_potassium_pump */
    auto AV_f_NaK = 1.0 / (1.0 + 0.1245 * exp((-0.1) * NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T)) + 0.0365 * AC_sigma * exp((-NV_Ith_S(y, 0)) * AC_F / (AC_R * AC_T)));
    auto AV_i_NaK = AC_i_NaK_max * AV_f_NaK / (1.0 + pow(AC_K_mNai / AC_Na_i, 1.5)) * AC_K_o / (AC_K_o + AC_K_mKo);
    
    /* time_independent_potassium_current_K1_gate */
    auto AV_K1_infinity = 1.0 / (2.0 + exp(1.62 * AC_F / (AC_R * AC_T) * (NV_Ith_S(y, 0) - AC_E_K)));
    
    /* transient_outward_potassium_current */
    auto AV_i_to = AC_g_to * NV_Ith_S(y, 6) * NV_Ith_S(y, 7) * (NV_Ith_S(y, 0) - AC_E_K);
    
    /* *remaining* */
    auto AV_i_Ca_max = AC_P_Ca / AC_C_sc * 4.0 * NV_Ith_S(y, 0) * pow(AC_F, 2.0) / (AC_R * AC_T) * (NV_Ith_S(y, 11) * exp(2.0 * NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T)) - 0.341 * AC_Ca_o) / (exp(2.0 * NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T)) - 1.0);
    auto AV_i_NaCa = AC_K_NaCa / ((pow(AC_K_mNa, 3.0) + pow(AC_Na_o, 3.0)) * (AC_K_mCa + AC_Ca_o) * (1.0 + AC_K_sat * exp((AC_eta - 1.0) * NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T)))) * (exp(AC_eta * NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T)) * pow(AC_Na_i, 3.0) * AC_Ca_o - exp((AC_eta - 1.0) * NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T)) * pow(AC_Na_o, 3.0) * NV_Ith_S(y, 11));
    auto AV_E_Ca = AC_R * AC_T / (2.0 * AC_F) * log(AC_Ca_o / NV_Ith_S(y, 11));
    auto AV_i_Kp = AC_g_Kp * AV_Kp_V * (NV_Ith_S(y, 0) - AC_E_K);
    auto AV_i_K1 = AC_g_K1 * AV_K1_infinity * AC_K_o / (AC_K_o + AC_K_mK1) * (NV_Ith_S(y, 0) - AC_E_K);
    auto AV_i_Ca = AV_i_Ca_max * NV_Ith_S(y, 8) * NV_Ith_S(y, 9) * NV_Ith_S(y, 10);
    auto AV_i_CaK = AC_P_CaK / AC_C_sc * NV_Ith_S(y, 8) * NV_Ith_S(y, 9) * NV_Ith_S(y, 10) / (1.0 + AV_i_Ca_max / AC_i_Ca_half) * 1000.0 * NV_Ith_S(y, 0) * pow(AC_F, 2.0) / (AC_R * AC_T) * (AC_K_i * exp(NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T)) - AC_K_o) / (exp(NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T)) - 1.0);
    auto AV_i_Ca_b = AC_g_Cab * (NV_Ith_S(y, 0) - AV_E_Ca);
    auto AV_i_Na = AC_g_Na * pow(NV_Ith_S(y, 1), 3.0) * NV_Ith_S(y, 2) * NV_Ith_S(y, 3) * (NV_Ith_S(y, 0) - AC_E_Na);
    auto AV_i_Na_b = AC_g_Nab * (NV_Ith_S(y, 0) - AC_E_Na);
    NV_Ith_S(ydot, 11) = AV_beta_i * (AV_J_rel + AV_J_leak - AV_J_up - AC_A_Cap * AC_C_sc / (2.0 * AC_F * AC_V_myo) * (AV_i_Ca + AV_i_Ca_b + AV_i_p_Ca - 2.0 * AV_i_NaCa));
    NV_Ith_S(ydot, 0) = (-(AV_i_Na + AV_i_Ca + AV_i_CaK + AV_i_Kr + AV_i_Ks + AV_i_to + AV_i_K1 + AV_i_Kp + AV_i_NaCa + AV_i_NaK + AV_i_p_Ca + AV_i_Na_b + AV_i_Ca_b));
    

    return 0;
}

/* Set initial values */
static void
default_initial_values(N_Vector y)
{
    NV_Ith_S(y, 0) = -94.7;
    NV_Ith_S(y, 1) =  2.46760000000000002e-04;
    NV_Ith_S(y, 2) = 0.99869;
    NV_Ith_S(y, 3) = 0.99887;
    NV_Ith_S(y, 4) = 0.229;
    NV_Ith_S(y, 5) = 0.0001;
    NV_Ith_S(y, 6) = 3.742e-05;
    NV_Ith_S(y, 7) = 1.0;
    NV_Ith_S(y, 8) = 0.983;
    NV_Ith_S(y, 9) = 0.0001;
    NV_Ith_S(y, 10) = 0.942;
    NV_Ith_S(y, 11) = 0.0472;
    NV_Ith_S(y, 12) = 320.0;

}