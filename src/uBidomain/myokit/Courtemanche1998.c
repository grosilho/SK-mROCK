/*
courtemanche_1998
Generated on 2022-09-22 13:30:38

*/

#define N_STATE 21

/* Declare intermediary, temporary and system variables */
static realtype AC_CMDN_max;
static realtype AC_CSQN_max;
static realtype AV_Ca_CMDN;
static realtype AV_Ca_CSQN;
static realtype AV_Ca_TRPN;
static realtype AC_Km_CMDN;
static realtype AC_Km_CSQN;
static realtype AC_Km_TRPN;
static realtype AC_TRPN_max;
static realtype AC_Ca_up_max;
static realtype AV_i_up_leak;
static realtype AV_Fn;
static realtype AC_K_rel;
static realtype AV_i_rel;
static realtype AC_tau_u;
static realtype AV_u_infinity;
static realtype AV_tau_v;
static realtype AV_v_infinity;
static realtype AV_tau_w;
static realtype AV_w_infinity;
static realtype AC_I_up_max;
static realtype AC_K_up;
static realtype AV_i_up;
static realtype AC_g_Ca_L;
static realtype AV_i_Ca_L;
static realtype AV_d_infinity;
static realtype AV_tau_d;
static realtype AV_f_Ca_infinity;
static realtype AC_tau_f_Ca;
static realtype AV_f_infinity;
static realtype AV_tau_f;
static realtype AC_I_NaCa_max;
static realtype AC_K_mCa;
static realtype AC_K_mNa;
static realtype AC_K_sat;
static realtype AC_Na_Ca_exchanger_current_gamma;
static realtype AV_i_NaCa;
static realtype AV_E_Ca;
static realtype AC_g_B_Ca;
static realtype AC_g_B_K;
static realtype AC_g_B_Na;
static realtype AV_i_B_Ca;
static realtype AV_i_B_K;
static realtype AV_i_B_Na;
static realtype AV_time;
static realtype AV_E_Na;
static realtype AC_g_Na;
static realtype AV_i_Na;
static realtype AV_alpha_h;
static realtype AV_beta_h;
static realtype AV_h_inf;
static realtype AV_tau_h;
static realtype AV_alpha_j;
static realtype AV_beta_j;
static realtype AV_j_inf;
static realtype AV_tau_j;
static realtype AV_alpha_m;
static realtype AV_beta_m;
static realtype AV_m_inf;
static realtype AV_tau_m;
static realtype AV_B1;
static realtype AV_B2;
static realtype AC_V_cell;
static realtype AC_V_i;
static realtype AC_V_rel;
static realtype AC_V_up;
static realtype AC_Cm;
static realtype AC_F;
static realtype AC_R;
static realtype AC_T;
static realtype AV_i_st;
static realtype AC_g_Kr;
static realtype AV_i_Kr;
static realtype AV_alpha_xr;
static realtype AV_beta_xr;
static realtype AV_tau_xr;
static realtype AV_xr_infinity;
static realtype AV_i_CaP;
static realtype AC_i_CaP_max;
static realtype AC_g_Ks;
static realtype AV_i_Ks;
static realtype AV_alpha_xs;
static realtype AV_beta_xs;
static realtype AV_tau_xs;
static realtype AV_xs_infinity;
static realtype AC_Km_K_o;
static realtype AC_Km_Na_i;
static realtype AV_f_NaK;
static realtype AV_i_NaK;
static realtype AC_i_NaK_max;
static realtype AC_sigma;
static realtype AC_Ca_o;
static realtype AC_K_o;
static realtype AC_Na_o;
static realtype AV_E_K;
static realtype AC_g_K1;
static realtype AV_i_K1;
static realtype AV_i_tr;
static realtype AC_tau_tr;
static realtype AC_K_Q10;
static realtype AC_g_to;
static realtype AV_i_to;
static realtype AV_alpha_oa;
static realtype AV_beta_oa;
static realtype AV_oa_infinity;
static realtype AV_tau_oa;
static realtype AV_alpha_oi;
static realtype AV_beta_oi;
static realtype AV_oi_infinity;
static realtype AV_tau_oi;
static realtype AV_g_Kur;
static realtype AV_i_Kur;
static realtype AV_alpha_ua;
static realtype AV_beta_ua;
static realtype AV_tau_ua;
static realtype AV_ua_infinity;
static realtype AV_alpha_ui;
static realtype AV_beta_ui;
static realtype AV_tau_ui;
static realtype AV_ui_infinity;

/* Set values of constants */
static void
updateConstants(void)
{
    /* Ca_buffers */
    AC_CMDN_max = 0.05;
    AC_CSQN_max = 10.0;
    AC_Km_CMDN = 0.00238;
    AC_Km_CSQN = 0.8;
    AC_Km_TRPN = 0.0005;
    AC_TRPN_max = 0.07;
    
    /* Ca_uptake_current_by_the_NSR */
    AC_I_up_max = 0.005;
    AC_K_up = 0.00092;
    
    /* L_type_Ca_channel_f_Ca_gate */
    AC_tau_f_Ca = 2.0;
    
    /* standard_ionic_concentrations */
    AC_Ca_o = 1.8;
    AC_K_o = 5.4;
    AC_Na_o = 140.0;
    
    /* transfer_current_from_NSR_to_JSR */
    AC_tau_tr = 180.0;
    
    /* Ca_leak_current_by_the_NSR */
    AC_Ca_up_max = 15.0;
    
    /* Ca_release_current_from_JSR */
    AC_K_rel = 30.0;
    
    /* Ca_release_current_from_JSR_u_gate */
    AC_tau_u = 8.0;
    
    /* L_type_Ca_channel */
    AC_g_Ca_L = 0.12375;
    
    /* Na_Ca_exchanger_current */
    AC_I_NaCa_max = 1600.0;
    AC_K_mCa = 1.38;
    AC_K_mNa = 87.5;
    AC_K_sat = 0.1;
    AC_Na_Ca_exchanger_current_gamma = 0.35;
    
    /* background_currents */
    AC_g_B_Ca = 0.001131;
    AC_g_B_K = 0.0;
    AC_g_B_Na =  6.74437500000000015e-04;
    
    /* fast_sodium_current */
    AC_g_Na = 7.8;
    
    /* intracellular_ion_concentrations */
    AC_V_cell = 20100.0;
    AC_V_i = AC_V_cell * 0.68;
    AC_V_rel = 0.0048 * AC_V_cell;
    AC_V_up = 0.0552 * AC_V_cell;
    
    /* membrane */
    AC_Cm = 100.0;
    AC_F = 96.4867;
    AC_R = 8.3143;
    AC_T = 310.0;
    
    /* rapid_delayed_rectifier_K_current */
    AC_g_Kr =  2.94117649999999994e-02;
    
    /* sarcolemmal_calcium_pump_current */
    AC_i_CaP_max = 0.275;
    
    /* slow_delayed_rectifier_K_current */
    AC_g_Ks =  1.29411759999999987e-01;
    
    /* sodium_potassium_pump */
    AC_Km_K_o = 1.5;
    AC_Km_Na_i = 10.0;
    AC_i_NaK_max =  5.99338739999999981e-01;
    AC_sigma = 1.0 / 7.0 * (exp(AC_Na_o / 67.3) - 1.0);
    
    /* time_independent_potassium_current */
    AC_g_K1 = 0.09;
    
    /* transient_outward_K_current */
    AC_K_Q10 = 3.0;
    AC_g_to = 0.1652;
    
}

/* Right-hand-side function of the model ODE */
static int rhs(realtype t, const N_Vector y, N_Vector ydot, void *f_data)
{
    /* Ca_buffers */
    auto AV_Ca_CMDN = AC_CMDN_max * NV_Ith_S(y, 17) / (NV_Ith_S(y, 17) + AC_Km_CMDN);
    auto AV_Ca_CSQN = AC_CSQN_max * NV_Ith_S(y, 19) / (NV_Ith_S(y, 19) + AC_Km_CSQN);
    auto AV_Ca_TRPN = AC_TRPN_max * NV_Ith_S(y, 17) / (NV_Ith_S(y, 17) + AC_Km_TRPN);
    
    /* Ca_release_current_from_JSR_w_gate */
    auto AV_tau_w = (abs(NV_Ith_S(y, 0) - 7.9) < 1e-10).select(6.0 * 0.2 / 1.3, 6.0 * (1.0 - exp((-(NV_Ith_S(y, 0) - 7.9)) / 5.0)) / ((1.0 + 0.3 * exp((-(NV_Ith_S(y, 0) - 7.9)) / 5.0)) * 1.0 * (NV_Ith_S(y, 0) - 7.9)));
    auto AV_w_infinity = 1.0 - pow(1.0 + exp((-(NV_Ith_S(y, 0) - 40.0)) / 17.0), (-1.0));
    NV_Ith_S(ydot, 15) = (AV_w_infinity - NV_Ith_S(y, 15)) / AV_tau_w;
    
    /* Ca_uptake_current_by_the_NSR */
    auto AV_i_up = AC_I_up_max / (1.0 + AC_K_up / NV_Ith_S(y, 17));
    
    /* L_type_Ca_channel_d_gate */
    auto AV_d_infinity = pow(1.0 + exp((NV_Ith_S(y, 0) + 10.0) / (-8.0)), (-1.0));
    auto AV_tau_d = (abs(NV_Ith_S(y, 0) + 10.0) < 1e-10).select(4.579 / (1.0 + exp((NV_Ith_S(y, 0) + 10.0) / (-6.24))), (1.0 - exp((NV_Ith_S(y, 0) + 10.0) / (-6.24))) / (0.035 * (NV_Ith_S(y, 0) + 10.0) * (1.0 + exp((NV_Ith_S(y, 0) + 10.0) / (-6.24)))));
    NV_Ith_S(ydot, 10) = (AV_d_infinity - NV_Ith_S(y, 10)) / AV_tau_d;
    
    /* L_type_Ca_channel_f_Ca_gate */
    auto AV_f_Ca_infinity = pow(1.0 + NV_Ith_S(y, 17) / 0.00035, (-1.0));
    NV_Ith_S(ydot, 12) = (AV_f_Ca_infinity - NV_Ith_S(y, 12)) / AC_tau_f_Ca;
    
    /* L_type_Ca_channel_f_gate */
    auto AV_f_infinity = exp((-(NV_Ith_S(y, 0) + 28.0)) / 6.9) / (1.0 + exp((-(NV_Ith_S(y, 0) + 28.0)) / 6.9));
    auto AV_tau_f = 9.0 * pow(0.0197 * exp((-pow(0.0337, 2.0)) * pow(NV_Ith_S(y, 0) + 10.0, 2.0)) + 0.02, (-1.0));
    NV_Ith_S(ydot, 11) = (AV_f_infinity - NV_Ith_S(y, 11)) / AV_tau_f;
    
    /* environment */
    auto AV_time = t;
    
    /* fast_sodium_current_h_gate */
    auto AV_alpha_h = (NV_Ith_S(y, 0) < (-40.0)).select(0.135 * exp((NV_Ith_S(y, 0) + 80.0) / (-6.8)), 0.0);
    auto AV_beta_h = (NV_Ith_S(y, 0) < (-40.0)).select(3.56 * exp(0.079 * NV_Ith_S(y, 0)) + 310000.0 * exp(0.35 * NV_Ith_S(y, 0)), 1.0 / (0.13 * (1.0 + exp((NV_Ith_S(y, 0) + 10.66) / (-11.1)))));
    auto AV_h_inf = AV_alpha_h / (AV_alpha_h + AV_beta_h);
    auto AV_tau_h = 1.0 / (AV_alpha_h + AV_beta_h);
    NV_Ith_S(ydot, 2) = (AV_h_inf - NV_Ith_S(y, 2)) / AV_tau_h;
    
    /* fast_sodium_current_j_gate */
    auto AV_alpha_j = (NV_Ith_S(y, 0) < (-40.0)).select(((-127140.0) * exp(0.2444 * NV_Ith_S(y, 0)) - 3.474e-05 * exp((-0.04391) * NV_Ith_S(y, 0))) * (NV_Ith_S(y, 0) + 37.78) / (1.0 + exp(0.311 * (NV_Ith_S(y, 0) + 79.23))), 0.0);
    auto AV_beta_j = (NV_Ith_S(y, 0) < (-40.0)).select(0.1212 * exp((-0.01052) * NV_Ith_S(y, 0)) / (1.0 + exp((-0.1378) * (NV_Ith_S(y, 0) + 40.14))), 0.3 * exp((-2.535e-07) * NV_Ith_S(y, 0)) / (1.0 + exp((-0.1) * (NV_Ith_S(y, 0) + 32.0))));
    auto AV_j_inf = AV_alpha_j / (AV_alpha_j + AV_beta_j);
    auto AV_tau_j = 1.0 / (AV_alpha_j + AV_beta_j);
    NV_Ith_S(ydot, 3) = (AV_j_inf - NV_Ith_S(y, 3)) / AV_tau_j;
    
    /* fast_sodium_current_m_gate */
    auto AV_alpha_m = (NV_Ith_S(y, 0) == (-47.13)).select(3.2, 0.32 * (NV_Ith_S(y, 0) + 47.13) / (1.0 - exp((-0.1) * (NV_Ith_S(y, 0) + 47.13))));
    auto AV_beta_m = 0.08 * exp((-NV_Ith_S(y, 0)) / 11.0);
    auto AV_m_inf = AV_alpha_m / (AV_alpha_m + AV_beta_m);
    auto AV_tau_m = 1.0 / (AV_alpha_m + AV_beta_m);
    NV_Ith_S(ydot, 1) = (AV_m_inf - NV_Ith_S(y, 1)) / AV_tau_m;
    
    /* rapid_delayed_rectifier_K_current_xr_gate */
    auto AV_alpha_xr = (abs(NV_Ith_S(y, 0) + 14.1) < 1e-10).select(0.0015, 0.0003 * (NV_Ith_S(y, 0) + 14.1) / (1.0 - exp((NV_Ith_S(y, 0) + 14.1) / (-5.0))));
    auto AV_beta_xr = (abs(NV_Ith_S(y, 0) - 3.3328) < 1e-10).select(3.78361180000000004e-04,  7.38980000000000030e-05 * (NV_Ith_S(y, 0) - 3.3328) / (exp((NV_Ith_S(y, 0) - 3.3328) / 5.1237) - 1.0));
    auto AV_xr_infinity = pow(1.0 + exp((NV_Ith_S(y, 0) + 14.1) / (-6.5)), (-1.0));
    auto AV_tau_xr = pow(AV_alpha_xr + AV_beta_xr, (-1.0));
    NV_Ith_S(ydot, 8) = (AV_xr_infinity - NV_Ith_S(y, 8)) / AV_tau_xr;
    
    /* slow_delayed_rectifier_K_current_xs_gate */
    auto AV_alpha_xs = (abs(NV_Ith_S(y, 0) - 19.9) < 1e-10).select(0.00068, 4e-05 * (NV_Ith_S(y, 0) - 19.9) / (1.0 - exp((NV_Ith_S(y, 0) - 19.9) / (-17.0))));
    auto AV_beta_xs = (abs(NV_Ith_S(y, 0) - 19.9) < 1e-10).select(0.000315, 3.5e-05 * (NV_Ith_S(y, 0) - 19.9) / (exp((NV_Ith_S(y, 0) - 19.9) / 9.0) - 1.0));
    auto AV_xs_infinity = pow(1.0 + exp((NV_Ith_S(y, 0) - 19.9) / (-12.7)), (-0.5));
    auto AV_tau_xs = 0.5 * pow(AV_alpha_xs + AV_beta_xs, (-1.0));
    NV_Ith_S(ydot, 9) = (AV_xs_infinity - NV_Ith_S(y, 9)) / AV_tau_xs;
    
    /* transfer_current_from_NSR_to_JSR */
    auto AV_i_tr = (NV_Ith_S(y, 20) - NV_Ith_S(y, 19)) / AC_tau_tr;
    
    /* Ca_leak_current_by_the_NSR */
    auto AV_i_up_leak = AC_I_up_max * NV_Ith_S(y, 20) / AC_Ca_up_max;
    
    /* Ca_release_current_from_JSR */
    auto AV_i_rel = AC_K_rel * pow(NV_Ith_S(y, 13), 2.0) * NV_Ith_S(y, 14) * NV_Ith_S(y, 15) * (NV_Ith_S(y, 19) - NV_Ith_S(y, 17));
    
    /* intracellular_ion_concentrations */
    NV_Ith_S(ydot, 19) = (AV_i_tr - AV_i_rel) * pow(1.0 + AC_CSQN_max * AC_Km_CSQN / pow(NV_Ith_S(y, 19) + AC_Km_CSQN, 2.0), (-1.0));
    auto AV_B2 = 1.0 + AC_TRPN_max * AC_Km_TRPN / pow(NV_Ith_S(y, 17) + AC_Km_TRPN, 2.0) + AC_CMDN_max * AC_Km_CMDN / pow(NV_Ith_S(y, 17) + AC_Km_CMDN, 2.0);
    NV_Ith_S(ydot, 20) = AV_i_up - (AV_i_up_leak + AV_i_tr * AC_V_rel / AC_V_up);
    
    /* sarcolemmal_calcium_pump_current */
    auto AV_i_CaP = AC_Cm * AC_i_CaP_max * NV_Ith_S(y, 17) / (0.0005 + NV_Ith_S(y, 17));
    
    /* sodium_potassium_pump */
    auto AV_f_NaK = pow(1.0 + 0.1245 * exp((-0.1) * AC_F * NV_Ith_S(y, 0) / (AC_R * AC_T)) + 0.0365 * AC_sigma * exp((-AC_F) * NV_Ith_S(y, 0) / (AC_R * AC_T)), (-1.0));
    auto AV_i_NaK = AC_Cm * AC_i_NaK_max * AV_f_NaK * 1.0 / (1.0 + pow(AC_Km_Na_i / NV_Ith_S(y, 16), 1.5)) * AC_K_o / (AC_K_o + AC_Km_K_o);
    
    /* time_independent_potassium_current */
    auto AV_E_K = AC_R * AC_T / AC_F * log(AC_K_o / NV_Ith_S(y, 18));
    auto AV_i_K1 = AC_Cm * AC_g_K1 * (NV_Ith_S(y, 0) - AV_E_K) / (1.0 + exp(0.07 * (NV_Ith_S(y, 0) + 80.0)));
    
    /* transient_outward_K_current */
    auto AV_i_to = AC_Cm * AC_g_to * pow(NV_Ith_S(y, 4), 3.0) * NV_Ith_S(y, 5) * (NV_Ith_S(y, 0) - AV_E_K);
    
    /* transient_outward_K_current_oa_gate */
    auto AV_alpha_oa = 0.65 * pow(exp((NV_Ith_S(y, 0) - (-10.0)) / (-8.5)) + exp((NV_Ith_S(y, 0) - (-10.0) - 40.0) / (-59.0)), (-1.0));
    auto AV_beta_oa = 0.65 * pow(2.5 + exp((NV_Ith_S(y, 0) - (-10.0) + 72.0) / 17.0), (-1.0));
    auto AV_oa_infinity = pow(1.0 + exp((NV_Ith_S(y, 0) - (-10.0) + 10.47) / (-17.54)), (-1.0));
    auto AV_tau_oa = pow(AV_alpha_oa + AV_beta_oa, (-1.0)) / AC_K_Q10;
    NV_Ith_S(ydot, 4) = (AV_oa_infinity - NV_Ith_S(y, 4)) / AV_tau_oa;
    
    /* transient_outward_K_current_oi_gate */
    auto AV_alpha_oi = pow(18.53 + 1.0 * exp((NV_Ith_S(y, 0) - (-10.0) + 103.7) / 10.95), (-1.0));
    auto AV_beta_oi = pow(35.56 + 1.0 * exp((NV_Ith_S(y, 0) - (-10.0) - 8.74) / (-7.44)), (-1.0));
    auto AV_oi_infinity = pow(1.0 + exp((NV_Ith_S(y, 0) - (-10.0) + 33.1) / 5.3), (-1.0));
    auto AV_tau_oi = pow(AV_alpha_oi + AV_beta_oi, (-1.0)) / AC_K_Q10;
    NV_Ith_S(ydot, 5) = (AV_oi_infinity - NV_Ith_S(y, 5)) / AV_tau_oi;
    
    /* ultrarapid_delayed_rectifier_K_current */
    auto AV_g_Kur = 0.005 + 0.05 / (1.0 + exp((NV_Ith_S(y, 0) - 15.0) / (-13.0)));
    auto AV_i_Kur = AC_Cm * AV_g_Kur * pow(NV_Ith_S(y, 6), 3.0) * NV_Ith_S(y, 7) * (NV_Ith_S(y, 0) - AV_E_K);
    
    /* ultrarapid_delayed_rectifier_K_current_ua_gate */
    auto AV_alpha_ua = 0.65 * pow(exp((NV_Ith_S(y, 0) - (-10.0)) / (-8.5)) + exp((NV_Ith_S(y, 0) - (-10.0) - 40.0) / (-59.0)), (-1.0));
    auto AV_beta_ua = 0.65 * pow(2.5 + exp((NV_Ith_S(y, 0) - (-10.0) + 72.0) / 17.0), (-1.0));
    auto AV_ua_infinity = pow(1.0 + exp((NV_Ith_S(y, 0) - (-10.0) + 20.3) / (-9.6)), (-1.0));
    auto AV_tau_ua = pow(AV_alpha_ua + AV_beta_ua, (-1.0)) / AC_K_Q10;
    NV_Ith_S(ydot, 6) = (AV_ua_infinity - NV_Ith_S(y, 6)) / AV_tau_ua;
    
    /* ultrarapid_delayed_rectifier_K_current_ui_gate */
    auto AV_alpha_ui = pow(21.0 + 1.0 * exp((NV_Ith_S(y, 0) - (-10.0) - 195.0) / (-28.0)), (-1.0));
    auto AV_beta_ui = 1.0 / exp((NV_Ith_S(y, 0) - (-10.0) - 168.0) / (-16.0));
    auto AV_ui_infinity = pow(1.0 + exp((NV_Ith_S(y, 0) - (-10.0) - 109.45) / 27.48), (-1.0));
    auto AV_tau_ui = pow(AV_alpha_ui + AV_beta_ui, (-1.0)) / AC_K_Q10;
    NV_Ith_S(ydot, 7) = (AV_ui_infinity - NV_Ith_S(y, 7)) / AV_tau_ui;
    
    /* *remaining* */
    auto AV_i_Ca_L = AC_Cm * AC_g_Ca_L * NV_Ith_S(y, 10) * NV_Ith_S(y, 11) * NV_Ith_S(y, 12) * (NV_Ith_S(y, 0) - 65.0);
    auto AV_i_NaCa = AC_Cm * AC_I_NaCa_max * (exp(AC_Na_Ca_exchanger_current_gamma * AC_F * NV_Ith_S(y, 0) / (AC_R * AC_T)) * pow(NV_Ith_S(y, 16), 3.0) * AC_Ca_o - exp((AC_Na_Ca_exchanger_current_gamma - 1.0) * AC_F * NV_Ith_S(y, 0) / (AC_R * AC_T)) * pow(AC_Na_o, 3.0) * NV_Ith_S(y, 17)) / ((pow(AC_K_mNa, 3.0) + pow(AC_Na_o, 3.0)) * (AC_K_mCa + AC_Ca_o) * (1.0 + AC_K_sat * exp((AC_Na_Ca_exchanger_current_gamma - 1.0) * NV_Ith_S(y, 0) * AC_F / (AC_R * AC_T))));
    auto AV_E_Ca = AC_R * AC_T / (2.0 * AC_F) * log(AC_Ca_o / NV_Ith_S(y, 17));
    auto AV_i_B_K = AC_Cm * AC_g_B_K * (NV_Ith_S(y, 0) - AV_E_K);
    auto AV_E_Na = AC_R * AC_T / AC_F * log(AC_Na_o / NV_Ith_S(y, 16));
    auto AV_i_Kr = AC_Cm * AC_g_Kr * NV_Ith_S(y, 8) * (NV_Ith_S(y, 0) - AV_E_K) / (1.0 + exp((NV_Ith_S(y, 0) + 15.0) / 22.4));
    auto AV_i_Ks = AC_Cm * AC_g_Ks * pow(NV_Ith_S(y, 9), 2.0) * (NV_Ith_S(y, 0) - AV_E_K);
    auto AV_Fn = 1000.0 * (1e-15 * AC_V_rel * AV_i_rel - 1e-15 / (2.0 * AC_F) * (0.5 * AV_i_Ca_L - 0.2 * AV_i_NaCa));
    auto AV_i_B_Ca = AC_Cm * AC_g_B_Ca * (NV_Ith_S(y, 0) - AV_E_Ca);
    auto AV_i_B_Na = AC_Cm * AC_g_B_Na * (NV_Ith_S(y, 0) - AV_E_Na);
    auto AV_i_Na = AC_Cm * AC_g_Na * pow(NV_Ith_S(y, 1), 3.0) * NV_Ith_S(y, 2) * NV_Ith_S(y, 3) * (NV_Ith_S(y, 0) - AV_E_Na);
    NV_Ith_S(ydot, 18) = (2.0 * AV_i_NaK - (AV_i_K1 + AV_i_to + AV_i_Kur + AV_i_Kr + AV_i_Ks + AV_i_B_K)) / (AC_V_i * AC_F);
    auto AV_u_infinity = pow(1.0 + exp((-(AV_Fn -  3.41749999999999983e-13)) / 1.367e-15), (-1.0));
    auto AV_tau_v = 1.91 + 2.09 * pow(1.0 + exp((-(AV_Fn -  3.41749999999999983e-13)) / 1.367e-15), (-1.0));
    auto AV_v_infinity = 1.0 - pow(1.0 + exp((-(AV_Fn - 6.835e-14)) / 1.367e-15), (-1.0));
    NV_Ith_S(ydot, 16) = ((-3.0) * AV_i_NaK - (3.0 * AV_i_NaCa + AV_i_B_Na + AV_i_Na)) / (AC_V_i * AC_F);
    auto AV_B1 = (2.0 * AV_i_NaCa - (AV_i_CaP + AV_i_Ca_L + AV_i_B_Ca)) / (2.0 * AC_V_i * AC_F) + (AC_V_up * (AV_i_up_leak - AV_i_up) + AV_i_rel * AC_V_rel) / AC_V_i;
    NV_Ith_S(ydot, 0) = (-(AV_i_Na + AV_i_K1 + AV_i_to + AV_i_Kur + AV_i_Kr + AV_i_Ks + AV_i_B_Na + AV_i_B_Ca + AV_i_NaK + AV_i_CaP + AV_i_NaCa + AV_i_Ca_L)) / AC_Cm;
    NV_Ith_S(ydot, 13) = (AV_u_infinity - NV_Ith_S(y, 13)) / AC_tau_u;
    NV_Ith_S(ydot, 14) = (AV_v_infinity - NV_Ith_S(y, 14)) / AV_tau_v;
    NV_Ith_S(ydot, 17) = AV_B1 / AV_B2;
    
    return 0;
}

/* Set initial values */
static void
default_initial_values(N_Vector y)
{
    NV_Ith_S(y, 0) = -81.18;
    NV_Ith_S(y, 1) = 0.002908;
    NV_Ith_S(y, 2) = 0.9649;
    NV_Ith_S(y, 3) = 0.9775;
    NV_Ith_S(y, 4) = 0.03043;
    NV_Ith_S(y, 5) = 0.9992;
    NV_Ith_S(y, 6) = 0.004966;
    NV_Ith_S(y, 7) = 0.9986;
    NV_Ith_S(y, 8) = 3.296e-05;
    NV_Ith_S(y, 9) = 0.01869;
    NV_Ith_S(y, 10) = 0.0001367;
    NV_Ith_S(y, 11) = 0.9996;
    NV_Ith_S(y, 12) = 0.7755;
    NV_Ith_S(y, 13) = 2.35e-112;
    NV_Ith_S(y, 14) = 1.0;
    NV_Ith_S(y, 15) = 0.9992;
    NV_Ith_S(y, 16) = 11.17;
    NV_Ith_S(y, 17) = 0.0001013;
    NV_Ith_S(y, 18) = 139.0;
    NV_Ith_S(y, 19) = 1.488;
    NV_Ith_S(y, 20) = 1.488;

}