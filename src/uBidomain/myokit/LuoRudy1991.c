/*
 * Luo-Rudy model (1991)
 * Generated on 2022-09-22 11:00:59
 * 
 * The well-known Luo Rudy model describes the mammalian ventricular 
 * action potential. The myokit implementation given below is based 
 * on an updated version by Leonid Livshitz and Yoram Rudy (2006).
*/

#define N_STATE 8

/* Declare intermediary, temporary and system variables */
static realtype AV_i_ion;
static realtype AV_i_stim;
static realtype AC_stim_amplitude;
static realtype AC_ENa;
static realtype AV_a;
static realtype AV_ina_m_alpha;
static realtype AV_ina_m_beta;
static realtype AV_ina_h_alpha;
static realtype AV_ina_h_beta;
static realtype AV_ina_j_alpha;
static realtype AV_ina_j_beta;
static realtype AC_gNa;
static realtype AV_INa;
static realtype AC_PNa_K;
static realtype AC_gK;
static realtype AC_ik_E;
static realtype AV_xi;
static realtype AV_ik_x_alpha;
static realtype AV_ik_x_beta;
static realtype AV_IK;
static realtype AC_gKp;
static realtype AV_IKp;
static realtype AV_ica_E;
static realtype AV_ica_d_alpha;
static realtype AV_ica_d_beta;
static realtype AV_ica_f_alpha;
static realtype AV_ica_f_beta;
static realtype AC_gCa;
static realtype AV_ICa;
static realtype AC_ik1_E;
static realtype AV_gK1;
static realtype AV_ik1_gK1_alpha;
static realtype AV_ik1_gK1_beta;
static realtype AV_IK1;
static realtype AC_gb;
static realtype AV_Ib;
static realtype AC_K_o;
static realtype AC_K_i;
static realtype AC_Na_o;
static realtype AC_Na_i;
static realtype AC_Ca_o;
static realtype AC_RTF;
static realtype AC_R;
static realtype AC_T;
static realtype AC_F;

/* Set values of constants */
static void
updateConstants(void)
{
    /* cell */
    AC_Ca_o = 1.8;
    AC_K_i = 145.0;
    AC_K_o = 5.4;
    AC_Na_i = 10.0;
    AC_Na_o = 140.0;
    AC_F = 96500.0;
    AC_R = 8314.0;
    AC_T = 310.0;
    AC_RTF = AC_R * AC_T / AC_F;
    
    /* ib */
    AC_gb = 0.03921;
    
    /* ikp */
    AC_gKp = 0.0183;
    
    /* ica */
    AC_gCa = 0.09;
    
    /* ik */
    AC_PNa_K = 0.01833;
    AC_gK = 0.282 * sqrt(AC_K_o / 5.4);
    AC_ik_E = AC_RTF * log((AC_K_o + AC_PNa_K * AC_Na_o) / (AC_K_i + AC_PNa_K * AC_Na_i));
    
    /* ik1 */
    AC_ik1_E = AC_RTF * log(AC_K_o / AC_K_i);
    
    /* ina */
    AC_ENa = AC_RTF * log(AC_Na_o / AC_Na_i);
    AC_gNa = 16.0;
    
}

/* Right-hand-side function of the model ODE */
static int rhs(realtype t, const N_Vector y, N_Vector ydot, void *f_data)
{
    /* ib */
    auto AV_Ib = AC_gb * (NV_Ith_S(y, 0) + 59.87);
    
    /* ikp */
    auto AV_IKp = AC_gKp * (NV_Ith_S(y, 0) + 87.8789) / (1.0 + exp((7.488 - NV_Ith_S(y, 0)) / 5.98));
    
    /* ica */
    auto AV_ica_E = 7.7 - 13.0287 * log(NV_Ith_S(y, 7) / AC_Ca_o);
    auto AV_ica_d_alpha = 0.095 * exp((-0.01) * (NV_Ith_S(y, 0) - 5.0)) / (1.0 + exp((-0.072) * (NV_Ith_S(y, 0) - 5.0)));
    auto AV_ica_d_beta = 0.07 * exp((-0.017) * (NV_Ith_S(y, 0) + 44.0)) / (1.0 + exp(0.05 * (NV_Ith_S(y, 0) + 44.0)));
    NV_Ith_S(ydot, 4) = AV_ica_d_alpha * (1.0 - NV_Ith_S(y, 4)) - AV_ica_d_beta * NV_Ith_S(y, 4);
    auto AV_ica_f_alpha = 0.012 * exp((-0.008) * (NV_Ith_S(y, 0) + 28.0)) / (1.0 + exp(0.15 * (NV_Ith_S(y, 0) + 28.0)));
    auto AV_ica_f_beta = 0.0065 * exp((-0.02) * (NV_Ith_S(y, 0) + 30.0)) / (1.0 + exp((-0.2) * (NV_Ith_S(y, 0) + 30.0)));
    NV_Ith_S(ydot, 5) = AV_ica_f_alpha * (1.0 - NV_Ith_S(y, 5)) - AV_ica_f_beta * NV_Ith_S(y, 5);
    auto AV_ICa = AC_gCa * NV_Ith_S(y, 4) * NV_Ith_S(y, 5) * (NV_Ith_S(y, 0) - AV_ica_E);
    NV_Ith_S(ydot, 7) = (-0.0001) * AV_ICa + 0.07 * (0.0001 - NV_Ith_S(y, 7));
    
    /* ik */
    auto AV_xi = (NV_Ith_S(y, 0) < (-100.0)).select(1.0, 
            (NV_Ith_S(y, 0) == (-77.0)).select(2.837 * 0.04 / exp(0.04 * (NV_Ith_S(y, 0) + 35.0)), 2.837 * (exp(0.04 * (NV_Ith_S(y, 0) + 77.0)) - 1.0) / ((NV_Ith_S(y, 0) + 77.0) * exp(0.04 * (NV_Ith_S(y, 0) + 35.0)))));
    auto AV_ik_x_alpha = 0.0005 * exp(0.083 * (NV_Ith_S(y, 0) + 50.0)) / (1.0 + exp(0.057 * (NV_Ith_S(y, 0) + 50.0)));
    auto AV_ik_x_beta = 0.0013 * exp((-0.06) * (NV_Ith_S(y, 0) + 20.0)) / (1.0 + exp((-0.04) * (NV_Ith_S(y, 0) + 20.0)));
    NV_Ith_S(ydot, 6) = AV_ik_x_alpha * (1.0 - NV_Ith_S(y, 6)) - AV_ik_x_beta * NV_Ith_S(y, 6);
    auto AV_IK = AC_gK * AV_xi * NV_Ith_S(y, 6) * (NV_Ith_S(y, 0) - AC_ik_E);
    
    /* ik1 */
    auto AV_ik1_gK1_alpha = 1.02 / (1.0 + exp(0.2385 * (NV_Ith_S(y, 0) - AC_ik1_E - 59.215)));
    auto AV_ik1_gK1_beta = (0.49124 * exp(0.08032 * (NV_Ith_S(y, 0) - AC_ik1_E + 5.476)) + exp(0.06175 * (NV_Ith_S(y, 0) - AC_ik1_E - 594.31))) / (1.0 + exp((-0.5143) * (NV_Ith_S(y, 0) - AC_ik1_E + 4.753)));
    auto AV_gK1 = 0.6047 * sqrt(AC_K_o / 5.4) * AV_ik1_gK1_alpha / (AV_ik1_gK1_alpha + AV_ik1_gK1_beta);
    auto AV_IK1 = AV_gK1 * (NV_Ith_S(y, 0) - AC_ik1_E);
    
    /* ina */
    auto AV_a = 1.0 - 1.0 / (1.0 + exp((-(NV_Ith_S(y, 0) + 40.0)) / 0.24));
    auto AV_ina_m_alpha = (NV_Ith_S(y, 0) == (-47.13)).select((-3.2), 0.32 * (NV_Ith_S(y, 0) + 47.13) / (1.0 - exp((-0.1) * (NV_Ith_S(y, 0) + 47.13))));
    auto AV_ina_m_beta = 0.08 * exp((-NV_Ith_S(y, 0)) / 11.0);
    NV_Ith_S(ydot, 1) = AV_ina_m_alpha * (1.0 - NV_Ith_S(y, 1)) - AV_ina_m_beta * NV_Ith_S(y, 1);
    auto AV_INa = AC_gNa * pow(NV_Ith_S(y, 1), 3.0) * NV_Ith_S(y, 2) * NV_Ith_S(y, 3) * (NV_Ith_S(y, 0) - AC_ENa);
    auto AV_ina_h_alpha = AV_a * 0.135 * exp((80.0 + NV_Ith_S(y, 0)) / (-6.8));
    auto AV_ina_h_beta = AV_a * (3.56 * exp(0.079 * NV_Ith_S(y, 0)) + 310000.0 * exp(0.35 * NV_Ith_S(y, 0))) + (1.0 - AV_a) / (0.13 * (1.0 + exp((NV_Ith_S(y, 0) + 10.66) / (-11.1))));
    NV_Ith_S(ydot, 2) = AV_ina_h_alpha * (1.0 - NV_Ith_S(y, 2)) - AV_ina_h_beta * NV_Ith_S(y, 2);
    auto AV_ina_j_alpha = AV_a * ((-127140.0) * exp(0.2444 * NV_Ith_S(y, 0)) - 3.474e-05 * exp((-0.04391) * NV_Ith_S(y, 0))) * (NV_Ith_S(y, 0) + 37.78) / (1.0 + exp(0.311 * (NV_Ith_S(y, 0) + 79.23)));
    auto AV_ina_j_beta = AV_a * (0.1212 * exp((-0.01052) * NV_Ith_S(y, 0)) / (1.0 + exp((-0.1378) * (NV_Ith_S(y, 0) + 40.14)))) + (1.0 - AV_a) * (0.3 * exp((-2.535e-07) * NV_Ith_S(y, 0)) / (1.0 + exp((-0.1) * (NV_Ith_S(y, 0) + 32.0))));
    NV_Ith_S(ydot, 3) = AV_ina_j_alpha * (1.0 - NV_Ith_S(y, 3)) - AV_ina_j_beta * NV_Ith_S(y, 3);
    
    /* membrane */
    NV_Ith_S(ydot, 0) = -(AV_INa + AV_IK + AV_Ib + AV_IKp + AV_IK1 + AV_ICa);
    
    return 0;
}

/* Set initial values */
static void
default_initial_values(N_Vector y)
{
    NV_Ith_S(y, 0) = -84.5286;
    NV_Ith_S(y, 1) = 0.0017;
    NV_Ith_S(y, 2) = 0.9832;
    NV_Ith_S(y, 3) = 0.995484;
    NV_Ith_S(y, 4) = 3e-06;
    NV_Ith_S(y, 5) = 1.0;
    NV_Ith_S(y, 6) = 0.0057;
    NV_Ith_S(y, 7) = 0.0002;

}