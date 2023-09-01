/*
HodgkinHuxley1952d_Modern
Generated on 2022-09-21 15:48:49
*/

#define N_STATE 4

/* Declare intermediary, temporary and system variables */
static realtype AC_C;
static realtype AC_A;
static realtype AV_i_stim;
static realtype AC_Ek;
static realtype AC_potassium_g_max;
static realtype AV_potassium_i;
static realtype AV_potassium_n_a;
static realtype AV_potassium_n_b;
static realtype AC_ENa;
static realtype AC_sodium_g_max;
static realtype AV_sodium_i;
static realtype AV_sodium_m_a;
static realtype AV_sodium_m_b;
static realtype AV_sodium_h_a;
static realtype AV_sodium_h_b;
static realtype AC_Eleak;
static realtype AC_leak_g_max;
static realtype AV_leak_i;

/* Set values of constants */
static void
updateConstants(void)
{
    /* leak */
    AC_Eleak = (-64.387);
    AC_leak_g_max = 0.3;
    
    /* potassium */
    AC_Ek = (-87.0);
    AC_potassium_g_max = 36.0;
    
    /* sodium */
    AC_ENa = 40.0;
    AC_sodium_g_max = 120.0;
    
    /* membrane */
    AC_C = 1.0;    
}

/* Right-hand-side function of the model ODE */
static int rhs(realtype t, const N_Vector y, N_Vector ydot, void *f_data)
{
    
    /* leak */
    auto AV_leak_i = AC_leak_g_max * (NV_Ith_S(y, 0) - AC_Eleak);
    
    /* potassium */
    auto AV_potassium_n_a = 0.01 * ((-NV_Ith_S(y, 0)) - 65.0) / (exp(((-NV_Ith_S(y, 0)) - 65.0) / 10.0) - 1.0);
    auto AV_potassium_n_b = 0.125 * exp(((-NV_Ith_S(y, 0)) - 75.0) / 80.0);
    NV_Ith_S(ydot, 1) = AV_potassium_n_a * (1.0 - NV_Ith_S(y, 1)) - AV_potassium_n_b * NV_Ith_S(y, 1);
    auto AV_potassium_i = AC_potassium_g_max * pow(NV_Ith_S(y, 1), 4.0) * (NV_Ith_S(y, 0) - AC_Ek);
    
    /* sodium */
    auto AV_sodium_h_a = 0.07 * exp(((-NV_Ith_S(y, 0)) - 75.0) / 20.0);
    auto AV_sodium_h_b = 1.0 / (exp(((-NV_Ith_S(y, 0)) - 45.0) / 10.0) + 1.0);
    NV_Ith_S(ydot, 3) = AV_sodium_h_a * (1.0 - NV_Ith_S(y, 3)) - AV_sodium_h_b * NV_Ith_S(y, 3);
    auto AV_sodium_m_a = 0.1 * ((-NV_Ith_S(y, 0)) - 50.0) / (exp(((-NV_Ith_S(y, 0)) - 50.0) / 10.0) - 1.0);
    auto AV_sodium_m_b = 4.0 * exp(((-NV_Ith_S(y, 0)) - 75.0) / 18.0);
    NV_Ith_S(ydot, 2) = AV_sodium_m_a * (1.0 - NV_Ith_S(y, 2)) - AV_sodium_m_b * NV_Ith_S(y, 2);
    auto AV_sodium_i = AC_sodium_g_max * pow(NV_Ith_S(y, 2), 3.0) * NV_Ith_S(y, 3) * (NV_Ith_S(y, 0) - AC_ENa);
    
    /* membrane */
    NV_Ith_S(ydot, 0) = -(AV_sodium_i + AV_potassium_i + AV_leak_i);
    

    return 0;
}

/* Set initial values */
static void
default_initial_values(N_Vector y)
{
    NV_Ith_S(y, 0) = -75.0;
    NV_Ith_S(y, 1) = 0.317;
    NV_Ith_S(y, 2) = 0.05;
    NV_Ith_S(y, 3) = 0.595;
}