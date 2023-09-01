/*
 * ohara-2011
 * Generated on 2022-09-21 16:15:31
 * The 2011 O'Hara-Rudy dynamical model simulates epicardial, 
 * endocardial and mid-myocardial ventricular cells based on human data.
*/

#define N_STATE 41

/* Declare intermediary, temporary and system variables */
static realtype AC_amplitude;
static realtype AC_mode;
static realtype AC_L;
static realtype AC_cell_r;
static realtype AC_vcell;
static realtype AC_Ageo;
static realtype AC_Acap;
static realtype AC_AF;
static realtype AC_vmyo;
static realtype AC_vnsr;
static realtype AC_vjsr;
static realtype AC_vss;
static realtype AC_R;
static realtype AC_T;
static realtype AC_F;
static realtype AC_RTF;
static realtype AC_FRT;
static realtype AC_FFRT;
static realtype AC_Nao;
static realtype AC_Cao;
static realtype AC_Ko;
static realtype AV_ENa;
static realtype AV_EK;
static realtype AC_PNaK;
static realtype AV_EKs;
static realtype AV_tm;
static realtype AV_ina_sm;
static realtype AV_ina_sh;
static realtype AV_thf;
static realtype AV_ths;
static realtype AC_Ahf;
static realtype AC_Ahs;
static realtype AV_ina_h;
static realtype AV_tj;
static realtype AV_sj;
static realtype AV_thsp;
static realtype AV_shsp;
static realtype AV_ina_hp;
static realtype AV_tjp;
static realtype AC_GNa;
static realtype AV_INa;
static realtype AV_inal_sm;
static realtype AC_th;
static realtype AV_inal_sh;
static realtype AC_thp;
static realtype AV_shp;
static realtype AC_GNaL;
static realtype AC_f_gnal;
static realtype AV_INaL;
static realtype AV_ta;
static realtype AV_one;
static realtype AV_two;
static realtype AV_sa;
static realtype AV_si;
static realtype AV_delta_epi;
static realtype AV_tif;
static realtype AV_tis;
static realtype AV_Aif;
static realtype AV_Ais;
static realtype AV_i;
static realtype AV_assp;
static realtype AV_dti_develop;
static realtype AV_dti_recover;
static realtype AV_tifp;
static realtype AV_tisp;
static realtype AV_ip;
static realtype AC_Gto;
static realtype AV_Ito;
static realtype AV_vf;
static realtype AV_vff;
static realtype AV_sd;
static realtype AV_td;
static realtype AV_sf;
static realtype AV_tff;
static realtype AV_tfs;
static realtype AC_Aff;
static realtype AC_Afs;
static realtype AV_ical_f;
static realtype AV_sfca;
static realtype AV_tfcaf;
static realtype AV_tfcas;
static realtype AV_Afcaf;
static realtype AV_Afcas;
static realtype AV_fca;
static realtype AC_tjca;
static realtype AV_tffp;
static realtype AV_fp;
static realtype AV_tfcafp;
static realtype AV_fcap;
static realtype AV_anca;
static realtype AC_Kmn;
static realtype AC_k2n;
static realtype AV_km2n;
static realtype AV_PhiCaL;
static realtype AV_PhiCaNa;
static realtype AV_PhiCaK;
static realtype AC_PCa;
static realtype AC_ical_PCa_base;
static realtype AC_PCap;
static realtype AC_PCaNa;
static realtype AC_PCaK;
static realtype AC_PCaNap;
static realtype AC_PCaKp;
static realtype AV_g;
static realtype AV_gp;
static realtype AV_ICaL;
static realtype AV_ICaNa;
static realtype AV_ICaK;
static realtype AV_ikr_sx;
static realtype AV_txf;
static realtype AV_txs;
static realtype AV_Axf;
static realtype AV_Axs;
static realtype AV_ikr_x;
static realtype AV_ikr_r;
static realtype AC_GKr;
static realtype AC_ikr_GKr_base;
static realtype AV_IKr;
static realtype AV_iks_sx;
static realtype AV_tx1;
static realtype AV_tx2;
static realtype AV_KsCa;
static realtype AC_GKs;
static realtype AC_iks_GKs_base;
static realtype AV_IKs;
static realtype AV_ik1_sx;
static realtype AV_tx;
static realtype AV_ik1_r;
static realtype AC_GK1;
static realtype AC_ik1_GK1_base;
static realtype AV_IK1;
static realtype AC_kna1;
static realtype AC_kna2;
static realtype AC_kna3;
static realtype AC_kasymm;
static realtype AC_wna;
static realtype AC_wca;
static realtype AC_wnaca;
static realtype AC_kcaon;
static realtype AC_kcaoff;
static realtype AC_qna;
static realtype AC_qca;
static realtype AV_hca;
static realtype AV_hna;
static realtype AV_inaca_h1;
static realtype AV_inaca_h2;
static realtype AV_inaca_h3;
static realtype AV_inaca_h4;
static realtype AV_inaca_h5;
static realtype AV_inaca_h6;
static realtype AV_inaca_h7;
static realtype AV_inaca_h8;
static realtype AV_inaca_h9;
static realtype AC_inaca_h10;
static realtype AC_inaca_h11;
static realtype AC_inaca_h12;
static realtype AC_inaca_k1;
static realtype AC_inaca_k2;
static realtype AV_inaca_k3p;
static realtype AV_inaca_k3pp;
static realtype AV_inaca_k3;
static realtype AV_inaca_k4p;
static realtype AV_inaca_k4pp;
static realtype AV_inaca_k4;
static realtype AC_inaca_k5;
static realtype AV_inaca_k6;
static realtype AV_inaca_k7;
static realtype AV_inaca_k8;
static realtype AV_inaca_x1;
static realtype AV_inaca_x2;
static realtype AV_inaca_x3;
static realtype AV_inaca_x4;
static realtype AV_inaca_E1;
static realtype AV_inaca_E2;
static realtype AV_inaca_E3;
static realtype AV_inaca_E4;
static realtype AC_inaca_KmCaAct;
static realtype AV_inaca_allo;
static realtype AV_inaca_JncxNa;
static realtype AV_inaca_JncxCa;
static realtype AC_Gncx;
static realtype AC_inaca_Gncx_base;
static realtype AV_INaCa;
static realtype AV_inacass_h1;
static realtype AV_inacass_h2;
static realtype AV_inacass_h3;
static realtype AV_inacass_h4;
static realtype AV_inacass_h5;
static realtype AV_inacass_h6;
static realtype AV_inacass_h7;
static realtype AV_inacass_h8;
static realtype AV_inacass_h9;
static realtype AC_inacass_h10;
static realtype AC_inacass_h11;
static realtype AC_inacass_h12;
static realtype AC_inacass_k1;
static realtype AC_inacass_k2;
static realtype AV_inacass_k3p;
static realtype AV_inacass_k3pp;
static realtype AV_inacass_k3;
static realtype AV_inacass_k4p;
static realtype AV_inacass_k4pp;
static realtype AV_inacass_k4;
static realtype AC_inacass_k5;
static realtype AV_inacass_k6;
static realtype AV_inacass_k7;
static realtype AV_inacass_k8;
static realtype AV_inacass_x1;
static realtype AV_inacass_x2;
static realtype AV_inacass_x3;
static realtype AV_inacass_x4;
static realtype AV_inacass_E1;
static realtype AV_inacass_E2;
static realtype AV_inacass_E3;
static realtype AV_inacass_E4;
static realtype AC_inacass_KmCaAct;
static realtype AV_inacass_allo;
static realtype AV_inacass_JncxNa;
static realtype AV_inacass_JncxCa;
static realtype AV_INaCa_ss;
static realtype AC_k1p;
static realtype AC_k1m;
static realtype AC_k2p;
static realtype AC_k2m;
static realtype AC_inak_k3p;
static realtype AC_k3m;
static realtype AC_inak_k4p;
static realtype AC_k4m;
static realtype AC_Knai0;
static realtype AC_Knao0;
static realtype AC_delta;
static realtype AV_Knai;
static realtype AV_Knao;
static realtype AC_Kki;
static realtype AC_Kko;
static realtype AC_MgADP;
static realtype AC_MgATP;
static realtype AC_Kmgatp;
static realtype AC_H;
static realtype AC_eP;
static realtype AC_Khp;
static realtype AC_Knap;
static realtype AC_Kxkur;
static realtype AV_P;
static realtype AV_a1;
static realtype AC_b1;
static realtype AC_a2;
static realtype AV_b2;
static realtype AV_a3;
static realtype AV_b3;
static realtype AC_a4;
static realtype AV_b4;
static realtype AV_inak_x1;
static realtype AV_inak_x2;
static realtype AV_inak_x3;
static realtype AV_inak_x4;
static realtype AV_inak_E1;
static realtype AV_inak_E2;
static realtype AV_inak_E3;
static realtype AV_inak_E4;
static realtype AV_JnakNa;
static realtype AV_JnakK;
static realtype AC_Pnak;
static realtype AC_inak_Pnak_base;
static realtype AV_INaK;
static realtype AV_xkb;
static realtype AC_GKb;
static realtype AV_IKb;
static realtype AC_PNab;
static realtype AV_INab;
static realtype AV_evf;
static realtype AC_PCab;
static realtype AV_ICab;
static realtype AV_evf2;
static realtype AC_GpCa;
static realtype AV_IpCa;
static realtype AC_bt;
static realtype AC_a_rel;
static realtype AV_Jrel_inf;
static realtype AV_ryr_Jrel_inf_base;
static realtype AV_tau_rel;
static realtype AV_ryr_Jrelnp_value;
static realtype AC_btp;
static realtype AC_a_relp;
static realtype AV_Jrel_infp;
static realtype AV_ryr_Jrel_infp_base;
static realtype AV_tau_relp;
static realtype AV_ryr_Jrelp_value;
static realtype AV_Jrel;
static realtype AC_serca_f;
static realtype AV_Jupnp;
static realtype AV_Jupp;
static realtype AV_Jleak;
static realtype AV_Jup;
static realtype AV_Jtr;
static realtype AV_JdiffNa;
static realtype AV_JdiffK;
static realtype AV_Jdiff;
static realtype AV_INa_tot;
static realtype AV_INa_ss_tot;
static realtype AV_IK_tot;
static realtype AV_IK_ss_tot;
static realtype AC_cmdnmax;
static realtype AC_calcium_cmdnmax_base;
static realtype AC_kmcmdn;
static realtype AC_trpnmax;
static realtype AC_kmtrpn;
static realtype AC_BSRmax;
static realtype AC_KmBSR;
static realtype AC_BSLmax;
static realtype AC_KmBSL;
static realtype AC_csqnmax;
static realtype AC_kmcsqn;
static realtype AV_ICa_tot;
static realtype AV_calcium_Cai_buff;
static realtype AV_calcium_Cai_a;
static realtype AV_calcium_Cai_b;
static realtype AV_ICa_ss_tot;
static realtype AV_calcium_Ca_ss_buff;
static realtype AV_calcium_Ca_ss_a;
static realtype AV_calcium_Ca_ss_b;
static realtype AV_calcium_Ca_jsr_buff;
static realtype AV_calcium_Ca_jsr_a;
static realtype AC_KmCaMK;
static realtype AC_aCaMK;
static realtype AC_bCaMK;
static realtype AC_CaMKo;
static realtype AC_KmCaM;
static realtype AV_CaMKb;
static realtype AV_CaMKa;
static realtype AV_camk_f;

/* Set values of constants */
static void
updateConstants(void)
{
    /* camk */
    AC_CaMKo = 0.05;
    AC_KmCaM = 0.0015;
    AC_KmCaMK = 0.15;
    AC_aCaMK = 0.05;
    AC_bCaMK = 0.00068;
    
    /* extra */
    AC_Cao = 1.8;
    AC_Ko = 5.4;
    AC_Nao = 140.0;
    
    /* ipca */
    AC_GpCa = 0.0005;
    
    /* phys */
    AC_F = 96485.0;
    AC_R = 8314.0;
    AC_T = 310.0;
    AC_FFRT = AC_F * AC_F / (AC_R * AC_T);
    AC_FRT = AC_F / (AC_R * AC_T);
    AC_RTF = AC_R * AC_T / AC_F;
    
    /* cell */
    AC_L = 0.01;
    AC_mode = 0.0;
    AC_cell_r = 0.0011;
    AC_Ageo = 2.0 * 3.14 * AC_cell_r * AC_cell_r + 2.0 * 3.14 * AC_cell_r * AC_L;
    AC_vcell = 1000.0 * 3.14 * AC_cell_r * AC_cell_r * AC_L;
    AC_Acap = 2.0 * AC_Ageo;
    AC_vjsr = 0.0048 * AC_vcell;
    AC_vmyo = 0.68 * AC_vcell;
    AC_vnsr = 0.0552 * AC_vcell;
    AC_vss = 0.02 * AC_vcell;
    AC_AF = AC_Acap / AC_F;
    
    /* icab */
    AC_PCab = 2.5e-08;
    
    /* inab */
    AC_PNab = 3.75e-10;
    
    /* nernst */
    AC_PNaK = 0.01833;
    
    /* stimulus */
    AC_amplitude = (-80.0);
    
    /* ical */
    AC_Aff = 0.6;
    AC_tjca = 75.0;
    AC_Afs = 1.0 - AC_Aff;
    AC_ical_PCa_base = 0.0001;
    AC_PCa = ((AC_mode == 0.0) ? AC_ical_PCa_base : ((AC_mode == 1.0) ? 1.2 * AC_ical_PCa_base : 2.5 * AC_ical_PCa_base));
    AC_Kmn = 0.002;
    AC_k2n = 1000.0;
    AC_PCaK = 0.0003574 * AC_PCa;
    AC_PCaNa = 0.00125 * AC_PCa;
    AC_PCap = 1.1 * AC_PCa;
    AC_PCaKp = 0.0003574 * AC_PCap;
    AC_PCaNap = 0.00125 * AC_PCap;
    
    /* ik1 */
    AC_ik1_GK1_base = 0.1908;
    AC_GK1 = ((AC_mode == 0.0) ? AC_ik1_GK1_base : ((AC_mode == 1.0) ? 1.2 * AC_ik1_GK1_base : 1.3 * AC_ik1_GK1_base));
    
    /* ikb */
    AC_GKb = ((AC_mode == 1.0) ? 0.0018 : 0.003);
    
    /* ikr */
    AC_ikr_GKr_base = 0.046;
    AC_GKr = ((AC_mode == 0.0) ? AC_ikr_GKr_base : ((AC_mode == 1.0) ? 1.3 * AC_ikr_GKr_base : 0.8 * AC_ikr_GKr_base));
    
    /* iks */
    AC_iks_GKs_base = 0.0034;
    AC_GKs = ((AC_mode == 1.0) ? 1.4 * AC_iks_GKs_base : AC_iks_GKs_base);
    
    /* ina */
    AC_Ahf = 0.99;
    AC_GNa = 75.0;
    AC_Ahs = 1.0 - AC_Ahf;
    
    /* inaca */
    AC_inaca_KmCaAct = 0.00015;
    AC_kasymm = 12.5;
    AC_kcaoff = 5000.0;
    AC_kcaon = 1500000.0;
    AC_kna1 = 15.0;
    AC_kna2 = 5.0;
    AC_kna3 = 88.12;
    AC_qca = 0.167;
    AC_qna = 0.5224;
    AC_wca = 60000.0;
    AC_wna = 60000.0;
    AC_wnaca = 5000.0;
    AC_inaca_Gncx_base = 0.0008;
    AC_Gncx = ((AC_mode == 0.0) ? AC_inaca_Gncx_base : ((AC_mode == 1.0) ? 1.1 * AC_inaca_Gncx_base : 1.4 * AC_inaca_Gncx_base));
    AC_inaca_h10 = AC_kasymm + 1.0 + AC_Nao / AC_kna1 * (1.0 + AC_Nao / AC_kna2);
    AC_inaca_k2 = AC_kcaoff;
    AC_inaca_k5 = AC_kcaoff;
    AC_inaca_h11 = AC_Nao * AC_Nao / (AC_inaca_h10 * AC_kna1 * AC_kna2);
    AC_inaca_h12 = 1.0 / AC_inaca_h10;
    AC_inaca_k1 = AC_inaca_h12 * AC_Cao * AC_kcaon;
    
    /* inak */
    AC_H = 1e-07;
    AC_Khp = 1.698e-07;
    AC_Kki = 0.5;
    AC_Kko = 0.3582;
    AC_Kmgatp = 1.698e-07;
    AC_Knai0 = 9.073;
    AC_Knao0 = 27.78;
    AC_Knap = 224.0;
    AC_Kxkur = 292.0;
    AC_MgADP = 0.05;
    AC_MgATP = 9.8;
    AC_delta = (-0.155);
    AC_eP = 4.2;
    AC_k1m = 182.4;
    AC_k1p = 949.5;
    AC_k2m = 39.4;
    AC_k2p = 687.2;
    AC_k3m = 79300.0;
    AC_inak_k3p = 1899.0;
    AC_k4m = 40.0;
    AC_inak_k4p = 639.0;
    AC_inak_Pnak_base = 30.0;
    AC_Pnak = ((AC_mode == 0.0) ? AC_inak_Pnak_base : ((AC_mode == 1.0) ? 0.9 * AC_inak_Pnak_base : 0.7 * AC_inak_Pnak_base));
    AC_a2 = AC_k2p;
    AC_a4 = AC_inak_k4p * AC_MgATP / AC_Kmgatp / (1.0 + AC_MgATP / AC_Kmgatp);
    AC_b1 = AC_k1m * AC_MgADP;
    
    /* ito */
    AC_Gto = ((AC_mode == 0.0) ? 0.02 : 0.08);
    
    /* serca */
    AC_serca_f = ((AC_mode == 1.0) ? 1.3 : 1.0);
    
    /* inacass */
    AC_inacass_KmCaAct = 0.00015;
    AC_inacass_h10 = AC_kasymm + 1.0 + AC_Nao / AC_kna1 * (1.0 + AC_Nao / AC_kna2);
    AC_inacass_k2 = AC_kcaoff;
    AC_inacass_k5 = AC_kcaoff;
    AC_inacass_h11 = AC_Nao * AC_Nao / (AC_inacass_h10 * AC_kna1 * AC_kna2);
    AC_inacass_h12 = 1.0 / AC_inacass_h10;
    AC_inacass_k1 = AC_inacass_h12 * AC_Cao * AC_kcaon;
    
    /* inal */
    AC_GNaL = 0.0075;
    AC_f_gnal = ((AC_mode == 1.0) ? 0.6 : 1.0);
    AC_th = 200.0;
    AC_thp = 3.0 * AC_th;
    
    /* ryr */
    AC_bt = 4.75;
    AC_a_rel = 0.5 * AC_bt;
    AC_btp = 1.25 * AC_bt;
    AC_a_relp = 0.5 * AC_btp;
    
    /* calcium */
    AC_BSLmax = 1.124;
    AC_BSRmax = 0.047;
    AC_KmBSL = 0.0087;
    AC_KmBSR = 0.00087;
    AC_csqnmax = 10.0;
    AC_kmcmdn = 0.00238;
    AC_kmcsqn = 0.8;
    AC_kmtrpn = 0.0005;
    AC_trpnmax = 0.07;
    AC_calcium_cmdnmax_base = 0.05;
    AC_cmdnmax = ((AC_mode == 1.0) ? 1.3 * AC_calcium_cmdnmax_base : AC_calcium_cmdnmax_base);
    
}

/* Right-hand-side function of the model ODE */
static int rhs(realtype t, const N_Vector y, N_Vector ydot, void *f_data)
{
    cout<<" begin rhs "<<endl;
    /* camk */
    auto AV_CaMKb = AC_CaMKo * (1.0 - NV_Ith_S(y, 40)) / (1.0 + AC_KmCaM / NV_Ith_S(y, 6));
    auto AV_CaMKa = AV_CaMKb + NV_Ith_S(y, 40);
    auto AV_camk_f = 1.0 / (1.0 + AC_KmCaMK / AV_CaMKa);
    NV_Ith_S(ydot, 40) = AC_aCaMK * AV_CaMKb * AV_CaMKa - AC_bCaMK * NV_Ith_S(y, 40);
    
    /* diff */
    auto AV_Jdiff = (NV_Ith_S(y, 6) - NV_Ith_S(y, 5)) / 0.2;
    auto AV_JdiffK = (NV_Ith_S(y, 4) - NV_Ith_S(y, 3)) / 2.0;
    auto AV_JdiffNa = (NV_Ith_S(y, 2) - NV_Ith_S(y, 1)) / 2.0;
    
    /* ipca */
    auto AV_IpCa = AC_GpCa * NV_Ith_S(y, 5) / (0.0005 + NV_Ith_S(y, 5));
    
    /* icab */
    auto AV_evf2 = exp(2.0 * NV_Ith_S(y, 0) * AC_FRT);
    auto AV_ICab = AC_PCab * 4.0 * NV_Ith_S(y, 0) * AC_FFRT * (NV_Ith_S(y, 5) * AV_evf2 - 0.341 * AC_Cao) / (AV_evf2 - 1.0);
    
    /* inab */
    auto AV_evf = exp(NV_Ith_S(y, 0) * AC_FRT);
    auto AV_INab = AC_PNab * NV_Ith_S(y, 0) * AC_FFRT * (NV_Ith_S(y, 1) * AV_evf - AC_Nao) / (AV_evf - 1.0);
    
    /* nernst */
    auto AV_EK = AC_RTF * log(AC_Ko / NV_Ith_S(y, 3));
    auto AV_ENa = AC_RTF * log(AC_Nao / NV_Ith_S(y, 1));
    auto AV_EKs = AC_RTF * log((AC_Ko + AC_PNaK * AC_Nao) / (NV_Ith_S(y, 3) + AC_PNaK * NV_Ith_S(y, 1)));
        
    /* ical */
    auto AV_Afcaf = 0.3 + 0.6 / (1.0 + exp((NV_Ith_S(y, 0) - 10.0) / 10.0));
    auto AV_sd = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 3.94) / (-4.23)));
    auto AV_sf = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 19.58) / 3.696));
    auto AV_td = 0.6 + 1.0 / (exp((-0.05) * (NV_Ith_S(y, 0) + 6.0)) + exp(0.09 * (NV_Ith_S(y, 0) + 14.0)));
    auto AV_tfcaf = 7.0 + 1.0 / (0.04 * exp((NV_Ith_S(y, 0) - 4.0) / (-7.0)) + 0.04 * exp((NV_Ith_S(y, 0) - 4.0) / 7.0));
    auto AV_tfcas = 100.0 + 1.0 / (0.00012 * exp(NV_Ith_S(y, 0) / (-3.0)) + 0.00012 * exp(NV_Ith_S(y, 0) / 7.0));
    auto AV_tff = 7.0 + 1.0 / (0.0045 * exp((NV_Ith_S(y, 0) + 20.0) / (-10.0)) + 0.0045 * exp((NV_Ith_S(y, 0) + 20.0) / 10.0));
    auto AV_tfs = 1000.0 + 1.0 / (3.5e-05 * exp((NV_Ith_S(y, 0) + 5.0) / (-4.0)) + 3.5e-05 * exp((NV_Ith_S(y, 0) + 5.0) / 6.0));
    auto AV_vf = NV_Ith_S(y, 0) * AC_FRT;
    auto AV_vff = NV_Ith_S(y, 0) * AC_FFRT;
    NV_Ith_S(ydot, 24) = (AV_sd - NV_Ith_S(y, 24)) / AV_td;
    NV_Ith_S(ydot, 25) = (AV_sf - NV_Ith_S(y, 25)) / AV_tff;
    NV_Ith_S(ydot, 26) = (AV_sf - NV_Ith_S(y, 26)) / AV_tfs;
    auto AV_Afcas = 1.0 - AV_Afcaf;
    auto AV_PhiCaK = 1.0 * AV_vff * (0.75 * NV_Ith_S(y, 4) * exp(1.0 * AV_vf) - 0.75 * AC_Ko) / (exp(1.0 * AV_vf) - 1.0);
    auto AV_PhiCaL = 4.0 * AV_vff * (NV_Ith_S(y, 6) * exp(2.0 * AV_vf) - 0.341 * AC_Cao) / (exp(2.0 * AV_vf) - 1.0);
    auto AV_PhiCaNa = 1.0 * AV_vff * (0.75 * NV_Ith_S(y, 2) * exp(1.0 * AV_vf) - 0.75 * AC_Nao) / (exp(1.0 * AV_vf) - 1.0);
    auto AV_sfca = AV_sf;
    auto AV_tfcafp = 2.5 * AV_tfcaf;
    auto AV_tffp = 2.5 * AV_tff;
    NV_Ith_S(ydot, 27) = (AV_sfca - NV_Ith_S(y, 27)) / AV_tfcaf;
    NV_Ith_S(ydot, 32) = (AV_sfca - NV_Ith_S(y, 32)) / AV_tfcafp;
    NV_Ith_S(ydot, 28) = (AV_sfca - NV_Ith_S(y, 28)) / AV_tfcas;
    NV_Ith_S(ydot, 31) = (AV_sf - NV_Ith_S(y, 31)) / AV_tffp;
    NV_Ith_S(ydot, 29) = (AV_sfca - NV_Ith_S(y, 29)) / AC_tjca;
    auto AV_km2n = NV_Ith_S(y, 29) * 1.0;
    auto AV_anca = 1.0 / (AC_k2n / AV_km2n + pow(1.0 + AC_Kmn / NV_Ith_S(y, 6), 4.0));
    NV_Ith_S(ydot, 30) = AV_anca * AC_k2n - NV_Ith_S(y, 30) * AV_km2n;
    auto AV_ical_f = AC_Aff * NV_Ith_S(y, 25) + AC_Afs * NV_Ith_S(y, 26);
    auto AV_fca = AV_Afcaf * NV_Ith_S(y, 27) + AV_Afcas * NV_Ith_S(y, 28);
    auto AV_fcap = AV_Afcaf * NV_Ith_S(y, 32) + AV_Afcas * NV_Ith_S(y, 28);
    auto AV_fp = AC_Aff * NV_Ith_S(y, 31) + AC_Afs * NV_Ith_S(y, 26);
    auto AV_g = NV_Ith_S(y, 24) * (AV_ical_f * (1.0 - NV_Ith_S(y, 30)) + NV_Ith_S(y, 29) * AV_fca * NV_Ith_S(y, 30));
    auto AV_gp = NV_Ith_S(y, 24) * (AV_fp * (1.0 - NV_Ith_S(y, 30)) + NV_Ith_S(y, 29) * AV_fcap * NV_Ith_S(y, 30));
    auto AV_ICaK = (1.0 - AV_camk_f) * AC_PCaK * AV_PhiCaK * AV_g + AV_camk_f * AC_PCaKp * AV_PhiCaK * AV_gp;
    auto AV_ICaL = (1.0 - AV_camk_f) * AC_PCa * AV_PhiCaL * AV_g + AV_camk_f * AC_PCap * AV_PhiCaL * AV_gp;
    auto AV_ICaNa = (1.0 - AV_camk_f) * AC_PCaNa * AV_PhiCaNa * AV_g + AV_camk_f * AC_PCaNap * AV_PhiCaNa * AV_gp;
    
    /* ik1 */
    auto AV_ik1_r = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 105.8 - 2.6 * AC_Ko) / 9.493));
    auto AV_ik1_sx = 1.0 / (1.0 + exp((-(NV_Ith_S(y, 0) + 2.5538 * AC_Ko + 144.59)) / (1.5692 * AC_Ko + 3.8115)));
    auto AV_tx = 122.2 / (exp((NV_Ith_S(y, 0) + 127.2) / (-20.36)) + exp((NV_Ith_S(y, 0) + 236.8) / 69.33));
    NV_Ith_S(ydot, 37) = (AV_ik1_sx - NV_Ith_S(y, 37)) / AV_tx;
    auto AV_IK1 = AC_GK1 * sqrt(AC_Ko) * AV_ik1_r * NV_Ith_S(y, 37) * (NV_Ith_S(y, 0) - AV_EK);
    
    /* ikb */
    auto AV_xkb = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) - 14.48) / (-18.34)));
    auto AV_IKb = AC_GKb * AV_xkb * (NV_Ith_S(y, 0) - AV_EK);
    
    /* ikr */
    auto AV_Axf = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 54.81) / 38.21));
    auto AV_ikr_r = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 55.0) / 75.0)) * 1.0 / (1.0 + exp((NV_Ith_S(y, 0) - 10.0) / 30.0));
    auto AV_ikr_sx = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 8.337) / (-6.789)));
    auto AV_txf = 12.98 + 1.0 / (0.3652 * exp((NV_Ith_S(y, 0) - 31.66) / 3.869) + 4.123e-05 * exp((NV_Ith_S(y, 0) - 47.78) / (-20.38)));
    auto AV_txs = 1.865 + 1.0 / (0.06629 * exp((NV_Ith_S(y, 0) - 34.7) / 7.355) + 1.128e-05 * exp((NV_Ith_S(y, 0) - 29.74) / (-25.94)));
    NV_Ith_S(ydot, 33) = (AV_ikr_sx - NV_Ith_S(y, 33)) / AV_txf;
    NV_Ith_S(ydot, 34) = (AV_ikr_sx - NV_Ith_S(y, 34)) / AV_txs;
    auto AV_Axs = 1.0 - AV_Axf;
    auto AV_ikr_x = AV_Axf * NV_Ith_S(y, 33) + AV_Axs * NV_Ith_S(y, 34);
    auto AV_IKr = AC_GKr * sqrt(AC_Ko / 5.4) * AV_ikr_x * AV_ikr_r * (NV_Ith_S(y, 0) - AV_EK);
    
    /* iks */
    auto AV_KsCa = 1.0 + 0.6 / (1.0 + pow(3.8e-05 / NV_Ith_S(y, 5), 1.4));
    auto AV_iks_sx = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 11.6) / (-8.932)));
    auto AV_tx1 = 817.3 + 1.0 / (0.0002326 * exp((NV_Ith_S(y, 0) + 48.28) / 17.8) + 0.001292 * exp((NV_Ith_S(y, 0) + 210.0) / (-230.0)));
    auto AV_tx2 = 1.0 / (0.01 * exp((NV_Ith_S(y, 0) - 50.0) / 20.0) + 0.0193 * exp((NV_Ith_S(y, 0) + 66.54) / (-31.0)));
    NV_Ith_S(ydot, 35) = (AV_iks_sx - NV_Ith_S(y, 35)) / AV_tx1;
    NV_Ith_S(ydot, 36) = (AV_iks_sx - NV_Ith_S(y, 36)) / AV_tx2;
    auto AV_IKs = AC_GKs * AV_KsCa * NV_Ith_S(y, 35) * NV_Ith_S(y, 36) * (NV_Ith_S(y, 0) - AV_EKs);
    
    /* ina */
    auto AV_ina_sh = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 82.9) / 6.086));
    auto AV_shsp = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 89.1) / 6.086));
    auto AV_ina_sm = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 39.57) / (-9.871)));
    auto AV_thf = 1.0 / (1.432e-05 * exp((NV_Ith_S(y, 0) + 1.196) / (-6.285)) + 6.149 * exp((NV_Ith_S(y, 0) + 0.5096) / 20.27));
    auto AV_ths = 1.0 / (0.009794 * exp((NV_Ith_S(y, 0) + 17.95) / (-28.05)) + 0.3343 * exp((NV_Ith_S(y, 0) + 5.73) / 56.66));
    auto AV_tj = 2.038 + 1.0 / (0.02136 * exp((NV_Ith_S(y, 0) + 100.6) / (-8.281)) + 0.3052 * exp((NV_Ith_S(y, 0) + 0.9941) / 38.45));
    auto AV_tm = 1.0 / (6.765 * exp((NV_Ith_S(y, 0) + 11.64) / 34.77) + 8.552 * exp((-(NV_Ith_S(y, 0) + 77.42)) / 5.955));
    NV_Ith_S(ydot, 10) = (AV_ina_sh - NV_Ith_S(y, 10)) / AV_thf;
    NV_Ith_S(ydot, 11) = (AV_ina_sh - NV_Ith_S(y, 11)) / AV_ths;
    NV_Ith_S(ydot, 9) = (AV_ina_sm - NV_Ith_S(y, 9)) / AV_tm;
    auto AV_sj = AV_ina_sh;
    auto AV_thsp = 3.0 * AV_ths;
    auto AV_tjp = 1.46 * AV_tj;
    NV_Ith_S(ydot, 13) = (AV_shsp - NV_Ith_S(y, 13)) / AV_thsp;
    NV_Ith_S(ydot, 12) = (AV_sj - NV_Ith_S(y, 12)) / AV_tj;
    NV_Ith_S(ydot, 14) = (AV_sj - NV_Ith_S(y, 14)) / AV_tjp;
    auto AV_ina_h = AC_Ahf * NV_Ith_S(y, 10) + AC_Ahs * NV_Ith_S(y, 11);
    auto AV_ina_hp = AC_Ahf * NV_Ith_S(y, 10) + AC_Ahs * NV_Ith_S(y, 13);
    auto AV_INa = AC_GNa * (NV_Ith_S(y, 0) - AV_ENa) * pow(NV_Ith_S(y, 9), 3.0) * ((1.0 - AV_camk_f) * AV_ina_h * NV_Ith_S(y, 12) + AV_camk_f * AV_ina_hp * NV_Ith_S(y, 14));
    
    /* inaca */
    auto AV_inaca_allo = 1.0 / (1.0 + pow(AC_inaca_KmCaAct / NV_Ith_S(y, 5), 2.0));
    auto AV_inaca_h4 = 1.0 + NV_Ith_S(y, 1) / AC_kna1 * (1.0 + NV_Ith_S(y, 1) / AC_kna2);
    auto AV_hca = exp(AC_qca * NV_Ith_S(y, 0) * AC_FRT);
    auto AV_hna = exp(AC_qna * NV_Ith_S(y, 0) * AC_FRT);
    auto AV_inaca_h1 = 1.0 + NV_Ith_S(y, 1) / AC_kna3 * (1.0 + AV_hna);
    auto AV_inaca_h5 = NV_Ith_S(y, 1) * NV_Ith_S(y, 1) / (AV_inaca_h4 * AC_kna1 * AC_kna2);
    auto AV_inaca_h6 = 1.0 / AV_inaca_h4;
    auto AV_inaca_h7 = 1.0 + AC_Nao / AC_kna3 * (1.0 + 1.0 / AV_hna);
    auto AV_inaca_h2 = NV_Ith_S(y, 1) * AV_hna / (AC_kna3 * AV_inaca_h1);
    auto AV_inaca_h3 = 1.0 / AV_inaca_h1;
    auto AV_inaca_h8 = AC_Nao / (AC_kna3 * AV_hna * AV_inaca_h7);
    auto AV_inaca_h9 = 1.0 / AV_inaca_h7;
    auto AV_inaca_k6 = AV_inaca_h6 * NV_Ith_S(y, 5) * AC_kcaon;
    auto AV_inaca_k3p = AV_inaca_h9 * AC_wca;
    auto AV_inaca_k3pp = AV_inaca_h8 * AC_wnaca;
    auto AV_inaca_k4p = AV_inaca_h3 * AC_wca / AV_hca;
    auto AV_inaca_k4pp = AV_inaca_h2 * AC_wnaca;
    auto AV_inaca_k7 = AV_inaca_h5 * AV_inaca_h2 * AC_wna;
    auto AV_inaca_k8 = AV_inaca_h8 * AC_inaca_h11 * AC_wna;
    auto AV_inaca_k3 = AV_inaca_k3p + AV_inaca_k3pp;
    auto AV_inaca_k4 = AV_inaca_k4p + AV_inaca_k4pp;
    auto AV_inaca_x1 = AC_inaca_k2 * AV_inaca_k4 * (AV_inaca_k7 + AV_inaca_k6) + AC_inaca_k5 * AV_inaca_k7 * (AC_inaca_k2 + AV_inaca_k3);
    auto AV_inaca_x2 = AC_inaca_k1 * AV_inaca_k7 * (AV_inaca_k4 + AC_inaca_k5) + AV_inaca_k4 * AV_inaca_k6 * (AC_inaca_k1 + AV_inaca_k8);
    auto AV_inaca_x3 = AC_inaca_k1 * AV_inaca_k3 * (AV_inaca_k7 + AV_inaca_k6) + AV_inaca_k8 * AV_inaca_k6 * (AC_inaca_k2 + AV_inaca_k3);
    auto AV_inaca_x4 = AC_inaca_k2 * AV_inaca_k8 * (AV_inaca_k4 + AC_inaca_k5) + AV_inaca_k3 * AC_inaca_k5 * (AC_inaca_k1 + AV_inaca_k8);
    auto AV_inaca_E1 = AV_inaca_x1 / (AV_inaca_x1 + AV_inaca_x2 + AV_inaca_x3 + AV_inaca_x4);
    auto AV_inaca_E2 = AV_inaca_x2 / (AV_inaca_x1 + AV_inaca_x2 + AV_inaca_x3 + AV_inaca_x4);
    auto AV_inaca_E3 = AV_inaca_x3 / (AV_inaca_x1 + AV_inaca_x2 + AV_inaca_x3 + AV_inaca_x4);
    auto AV_inaca_E4 = AV_inaca_x4 / (AV_inaca_x1 + AV_inaca_x2 + AV_inaca_x3 + AV_inaca_x4);
    auto AV_inaca_JncxCa = AV_inaca_E2 * AC_inaca_k2 - AV_inaca_E1 * AC_inaca_k1;
    auto AV_inaca_JncxNa = 3.0 * (AV_inaca_E4 * AV_inaca_k7 - AV_inaca_E1 * AV_inaca_k8) + AV_inaca_E3 * AV_inaca_k4pp - AV_inaca_E2 * AV_inaca_k3pp;
    auto AV_INaCa = 0.8 * AC_Gncx * AV_inaca_allo * (AV_inaca_JncxNa + 2.0 * AV_inaca_JncxCa);
    
    /* inak */
    auto AV_Knai = AC_Knai0 * exp(AC_delta * NV_Ith_S(y, 0) * AC_FRT / 3.0);
    auto AV_Knao = AC_Knao0 * exp((1.0 - AC_delta) * NV_Ith_S(y, 0) * AC_FRT / 3.0);
    auto AV_P = AC_eP / (1.0 + AC_H / AC_Khp + NV_Ith_S(y, 1) / AC_Knap + NV_Ith_S(y, 3) / AC_Kxkur);
    auto AV_a1 = AC_k1p * pow(NV_Ith_S(y, 1) / AV_Knai, 3.0) / (pow(1.0 + NV_Ith_S(y, 1) / AV_Knai, 3.0) + pow(1.0 + NV_Ith_S(y, 3) / AC_Kki, 2.0) - 1.0);
    auto AV_a3 = AC_inak_k3p * pow(AC_Ko / AC_Kko, 2.0) / (pow(1.0 + AC_Nao / AV_Knao, 3.0) + pow(1.0 + AC_Ko / AC_Kko, 2.0) - 1.0);
    auto AV_b2 = AC_k2m * pow(AC_Nao / AV_Knao, 3.0) / (pow(1.0 + AC_Nao / AV_Knao, 3.0) + pow(1.0 + AC_Ko / AC_Kko, 2.0) - 1.0);
    auto AV_b3 = AC_k3m * AV_P * AC_H / (1.0 + AC_MgATP / AC_Kmgatp);
    auto AV_b4 = AC_k4m * pow(NV_Ith_S(y, 3) / AC_Kki, 2.0) / (pow(1.0 + NV_Ith_S(y, 1) / AV_Knai, 3.0) + pow(1.0 + NV_Ith_S(y, 3) / AC_Kki, 2.0) - 1.0);
    auto AV_inak_x1 = AC_a4 * AV_a1 * AC_a2 + AV_b2 * AV_b4 * AV_b3 + AC_a2 * AV_b4 * AV_b3 + AV_b3 * AV_a1 * AC_a2;
    auto AV_inak_x2 = AV_b2 * AC_b1 * AV_b4 + AV_a1 * AC_a2 * AV_a3 + AV_a3 * AC_b1 * AV_b4 + AC_a2 * AV_a3 * AV_b4;
    auto AV_inak_x3 = AC_a2 * AV_a3 * AC_a4 + AV_b3 * AV_b2 * AC_b1 + AV_b2 * AC_b1 * AC_a4 + AV_a3 * AC_a4 * AC_b1;
    auto AV_inak_x4 = AV_b4 * AV_b3 * AV_b2 + AV_a3 * AC_a4 * AV_a1 + AV_b2 * AC_a4 * AV_a1 + AV_b3 * AV_b2 * AV_a1;
    auto AV_inak_E1 = AV_inak_x1 / (AV_inak_x1 + AV_inak_x2 + AV_inak_x3 + AV_inak_x4);
    auto AV_inak_E2 = AV_inak_x2 / (AV_inak_x1 + AV_inak_x2 + AV_inak_x3 + AV_inak_x4);
    auto AV_inak_E3 = AV_inak_x3 / (AV_inak_x1 + AV_inak_x2 + AV_inak_x3 + AV_inak_x4);
    auto AV_inak_E4 = AV_inak_x4 / (AV_inak_x1 + AV_inak_x2 + AV_inak_x3 + AV_inak_x4);
    auto AV_JnakK = 2.0 * (AV_inak_E4 * AC_b1 - AV_inak_E3 * AV_a1);
    auto AV_JnakNa = 3.0 * (AV_inak_E1 * AV_a3 - AV_inak_E2 * AV_b3);
    auto AV_INaK = AC_Pnak * (AV_JnakNa + AV_JnakK);
    
    /* ito */
    auto AV_Aif = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) - 213.6) / 151.2));
//    auto AV_delta_epi = ((AC_mode == 1.0) ? 1.0 - 0.95 / (1.0 + exp((NV_Ith_S(y, 0) + 70.0) / 5.0)) : 1.0);
    Array AV_delta_epi = Array::Ones(NV_Ith_S(y, 0).size());
    if(AC_mode == 1.0)
        AV_delta_epi += - 0.95 / (1.0 + exp((NV_Ith_S(y, 0) + 70.0) / 5.0));
                
    auto AV_dti_develop = 1.354 + 0.0001 / (exp((NV_Ith_S(y, 0) - 167.4) / 15.89) + exp((NV_Ith_S(y, 0) - 12.23) / (-0.2154)));
    auto AV_dti_recover = 1.0 - 0.5 / (1.0 + exp((NV_Ith_S(y, 0) + 70.0) / 20.0));
    auto AV_sa = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) - 14.34) / (-14.82)));
    auto AV_si = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 43.94) / 5.711));
    auto AV_Ais = 1.0 - AV_Aif;
    auto AV_one = 1.0 / (1.2089 * (1.0 + exp((NV_Ith_S(y, 0) - 18.4099) / (-29.3814))));
    auto AV_two = 3.5 / (1.0 + exp((NV_Ith_S(y, 0) + 100.0) / 29.3814));
    auto AV_ta = 1.0515 / (AV_one + AV_two);
    auto AV_tif = (4.562 + 1.0 / (0.3933 * exp((NV_Ith_S(y, 0) + 100.0) / (-100.0)) + 0.08004 * exp((NV_Ith_S(y, 0) + 50.0) / 16.59))) * AV_delta_epi;
    auto AV_tis = (23.62 + 1.0 / (0.001416 * exp((NV_Ith_S(y, 0) + 96.52) / (-59.05)) + 1.78e-08 * exp((NV_Ith_S(y, 0) + 114.1) / 8.079))) * AV_delta_epi;
    NV_Ith_S(ydot, 18) = (AV_sa - NV_Ith_S(y, 18)) / AV_ta;
    auto AV_assp = 1.0 / (1.0 + exp((-(NV_Ith_S(y, 0) - 24.34)) / 14.82));
    NV_Ith_S(ydot, 21) = (AV_assp - NV_Ith_S(y, 21)) / AV_ta;
    NV_Ith_S(ydot, 19) = (AV_si - NV_Ith_S(y, 19)) / AV_tif;
    NV_Ith_S(ydot, 20) = (AV_si - NV_Ith_S(y, 20)) / AV_tis;
    auto AV_i = AV_Aif * NV_Ith_S(y, 19) + AV_Ais * NV_Ith_S(y, 20);
    auto AV_ip = AV_Aif * NV_Ith_S(y, 22) + AV_Ais * NV_Ith_S(y, 23);
    auto AV_tifp = AV_dti_develop * AV_dti_recover * AV_tif;
    auto AV_tisp = AV_dti_develop * AV_dti_recover * AV_tis;
    NV_Ith_S(ydot, 22) = (AV_si - NV_Ith_S(y, 22)) / AV_tifp;
    NV_Ith_S(ydot, 23) = (AV_si - NV_Ith_S(y, 23)) / AV_tisp;
    auto AV_Ito = AC_Gto * (NV_Ith_S(y, 0) - AV_EK) * ((1.0 - AV_camk_f) * NV_Ith_S(y, 18) * AV_i + AV_camk_f * NV_Ith_S(y, 21) * AV_ip);
    
    /* serca */
    auto AV_Jleak = 0.0039375 * NV_Ith_S(y, 7) / 15.0;
    auto AV_Jtr = (NV_Ith_S(y, 7) - NV_Ith_S(y, 8)) / 100.0;
    auto AV_Jupnp = AC_serca_f * (0.004375 * NV_Ith_S(y, 5) / (NV_Ith_S(y, 5) + 0.00092));
    auto AV_Jupp = AC_serca_f * (2.75 * 0.004375 * NV_Ith_S(y, 5) / (NV_Ith_S(y, 5) + 0.00092 - 0.00017));
    auto AV_Jup = (1.0 - AV_camk_f) * AV_Jupnp + AV_camk_f * AV_Jupp - AV_Jleak;
    
    /* inacass */
    auto AV_inacass_h1 = 1.0 + NV_Ith_S(y, 2) / AC_kna3 * (1.0 + AV_hna);
    auto AV_inacass_h4 = 1.0 + NV_Ith_S(y, 2) / AC_kna1 * (1.0 + NV_Ith_S(y, 2) / AC_kna2);
    auto AV_inacass_h7 = 1.0 + AC_Nao / AC_kna3 * (1.0 + 1.0 / AV_hna);
    auto AV_inacass_allo = 1.0 / (1.0 + pow(AC_inacass_KmCaAct / NV_Ith_S(y, 6), 2.0));
    auto AV_inacass_h2 = NV_Ith_S(y, 2) * AV_hna / (AC_kna3 * AV_inacass_h1);
    auto AV_inacass_h3 = 1.0 / AV_inacass_h1;
    auto AV_inacass_h5 = NV_Ith_S(y, 2) * NV_Ith_S(y, 2) / (AV_inacass_h4 * AC_kna1 * AC_kna2);
    auto AV_inacass_h6 = 1.0 / AV_inacass_h4;
    auto AV_inacass_h8 = AC_Nao / (AC_kna3 * AV_hna * AV_inacass_h7);
    auto AV_inacass_h9 = 1.0 / AV_inacass_h7;
    auto AV_inacass_k3p = AV_inacass_h9 * AC_wca;
    auto AV_inacass_k3pp = AV_inacass_h8 * AC_wnaca;
    auto AV_inacass_k4p = AV_inacass_h3 * AC_wca / AV_hca;
    auto AV_inacass_k4pp = AV_inacass_h2 * AC_wnaca;
    auto AV_inacass_k6 = AV_inacass_h6 * NV_Ith_S(y, 6) * AC_kcaon;
    auto AV_inacass_k7 = AV_inacass_h5 * AV_inacass_h2 * AC_wna;
    auto AV_inacass_k8 = AV_inacass_h8 * AC_inacass_h11 * AC_wna;
    auto AV_inacass_k3 = AV_inacass_k3p + AV_inacass_k3pp;
    auto AV_inacass_k4 = AV_inacass_k4p + AV_inacass_k4pp;
    auto AV_inacass_x1 = AC_inacass_k2 * AV_inacass_k4 * (AV_inacass_k7 + AV_inacass_k6) + AC_inacass_k5 * AV_inacass_k7 * (AC_inacass_k2 + AV_inacass_k3);
    auto AV_inacass_x2 = AC_inacass_k1 * AV_inacass_k7 * (AV_inacass_k4 + AC_inacass_k5) + AV_inacass_k4 * AV_inacass_k6 * (AC_inacass_k1 + AV_inacass_k8);
    auto AV_inacass_x3 = AC_inacass_k1 * AV_inacass_k3 * (AV_inacass_k7 + AV_inacass_k6) + AV_inacass_k8 * AV_inacass_k6 * (AC_inacass_k2 + AV_inacass_k3);
    auto AV_inacass_x4 = AC_inacass_k2 * AV_inacass_k8 * (AV_inacass_k4 + AC_inacass_k5) + AV_inacass_k3 * AC_inacass_k5 * (AC_inacass_k1 + AV_inacass_k8);
    auto AV_inacass_E1 = AV_inacass_x1 / (AV_inacass_x1 + AV_inacass_x2 + AV_inacass_x3 + AV_inacass_x4);
    auto AV_inacass_E2 = AV_inacass_x2 / (AV_inacass_x1 + AV_inacass_x2 + AV_inacass_x3 + AV_inacass_x4);
    auto AV_inacass_E3 = AV_inacass_x3 / (AV_inacass_x1 + AV_inacass_x2 + AV_inacass_x3 + AV_inacass_x4);
    auto AV_inacass_E4 = AV_inacass_x4 / (AV_inacass_x1 + AV_inacass_x2 + AV_inacass_x3 + AV_inacass_x4);
    auto AV_inacass_JncxCa = AV_inacass_E2 * AC_inacass_k2 - AV_inacass_E1 * AC_inacass_k1;
    auto AV_inacass_JncxNa = 3.0 * (AV_inacass_E4 * AV_inacass_k7 - AV_inacass_E1 * AV_inacass_k8) + AV_inacass_E3 * AV_inacass_k4pp - AV_inacass_E2 * AV_inacass_k3pp;
    auto AV_INaCa_ss = 0.2 * AC_Gncx * AV_inacass_allo * (AV_inacass_JncxNa + 2.0 * AV_inacass_JncxCa);
    
    /* inal */
    auto AV_inal_sh = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 87.61) / 7.488));
    auto AV_shp = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 93.81) / 7.488));
    auto AV_inal_sm = 1.0 / (1.0 + exp((NV_Ith_S(y, 0) + 42.85) / (-5.264)));
    NV_Ith_S(ydot, 16) = (AV_inal_sh - NV_Ith_S(y, 16)) / AC_th;
    NV_Ith_S(ydot, 15) = (AV_inal_sm - NV_Ith_S(y, 15)) / AV_tm;
    auto AV_INaL = AC_f_gnal * AC_GNaL * (NV_Ith_S(y, 0) - AV_ENa) * NV_Ith_S(y, 15) * ((1.0 - AV_camk_f) * NV_Ith_S(y, 16) + AV_camk_f * NV_Ith_S(y, 17));
    NV_Ith_S(ydot, 17) = (AV_shp - NV_Ith_S(y, 17)) / AC_thp;
    
    /* potassium */
    auto AV_IK_ss_tot = AV_ICaK;
    auto AV_IK_tot = AV_Ito + AV_IKr + AV_IKs + AV_IK1 + AV_IKb - 2.0 * AV_INaK;
    NV_Ith_S(ydot, 4) = (-AV_IK_ss_tot) * AC_AF / AC_vss - AV_JdiffK;
    NV_Ith_S(ydot, 3) = (-(AV_IK_tot)) * AC_AF / AC_vmyo + AV_JdiffK * AC_vss / AC_vmyo;
    
    /* ryr */
    auto AV_Jrel = (1.0 - AV_camk_f) * NV_Ith_S(y, 38) + AV_camk_f * NV_Ith_S(y, 39);
    auto AV_ryr_Jrel_inf_base = AC_a_rel * (-AV_ICaL) / (1.0 + pow(1.5 / NV_Ith_S(y, 8), 8.0));
//    auto AV_Jrel_inf = ((AC_mode == 2.0) ? 1.7 * AV_ryr_Jrel_inf_base : AV_ryr_Jrel_inf_base);
    Array AV_Jrel_inf = AV_ryr_Jrel_inf_base;
    if(AC_mode == 2.0)
        AV_Jrel_inf *= 1.7;
    auto AV_ryr_Jrelnp_value = AC_bt / (1.0 + 0.0123 / NV_Ith_S(y, 8));
    auto AV_tau_rel = (AV_ryr_Jrelnp_value < 0.001).select(0.001, AV_ryr_Jrelnp_value);
    NV_Ith_S(ydot, 38) = (AV_Jrel_inf - NV_Ith_S(y, 38)) / AV_tau_rel;
    auto AV_ryr_Jrel_infp_base = AC_a_relp * (-AV_ICaL) / (1.0 + pow(1.5 / NV_Ith_S(y, 8), 8.0));
//    auto AV_Jrel_infp = ((AC_mode == 2.0) ? 1.7 * AV_ryr_Jrel_infp_base : AV_ryr_Jrel_infp_base);
    Array AV_Jrel_infp = AV_ryr_Jrel_infp_base;
    if(AC_mode == 2.0)
        AV_Jrel_infp *= 1.7;
    auto AV_ryr_Jrelp_value = AC_btp / (1.0 + 0.0123 / NV_Ith_S(y, 8));
    auto AV_tau_relp = (AV_ryr_Jrelp_value < 0.001).select(0.001, AV_ryr_Jrelp_value);
    NV_Ith_S(ydot, 39) = (AV_Jrel_infp - NV_Ith_S(y, 39)) / AV_tau_relp;
    
    /* calcium */
    auto AV_ICa_ss_tot = AV_ICaL - 2.0 * AV_INaCa_ss;
    auto AV_ICa_tot = AV_IpCa + AV_ICab - 2.0 * AV_INaCa;
    NV_Ith_S(ydot, 7) = AV_Jup - AV_Jtr * AC_vjsr / AC_vnsr;
    auto AV_calcium_Ca_jsr_a = AC_kmcsqn + NV_Ith_S(y, 8);
    auto AV_calcium_Ca_jsr_buff = 1.0 / (1.0 + AC_csqnmax * AC_kmcsqn / (AV_calcium_Ca_jsr_a * AV_calcium_Ca_jsr_a));
    NV_Ith_S(ydot, 8) = AV_calcium_Ca_jsr_buff * (AV_Jtr - AV_Jrel);
    auto AV_calcium_Ca_ss_a = AC_KmBSR + NV_Ith_S(y, 6);
    auto AV_calcium_Ca_ss_b = AC_KmBSL + NV_Ith_S(y, 6);
    auto AV_calcium_Ca_ss_buff = 1.0 / (1.0 + AC_BSRmax * AC_KmBSR / (AV_calcium_Ca_ss_a * AV_calcium_Ca_ss_a) + AC_BSLmax * AC_KmBSL / (AV_calcium_Ca_ss_b * AV_calcium_Ca_ss_b));
    NV_Ith_S(ydot, 6) = AV_calcium_Ca_ss_buff * ((-AV_ICa_ss_tot) * AC_AF / (2.0 * AC_vss) + AV_Jrel * AC_vjsr / AC_vss - AV_Jdiff);
    auto AV_calcium_Cai_a = AC_kmcmdn + NV_Ith_S(y, 5);
    auto AV_calcium_Cai_b = AC_kmtrpn + NV_Ith_S(y, 5);
    auto AV_calcium_Cai_buff = 1.0 / (1.0 + AC_cmdnmax * AC_kmcmdn / (AV_calcium_Cai_a * AV_calcium_Cai_a) + AC_trpnmax * AC_kmtrpn / (AV_calcium_Cai_b * AV_calcium_Cai_b));
    NV_Ith_S(ydot, 5) = AV_calcium_Cai_buff * ((-AV_ICa_tot) * AC_AF / (2.0 * AC_vmyo) - AV_Jup * AC_vnsr / AC_vmyo + AV_Jdiff * AC_vss / AC_vmyo);
    
    /* sodium */
    auto AV_INa_ss_tot = AV_ICaNa + 3.0 * AV_INaCa_ss;
    auto AV_INa_tot = AV_INa + AV_INaL + AV_INab + 3.0 * AV_INaCa + 3.0 * AV_INaK;
    NV_Ith_S(ydot, 2) = (-AV_INa_ss_tot) * AC_AF / AC_vss - AV_JdiffNa;
    NV_Ith_S(ydot, 1) = (-AV_INa_tot) * AC_AF / AC_vmyo + AV_JdiffNa * AC_vss / AC_vmyo;
    
    /* membrane */
    NV_Ith_S(ydot, 0) = -(AV_INa_tot + AV_INa_ss_tot + AV_ICa_tot + AV_ICa_ss_tot + AV_IK_tot + AV_IK_ss_tot);
    
    return 0;
}

/* Set initial values */
static void
default_initial_values(N_Vector y)
{
    NV_Ith_S(y, 0) = -87.0;
    NV_Ith_S(y, 1) = 7.0;
    NV_Ith_S(y, 2) = 7.0;
    NV_Ith_S(y, 3) = 145.0;
    NV_Ith_S(y, 4) = 145.0;
    NV_Ith_S(y, 5) = 0.0001;
    NV_Ith_S(y, 6) = 0.0001;
    NV_Ith_S(y, 7) = 1.2;
    NV_Ith_S(y, 8) = 1.2;
    NV_Ith_S(y, 9) = 0.0;
    NV_Ith_S(y, 10) = 1.0;
    NV_Ith_S(y, 11) = 1.0;
    NV_Ith_S(y, 12) = 1.0;
    NV_Ith_S(y, 13) = 1.0;
    NV_Ith_S(y, 14) = 1.0;
    NV_Ith_S(y, 15) = 0.0;
    NV_Ith_S(y, 16) = 1.0;
    NV_Ith_S(y, 17) = 1.0;
    NV_Ith_S(y, 18) = 0.0;
    NV_Ith_S(y, 19) = 1.0;
    NV_Ith_S(y, 20) = 1.0;
    NV_Ith_S(y, 21) = 0.0;
    NV_Ith_S(y, 22) = 1.0;
    NV_Ith_S(y, 23) = 1.0;
    NV_Ith_S(y, 24) = 0.0;
    NV_Ith_S(y, 25) = 1.0;
    NV_Ith_S(y, 26) = 1.0;
    NV_Ith_S(y, 27) = 1.0;
    NV_Ith_S(y, 28) = 1.0;
    NV_Ith_S(y, 29) = 1.0;
    NV_Ith_S(y, 30) = 0.0;
    NV_Ith_S(y, 31) = 1.0;
    NV_Ith_S(y, 32) = 1.0;
    NV_Ith_S(y, 33) = 0.0;
    NV_Ith_S(y, 34) = 0.0;
    NV_Ith_S(y, 35) = 0.0;
    NV_Ith_S(y, 36) = 0.0;
    NV_Ith_S(y, 37) = 1.0;
    NV_Ith_S(y, 38) = 0.0;
    NV_Ith_S(y, 39) = 0.0;
    NV_Ith_S(y, 40) = 0.0;
}