      ! generic operator names:
      character(len=18), parameter ::
     &     op_unity      = 'I  ',
     &     op_ham        = 'H  ',
     &     op_fock       = 'F  ',
     &     op_hhat       = 'Hhat ',
     &     op_hbar       = 'Hbar ',
     &     op_top        = 'T ',
     &     op_tpt        = 'T(pt) ',
     &     op_tbar       = 'TBAR ',
     &     op_omg        = 'OMG ',
     &     op_omg_l      = 'OMG_L ',
     &     op_omg_t      = 'OMG_T ',
     &     op_dia        = 'DIA ',
     &     op_dia_ip     = 'DIA(-1) ',
     &     op_dia_ea     = 'DIA(+1) ',
     &     op_dia_pt     = 'DIA_PT '
      character(len=18), parameter ::
     &     op_r          = 'R ',
     &     op_l          = 'L ',
     &     op_rip        = 'R(-1) ',
     &     op_lip        = 'L(-1) ',
     &     op_rea        = 'R(+1) ',
     &     op_lea        = 'L(+1) ',
     &     op_cclg       = 'L(CC) ',
     &     op_ecclg      = 'L(ECC) ',
     &     op_ccen       = 'E(CC) ',
     &     op_mplg       = 'L(MP) ',
     &     op_mpen       = 'E(MP) ',
     &     op_dept       = 'DE(PT) ',
     &     op_dlpt       = 'DL(PT) ',
     &     op_ccr12s0    = 'S0(R12) ',
     &     op_tbar_a     = 'TBAR.A ',
     &     op_etapt      = 'ETA(PT) ',
     &     op_eta        = 'ETA ',
     &     op_1dens      = '1DEN ',
     &     op_h0_r       = 'H0.R ',
     &     op_a_r        = 'A.R ',
     &     op_l_a        = 'L.A ',
     &     op_s_t        = 'S.T ',
     &     op_s_r        = 'S.R ',
     &     op_tbar_s     = 'TBAR.S ',
     &     op_l_s        = 'L.S ',
     &     op_l_s_r      = 'L.S.R ',
     &     op_l_a_r      = 'L.A.R ',
     &     op_a_rip      = 'A.R(-1) ',
     &     op_lip_a      = 'L(-1).A ',
     &     op_a_rea      = 'A.R(+1) ',
     &     op_lea_a      = 'L(+1).A '
      ! R12 extensions
      character(len=18), parameter ::
     &     op_fpk        = 'F+K ',
     &     op_mpr12lg    = 'L(MPR12) ',
     &     op_mpr12en    = 'E(MPR12) ',
     &     op_ccr12lg    = 'L(CCR12) ',
     &     op_ccr12en    = 'E(CCR12) ',
     &     op_r12        = 'R12 ',
     &     op_g_x        = 'G ',
     &     op_g_z        = 'G_Z ',
     &     op_rba        = 'R12+ ',
c     &     op_c12       = 'C ',
c     &     op_cba       = 'CBAR ',
     &     op_c12        = 'T'' ',
     &     op_cba        = 'TBAR'' ',
     &     op_omgr12     = 'OMG-R12 ',
     &     op_diar12     = 'DIA-R12 ',
     &     op_rint       = 'R12-INT ',
     &     op_rintx      = 'R12-INTX ',
     &     op_rintbar    = 'R12BAR ',
     &     op_rdagbar    = 'R12+BAR ',
     &     op_rintbreve  = 'R12BRV ',
     &     op_rinttilde  = 'R12TLD ',
     &     op_rintc      = 'R12C-I ',
     &     op_rinba      = 'R12+-INT ',
     &     op_rttr       = 'RTTR_COM ',
     &     op_ttr        = 'TTR_COM ',
     &     op_ff         = 'R.R ',
     &     op_ffg        = 'R.R.G ',
     &     op_ffbar      = 'R.RBAR ',
     &     op_gr         = 'G.R12 ',
     &     op_gr_x       = 'G.R12-X ',
     &     op_rg         = 'R12.G ',
     &     op_v_inter    = 'V ',
     &     op_v_x        = 'V-X ',
     &     op_vp_inter   = 'U ',
     &     op_v3_inter   = 'V3 ',
     &     op_v4_inter   = 'V4 ',
     &     op_v0_inter   = 'V0 ',
     &     op_vbar_inter = 'V+ ',
     &     op_b_inter    = 'B ',
     &     op_bh_inter   = 'Bh ',
     &     op_c_inter    = 'C-INT ',
     &     op_cbar_inter = 'C+-INT ',
     &     op_p_inter    = 'P-INT ',
     &     op_p3f_inter  = 'P3F ',
     &     op_p3g_inter  = 'P3G ',
     &     op_z_inter    = 'Z-INT ',
     &     op_z4_inter   = 'Z4-INT ',
     &     op_k4_inter   = 'K4-INT ',
     &     op_b_inv      = 'B_INV ',
     &     op_x_inter    = 'X ',
     &     op_xh_inter   = 'Xh ',
     &     op_xp_inter   = 'X'' ',
     &     op_x1_inter   = 'X1 ',
     &     op_x_inv      = 'X_INV ',
     &     op_exchange   = 'K ',
     &     op_hartree    = 'F+K ',
c dbg
     &     op_x_test     = 'V-TEST ',
     &     op_v_test     = 'X-TEST ',
     &     op_p_test     = 'P-TEST ',
     &     op_z_test     = 'Z-TEST ',
c dbg
     &     op_cex        = 'T12'' ',
     &     op_cex_pt     = 'T12''(pt) ',
     &     op_cexbar     = 'T12BAR'' ',
     &     op_omgcex     = 'OMG_T12'' ',
     &     op_cexx       = 'T12" ',
     &     op_cexx_pt    = 'T12"(pt) ',
     &     op_cexxbar    = 'T12BAR" ',
     &     op_omgcexx    = 'OMG_T12" '
      character(len=18), parameter ::
     &     op_rp         = 'R'' ',
     &     op_a_rp       = '(A.R)'' ',
     &     op_s_rp       = '(S.R)'' ',
     &     op_s_c        = '(S.T)'' ',
     &     op_evs        = 'V_S ',
     &     op_s_evs      = 'S.V_S ',
     &     op_lp         = 'L'' ',
     &     op_l_ap       = '(L.A)'' ',
     &     op_l_sp       = '(L.S)'' ',
     &     op_tbar_ap    = '(TBR.A)'' ',
     &     op_tbar_sp    = '(TBR.S)'' '
      ! generic ME-list names
      character(len=18), parameter ::
     &     mel_cclg0      = 'L0(CC) ',
     &     mel_ecclg0     = 'L0(ECC) ',
     &     mel_ccen0      = 'E0(CC) ',
     &     mel_dlpt0      = 'DL0(PT) ',
     &     mel_dept0      = 'DE0(PT) ',
     &     mel_dept0def   = 'DEF-DE0(PT) ',
     &     mel_ham        = 'H0 ',
     &     mel_hhat       = 'H0hat ',
     &     mel_hbar       = 'H0bar ',
     &     mel_top        = 'T0 ',
     &     mel_tpt        = 'T(pt)0 ',
     &     mel_tptdef     = 'DEF-T(pt)0 ',
     &     mel_h0_tpt     = 'H0.T(pt)0 ',
     &     mel_h0_tptdef  = 'DEF-H0.T(pt)0 ',
     &     mel_tbar       = 'L0 ',
     &     mel_tbar_a     = 'L0.A ',
     &     mel_eta        = 'ETA0 ',
     &     mel_etapt      = 'ETA(pt)0 ',
     &     mel_etaptdef   = 'DEF-ETA(pt)0 ',
     &     mel_omg        = 'O0 ',
     &     mel_omg_l      = 'O0_TBAR ',
     &     mel_omg_t      = 'O0_T ',
     &     mel_topdef     = 'DEF-T0 ',   ! list-definition target
     &     mel_omgdef     = 'DEF-O0 ',
     &     mel_hhatdef    = 'DEF-Hhat ',
     &     meldef_hbar    = 'DEF-Hbar ',
     &     mel_ccen0def   = 'DEF-E0(CC) ',
     &     mel_tbardef    = 'DEF-L0 ',
     &     mel_etadef     = 'DEF-ETA0 ',
     &     mel_tbar_adef  = 'DEF-L0.A ',
     &     mel_dia        = 'DIA ',
     &     mel_dia_pt     = 'DIA_PT ',
     &     mel_dia_ip     = 'DIA(-1) ',
     &     mel_dia_ea     = 'DIA(+1) ',
     &     mel_1dens      = '1DEN0 ',
     &     meldef_1dens   = 'DEF-1DEN0 ',
     &     meldef_ecclg0  = 'DEF-L0(ECC) ',
     &     meldef_omg_l   = 'DEF-O0_TBAR ',
     &     meldef_omg_t   = 'DEF-O0_T '
      character(len=18), parameter ::
     &     mel_rex        = 'RE0 ',
     &     meldef_rex     = 'DEF-RE0 ',
     &     mel_a_rex      = 'A.RE0 ',
     &     meldef_a_rex   = 'DEF-A.RE0 ',
     &     mel_s_rex      = 'S.RE0 ',
     &     meldef_s_rex   = 'DEF-S.RE0 ',
     &     mel_lex        = 'LE0 ',
     &     meldef_lex     = 'DEF-LE0 ',
     &     mel_lex_a      = 'LE0.A ',
     &     meldef_lex_a   = 'DEF-LE0.A ',
     &     mel_l_s_r      = 'LE.S.RE0 ',
     &     meldef_l_s_r   = 'DEF-LE0.S.RE0 ',
     &     mel_l_a_r      = 'LE.A.RE0 ',
     &     meldef_l_a_r   = 'DEF-LE0.A.RE0 ',
     &     mel_lex_s      = 'LE0.S ',
     &     meldef_lex_s   = 'DEF-LE0.S ',
     &     mel_rip        = 'RI0 ',
     &     meldef_rip     = 'DEF-RI0 ',
     &     mel_a_rip      = 'A.RI0 ',
     &     meldef_a_rip   = 'DEF-A.RI0 ',
     &     mel_lip        = 'LI0 ',
     &     meldef_lip     = 'DEF-LI0 ',
     &     mel_lip_a      = 'LI0.A ',
     &     meldef_lip_a   = 'DEF-LI0.A ',
     &     mel_rea        = 'RA0 ',
     &     meldef_rea     = 'DEF-RA0 ',
     &     mel_a_rea      = 'A.RA0 ',
     &     meldef_a_rea   = 'DEF-A.RA0 ',
     &     mel_lea        = 'LA0 ',
     &     meldef_lea     = 'DEF-LA0 ',
     &     mel_lea_a      = 'LA0.A ',
     &     meldef_lea_a   = 'DEF-LA0.A '
      ! ME-lists for R12
      character(len=18), parameter ::
     &     mel_rint         = 'R12-GEM ',
     &     mel_rintx        = 'R12-GEMX ',
     &     mel_gintx        = 'G-X ',
     &     mel_gintz        = 'G-Z ',
     &     mel_rintbar      = 'R12BAR-GEM ',
     &     mel_rdagbar      = 'R12BAR+-GEM ',
     &     mel_rinttilde    = 'R12TLD-GEM ',
     &     mel_rintbreve    = 'R12BRV-GEM ',
     &     mel_rintc        = 'R12C-GEM ',
     &     meldef_rintbar   = 'DEF-R12BAR-GEM ',
     &     meldef_rdagbar   = 'DEF-R12BAR+-GEM ',
     &     meldef_rinttilde = 'DEF-R12TILDE-GEM ',
     &     meldef_rintc    = 'DEF-R12C-GEM ',
     &     mel_rinba        = 'R12+-GEM ',
     &     mel_ttr          = 'TTR-INT ',
     &     mel_rttr         = 'RTTR-INT ',
     &     mel_ff           = 'RSQ-INT ',
     &     mel_ffg          = 'RSQ.G-INT ',
     &     mel_ffbar        = 'RSQBAR-INT ',
     &     mel_rg           = 'R.G-INT ',
     &     mel_gr           = 'G.R-INT ',
     &     mel_exchange     = 'K-INT ',
     &     mel_hartree      = 'F+K-INT ',
     &     meldef_hartree   = 'DEF-F+K-INT ',
     &     mel_v_inter      = 'V-INTER ',
     &     mel_v_def        =   'DEF-V-INTER ',
     &     mel_vp_inter     = 'U-INTER ',
     &     mel_vp_def       = 'DEF-U-INTER ',
     &     mel_v0_inter     = 'V0-INTER ',
     &     mel_v0_def       = 'DEF-V0-INTER ',
     &     mel_vbar_inter   = 'V+-INTER ',
     &     mel_vbar_def     = 'DEF-V+-INTER ',
     &     mel_v0bar_inter  = 'V0+-INTER ',
     &     mel_v0bar_def    = 'DEF-V0+-INTER ',
     &     mel_x_inter      = 'X-INTER ',
     &     mel_x_def        = 'DEF-X-INTER ',
     &     mel_xh_inter     = 'Xh-INTER ',
     &     mel_xh_def       = 'DEF-Xh-INTER ',
     &     mel_xp_inter     = 'X''-INTER ',
     &     mel_xp_def       = 'DEF-X''-INTER ',
     &     mel_x1_inter     = 'X1-INTER ',
     &     mel_x1_def       = 'DEF-X1-INTER ',
     &     mel_x_inv        = 'XINV ',
     &     mel_b_inter      = 'B-INTER ',
     &     mel_b_def        = 'DEF-B-INTER ',
     &     mel_bh_inter     = 'Bh-INTER ',
     &     mel_bh_def       = 'DEF-Bh-INTER ',
     &     mel_c_inter      = 'C-INTER ',
     &     mel_c_def        = 'DEF-C-INTER ',
     &     mel_cbar_inter   = 'C+-INTER ',
     &     mel_cbar_def     = 'DEF-C+-INTER ',
     &     mel_p_inter      = 'P-INTER ',
     &     mel_p_def        = 'DEF-P-INTER ',
     &     mel_z_inter      = 'Z-INTER ',
     &     mel_z_def        = 'DEF-Z-INTER ',
     &     mel_p3f_inter    = 'P3F-INTER ',
     &     mel_p3f_def      = 'DEF-P3F-INTER ',
     &     mel_p3g_inter    = 'P3G-INTER ',
     &     mel_p3g_def      = 'DEF-P3G-INTER ',
     &     mel_b_inv        = 'BINV ',
     &     mel_b_dia        = 'BDIA ',
c dbg
     &     mel_x_test_def  = 'DEF-X-TEST ',
     &     mel_x_test      = 'X-INTER-TEST ',
     &     mel_v_test_def  = 'DEF-V-TEST ',
     &     mel_v_test      = 'V-INTER-TEST ',
     &     mel_p_test_def  = 'DEF-P-TEST ',
     &     mel_p_test      = 'P-INTER-TEST ',
     &     mel_z_test_def  = 'DEF-Z-TEST ',
     &     mel_z_test      = 'Z-INTER-TEST ',
c dbg
     &     mel_diar12      = 'DIA(R12) '
      character(len=18), parameter ::
     &     mel_c12def      = 'DEF-T0''-R12 ',
     &     mel_c12         = 'T0''-R12 ',
     &     mel_cbardef     = 'DEF-T0BAR''-R12 ',
     &     mel_cbar        = 'T0BAR''-R12 ',
     &     mel_omgr12def   = 'DEF-OMG0-R12 ',
     &     mel_omgr12      = 'OMG0-R12 ',
     &     mel_mpr12lg0def = 'DEF-L0(MPR12) ',
     &     mel_mpr12lg0    = 'L0(MPR12) ',
     &     mel_mpr12en0def = 'DEF-E0(MPR12) ',
     &     mel_mpr12en0    = 'E0(MPR12) ',
     &     mel_ccr12lg0def = 'DEF-L0(CCR12) ',
     &     mel_ccr12lg0    = 'L0(CCR12) ',
     &     mel_ccr12en0def = 'DEF-E0(CCR12) ',
     &     mel_ccr12en0    = 'E0(CCR12) ',
     &     mel_cex_def     = 'DEF-T12''-R12 ',      
     &     mel_cex         = 'T12''-R12 ',
     &     mel_cexbar_def  = 'DEF-T12BAR''-R12 ',
     &     mel_cexbar    = 'T12BAR''-R12 ',
     &     mel_omgcexdef = 'DEF-OMG12''-R12 ',
     &     mel_omgcex     = 'OMG12''-R12 ',
     &     mel_cexx_def   = 'DEF-T12"-R12 ',      
     &     mel_cexx        = 'T12"-R12 ',
     &     mel_cexxbar_def = 'DEF-T12BAR"-R12 ',
     &     mel_cexxbar    = 'T12BAR"-R12 ',
     &     mel_omgcexxdef = 'DEF-OMG12"-R12 ',
     &     mel_omgcexx     = 'OMG12"-R12 '
      character(len=18), parameter ::
     &     mel_rpex        = 'RE0P ',
     &     meldef_rpex     = 'DEF-RE0P ',
     &     mel_a_rpex      = 'A.RE0P ',
     &     meldef_a_rpex  = 'DEF-A.RE0P ',
     &     mel_s_rpex      = 'S.RE0P ',
     &     meldef_s_rpex  = 'DEF-S.RE0P ',
     &     mel_lpex        = 'LE0P ',
     &     meldef_lpex     = 'DEF-LE0P ',
     &     mel_lex_ap      = 'LE0.AP ',
     &     meldef_lex_ap  = 'DEF-LE0.AP ',
     &     mel_lex_sp      = 'LE0.SP ',
     &     meldef_lex_sp  = 'DEF-LE0.SP ',
     &     mel_evs         = 'V0_S ',
     &     meldef_evs      = 'DEF-V0_S ',
     &     mel_s_evs   = 'S.V0_S ',
     &     meldef_s_evs   = 'DEF-S.V0_S '

	  
      