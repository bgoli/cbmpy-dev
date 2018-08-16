\ENCODING=ISO-8859-1
\Problem name: MODELID_1461534

Minimize
 MNAF: bin_R_EXinput_pi_e + bin_R_EXinput_nh4_e + bin_R_EXinput_uri_e
       + bin_R_EXinput_asp_L_e + bin_R_EXinput_glyclt_e + bin_R_EXinput_tre_e
       + bin_R_EXinput_glcn_e + bin_R_EXinput_k_e + bin_R_EXinput_glc_e
       + bin_R_EXinput_dgsn_e + bin_R_EXinput_gly_e + bin_R_EXinput_12ppd_S_e
       + bin_R_EXinput_thm_e + bin_R_EXinput_so4_e + bin_R_EXinput_gam_e
       + bin_R_EXinput_galct_D_e + bin_R_EXinput_3hcinnm_e
       + bin_R_EXinput_malttr_e + bin_R_EXinput_man6p_e + bin_R_EXinput_ile_L_e
       + bin_R_EXinput_lys_L_e + bin_R_EXinput_3hpppn_e + bin_R_EXinput_ura_e
       + bin_R_EXinput_xan_e + bin_R_EXinput_fe2_e + bin_R_EXinput_succ_e
       + bin_R_EXinput_acnam_e + bin_R_EXinput_val_L_e + bin_R_EXinput_ala_L_e
Subject To
 M_10fthf_c:      - R_AICART + R_MTHFC - R_FTHFD - R_GARFT  = 0
 M_12dgr_EC_c:    0.02 R_PAPA_EC - 0.02 R_DAGK_EC + 0.02 R_PEPT_EC  = 0
 M_12ppd_S_c:     R_12PPDt + R_LCARS  = 0
 M_12ppd_S_e:     - R_12PPDt - R_EXoutput_12ppd_S_e + R_EXinput_12ppd_S_e  = 0
 M_13dpg_c:       R_PGK + R_GAPD  = 0
 M_15dap_c:       R_LYSDC - R_CADVt  = 0
 M_15dap_e:       - R_EX_15dap_e + R_CADVt  = 0
 M_1pyr5c_c:      R_G5SADs - R_P5CR - R_P5CD + R_PROD2  = 0
 M_23ddhb_c:      R_ICHORT - R_DHBD  = 0
 M_23dhb_c:       - R_DHBSr + R_DHBD  = 0
 M_23dhba_c:      - 3 R_ENTCS + R_DHBSr  = 0
 M_23dhdp_c:      - R_DHDPRy + R_DHDPS  = 0
 M_23dhmb_c:      R_KARA1i - R_DHAD1  = 0
 M_23dhmp_c:      R_KARA2i - R_DHAD2  = 0
 M_23doguln_c:    - R_DOGULNR  = 0
 M_25aics_c:      - R_ADSL2r + R_PRASCS  = 0
 M_25dkglcn_c:    - R_DKGLCNR2y - R_DKGLCNR1 - R_DKGLCNR2x  = 0
 M_25drapp_c:     - R_DHPPDA2 + R_GTPCII2  = 0
 M_26dap_LL_c:    R_SDPDS - R_DAPE  = 0
 M_26dap_M_c:     R_DAPabc - R_DAPDC - R_UAAGDS + R_DAPE  = 0
 M_26dap_M_e:     - R_DAPabc - R_EX_26dap_M_e  = 0
 M_2ahbut_c:      R_ACHBS - R_KARA2i  = 0
 M_2aobut_c:      - R_GLYATi + R_THRD  = 0
 M_2cpr5p_c:      - R_IGPS + R_PRAIi  = 0
 M_2dda7p_c:      - R_DHQS + R_DDPA  = 0
 M_2ddg6p_c:      R_EDD - R_EDA + R_DDGLK  = 0
 M_2ddglcn_c:     R_DDGLCNt2r - R_DDGLK + R_MNNH + R_ALTRH  = 0
 M_2ddglcn_e:     - R_DDGLCNt2r - R_EX_2ddglcn_e  = 0
 M_2dh3dgal6p_c:  R_DDGALK - R_DDPGALA  = 0
 M_2dh3dgal_c:    - R_DDGALK + R_GALCTND  = 0
 M_2dhglcn_c:     - R_2DGLCNRy - R_2DGLCNRx  = 0
 M_2dhguln_c:     - R_2DGULRx - R_2DGULRy + R_DKGLCNR1  = 0
 M_2dhp_c:        - R_DPR + R_MOHMT  = 0
 M_2dmmq8_c:      R_DHNAOT - R_AMMQT8_2 + R_DMSOR2 - R_HYD3 + R_TMAOR2
                  - R_GLYCTO4 - R_G3PD7 + R_DMSOR2e + R_TMAOR2e - R_NADH9
                  - R_NADH8 + R_FRD3  = 0
 M_2dmmql8_c:     - R_DMSOR2 + R_HYD3 - R_TMAOR2 + R_GLYCTO4 + R_G3PD7
                  - R_DMSOR2e - R_TMAOR2e + R_NADH9 + R_NADH8 - R_FRD3  = 0
 M_2dr1p_c:       - R_PPM2 + R_DURIPP + R_PUNP4 + R_PUNP6 + R_PUNP2 + R_TMDPP
                   = 0
 M_2dr5p_c:       R_PPM2 - R_DRPA  = 0
 M_2h3oppan_c:    R_GLCRAL + R_GLXCL - R_TRSAR + R_HPYRI  = 0
 M_2ippm_c:       R_IPPMIa - R_IPPMIb  = 0
 M_2kmb_c:        R_DKMPPD + R_DKMPPD2 - R_UNK3  = 0
 M_2mahmp_c:      - R_TMPPP + R_PMPK  = 0
 M_2mcacn_c:      R_MCITD - R_MICITD  = 0
 M_2mcit_c:       R_MCITS - R_MCITD  = 0
 M_2me4p_c:       - R_MEPCT + R_DXPRIi  = 0
 M_2mecdp_c:      R_MECDPS - R_MECDPDH  = 0
 M_2obut_c:       - R_ACHBS + R_THRD_L - R_OBTFL  = 0
 M_2ohph_c:       - R_OHPHM + R_OPHHX  = 0
 M_2ombzl_c:      R_OMPHHX - R_OMBZLM  = 0
 M_2omhmbl_c:     R_OMMBLHX - R_DMQMT  = 0
 M_2ommbl_c:      - R_OMMBLHX + R_OMBZLM  = 0
 M_2omph_c:       - R_OMPHHX + R_OHPHM  = 0
 M_2oph_c:        R_OPHBDC - R_OPHHX  = 0
 M_2p4c2me_c:     R_CDPMEK - R_MECDPS  = 0
 M_2pg_c:         - R_ENO - R_PGM  = 0
 M_2pglyc_c:      - R_PGLYCP  = 0
 M_2shchc_c:      R_SHCHCS2 - R_SUCBZS  = 0
 M_34hpp_c:       R_PPND + R_TYRTA  = 0
 M_3c2hmp_c:      - R_IPMD - R_IPPMIa  = 0
 M_3c3hmp_c:      R_IPPMIb + R_IPPS  = 0
 M_3c4mop_c:      R_IPMD - R_OMCDC  = 0
 M_3dhguln_c:     R_DOGULNR  = 0
 M_3dhgulnp_c:    - R_KG6PDC  = 0
 M_3dhq_c:        - R_DHQD + R_DHQS  = 0
 M_3dhsk_c:       - R_SHK3Dr + R_DHQD  = 0
 M_3hcinnm_c:     - R_3HCINNMH + R_HCINNMt2r  = 0
 M_3hcinnm_e:     R_EXinput_3hcinnm_e - R_HCINNMt2r - R_EXoutput_3hcinnm_e  = 0
 M_3hmrsACP_c:    - R_U23GAAT - R_UAGAAT + R_KAS16  = 0
 M_3hpppn_c:      - R_3HPPPNH + R_HPPPNt2r  = 0
 M_3hpppn_e:      - R_EXoutput_3hpppn_e - R_HPPPNt2r + R_EXinput_3hpppn_e  = 0
 M_3ig3p_c:       R_IGPS - R_TRPS3 - R_TRPS1  = 0
 M_3mob_c:        R_VALTA - R_MOHMT - R_VPAMT - R_IPPS + R_DHAD1  = 0
 M_3mop_c:        R_ILETA + R_DHAD2  = 0
 M_3ophb_c:       R_HBZOPT - R_OPHBDC  = 0
 M_3pg_c:         R_GLYCK + R_PGM - R_PGK - R_PGCD  = 0
 M_3php_c:        - R_PSERT + R_PGCD  = 0
 M_3psme_c:       R_PSCVT - R_CHORS  = 0
 M_4abut_c:       R_GLUDC - R_ABTA + R_ABUTt2 - R_GLUABUTt7 + R_ABUTD  = 0
 M_4abut_e:       - R_ABUTt2 + R_GLUABUTt7 - R_EX_4abut_e  = 0
 M_4abutn_c:      R_PTRCTA - R_ABUTD  = 0
 M_4abz_c:        R_ADCL - R_DHPS2  = 0
 M_4adcho_c:      R_ADCS - R_ADCL  = 0
 M_4ahmmp_c:      - R_HMPK1  = 0
 M_4ampm_c:       R_AMPMS - R_PMPK + R_HMPK1  = 0
 M_4c2me_c:       - R_CDPMEK + R_MEPCT  = 0
 M_4h2opntn_c:    R_OP4ENH - R_HOPNTAL  = 0
 M_4hba_c:        R_THZPSN  = 0
 M_4hbz_c:        - R_HBZOPT + R_CHRPL  = 0
 M_4hthr_c:       R_4HTHRS  = 0
 M_4mhetz_c:      - R_HETZK  = 0
 M_4mop_c:        - R_LEUTAi + R_OMCDC  = 0
 M_4mpetz_c:      R_HETZK + R_THZPSN - R_TMPPP  = 0
 M_4pasp_c:       R_ASAD + R_ASPK  = 0
 M_4per_c:        - R_PERD + R_E4PD  = 0
 M_4ppan_c:       R_PNTK - R_PPNCL2  = 0
 M_4ppcys_c:      - R_PPCDC + R_PPNCL2  = 0
 M_4r5au_c:       R_RBFSb - R_RBFSa + R_PMDPHT  = 0
 M_5aizc_c:       - R_AIRC3 - R_PRASCS  = 0
 M_5aop_c:        - 2 R_PPBNGS + R_G1SATi  = 0
 M_5aprbu_c:      R_APRAUR - R_PMDPHT  = 0
 M_5apru_c:       R_DHPPDA2 - R_APRAUR  = 0
 M_5caiz_c:       R_AIRC3 + R_AIRC2  = 0
 M_5dglcn_c:      - R_IDOND - R_IDOND2 + R_DKGLCNR2y - R_5DGLCNR + R_DKGLCNR2x
                   = 0
 M_5dh4dglc_c:    - R_GLCRAL + R_GLCRD + R_GALCTD  = 0
 M_5mdr1p_c:      R_MTRK - R_MTRI  = 0
 M_5mdru1p_c:     - R_MDRPD + R_MTRI  = 0
 M_5mta_c:        R_SPMS - R_MTAN  = 0
 M_5mthf_c:       - 0.05 R_BiomassEcoli + R_MTHFR2 - R_METS  = 0
 M_5mtr_c:        - R_MTRK + R_MTAN  = 0
 M_5prdmbz_c:     R_NNDMBRT - R_RZ5PP  = 0
 M_6hmhpt_c:      R_DHNPA2 - R_HPPK2  = 0
 M_6hmhptpp_c:    - R_DHPS2 + R_HPPK2  = 0
 M_6pgc_c:        - R_EDD + R_GNK - R_GND + R_PGL  = 0
 M_6pgl_c:        - R_PGL + R_G6PDH2r  = 0
 M_8aonn_c:       R_AOXSr - R_AMAOTr  = 0
 M_ACP_c:         R_U23GAAT - R_AACPS5 - R_AACPS4 - R_AACPS2 - R_AACPS1
                  - R_ACOATA + R_ACPS1 + 5 R_C141SN + 4 R_C120SN + R_EDTXS4
                  + R_EDTXS3 + R_EDTXS2 + R_EDTXS1 + 6 R_C161SN + 6 R_C160SN
                  + 7 R_C181SN - R_AACPS3 + 5 R_C140SN + R_UAGAAT + R_KAS14
                  + R_KAS16 + 2 R_PASYN_EC - R_MCOATA  = 0
 M_N1aspmd_c:     R_SPMDAT1  = 0
 M_aacald_c:      R_TAUDO  = 0
 M_aacoa_c:       R_FAO4 + R_ACACCT + R_ACACT1r  = 0
 M_acACP_c:       R_ACOATA + R_MACPD - R_KAS14  = 0
 M_ac_c:          R_ALDD2x + R_NACODA + R_BUTCT + R_ACODA + R_UHGADA + R_ACACCT
                  + R_AGDC - R_ACKr + R_POX + R_ACt2r + R_CYSS - R_ACS + R_CITL
                   = 0
 M_ac_e:          - R_EX_ac_e - R_ACt2r  = 0
 M_acac_c:        - R_ACACCT + R_ACACt2  = 0
 M_acac_e:        - R_ACACt2 - R_EX_acac_e  = 0
 M_acald_c:       - R_ALDD2x - R_ACALDi + R_HOPNTAL + R_ETHAAL + R_DRPA
                  + R_ACALDt + R_THRAr  = 0
 M_acald_e:       - R_EX_acald_e - R_ACALDt  = 0
 M_accoa_c:       R_PDH - R_BUTCT - R_ACOATA - R_ACCOACr - R_ADHEr + R_GLYATi
                  - R_PTAr - R_ACACCT - 0.00005 R_BiomassEcoli - R_MALS
                  + R_ACALDi - R_ACGS - R_TDPADGAT - 2 R_ACACT1r + R_PFL
                  - R_SERAT + R_ACS - R_IPPS - R_CS - R_KAS15 - R_G1PACT
                  + 7 R_FAO1 + 8 R_FAO2 - R_SPMDAT1 - R_SPMDAT2 + 9 R_FAO3  = 0
 M_acg5p_c:       R_AGPR + R_ACGK  = 0
 M_acg5sa_c:      - R_NACODA - R_AGPR + R_ACOTA  = 0
 M_acgam1p_c:     R_G1PACT - R_UAGDP  = 0
 M_acgam6p_c:     - R_AGDC + R_ACGApts + R_AMANAPE  = 0
 M_acgam_e:       - R_ACGApts - R_EX_acgam_e  = 0
 M_acglu_c:       - R_ACGK + R_ACGS  = 0
 M_acmana_c:      R_ACNML - R_AMANK  = 0
 M_acmana_e:      - R_EX_acmana_e - R_ACMANApts  = 0
 M_acmanap_c:     R_AMANK - R_AMANAPE + R_ACMANApts  = 0
 M_acnam_c:       - R_ACNML + R_ACNAMt2  = 0
 M_acnam_e:       - R_EXoutput_acnam_e - R_ACNAMt2 + R_EXinput_acnam_e  = 0
 M_acon_T_c:      - R_ACONMT  = 0
 M_aconm_c:       R_ACONMT  = 0
 M_acorn_c:       - R_ACODA - R_ACOTA  = 0
 M_acser_c:       - R_CYSS + R_SERAT  = 0
 M_actACP_c:      - R_C141SN - R_C120SN - R_C161SN - R_C160SN - R_C181SN
                  - R_C140SN + R_KAS14 + R_KAS15  = 0
 M_actp_c:        R_PTAr + R_ACKr  = 0
 M_ade_c:         R_ADEt2r + R_AMPN + R_PUNP1 + R_PUNP2 - R_ADD + R_MTAN
                  - R_ADPT + R_AHCYSNS  = 0
 M_ade_e:         - R_ADEt2r - R_EX_ade_e  = 0
 M_adn_c:         - R_ADNK1 - R_PUNP1 + R_ADNt2r + R_NTD7 + R_AHC - R_ADA  = 0
 M_adn_e:         - R_ADNt2r - R_EX_adn_e  = 0
 M_adocbi_c:      R_CBIAT_DELETE - R_ADOCBIK  = 0
 M_adocbip_c:     - R_ACBIPGT + R_ADOCBIK  = 0
 M_adocbl_c:      R_ADOCBLS + R_CBLAT_DELETE  = 0
 M_adp_c:         R_CTPS2 + R_PYDXNK + R_GLYCK + R_CYTK2 + R_INSK + R_ORNabc
                  + R_TMDK1 + R_FCLK + R_TDSK + R_TMKr + R_DGK1 + R_PRAGSr
                  + R_ARBabc + R_ATPM + R_UMPK + R_VALabc + R_CBMK + R_GALabc
                  + R_CYSabc + R_SHKK + R_CBL1abc + R_MALTHXabc + R_GLU5K
                  + R_DPCOAK + R_FRUK + R_PNTK + R_ADNK1 + R_RBK_L1 + R_ACCOACr
                  + R_HSK + R_ACCOAL + R_DAPabc + R_NDPK8 + R_NDPK7 + R_NDPK6
                  + R_NDPK5 + R_NDPK4 + R_NDPK3 + R_NDPK2 + R_NDPK1 + R_DURIK1
                  + R_GLNabc + R_TAURabc + R_PFK + R_GNK + R_GSNK + R_DXYLK
                  + R_NAMNPP + R_TMPKr + R_TSULabc + R_HETZK + R_GLUCYS + R_PGK
                  + R_PTRCabc + R_Kabc + 45.5608 R_BiomassEcoli + R_FE2abc
                  + R_GLUabc + R_MALTabc + R_METDabc + R_PROabc + R_SUCCabc
                  + R_DDGLK + R_UAMAS + R_CDPMEK + R_UAMAGS + 3 R_LPSSYN_EC
                  + R_ASNabc + R_AMANK - R_RNDR1 - R_PYK + R_MALTTRabc + R_ACGK
                  + R_PYDXK + R_HEX1 + R_DTMPK + R_GLCS1 + R_MALTPTabc
                  + R_SPMDabc + R_CYTK1 + R_RBFK + R_DAGK_EC + R_MALTTTRabc
                  + R_DDGALK + R_AIRC2 + R_ACKr + R_METabc + R_XYLK + R_XYLabc
                  + R_PIabc + R_ADSK + R_MTRK + R_DHFS + R_RMK + R_URIDK2r
                  + R_UAAGDS + R_AP5AH + R_ARGabc + R_GMHEPK + R_ALAabc
                  + R_NADK - R_ATPS4r + R_SULabc + R_GTHS + R_GK1 + R_GSPMDS
                  + 2 R_AP4AH + R_PMPK + R_ADOCBIK + R_GLNS + R_UGMDDS + R_RBK
                  + R_SUCOAS + R_GLYC3Pabc + R_LYSabc + 2 R_CBPS + R_PYDAMK
                  + 2 R_ADK1 + R_ADK3 + R_ADK4 + R_PRFGS + R_HMPK1 + R_THMabc
                  + R_RIBabc + R_PFK_2 + R_HISabc + R_ILEabc + R_GLYBabc
                  + R_GLYK + R_GALKr + R_CHLabc + R_THRabc + R_GART + R_PRASCS
                  + R_DADK + R_ALAALAr - R_PPAKr + R_DBTSr + R_PPCK + R_LEUabc
                  + R_ASPabc + R_PRAIS + R_ASPK  = 0
 M_adpglc_c:      R_GLGC - R_GLCS1  = 0
 M_adphep_DD_c:   R_GMHEPAT - R_AGMHE  = 0
 M_adphep_LD_c:   - 3 R_LPSSYN_EC + R_AGMHE  = 0
 M_agdpcbi_c:     - R_ADOCBLS + R_ACBIPGT  = 0
 M_agm_c:         R_ARGDC - R_AGMT  = 0
 M_agpc_EC_c:     0.02 R_PLIPA3 - 0.02 R_LPLIPA5 - 0.02 R_LPLIPA3  = 0
 M_agpe_EC_c:     0.02 R_PLIPA2 - 0.02 R_LPLIPA2 - 0.02 R_LPLIPA4  = 0
 M_agpg_EC_c:     - 0.02 R_LPLIPA6 + 0.02 R_PLIPA1 - 0.02 R_LPLIPA1  = 0
 M_ahcys_c:       R_AMMQT8_2 + R_OHPHM + R_OMBZLM + R_DMQMT + R_ACONMT
                  + 2 R_UPP3MT - R_AHC - R_AHCYSNS  = 0
 M_ahdt_c:        R_GTPCI - R_DNTPPA  = 0
 M_aicar_c:       R_ADSL2r - R_AICART + R_IG3PS  = 0
 M_air_c:         - R_AMPMS - R_AIRC2 + R_PRAIS  = 0
 M_akg_c:         R_GLUDy + R_HSTPT - R_ASPTA - R_AKGDH + R_LEUTAi - R_TYRTA
                  + R_AKGt2r - R_ABTA - R_VALTA - R_PHETA1 - R_ILETA - R_TAUDO
                  + R_PSERT - R_ACOTA - R_PTRCTA - R_ALATA_L + R_ICDHyr
                  - R_SOTA + R_OHPBAT - R_SDPTA - R_OXGDC2 + R_UNK3 + R_TDPAGTA
                  - R_GLUSy - R_ORNTA  = 0
 M_akg_e:         - R_AKGt2r - R_EX_akg_e  = 0
 M_ala_B_c:       - R_PANTS + R_ASP1DC  = 0
 M_ala_D_c:       R_ALAR - R_DAAD - R_ALATA_D2 + R_DALAt2r - 2 R_ALAALAr  = 0
 M_ala_D_e:       - R_DALAt2r - R_EX_ala_D_e  = 0
 M_ala_L_c:       - R_AOXSr - R_ALAR - 0.488 R_BiomassEcoli - R_UAMAS
                  - R_ALATA_L2 + R_ALAt2r + R_THZPSN - R_ALATA_L + R_ALAabc
                  - R_VPAMT + R_BTS2  = 0
 M_ala_L_e:       - R_EXoutput_ala_L_e - R_ALAt2r + R_EXinput_ala_L_e
                  - R_ALAabc  = 0
 M_alaala_c:      - R_UGMDDS + R_ALAALAr  = 0
 M_alac_S_c:      - R_KARA1i + R_ACLS  = 0
 M_alltn_c:       R_ALLTNt2r - R_ALLTN  = 0
 M_alltn_e:       - R_ALLTNt2r - R_EX_alltn_e  = 0
 M_alltt_c:       - R_ALLTAH + R_ALLTN  = 0
 M_altrn_c:       - R_TAGURr - R_ALTRH  = 0
 M_amet_c:        - R_AMMQT8_2 - R_OHPHM - R_AMAOTr - R_OMBZLM - R_DMQMT
                  + R_METAT - R_ACONMT - 2 R_UPP3MT - R_ADMDCr  = 0
 M_ametam_c:      - R_SPMS + R_ADMDCr  = 0
 M_amob_c:        R_AMAOTr  = 0
 M_amp_c:         R_PRPPS + R_PANTS + R_AACPS5 + R_AACPS4 + R_AACPS2 + R_AACPS1
                  + R_NADDP + R_GLUTRS + R_ADNK1 + R_ADSL1r + 6 R_ENTCS
                  + R_BPNT - 0.001 R_BiomassEcoli + R_SUCBZL + R_NADS1 - R_AMPN
                  + R_NTPP6 + R_PPS + R_THZPSN + R_AACPS3 - R_NTD7 + R_PACCOAL
                  - R_ADK1 - R_ADK3 - R_ADK4 + R_ACS + R_ASNS1 + R_ADPT
                  + R_GMPS2 + R_FAO1 + R_FAO2 + R_ARGSS + R_HPPK2 + R_ASNS2
                  + R_FAO3 + R_SELNPS  = 0
 M_amp_e:         - R_EX_amp_e + R_NADDPe  = 0
 M_anth_c:        - R_ANPRT + R_ANS  = 0
 M_ap4a_c:        - R_AP4AH  = 0
 M_ap5a_c:        - R_AP5AH  = 0
 M_apg_EC_c:      0.02 R_LPLIPA6 + 0.02 R_LPLIPA5 + 0.02 R_LPLIPA4  = 0
 M_apoACP_c:      - R_ACPS1  = 0
 M_aps_c:         R_SADT2 - R_ADSK  = 0
 M_ara5p_c:       - R_KDOPS + R_A5PISO  = 0
 M_arab_L_c:      R_ARBabc + R_ARBt2r - R_ARAI  = 0
 M_arab_L_e:      - R_ARBabc - R_ARBt2r - R_EX_arab_L_e  = 0
 M_arbt6p_c:      - R_AB6PGH  = 0
 M_arg_L_c:       - R_ARGDC - R_AST - 0.281 R_BiomassEcoli + R_ARGSL + R_ARGabc
                  + R_ARGORNt7  = 0
 M_arg_L_e:       - R_EX_arg_L_e - R_ARGabc - R_ARGORNt7  = 0
 M_argsuc_c:      - R_ARGSL + R_ARGSS  = 0
 M_asn_L_c:       - R_ASNN - 0.229 R_BiomassEcoli + R_ASNabc + R_ASNt2r
                  + R_ASNS1 + R_ASNS2  = 0
 M_asn_L_e:       - R_EX_asn_L_e - R_ASNabc - R_ASNt2r  = 0
 M_asp_L_c:       R_ASNN - R_ASP1DC - R_ASPTA + R_ASPt2 - R_ASPO4 - R_ASPO6
                  - 0.229 R_BiomassEcoli - R_ASPO3 - R_ADSS + R_ASPt2_3
                  + R_ASPt2_2 - R_ASPO5 - R_ASNS1 - R_ASPCT - R_PRASCS
                  - R_ARGSS - R_ASNS2 + R_ASPabc - R_ASPT - R_ASPK  = 0
 M_asp_L_e:       - R_ASPt2 - R_EXoutput_asp_L_e - R_ASPt2_3 - R_ASPt2_2
                  + R_EXinput_asp_L_e - R_ASPabc  = 0
 M_aspsa_c:       - R_ASAD - R_DHDPS + R_HSDy  = 0
 M_atp_c:         - R_CTPS2 - R_PRPPS - R_PYDXNK - R_PANTS - R_GLYCK - R_CYTK2
                  - R_INSK - R_ORNabc - R_TMDK1 - R_FMNAT - R_FCLK - R_TDSK
                  - R_TMKr - R_DGK1 - R_PRAGSr - R_ARBabc - R_ATPM - R_AACPS5
                  - R_AACPS4 - R_AACPS2 - R_AACPS1 - R_UMPK - R_VALabc - R_CBMK
                  - R_GALabc - R_CYSabc - R_SHKK - R_GLUTRS - R_CBL1abc
                  - R_MALTHXabc - R_GLU5K - R_DPCOAK - R_FRUK - R_PNTK
                  - R_ADNK1 - R_RBK_L1 - R_ACCOACr - R_HSK - R_ACCOAL
                  - R_DAPabc - R_NDPK8 - R_NDPK7 - R_NDPK6 - R_NDPK5 - R_NDPK4
                  - R_NDPK3 - R_NDPK2 - R_NDPK1 - R_DURIK1 - R_GLNabc
                  - R_TAURabc - R_PFK - R_GNK - R_GMHEPAT - R_NMNAT - R_GSNK
                  - R_DXYLK - R_NAMNPP - R_TMPKr - R_TSULabc - R_GLGC - R_HETZK
                  - R_GLUCYS - R_PGK - R_PTRCabc - R_Kabc
                  - 45.7318 R_BiomassEcoli - R_FE2abc - R_GLUabc - R_SUCBZL
                  - R_CBLAT_DELETE - R_MALTabc - R_METDabc - R_PROabc - R_NADS1
                  - R_SERASr - R_SUCCabc - R_DDGLK - R_UAMAS - R_DHBSr
                  - R_SADT2 - R_CDPMEK - R_UAMAGS - R_RNTR1 - R_ASNabc
                  - R_AMANK - R_METAT - R_NTPP6 + R_PYK - R_MALTTRabc - R_ACGK
                  - R_PYDXK - R_ATPPRT - R_PPS - R_HEX1 - R_DTMPK
                  - R_CBIAT_DELETE - R_MALTPTabc - R_THZPSN - R_SPMDabc
                  - R_CYTK1 - R_RBFK - R_DAGK_EC - R_MALTTTRabc - R_DDGALK
                  - R_AIRC2 - R_ACKr - R_METabc - R_AACPS3 - R_XYLK - R_XYLabc
                  - R_PIabc - R_ADSK - R_MTRK - R_DHFS - R_RMK - R_URIDK2r
                  - R_UAAGDS + R_AP5AH - R_ARGabc - R_GMHEPK - R_ALAabc
                  - R_NADK + R_ATPS4r - R_SULabc - R_PACCOAL - R_GTHS - R_GK1
                  - R_GSPMDS - R_PMPK - R_ADOCBIK - R_GLNS - R_UGMDDS - R_RBK
                  - R_SUCOAS - R_GLYC3Pabc - R_LYSabc - 2 R_CBPS - R_PYDAMK
                  - R_ADK1 - R_ACS - R_PRFGS - R_ASNS1 - R_HMPK1 - R_THMabc
                  - R_RIBabc - R_PFK_2 - R_HISabc - R_ILEabc - R_GLYBabc
                  - R_GLYK - R_GALKr - R_CHLabc - R_GMPS2 - R_NNAT - R_THRabc
                  - R_GART - R_PRASCS - R_DADK - R_ALAALAr - R_FAO1 - R_FAO2
                  + R_PPAKr - R_ARGSS - R_HPPK2 - R_DBTSr - R_PPCK - R_LEUabc
                  - R_PTPATi - R_ASNS2 - R_ASPabc - R_ADNCYC - R_FAO3
                  - R_SELNPS - R_PRAIS - R_ASPK  = 0
 M_bbtcoa_c:      - R_CRNBTCT  = 0
 M_betald_c:      - R_BETALDHy - R_BETALDHx  = 0
 M_btcoa_c:       R_BUTCT - R_FAO4  = 0
 M_btn_c:         R_BSORy + R_BTS2 + R_BSORx  = 0
 M_btnso_c:       - R_BSORy - R_BSORx  = 0
 M_but_c:         - R_BUTCT + R_BUTt2r  = 0
 M_but_e:         - R_BUTt2r - R_EX_but_e  = 0
 M_camp_c:        R_ADNCYC  = 0
 M_cbasp_c:       R_DHORTS + R_ASPCT  = 0
 M_cbi_c:         - R_CBIAT_DELETE  = 0
 M_cbl1_c:        R_CBL1abc - R_CBLAT_DELETE  = 0
 M_cbl1_e:        - R_CBL1abc - R_EX_cbl1_e  = 0
 M_cbp_c:         R_CBMK - R_OCBT + R_CBPS - R_ASPCT  = 0
 M_cdp_c:         - R_NDPK3 + 2 R_LPSSYN_EC + R_CYTK1 - R_RNDR3  = 0
 M_cdpdag_EC_c:   0.02 R_DASYN_EC - 0.02 R_PSSA_EC - 0.02 R_PGSA_EC
                  - 0.02 R_CDAPPA_EC  = 0
 M_cdpea_c:       - 2 R_LPSSYN_EC + R_PEPT_EC  = 0
 M_cechddd_c:     R_PPPNDO - R_DHPPD  = 0
 M_cenchddd_c:    - R_DHCIND + R_CINNDO  = 0
 M_chol_c:        R_CHLt2r + R_GPDDA1 + R_CHLabc  = 0
 M_chol_e:        - R_EX_chol_e - R_CHLt2r - R_CHLabc  = 0
 M_chor_c:        - R_ADCS - R_CHRPL - R_CHORM + R_CHORS - R_ANS - R_ICHORSi
                   = 0
 M_cinnm_c:       - R_CINNDO  = 0
 M_cit_c:         - R_ACONT + R_CITt7 + R_CS - R_CITL  = 0
 M_cit_e:         - R_EX_cit_e - R_CITt7  = 0
 M_citr_L_c:      R_OCBT - R_ARGSS  = 0
 M_ckdo_c:        - R_MOAT - 3 R_LPSSYN_EC + R_KDOCT2 - R_MOAT2  = 0
 M_clpn_EC_c:     0.02 R_CLPNS_EC - 0.000129 R_BiomassEcoli  = 0
 M_cmp_c:         R_MOAT + R_PSSA_EC + R_PPNCL2 + 3 R_LPSSYN_EC + R_MECDPS
                  + R_NTPP4 + R_PGSA_EC - R_CYTK1 + R_CYTDK2 - R_CMPN
                  - R_PEPT_EC - R_NTD4 + R_MOAT2 + R_CDAPPA_EC  = 0
 M_co2_c:         R_LYSDC + R_AOXSr + R_GLYCL + R_DHNAOT + R_PPND + R_PDH
                  + R_GLUDC + R_ORNDC + R_PPCDC + R_ASP1DC - R_CBMK + R_ACHBS
                  + R_PSD_EC + 5 R_C141SN + R_AKGDH + R_GND + 4 R_C120SN
                  + R_OMCDC + R_IGPS + R_ARGDC + R_UGLYCH + R_TAUDO + R_PPNDH
                  + 6 R_C161SN + 2 R_CPPPGO - R_HCO3E - R_PPC + R_MACPD
                  + R_ACLS + 6 R_C160SN + R_THZPSN + 7 R_C181SN + 2 R_CYNTAH
                  + R_DXPS + R_POX + R_DAPDC + R_SADH + R_FDH3 + R_FDH2
                  + R_ICDHyr + 5 R_C140SN + R_PDX5PS + R_FHL + R_OXGDC2
                  + R_OPHBDC + R_KG6PDC + R_ME1 + R_ME2 + 4 R_UPPDC1 + R_OMPDC
                  + R_KAS14 + R_KAS16 + R_GLXCL + R_NNDPR + R_KAS15 + R_ADMDCr
                  - R_DBTSr + R_PPCK + R_CO2t + R_MMCD  = 0
 M_co2_e:         - R_EX_co2_e - R_CO2t  = 0
 M_coa_c:         R_AOXSr - R_PDH + R_ACOATA - R_ACPS1 + R_DPCOAK - R_ACCOAL
                  + R_HSST - R_AKGDH + R_ADHEr - R_GLYATi + R_PTAr + R_PTA2
                  + R_AST - 0.000006 R_BiomassEcoli + R_MCITS - R_SUCBZL
                  + R_MALS - R_ACALDi + R_THDPS + R_ACGS + R_TDPADGAT
                  + R_ACACT1r - R_PFL - R_PACCOAL - R_SUCOAS + R_SERAT + R_NPHS
                  - R_OBTFL - R_ACS + R_IPPS + R_CS + R_KAS15 + R_G1PACT
                  - 7 R_FAO1 - 8 R_FAO2 + R_SPMDAT1 + R_SPMDAT2 - 9 R_FAO3
                  + R_MCOATA  = 0
 M_cpppg3_c:      - R_CPPPGO + R_UPPDC1  = 0
 M_crn_c:         R_CRNt7 - R_CRNCBCT - R_CRNBTCT  = 0
 M_crn_e:         - R_EX_crn_e - R_CRNt7  = 0
 M_crncoa_c:      R_CRNCBCT - R_CRNCDH + R_CRNBTCT  = 0
 M_csn_c:         R_CSNt2 + R_CMPN - R_CSND  = 0
 M_csn_e:         - R_CSNt2 - R_EX_csn_e  = 0
 M_ctbt_c:        R_CRNCBCT  = 0
 M_ctbtcoa_c:     - R_CRNCBCT + R_CRNCDH  = 0
 M_ctp_c:         R_CTPS2 - R_DASYN_EC + R_NDPK3 - R_RNTR3 - R_PPNCL2
                  - 0.126 R_BiomassEcoli - R_NTPP4 - R_MEPCT - R_KDOCT2  = 0
 M_cyan_c:        - R_CYANST  = 0
 M_cynt_c:        R_CYNTt2 - R_CYNTAH  = 0
 M_cynt_e:        - R_CYNTt2 - R_EX_cynt_e  = 0
 M_cys_L_c:       R_CYSabc - R_PPNCL2 - R_GLUCYS - 0.087 R_BiomassEcoli
                  - R_THZPSN - R_CYSDS + R_CYSS - R_BTS2 - R_SHSL1  = 0
 M_cys_L_e:       - R_CYSabc - R_EX_cys_L_e  = 0
 M_cyst_L_c:      R_SHSL1 - R_CYSTL  = 0
 M_cytd_c:        R_CYTDt2r - R_CYTDK2 + R_NTD4 - R_CYTD  = 0
 M_cytd_e:        - R_CYTDt2r - R_EX_cytd_e  = 0
 M_dad_2_c:       - R_PUNP2 + R_DADNt2 + R_NTD6 - R_DADA  = 0
 M_dad_2_e:       - R_EX_dad_2_e - R_DADNt2  = 0
 M_dadp_c:        - R_NDPK8 + R_RNDR1 + R_DADK  = 0
 M_damp_c:        R_NTPP5 - R_NTD6 - R_DADK  = 0
 M_dann_c:        R_AMAOTr - R_DBTSr  = 0
 M_datp_c:        R_NDPK8 - 0.0247 R_BiomassEcoli + R_RNTR1 - R_NTPP5  = 0
 M_db4p_c:        R_DB4PS - R_RBFSa  = 0
 M_dcamp_c:       - R_ADSL1r + R_ADSS  = 0
 M_dcdp_c:        R_CYTK2 - R_NDPK7 + R_RNDR3  = 0
 M_dcmp_c:        - R_CYTK2 - R_NTD3 + R_NTPP3  = 0
 M_dctp_c:        R_NDPK7 + R_RNTR3 - 0.0254 R_BiomassEcoli - R_NTPP3 - R_DCTPD
                   = 0
 M_dcyt_c:        R_NTD3 + R_DCYTt2 - R_DCYTD  = 0
 M_dcyt_e:        - R_EX_dcyt_e - R_DCYTt2  = 0
 M_ddcaACP_c:     R_C120SN - R_EDTXS1 - R_KAS16  = 0
 M_dgdp_c:        R_DGK1 - R_NDPK5 + R_RNDR2  = 0
 M_dgmp_c:        - R_DGK1 - R_NTD8 + R_NTPP1  = 0
 M_dgsn_c:        R_NTD8 + R_DGSNt2 - R_PUNP4 + R_NTPTP1  = 0
 M_dgsn_e:        R_EXinput_dgsn_e - R_EXoutput_dgsn_e - R_DGSNt2  = 0
 M_dgtp_c:        R_NDPK5 + R_RNTR2 - 0.0254 R_BiomassEcoli - R_NTPP1
                  - R_NTPTP1  = 0
 M_dha_c:         - R_DHAPT + R_DHAt + R_F6PA + R_GLYCDx  = 0
 M_dha_e:         - R_DHAt - R_EX_dha_e  = 0
 M_dhap_c:        R_RMPA + R_DHAPT + R_FCLPA - R_TPI + R_G3PD7 + R_G3PD5
                  + R_G3PD2 - R_MGSA - R_QULNS + R_TGBPA + R_FBA + R_G3PD6  = 0
 M_dhcinnm_c:     R_DHCIND - R_DHCINDO + R_3HCINNMH  = 0
 M_dhf_c:         R_TMDS + R_DHFS - R_DHFR  = 0
 M_dhna_c:        - R_DHNAOT + R_NPHS  = 0
 M_dhnpt_c:       R_DNMPPA - R_DHNPA2  = 0
 M_dhor_S_c:      - R_DHORD5 - R_DHORTS - R_DHORD2  = 0
 M_dhpmp_c:       - R_DNMPPA + R_DNTPPA  = 0
 M_dhpppn_c:      R_3HPPPNH - R_HPPPNDO + R_DHPPD  = 0
 M_dhpt_c:        R_DHPS2 - R_DHFS  = 0
 M_dhptd_c:       R_RHCCE - R_DHPTDCs  = 0
 M_din_c:         - R_PUNP6 + R_DINSt2 + R_DADA  = 0
 M_din_e:         - R_EX_din_e - R_DINSt2  = 0
 M_dkmpp_c:       - R_DKMPPD - R_DKMPPD2 + R_MDRPD  = 0
 M_dmbzid_c:      - R_NNDMBRT  = 0
 M_dmlz_c:        - 2 R_RBFSb + R_RBFSa  = 0
 M_dmpp_c:        R_DMPPS + R_IPDDIi - R_DMATT  = 0
 M_dms_c:         R_DMSOR1 + R_DMSOR2  = 0
 M_dms_e:         R_DMSOR1e + R_DMSOR2e - R_EX_dms_e  = 0
 M_dmso_c:        - R_DMSOR1 - R_DMSOR2  = 0
 M_dmso_e:        - R_DMSOR1e - R_EX_dmso_e - R_DMSOR2e  = 0
 M_dnad_c:        - R_NADS1 + R_NNAT  = 0
 M_dpcoa_c:       - R_DPCOAK + R_PTPATi  = 0
 M_dscl_c:        - R_SHCHD2 + R_UPP3MT  = 0
 M_dtbt_c:        - R_BTS2 + R_DBTSr  = 0
 M_dtdp4aaddg_c:  - R_AADDGT + R_TDPADGAT  = 0
 M_dtdp4addg_c:   - R_TDPADGAT + R_TDPAGTA  = 0
 M_dtdp4d6dg_c:   R_TDPGDH - R_TDPDRE - R_TDPAGTA  = 0
 M_dtdp4d6dm_c:   - R_TDPDRR + R_TDPDRE  = 0
 M_dtdp_c:        R_AADDGT - R_NDPK4 + R_DTMPK  = 0
 M_dtdpglu_c:     - R_TDPGDH + R_G1PTT  = 0
 M_dtdprmn_c:     R_TDPDRR  = 0
 M_dtmp_c:        R_TMDS + R_TMDK1 - R_NTD5 + R_NTPP7 - R_DTMPK  = 0
 M_dttp_c:        R_NDPK4 - 0.0247 R_BiomassEcoli - R_NTPP7 - R_G1PTT  = 0
 M_dudp_c:        R_RNDR4 - R_NDPK6 + R_URIDK2r  = 0
 M_dump_c:        - R_TMDS - R_NTD1 + R_DUTPDP + R_DURIK1 - R_URIDK2r  = 0
 M_duri_c:        R_NTD1 + R_DURIt2 - R_DURIK1 - R_DURIPP + R_DCYTD  = 0
 M_duri_e:        - R_DURIt2 - R_EX_duri_e  = 0
 M_dutp_c:        R_NDPK6 - R_DUTPDP + R_RNTR4 + R_DCTPD  = 0
 M_dxyl5p_c:      R_DXYLK - R_THZPSN + R_DXPS - R_PDX5PS - R_DXPRIi  = 0
 M_dxyl_c:        - R_DXYLK  = 0
 M_e4p_c:         R_TALA - R_DDPA - R_TKT2 - R_E4PD  = 0
 M_eca_EC_c:      R_ECAP_EC  = 0
 M_eig3p_c:       - R_IGPDH + R_IG3PS  = 0
 M_enter_c:       R_ENTCS  = 0
 M_etha_c:        R_GPDDA2 - R_ETHAAL  = 0
 M_etoh_c:        R_ADHEr + R_ETOHt2r  = 0
 M_etoh_e:        - R_ETOHt2r - R_EX_etoh_e  = 0
 M_f1p_c:         - R_FRUK + R_FRUpts  = 0
 M_f6p_c:         R_FRUpts2 + R_SBTPD - R_GF6PTA + R_M1PD + R_G6PDA - R_PFK
                  + R_PGI - R_F6PA + R_MAN6PI + R_TALA + R_FBP + R_TKT2  = 0
 M_fad_c:         R_FMNAT - R_FAO4 - R_DAAD - 0.00001 R_BiomassEcoli + R_SUCD4
                  - R_PROD2 - R_SUCD1i - 6 R_FAO1 - 7 R_FAO2 - 8 R_FAO3  = 0
 M_fadh2_c:       R_FAO4 + R_DAAD - R_SUCD4 + R_PROD2 + R_SUCD1i + 6 R_FAO1
                  + 7 R_FAO2 + 8 R_FAO3  = 0
 M_fc1p_c:        R_FCLK - R_FCLPA  = 0
 M_fcl_L_c:       - R_FCLK + R_FCI  = 0
 M_fdp_c:         R_FRUK + R_PFK - R_FBP - R_FBA  = 0
 M_fe2_c:         - R_FCLT + R_FE2abc - R_SHCHF  = 0
 M_fe2_e:         R_EXinput_fe2_e - R_FE2abc - R_EXoutput_fe2_e  = 0
 M_fgam_c:        - R_PRFGS + R_GART + R_GARFT  = 0
 M_fmn_c:         - R_FMNAT + R_RBFK  = 0
 M_for_c:         R_DKMPPD + R_GTPCI + R_DKMPPD2 + R_FORt + R_DB4PS + 2 R_AMPMS
                  - R_FDH3 - R_FDH2 + R_PFL - R_FHL + R_OBTFL + R_FTHFD
                  + R_GTPCII2 - R_GART  = 0
 M_for_e:         - R_FORt - R_EX_for_e  = 0
 M_fpram_c:       R_PRFGS - R_PRAIS  = 0
 M_fprica_c:      R_IMPC + R_AICART  = 0
 M_frdp_c:        R_GRTT - R_OCTDPS - R_UDCPDPS - R_HEMEOS  = 0
 M_fru_c:         R_FFSD - R_XYLI2i  = 0
 M_fru_e:         - R_FRUpts2 - R_FRUpts - R_EX_fru_e  = 0
 M_fruur_c:       R_MANAO + R_GUI1  = 0
 M_fuc1p_L_c:     R_FUCPt6_2  = 0
 M_fuc1p_L_e:     - R_EX_fuc1p_L_e - R_FUCPt6_2  = 0
 M_fuc_L_c:       - R_FCI + R_FUCt  = 0
 M_fuc_L_e:       - R_EX_fuc_L_e - R_FUCt  = 0
 M_fum_c:         R_ADSL2r + R_ADSL1r + R_SUCFUMt + R_ARGSL - R_FUM + R_FUMt2_3
                  + R_FUMt2_2 + R_HKNTDH + R_SUCD1i - R_ASPO5 - R_FRD3 + R_ASPT
                  - R_FRD2  = 0
 M_fum_e:         - R_SUCFUMt - R_FUMt2_3 - R_FUMt2_2 - R_EX_fum_e  = 0
 M_g1p_c:         R_UGLT + R_GLCP + R_MLTP3 - R_GLGC + R_MLTP1 + R_MLTP2
                  - R_G1PP - R_GALU - R_PGMT - R_G1PTT  = 0
 M_g3p_c:         R_EDA + R_TPI + R_TRPS3 + R_TRPS1 + R_F6PA - R_TALA - R_GAPD
                  - R_DXPS + R_TKT1 + R_DRPA + R_TGBPA + R_TKT2 + R_DDPGALA
                  + R_FBA  = 0
 M_g3pc_c:        - R_GPDDA1 + R_LPLIPA5 + R_LPLIPA3  = 0
 M_g3pe_c:        - R_GPDDA2 + R_LPLIPA2 + R_LPLIPA4  = 0
 M_g3pg_c:        R_LPLIPA6 - R_GPDDA4 + R_LPLIPA1  = 0
 M_g3pi_c:        - R_GPDDA5  = 0
 M_g3ps_c:        - R_GPDDA3  = 0
 M_g6p_c:         R_GLCpts + R_FFSD - R_PGI + R_HEX1 - R_TRE6PS + R_TRE6PH
                  + R_G6Pt6_2 + R_PGMT + R_AB6PGH - R_G6PDH2r  = 0
 M_g6p_e:         - R_G6Pt6_2 - R_EX_g6p_e  = 0
 M_gal1p_c:       - R_UGLT + R_GALKr  = 0
 M_gal_c:         R_GALabc + R_GALS3 + R_LACZ + R_GALt2 - R_GALKr  = 0
 M_gal_e:         - R_GALabc - R_GALt2 - R_EX_gal_e  = 0
 M_galct_D_c:     R_GALCTt2r - R_GALCTD  = 0
 M_galct_D_e:     - R_GALCTt2r - R_EXoutput_galct_D_e + R_EXinput_galct_D_e
                   = 0
 M_galctn_D_c:    R_GALCTNt2r - R_GALCTND  = 0
 M_galctn_D_e:    - R_EX_galctn_D_e - R_GALCTNt2r  = 0
 M_galt1p_c:      R_GALTpts - R_GLTPD  = 0
 M_galt_e:        - R_GALTpts - R_EX_galt_e  = 0
 M_galur_c:       - R_GUI2 + R_GALURt2r  = 0
 M_galur_e:       - R_GALURt2r - R_EX_galur_e  = 0
 M_gam1p_c:       - R_PGAMT - R_G1PACT  = 0
 M_gam6p_c:       R_GF6PTA - R_G6PDA + R_AGDC + R_GAMpts + R_PGAMT  = 0
 M_gam_e:         - R_EXoutput_gam_e + R_EXinput_gam_e - R_GAMpts  = 0
 M_gar_c:         R_PRAGSr - R_GART - R_GARFT  = 0
 M_gbbtn_c:       - R_CRNt7 + R_CRNBTCT  = 0
 M_gbbtn_e:       R_CRNt7 - R_EX_gbbtn_e  = 0
 M_gcald_c:       - R_GCALDD + R_DHNPA2  = 0
 M_gdp_c:         - R_NDPK1 + 2 R_GP4GH + R_SADT2 + R_ADSS + R_CYTDK2
                  - R_MAN1PT2 + R_GK1 - R_RNDR2 + R_ADK3 + R_URIK2  = 0
 M_gdpddman_c:    - R_GDMANE + R_GMAND  = 0
 M_gdpfuc_c:      R_GOFUCR  = 0
 M_gdpmann_c:     R_MAN1PT2 - R_GMAND  = 0
 M_gdpofuc_c:     R_GDMANE - R_GOFUCR  = 0
 M_glc_D_c:       2 R_TREH + R_GALS3 + R_LACZ + R_MLTG1 + R_XYLI2i + R_GLCt2
                  - R_HEX1 + R_TRE6PH + R_G1PP + R_MLTG2 + R_MLTG3 + R_MLTG4
                  + R_MLTG5 + R_AMALT1 + R_AMALT3 + R_AMALT4 + R_AMALT2  = 0
 M_glc_D_e:       - R_GLCpts - R_EXoutput_glc_e + 2 R_TREHe - R_GLCt2
                  + R_EXinput_glc_e - R_GLCDe  = 0
 M_glcn_c:        - R_GNK + R_2DGLCNRy + R_2DGLCNRx + R_GLCNt2r + R_5DGLCNR
                   = 0
 M_glcn_e:        R_EXinput_glcn_e - R_GLCNt2r + R_GLCDe - R_EXoutput_glcn_e
                   = 0
 M_glcr_c:        R_GLCRt2r - R_GLCRD  = 0
 M_glcr_e:        - R_GLCRt2r - R_EX_glcr_e  = 0
 M_glcur_c:       - R_GUI1 + R_GLCURt2r  = 0
 M_glcur_e:       - R_GLCURt2r - R_EX_glcur_e  = 0
 M_gln_L_c:       - R_CTPS2 - R_ADCS - R_GF6PTA + R_GLNabc - R_GLUN
                  - 0.25 R_BiomassEcoli + R_GLNS - R_CBPS - R_PRFGS - R_ASNS1
                  - R_GLUSy - R_IG3PS - R_GMPS2 - R_ANS - R_GLUPRT  = 0
 M_gln_L_e:       - R_GLNabc - R_EX_gln_L_e  = 0
 M_glu1sa_c:      R_GLUTRR - R_G1SATi  = 0
 M_glu5p_c:       - R_G5SD + R_GLU5K  = 0
 M_glu5sa_c:      R_NACODA - R_G5SADs + R_G5SD + R_ORNTA  = 0
 M_glu_D_c:       - R_GLUR - R_UAMAGS  = 0
 M_glu_L_c:       R_CTPS2 + R_ADCS - R_GLUDC - R_GLUDy + R_GLUt2r - R_GLUTRS
                  - R_HSTPT - R_GLU5K + R_GF6PTA + R_ASPTA - R_LEUTAi + R_TYRTA
                  + R_GLUR + R_GLUN + R_ABTA + R_VALTA + R_PHETA1 - R_GLUCYS
                  + R_ILETA - 0.25 R_BiomassEcoli + R_GLUabc - R_PSERT
                  + R_GLUABUTt7 - R_ACGS + R_ACOTA + R_GLUt4 + R_PTRCTA
                  + R_ALATA_L + R_P5CD - R_DHFS + R_SOTA - R_OHPBAT + R_SDPTA
                  + R_SGDS - R_UNK3 - R_GLNS + R_CBPS - R_TDPAGTA + R_PRFGS
                  + R_ASNS1 + 2 R_GLUSy + R_IG3PS + R_ORNTA + R_GMPS2 + R_ANS
                  + R_GLUPRT  = 0
 M_glu_L_e:       - R_GLUt2r - R_GLUabc - R_GLUABUTt7 - R_GLUt4 - R_EX_glu_L_e
                   = 0
 M_glucys_c:      R_GLUCYS - R_GTHS  = 0
 M_glutrna_c:     R_GLUTRS - R_GLUTRR  = 0
 M_glx_c:         - R_GLYCLTDx - R_GLYCLTDy + R_GLYCTO2 + R_GLYCTO3 + R_ICL
                  + R_UGLYCH - R_MALS + R_GLYCTO4 - 2 R_GLXCL  = 0
 M_gly_c:         - R_GLYCL - R_PRAGSr + R_GHMT2 + R_GLYATi
                  - 0.582 R_BiomassEcoli + R_GLYt2r - R_GTHS + R_THRAr  = 0
 M_gly_e:         - R_EXoutput_gly_e + R_EXinput_gly_e - R_GLYt2r  = 0
 M_glyald_c:      R_GLYALDt - R_ALCD19  = 0
 M_glyald_e:      - R_EX_glyald_e - R_GLYALDt  = 0
 M_glyb_c:        R_GLYBt2r + R_BETALDHy + R_BETALDHx + R_GLYBabc  = 0
 M_glyb_e:        - R_GLYBt2r - R_EX_glyb_e - R_GLYBabc  = 0
 M_glyc3p_c:      R_GPDDA5 + R_GPDDA4 + R_GPDDA3 + R_GPDDA2 + R_GPDDA1
                  - R_PGSA_EC - R_G3PD7 - R_G3PD5 - R_G3PD2 + R_GLYC3Pt6
                  + R_GLYC3Pabc + R_GLYK - R_PASYN_EC - R_G3PD6  = 0
 M_glyc3p_e:      - R_EX_glyc3p_e - R_GLYC3Pt6 - R_GLYC3Pabc  = 0
 M_glyc_R_c:      - R_GLYCK + R_TRSAR + R_HPYRRx + R_HPYRRy  = 0
 M_glyc_c:        - R_GLYCt + R_CLPNS_EC + R_GPDDA4 - R_GLYCDx + R_ALCD19
                  - R_GLYK  = 0
 M_glyc_e:        R_GLYCt - R_EX_glyc_e  = 0
 M_glyclt_c:      R_GLYCLTDx + R_GLYCLTDy - R_GLYCTO2 - R_GLYCTO3 + R_GCALDD
                  + R_GLYCLTt2r - R_GLYCTO4 + R_PGLYCP  = 0
 M_glyclt_e:      R_EXinput_glyclt_e - R_GLYCLTt2r - R_EXoutput_glyclt_e  = 0
 M_glycogen_c:    - R_GLCP - 0.154 R_BiomassEcoli + R_GLCS1  = 0
 M_gmhep17bp_c:   - R_GMHEPPA + R_GMHEPK  = 0
 M_gmhep1p_c:     - R_GMHEPAT + R_GMHEPPA  = 0
 M_gmhep7p_c:     R_S7PI - R_GMHEPK  = 0
 M_gmp_c:         R_GUAPRT - R_NTD9 + R_GSNK + R_ADOCBLS - R_GMPR + R_NTPP2
                  - R_GK1 + R_GMPS2  = 0
 M_gp4g_c:        - R_GP4GH  = 0
 M_grdp_c:        - R_GRTT + R_DMATT  = 0
 M_gsn_c:         R_GSNt2 + R_NTD9 - R_GSNK + R_NTPTP2 - R_PUNP3  = 0
 M_gsn_e:         - R_GSNt2 - R_EX_gsn_e  = 0
 M_gthox_c:       - R_GTHOr  = 0
 M_gthrd_c:       - R_LGTHL + R_GLYOX + 2 R_GTHOr + R_GTHS + R_GSPMDA
                  - R_GSPMDS  = 0
 M_gtp_c:         - R_GTPCI + R_NDPK1 - R_RNTR2 - R_ACBIPGT - R_NTPTP2
                  - 0.203 R_BiomassEcoli - R_SADT2 - R_NTPP2 - R_ADSS
                  - R_CYTDK2 - R_ADK3 - R_URIK2 - R_GTPCII2  = 0
 M_gtspmd_c:      - R_GSPMDA + R_GSPMDS  = 0
 M_gua_c:         - R_GUAPRT + R_GUAt + R_GUAt2 + R_PUNP4 + R_PUNP3 - R_GUAD
                   = 0
 M_gua_e:         - R_GUAt - R_GUAt2 - R_EX_gua_e  = 0
 M_h2_c:          - R_HYD2 - R_HYD3 - R_HYD1 + R_FHL  = 0
 M_h2mb4p_c:      - R_DMPPS - R_IPDPS + R_MECDPDH  = 0
 M_h2o2_c:        - 2 R_CAT + R_SPODM + R_PDX5PO + R_ASPO6 + R_PEAMNO
                  + R_PYAM5PO  = 0
 M_h2o_c:         - R_PLIPA3 - R_TREH - R_CTPS2 - R_ALDD2x - R_NACODA - R_ASNN
                  - R_PYDXPP - R_LCAD - R_ORNabc + 2 R_CAT + R_EDD - R_ARBabc
                  - R_SDPDS - R_ATPM - R_ALLTAH - R_NADDP - R_NTD9 - R_NTD8
                  - R_ACODA - R_NTD1 - R_NTD3 - R_NTD2 - R_NTD5 - R_GLUDy
                  + R_RNDR4 - R_PRAMPC + R_G5SADs - R_FFSD - R_VALabc
                  - R_GALabc + R_DMPPS - R_CYSabc - R_CBL1abc - R_DHPPDA2
                  - R_MALTHXabc - R_DNMPPA - R_NNAM - R_USHD + 6 R_C141SN
                  - R_PLIPA1 - R_PLIPA2 - R_ICHORT - R_DKMPPD - R_UACMAMO
                  + R_DMSOR1 + R_DMSOR2 - R_PGPP_EC + R_3HPPPNH + R_GHMT2
                  - R_UDCPDP - R_FAO4 - R_GTPCI - R_DAPabc - R_DUTPDP
                  - R_GLNabc - R_GALS3 - R_TAURabc - R_LACZ - R_G6PDA + R_ENO
                  - R_UHGADA + R_RNTR2 + R_RNTR3 - R_GLUN + 5 R_C120SN
                  - R_DNTPPA - 3 R_DKMPPD2 - R_NAMNPP - R_TRPAS2 - R_DAAD
                  + R_IPDPS - R_NTPTP2 - R_TSULabc + R_H2Ot - R_SSALy - R_SSALx
                  - R_MLTG1 - R_BPNT + R_3HCINNMH - R_GCALDD + R_IGPS + R_BSORy
                  - R_AGDC - R_PGL - R_PTRCabc - R_UGLYCH + R_TRPS2 + R_TRPS1
                  - R_Kabc - 45.5608 R_BiomassEcoli - R_FE2abc - R_MCITS
                  - R_UDPGD + R_MCITD - R_GLUabc - R_NTPP1 - R_KDOPS + R_SUCBZS
                  - R_GMHEPPA - R_THRS - R_NADPPPS - R_CBLAT_DELETE - R_MALTabc
                  - R_METDabc - R_NMNN - R_PROabc - R_AGMT - R_SUCCabc - R_MALS
                  + R_TMAOR2 - R_IMPC - R_AMPN + R_PPNDH - R_GP4GH - R_AMPMS
                  - R_OP4ENH - R_PAPA_EC - R_SADT2 - R_IMPD + R_MNNH + R_RNTR1
                  - R_ASNabc + R_TARTD - R_GPDDA5 - R_GPDDA4 - R_GPDDA3
                  - R_GPDDA2 - R_GPDDA1 - R_METAT - R_NTPP8 - R_ALLTN - R_NTPP6
                  - R_NTPP7 - R_NTPP4 - R_NTPP5 - R_NTPP2 - R_NTPP3 + R_RNDR1
                  + 7 R_C161SN - R_THDPS - R_MALTTRabc + 2 R_CPPPGO - R_GLYOX
                  - R_DHORTS - R_HCO3E + 2 R_RBFSa - R_PPS - R_NTPTP1 - R_FUM
                  - R_4HTHRS - R_PPC - R_PPA - R_PEAMNO - R_TRE6PP + R_DHQD
                  - R_CBIAT_DELETE - R_MOHMT - R_TRE6PH - R_G1PP - R_MALTPTabc
                  + R_MDRPD + 7 R_C160SN + R_THZPSN - R_SPMDabc + 2 R_PPBNGS
                  - R_HMBS - R_CMPN - R_MALTTTRabc + 8 R_C181SN - R_PYAM5PO
                  - 2 R_P5CD - R_SGSAD - R_NMNt7 + R_IGPDH - R_POX - R_LPLIPA2
                  - R_LPLIPA3 - R_LPLIPA1 + R_CYTBD - R_METabc - R_PMDPHT
                  - 2 R_SADH - R_XYLabc + R_NO3R2 + R_CYTBO3 + R_NO3R1 - R_FBP
                  + R_RNTR4 - R_ALDD19x - R_PIabc - R_NTD10 + R_CRNCDH
                  + R_DHAD2 - R_BETALDHy - R_BETALDHx + R_TDPGDH + R_MECDPDH
                  - R_NTD4 - R_NTD7 - R_HKNDDH - R_PGLYCP - R_PSP_L - R_AP5AH
                  + 6 R_C140SN + 2 R_PDX5PS - R_CSND - R_RZ5PP - R_ARGabc
                  - R_DDPA + 2 R_QULNS - R_SGDS - R_ALAabc + R_ATPS4r
                  - R_SULabc - R_CYSDS - R_CYTD - R_ADD - R_MTHFC - R_GSPMDA
                  - R_MTAN - R_AP4AH - R_HKNTDH + R_TMAOR1 - R_MLTG2 - R_MLTG3
                  - R_MLTG4 - R_MLTG5 + R_ALTRH + R_DHPTDCs + R_IPPMIa
                  - R_IPPMIb - R_NTD11 + R_UPP3S - R_NTD6 - R_ABUTD
                  - R_GLYC3Pabc + 3 R_PPPGO - R_LYSabc + 2 R_NTRIR2x - R_MI1PP
                  - R_CBPS + R_RNDR2 + R_RNDR3 - R_AB6PGH - R_HYPOE - R_IPPS
                  + 2 R_DHDPS - R_PRFGS - R_ASNS1 - R_THMabc - R_AHC - R_RIBabc
                  - R_CS - R_HISabc - R_ILEabc - R_DCTPD - R_GLYBabc - R_HISTP
                  - R_MICITD - R_PDXPP - R_HISTD + R_GLCRD - R_NMNDA - R_FTHFD
                  - R_CHLabc - R_GMPS2 - R_THRabc + R_GALCTND - 3 R_GTPCII2
                  - R_HEMEOS + R_BSORx - R_CDAPPA_EC + R_GMAND + R_DHAD1
                  - R_DADA - R_GLUPRT - R_E4PD - 6 R_FAO1 - 7 R_FAO2
                  - R_AHCYSNS - R_LEUabc + R_GALCTD - R_ASPabc - R_CYSTL
                  - R_ADA - 3 R_SULR - R_GUAD - R_KDOPP - 8 R_FAO3 - R_SELNPS
                  - R_PRATPP - R_DCYTD  = 0
 M_h2o_e:         R_EXinput_h2o_e - R_TREHe - R_H2Ot + R_DMSOR1e + R_TMAOR1e
                  + R_DMSOR2e + R_TMAOR2e - R_EXoutput_h2o_e - R_NADDPe
                  - R_GLCDe  = 0
 M_h2s_c:         R_CYSDS - R_CYSS - R_SULR  = 0
 M_h_c:           R_PLIPA3 + 2 R_CTPS2 - R_LYSDC + 2 R_ALDD2x + R_PRPPS
                  - R_DASYN_EC + R_IPMD + R_PROt2r - R_IDOND + R_PYDXNK
                  + R_PANTS - R_AOXSr + R_GLYCK + R_ADCL + R_MOAT + R_RMNt
                  + R_SPMS + R_GLYBt2r + R_DHPS2 + R_SERt2r + R_DHNAOT
                  + 2 R_LCAD + R_INSK + 2 R_FCLT + R_ORNabc + R_TMDK1 + R_Kt2r
                  - R_FMNAT + R_FCLK + R_TDSK + 2 R_SUCCt2_2 + R_UAGPT3
                  + R_TMKr - R_SHK3Dr + R_U23GAAT + R_PRAGSr + R_ALLTNt2r
                  + R_AMMQT8_2 - 2 R_SPODM + R_ARBabc + R_DHCIND + R_ATPM
                  + R_ACMAMUT + R_SBTPD + R_GSNt2 + 2 R_PAPSR + 2 R_NADDP
                  - R_GLUDC - R_GLYCLTDx - R_GLYCLTDy + R_PSSA_EC + R_NAt3_1
                  + 2 R_NAt3_2 + R_GLUDy - R_ORNDC + R_AADDGT + R_OHPHM
                  - R_PPCDC + R_LYSt2r + R_G5SADs + R_DURIt2 + R_GLUt2r
                  + R_VALabc - R_ASP1DC + 2 R_CBMK + R_GALabc - R_ACHBS
                  - R_DMPPS + R_CYSabc + R_SHKK - R_GLUTRR + R_CBL1abc - R_G5SD
                  - R_DHPPDA2 + R_ACPS1 + R_MALTHXabc + 2 R_USHD - R_PSD_EC
                  - 16 R_C141SN + R_PLIPA1 + R_CYTDt2r + R_PLIPA2 + 2 R_DKMPPD
                  + R_DPCOAK + 3 R_UACMAMO + R_FRUK + R_PNTK + R_ADNK1
                  - R_IDOND2 + R_RBK_L1 - R_2DGULRx - R_2DGULRy + R_ACCOACr
                  - R_DHDPRy - R_3HPPPNH + R_UDCPDP + R_HSK + R_PHEt2r + R_FAO4
                  + R_DAPabc + R_DUTPDP + R_DURIK1 + R_ASPt2 + R_M1PD
                  + R_GLNabc + R_TAURabc + R_THMDt2r + R_PFK + R_GNK
                  + R_DDGLCNt2r - R_KARA1i + R_SHCHD2 - R_GMHEPAT - R_NMNAT
                  - R_LCARS + R_AKGt2r + R_MANAO + R_GSNK - 2 R_ADHEr
                  - 14 R_C120SN + R_DNTPPA + R_DXYLK + 6 R_DKMPPD2 + R_OCBT
                  + R_ADOCBLS - R_CINNDO - R_ACBIPGT + R_TRPt2r - 2 R_GMPR
                  + R_DGSNt2 - R_IPDPS + R_ILEt2r - 2 R_HYD2 - 2 R_HYD3
                  - 2 R_HYD1 + R_PPNCL2 + R_OMBZLM + R_PPTGS + 6 R_ENTCS
                  + R_TSULabc - R_GLGC - R_DPR - R_OMCDC + 3 R_SUCCt2_3
                  + 2 R_SSALy + 2 R_SSALx - R_3HCINNMH + 2 R_GCALDD + R_DMQMT
                  - R_IGPS + R_HETZK + R_MDH - R_DOGULNR - R_BSORy + R_LCTSt
                  - R_ARGDC + R_GLUCYS + R_NO2t2r + R_PGL + R_PTRCabc
                  - 2 R_UGLYCH + R_BUTt2r + R_AST + R_Kabc + R_GUAt2 + R_TAUDO
                  + R_GALt2 + 45.56035 R_BiomassEcoli + R_AGPR + R_ABUTt2
                  + R_CHLt2r - R_KARA2i - R_TRDR + R_FE2abc + R_MCITS + R_CSNt2
                  + 3 R_UDPGD - R_2DGLCNRy - R_2DGLCNRx + R_GLUabc + R_NTPP1
                  + R_LDH_D + R_MALTabc - R_APRAUR + R_METDabc + R_NMNN
                  + R_DALAt2r + R_PROabc + R_NADS1 - R_SERASr + R_DB4PS
                  + R_SUCCabc + R_MALS + R_ASAD - R_TMAOR2 + R_URAt2r
                  + R_ADEt2r + R_DDGLK - R_PPNDH + R_UAMAS + 2 R_GP4GH
                  + 4 R_AMPMS + R_NNDMBRT + R_ACALDi + R_DSERt2r + R_CDPMEK
                  + R_UAMAGS + R_GLCRt2r + R_IMPD + 10 R_LPSSYN_EC - R_MTHFR2
                  + R_ASNabc + R_GPDDA5 + R_GPDDA4 + R_GPDDA3 + R_GPDDA2
                  + R_GPDDA1 + R_GLCt2 + R_AMANK + R_GAPD + R_NTPP8 + R_ALLTN
                  + R_NTPP6 + R_NTPP7 + R_NTPP4 + R_NTPP5 + R_NTPP2 + R_NTPP3
                  - 19 R_C161SN - R_PYK + R_MALTTRabc - 2 R_CPPPGO + R_GLYOX
                  + R_DHORTS + R_PYDXK + R_HCO3E + R_ACGS + R_TDPADGAT
                  + R_PGSA_EC - R_GOFUCR + 2 R_PPS + R_HEX1 + R_PPC + R_PPA
                  - R_SUCCt2b - R_MACPD - R_TDPDRR + 2 R_ADSS + R_GLYCLTt2r
                  + 2 R_THD2 + 3 R_ASPt2_3 + 2 R_ASPt2_2 + R_TRE6PS + R_ALAt2r
                  - R_ACLS + R_INDOLEt2r + R_ASNt2r + R_GLCS1 + R_MALTPTabc
                  - 20 R_C160SN + R_THZPSN + R_SPMDabc + R_PPBNGS + R_CYTDK2
                  + R_RBFK + R_DAGK_EC + R_UPP3MT + R_DHBD + R_MALTTTRabc
                  - R_PEPT_EC + R_DDGALK - 22 R_C181SN - 2 R_P5CR + R_GALCTt2r
                  + R_CYNTt2 - R_GTHOr + R_G3PD2 + R_THRD + R_AIRC2
                  - 3 R_CYNTAH + R_TYRt2r + R_P5CD + R_ADNt2r + 2 R_SGSAD
                  + R_NMNt7 + R_OCDCAt2 - R_DXPS - R_GALU + R_CYANST
                  + R_LPLIPA2 + R_LPLIPA3 + R_LPLIPA1 - R_DAPDC + R_URIt2r
                  - 2 R_CYTBD + R_METabc + R_ETOHt2r + R_XYLK - 2 R_SADH
                  + R_XYLabc - 2 R_NO3R2 - 2.5 R_CYTBO3 - 2 R_NO3R1 + R_GLYt2r
                  + 3 R_FUMt2_3 + 2 R_FUMt2_2 + R_XANt2 + R_XTSNt2r
                  + 2 R_ALDD19x - 3 R_FDH3 - 3 R_FDH2 + R_PIabc + R_ADSK
                  + R_MTRK + 2 R_BETALDHy + 2 R_BETALDHx + R_LEUt2r - R_MECDPDH
                  - R_PPPNDO + R_RMK - R_DKGLCNR2y - R_MAN1PT2 - R_MEPCT
                  + 2 R_HKNDDH + R_PROD2 + R_UAAGDS + 2 R_AP5AH - 17 R_C140SN
                  + R_PDX5PS - R_CSND - R_UAPGR + R_DADNt2 - R_FHL + R_ARGabc
                  - R_TMPPP + R_QULNS - R_DKGLCNR1 + R_GMHEPK + R_TTDCAt2
                  + R_GALURt2r + R_ALAabc + R_NADK + 3 R_ATPS4r + R_SULabc
                  - R_CYTD + R_PERD + R_TAGURr + R_GLCURt2r - R_ADD + R_PPPNt2r
                  + R_ARBt2r + R_MTHFD + R_GTHS + R_PGCD + R_HDCAt2 + R_GSPMDS
                  - R_OXGDC2 - R_OPHBDC + 2 R_AP4AH + 2 R_HKNTDH + R_GALCTNt2r
                  - R_KG6PDC + 3 R_MALt2_3 + 2 R_MALt2_2 - R_TMAOR1 + R_ADOCBIK
                  + R_GLNS + R_GLYCDx + R_VALt2r + R_UGMDDS + 3 R_NAt3_15
                  + R_ACt2r - R_G1PTT + R_HISt2r + R_RBK + R_DINSt2 + R_INSt2r
                  - 4 R_UPPDC1 + R_DCYTt2 + 2 R_ABUTD + R_GLYC3Pabc - R_OMPDC
                  + R_LYSabc + R_ACNAMt2 - 5 R_NTRIR2x + R_GLCNt2r + R_IDONt2r
                  + R_CYSS + 2 R_CBPS + R_PYDAMK + R_MOAT2 + R_DHPPD - R_NADH10
                  - R_KAS14 - 2 R_KAS16 + R_PYRt2r - R_5DGLCNR + 3 R_SHCHF
                  + R_IPPS + R_XYLt2 + R_DHDPS + R_PRFGS + R_FUCt + R_ASNS1
                  + R_URIK2 - R_ALCD19 + R_HMPK1 + R_THMabc + R_ACACt2
                  - R_GLXCL + R_RIBabc + R_PFK_2 - R_DHFR + R_CS + R_HCINNMt2r
                  + R_HISabc - 2 R_NNDPR + R_ILEabc - R_DCTPD - R_GLUSy
                  - R_TRSAR + R_GLYBabc + R_GLYK + R_GALKr + R_IG3PS
                  + R_HPPPNt2r + R_THRt2r + 3 R_HISTD + 2 R_BTS2 - R_KAS15
                  + R_FTHFD + R_PIt2r + R_CHLabc + 2 R_GMPS2 - R_NNAT
                  - R_DKGLCNR2x + R_ASPCT + R_THRabc + 2 R_GTPCII2 + R_MELIBt2
                  + R_GART + R_GLTPD - R_HPYRRx - R_HPYRRy - R_BSORx - R_NADH5
                  + R_ANS - 3 R_NADH7 - 4.5 R_NADH6 - R_NADH9 - 3.8 R_NADH8
                  + R_PRASCS + 2 R_CDAPPA_EC - R_DADA + R_G1PACT + R_ALAALAr
                  + 2 R_E4PD + 6 R_FAO1 + R_CADVt + 7 R_FAO2 - R_UAGDP
                  + R_SHSL1 + R_LPADSS - R_ADMDCr + R_ARGSS + R_HPPK2
                  + 3 R_DBTSr + R_D_LACt2 + R_SPMDAT1 + R_SPMDAT2 + R_LEUabc
                  + R_ECAP_EC - R_PTPATi + R_L_LACt2r + R_ASNS2 + R_ASPabc
                  + R_PTRCt2r - R_ADA + R_HSDy + 5 R_SULR - R_GUAD + R_G6PDH2r
                  + 8 R_FAO3 - R_MMCD + R_PRATPP - R_DXPRIi + R_GARFT - R_DCYTD
                  + 2 R_PRAIS  = 0
 M_h_e:           - R_PROt2r - R_RMNt - R_GLYBt2r - R_SERt2r - R_Kt2r
                  - 2 R_SUCCt2_2 - R_ALLTNt2r - R_GSNt2 - R_NAt3_1 - 2 R_NAt3_2
                  - R_LYSt2r - R_DURIt2 - R_GLUt2r - R_CYTDt2r - R_PHEt2r
                  - R_ASPt2 - R_THMDt2r - R_DDGLCNt2r - R_AKGt2r - R_TRPt2r
                  - R_DGSNt2 - R_ILEt2r + 2 R_HYD2 + 2 R_HYD3 + 2 R_HYD1
                  - 3 R_SUCCt2_3 - R_LCTSt - R_NO2t2r - R_BUTt2r - R_GUAt2
                  - R_GALt2 - R_ABUTt2 - R_CHLt2r - R_CSNt2 - R_DALAt2r
                  - R_URAt2r - R_ADEt2r - R_DSERt2r - R_GLCRt2r - R_GLCt2
                  + R_SUCCt2b - R_GLYCLTt2r - 2 R_THD2 - 3 R_ASPt2_3
                  - 2 R_ASPt2_2 - R_ALAt2r - R_INDOLEt2r - R_ASNt2r - R_TMAOR1e
                  - R_GALCTt2r - R_CYNTt2 - R_TYRt2r - R_ADNt2r - R_OCDCAt2
                  - R_URIt2r + 2 R_CYTBD - R_ETOHt2r + 2 R_NO3R2 + 2.5 R_CYTBO3
                  + 2 R_NO3R1 - R_GLYt2r - 3 R_FUMt2_3 - 2 R_FUMt2_2 - R_XANt2
                  - R_XTSNt2r + 2 R_FDH3 + 2 R_FDH2 - R_TMAOR2e - R_LEUt2r
                  - R_DADNt2 - R_TTDCAt2 - R_GALURt2r - 4 R_ATPS4r - R_GLCURt2r
                  - R_PPPNt2r - R_ARBt2r - R_HDCAt2 - R_GALCTNt2r - 3 R_MALt2_3
                  - 2 R_MALt2_2 - R_VALt2r - 3 R_NAt3_15 - R_ACt2r - R_HISt2r
                  - R_DINSt2 - R_INSt2r - R_DCYTt2 - R_ACNAMt2 - R_GLCNt2r
                  - R_IDONt2r - R_PYRt2r - R_XYLt2 - R_FUCt - R_ACACt2
                  + 2 R_NADDPe - R_HCINNMt2r - R_HPPPNt2r - R_THRt2r - R_PIt2r
                  - R_EX_h_e - R_MELIBt2 + 2 R_NADH7 + 3.5 R_NADH6
                  + 2.8 R_NADH8 + R_GLCDe - R_CADVt - R_D_LACt2 - R_L_LACt2r
                  - R_PTRCt2r  = 0
 M_hco3_c:        - R_ACCOACr + R_HCO3E - R_AIRC2 - R_CYNTAH - R_CBPS  = 0
 M_hcys_L_c:      R_RHCCE - R_METS + R_AHC + R_CYSTL  = 0
 M_hdca_c:        0.36 R_PLIPA3 + 0.36 R_PLIPA1 + 0.36 R_PLIPA2
                  + 0.36 R_LPLIPA2 + 0.36 R_LPLIPA3 + 0.36 R_LPLIPA1 - R_AACPS3
                  + R_HDCAt2 - R_FAO2  = 0
 M_hdca_e:        - R_HDCAt2 - R_EX_hdca_e  = 0
 M_hdcea_c:       0.07 R_PLIPA3 - R_AACPS4 + 0.07 R_PLIPA1 + 0.07 R_PLIPA2
                  + 0.07 R_LPLIPA2 + 0.07 R_LPLIPA3 + 0.07 R_LPLIPA1  = 0
 M_hdeACP_c:      R_AACPS4 - R_EDTXS3 + R_C161SN - 0.14 R_PASYN_EC  = 0
 M_hemeO_c:       R_HEMEOS  = 0
 M_his_L_c:       - 0.09 R_BiomassEcoli + R_HISt2r + R_HISabc + R_HISTD  = 0
 M_his_L_e:       - R_HISt2r - R_HISabc - R_EX_his_L_e  = 0
 M_hisp_c:        R_HSTPT - R_HISTP  = 0
 M_histd_c:       R_HISTP - R_HISTD  = 0
 M_hkndd_c:       R_HPPPNDO - R_HKNDDH  = 0
 M_hkntd_c:       R_DHCINDO - R_HKNTDH  = 0
 M_hmbil_c:       R_HMBS - R_UPP3S  = 0
 M_hmfurn_c:      R_DHPTDCs  = 0
 M_hom_L_c:       - R_HSK - R_HSST - R_HSDy  = 0
 M_hpyr_c:        - R_HPYRRx - R_HPYRRy - R_HPYRI  = 0
 M_hqn_c:         R_AB6PGH  = 0
 M_hxan_c:        R_HYXNt + R_PUNP5 + R_PUNP6 + R_ADD - R_HXPRT  = 0
 M_hxan_e:        - R_HYXNt - R_EX_hxan_e  = 0
 M_iasp_c:        R_ASPO4 + R_ASPO6 + R_ASPO3 - R_QULNS + R_ASPO5  = 0
 M_ichor_c:       - R_SHCHCS2 - R_ICHORT + R_ICHORSi  = 0
 M_icit_c:        - R_ICL + R_ACONT - R_ICDHyr  = 0
 M_idon_L_c:      R_IDOND + R_IDOND2 + R_2DGULRx + R_2DGULRy + R_IDONt2r  = 0
 M_idon_L_e:      - R_EX_idon_L_e - R_IDONt2r  = 0
 M_idp_c:         R_ADK4  = 0
 M_ile_L_c:       R_ILEt2r - R_ILETA - 0.276 R_BiomassEcoli + R_ILEabc  = 0
 M_ile_L_e:       - R_ILEt2r + R_EXinput_ile_L_e - R_ILEabc
                  - R_EXoutput_ile_L_e  = 0
 M_imacp_c:       - R_HSTPT + R_IGPDH  = 0
 M_imp_c:         R_INSK + R_GMPR - R_IMPC - R_IMPD - R_ADSS - R_NTD11
                  + R_HXPRT  = 0
 M_indole_c:      R_TRPAS2 + R_TRPS3 - R_TRPS2 + R_INDOLEt2r  = 0
 M_indole_e:      - R_INDOLEt2r - R_EX_indole_e  = 0
 M_inost_c:       R_GPDDA5 + R_MI1PP  = 0
 M_ins_c:         - R_INSK - R_PUNP5 + R_NTD11 + R_INSt2r + R_ADA  = 0
 M_ins_e:         - R_INSt2r - R_EX_ins_e  = 0
 M_ipdp_c:        - R_GRTT + R_IPDPS - R_IPDDIi - 5 R_OCTDPS - 8 R_UDCPDPS
                  - R_DMATT  = 0
 M_itp_c:         - R_ADK4  = 0
 M_k_c:           R_Kt2r + R_Kabc  = 0
 M_k_e:           R_EXinput_k_e - R_Kt2r - R_Kabc - R_EXoutput_k_e  = 0
 M_kdo2lipid4L_c: - R_EDTXS2 + R_EDTXS1  = 0
 M_kdo2lipid4_c:  - R_EDTXS3 - R_EDTXS1 + R_MOAT2  = 0
 M_kdo2lipid4p_c: - R_EDTXS4 + R_EDTXS3  = 0
 M_kdo8p_c:       R_KDOPS - R_KDOPP  = 0
 M_kdo_c:         - R_KDOCT2 + R_KDOPP  = 0
 M_kdolipid4_c:   R_MOAT - R_MOAT2  = 0
 M_lac_D_c:       - R_LDH_D + R_GLYOX - R_LDH_D2 + R_D_LACt2  = 0
 M_lac_D_e:       - R_EX_lac_D_e - R_D_LACt2  = 0
 M_lac_L_c:       R_LCAD + R_L_LACt2r - R_L_LACD3 - R_L_LACD2  = 0
 M_lac_L_e:       - R_EX_lac_L_e - R_L_LACt2r  = 0
 M_lald_L_c:      - R_LCAD + R_RMPA + R_FCLPA - R_LCARS  = 0
 M_lcts_c:        - R_LACZ + R_LCTSt  = 0
 M_lcts_e:        - R_LCTSt - R_EX_lcts_e  = 0
 M_leu_L_c:       R_LEUTAi - 0.428 R_BiomassEcoli + R_LEUt2r + R_LEUabc  = 0
 M_leu_L_e:       - R_EX_leu_L_e - R_LEUt2r - R_LEUabc  = 0
 M_lgt_S_c:       R_LGTHL - R_GLYOX  = 0
 M_lipa_c:        R_EDTXS2 - R_LPSSYN_EC  = 0
 M_lipa_cold_c:   R_EDTXS4  = 0
 M_lipidA_c:      - R_MOAT + R_TDSK  = 0
 M_lipidAds_c:    - R_TDSK + R_LPADSS  = 0
 M_lipidX_c:      R_USHD - R_LPADSS  = 0
 M_lps_EC_c:      - 0.0084 R_BiomassEcoli + R_LPSSYN_EC  = 0
 M_lys_L_c:       - R_LYSDC + R_LYSt2r - 0.326 R_BiomassEcoli + R_DAPDC
                  + R_LYSabc + R_CADVt  = 0
 M_lys_L_e:       - R_EXoutput_lys_L_e - R_LYSt2r - R_LYSabc - R_CADVt
                  + R_EXinput_lys_L_e  = 0
 M_malACP_c:      - 5 R_C141SN - 4 R_C120SN - 6 R_C161SN - R_MACPD - 6 R_C160SN
                  - 7 R_C181SN - 5 R_C140SN - R_KAS14 - R_KAS16 - R_KAS15
                  + R_MCOATA  = 0
 M_mal_L_c:       - R_MDH + R_MALS + R_FUM - R_MDH2 - R_MDH3 + R_MALt2_3
                  + R_MALt2_2 - R_ME1 - R_ME2  = 0
 M_mal_L_e:       - R_MALt2_3 - R_MALt2_2 - R_EX_mal_L_e  = 0
 M_malcoa_c:      R_ACCOACr - R_MCOATA  = 0
 M_malt6p_c:      R_MALTpts  = 0
 M_malt_c:        R_MLTG1 + R_MALTabc - R_AMALT1 - R_AMALT3 - R_AMALT4
                  - R_AMALT2  = 0
 M_malt_e:        - R_MALTpts - R_MALTabc - R_EX_malt_e  = 0
 M_malthp_c:      - R_MLTP3 - R_MLTG5 + R_AMALT4  = 0
 M_malthx_c:      R_MALTHXabc + R_MLTP3 - R_MLTP2 - R_MLTG4 + R_MLTG5
                  + R_AMALT3 - R_AMALT4  = 0
 M_malthx_e:      - R_MALTHXabc - R_EX_malthx_e  = 0
 M_maltpt_c:      - R_MLTP1 + R_MLTP2 + R_MALTPTabc - R_MLTG3 + R_MLTG4
                  - R_AMALT3 + R_AMALT2  = 0
 M_maltpt_e:      - R_MALTPTabc - R_EX_maltpt_e  = 0
 M_malttr_c:      - R_MLTG1 + R_MALTTRabc + R_MLTG2 - R_AMALT1  = 0
 M_malttr_e:      R_EXinput_malttr_e - R_MALTTRabc - R_EXoutput_malttr_e  = 0
 M_maltttr_c:     R_MLTP1 + R_MALTTTRabc - R_MLTG2 + R_MLTG3 + R_AMALT1
                  - R_AMALT2  = 0
 M_maltttr_e:     - R_MALTTTRabc - R_EX_maltttr_e  = 0
 M_man1p_c:       - R_MAN1PT2 - R_PMANM  = 0
 M_man6p_c:       R_MAN6Pt6_2 + R_MANpts - R_MAN6PI + R_PMANM  = 0
 M_man6p_e:       - R_MAN6Pt6_2 - R_EXoutput_man6p_e + R_EXinput_man6p_e  = 0
 M_man_e:         - R_MANpts - R_EX_man_e  = 0
 M_mana_c:        - R_MANAO - R_MNNH  = 0
 M_melib_c:       - R_GALS3 + R_MELIBt2  = 0
 M_melib_e:       - R_EX_melib_e - R_MELIBt2  = 0
 M_met_D_c:       R_METDabc  = 0
 M_met_D_e:       - R_METDabc - R_EX_met_D_e  = 0
 M_met_L_c:       - 0.146 R_BiomassEcoli - R_METAT + R_METabc + R_METS + R_UNK3
                   = 0
 M_met_L_e:       - R_METabc - R_EX_met_L_e  = 0
 M_methf_c:       R_MTHFD - R_MTHFC  = 0
 M_mi1p_D_c:      - R_MI1PP  = 0
 M_micit_c:       - R_MCITL2 + R_MICITD  = 0
 M_mlthf_c:       - R_TMDS + R_GLYCL + R_GHMT2 - R_MTHFR2 - R_MOHMT - R_MTHFD
                   = 0
 M_mmcoa_R_c:     - R_MME + R_MMM2  = 0
 M_mmcoa_S_c:     R_MME - R_MMCD  = 0
 M_mnl1p_c:       R_MNLpts - R_M1PD  = 0
 M_mnl_e:         - R_MNLpts - R_EX_mnl_e  = 0
 M_mql8_c:        R_DHORD5 - R_DMSOR1 + R_GLYCTO3 + R_HYD2 + R_ASPO4
                  - R_DMSOR1e - R_TMAOR1e - R_NO3R2 + R_FDH3 + R_MDH3
                  - R_TMAOR1 + R_NADH10 + R_NADH7 + R_L_LACD3 + R_G3PD6
                  - R_FRD2  = 0
 M_mqn8_c:        R_AMMQT8_2 - R_DHORD5 + R_DMSOR1 - R_GLYCTO3 - R_HYD2
                  - R_ASPO4 + R_DMSOR1e + R_TMAOR1e + R_NO3R2 - R_FDH3 - R_MDH3
                  + R_TMAOR1 - R_NADH10 - R_NADH7 - R_L_LACD3 - R_G3PD6
                  + R_FRD2  = 0
 M_mthgxl_c:      - R_LGTHL + R_MGSA  = 0
 M_myrsACP_c:     R_AACPS1 - R_EDTXS4 - R_EDTXS2 + R_C140SN - 0.04 R_PASYN_EC
                   = 0
 M_n8aspmd_c:     R_SPMDAT2  = 0
 M_na1_c:         - R_NAt3_1 - R_NAt3_2 + R_PROt4 + R_SERt4 + R_GLUt4 + R_THRt4
                  - 2 R_NAt3_15 + R_PNTOt4  = 0
 M_na1_e:         R_NAt3_1 + R_NAt3_2 - R_PROt4 - R_EX_na1_e - R_SERt4
                  - R_GLUt4 - R_THRt4 + 2 R_NAt3_15 - R_PNTOt4  = 0
 M_nac_c:         R_NNAM - R_NAMNPP + R_NNDMBRT + R_NACUP  = 0
 M_nac_e:         - R_EX_nac_e - R_NACUP  = 0
 M_nad_c:         - R_ALDD2x - R_IPMD + R_IDOND - R_GLYCL - R_LCAD - R_DHCIND
                  - R_SBTPD - R_PPND - R_PDH - R_NADDP + R_GLYCLTDx + R_DMPPS
                  - 2 R_UACMAMO + R_2DGULRx + R_3HPPPNH - R_FAO4 - R_M1PD
                  - R_AKGDH - R_SHCHD2 + R_NMNAT + R_LCARS - R_MANAO
                  + 2 R_ADHEr + R_CINNDO + R_IPDPS - R_SSALx + R_3HCINNMH
                  - R_GCALDD - R_MDH + R_DOGULNR - 0.00215 R_BiomassEcoli
                  - R_NADTRHD - 2 R_UDPGD + R_2DGLCNRx + R_NADPPPS - R_LDH_D
                  + R_NADS1 - R_ACALDi - R_IMPD + R_MTHFR2 - R_GAPD + R_THD2
                  - R_DHBD - R_THRD - R_P5CD - R_SGSAD - R_ALDD19x - R_BETALDHx
                  + R_PPPNDO - R_PDX5PS - R_NADK - R_PERD - R_TAGURr - R_PGCD
                  - R_GLYCDx - R_ME1 - R_ABUTD + 3 R_NTRIR2x - R_DHPPD
                  + R_NADH10 + R_ALCD19 + R_TRSAR - 2 R_HISTD + R_DKGLCNR2x
                  - R_GLTPD + R_HPYRRx + R_BSORx + R_NADH5 + R_NADH7 + R_NADH6
                  + R_NADH9 + R_NADH8 - R_E4PD - 6 R_FAO1 - 7 R_FAO2 - 8 R_FAO3
                   = 0
 M_nad_e:         - R_NADDPe - R_EX_nad_e  = 0
 M_nadh_c:        R_ALDD2x + R_IPMD - R_IDOND + R_GLYCL + R_LCAD + R_DHCIND
                  + R_SBTPD + R_PPND + R_PDH - R_GLYCLTDx - R_DMPPS
                  + 2 R_UACMAMO - R_2DGULRx - R_3HPPPNH + R_FAO4 + R_M1PD
                  + R_AKGDH + R_SHCHD2 - R_LCARS + R_MANAO - 2 R_ADHEr
                  - R_CINNDO - R_IPDPS + R_SSALx - R_3HCINNMH + R_GCALDD
                  + R_MDH - R_DOGULNR - 0.00005 R_BiomassEcoli + R_NADTRHD
                  + 2 R_UDPGD - R_2DGLCNRx + R_LDH_D + R_ACALDi + R_IMPD
                  - R_MTHFR2 + R_GAPD - R_THD2 + R_DHBD + R_THRD + R_P5CD
                  + R_SGSAD + R_ALDD19x + R_BETALDHx - R_PPPNDO + R_PDX5PS
                  + R_PERD + R_TAGURr + R_PGCD + R_GLYCDx + R_ME1 + R_ABUTD
                  - 3 R_NTRIR2x + R_DHPPD - R_NADH10 - R_ALCD19 - R_TRSAR
                  + 2 R_HISTD - R_DKGLCNR2x + R_GLTPD - R_HPYRRx - R_BSORx
                  - R_NADH5 - R_NADH7 - R_NADH6 - R_NADH9 - R_NADH8 + R_E4PD
                  + 6 R_FAO1 + 7 R_FAO2 + 8 R_FAO3  = 0
 M_nadp_c:        R_SHK3Dr + R_GLYCLTDy - R_GLUDy + R_GLUTRR + R_G5SD
                  + 11 R_C141SN + R_IDOND2 + R_2DGULRy + R_DHDPRy + R_KARA1i
                  - R_GND + 10 R_C120SN + R_GMPR + R_DPR - R_SSALy + R_BSORy
                  - 0.00013 R_BiomassEcoli - R_AGPR + R_KARA2i + R_TRDR
                  + R_NADTRHD + R_2DGLCNRy - R_NADPPPS + R_APRAUR - R_ASAD
                  + 13 R_C161SN + R_GOFUCR + R_TDPDRR - R_THD2 + 14 R_C160SN
                  + 15 R_C181SN + R_P5CR + R_GTHOr - R_G3PD2 - R_BETALDHy
                  + R_DKGLCNR2y - R_ICDHyr + 12 R_C140SN + R_UAPGR + R_DKGLCNR1
                  + R_NADK - R_MTHFD - R_ME2 + R_KAS16 + R_5DGLCNR + R_DHFR
                  + R_GLUSy + R_HPYRRy - R_HSDy - 3 R_SULR - R_G6PDH2r
                  + R_DXPRIi  = 0
 M_nadph_c:       - R_SHK3Dr - R_GLYCLTDy + R_GLUDy - R_GLUTRR - R_G5SD
                  - 11 R_C141SN - R_IDOND2 - R_2DGULRy - R_DHDPRy - R_KARA1i
                  + R_GND - 10 R_C120SN - R_GMPR - R_DPR + R_SSALy - R_BSORy
                  - 0.0004 R_BiomassEcoli + R_AGPR - R_KARA2i - R_TRDR
                  - R_NADTRHD - R_2DGLCNRy - R_APRAUR + R_ASAD - 13 R_C161SN
                  - R_GOFUCR - R_TDPDRR + R_THD2 - 14 R_C160SN - 15 R_C181SN
                  - R_P5CR - R_GTHOr + R_G3PD2 + R_BETALDHy - R_DKGLCNR2y
                  + R_ICDHyr - 12 R_C140SN - R_UAPGR - R_DKGLCNR1 + R_MTHFD
                  + R_ME2 - R_KAS16 - R_5DGLCNR - R_DHFR - R_GLUSy - R_HPYRRy
                  + R_HSDy + 3 R_SULR + R_G6PDH2r - R_DXPRIi  = 0
 M_ncam_c:        - R_NNAM + R_NMNN + R_NMNt7  = 0
 M_nh4_c:         R_ASNN + R_GLYCL + R_SERD_D + R_SERD_L + R_GLUDy - R_CBMK
                  + R_DHPPDA2 + R_NNAM + R_G6PDA + R_GLUN + R_TRPAS2 + R_DAAD
                  + R_GMPR + 2 R_UGLYCH - R_NADS1 + R_PEAMNO + 4 R_HMBS
                  + R_NH4t + R_CYNTAH + R_PYAM5PO + 2 R_SADH + R_ETHAAL
                  + R_CSND + R_CYSDS + R_CYTD + R_ADD + R_THRD_L - R_GLNS
                  + R_NTRIR2x + R_DCTPD + R_NMNDA + R_DADA - R_ASNS2 + R_CYSTL
                  + R_ADA + R_GUAD + R_ASPT + R_DCYTD  = 0
 M_nh4_e:         - R_EXoutput_nh4_e + R_EXinput_nh4_e - R_NH4t  = 0
 M_nicrnt_c:      R_NAMNPP - R_NNDMBRT + R_NNDPR + R_NMNDA - R_NNAT  = 0
 M_nmn_c:         R_NADDP - R_NMNAT + R_NMNP - R_NMNN - R_NMNDA  = 0
 M_nmn_e:         - R_NMNP - R_EX_nmn_e - R_NMNt7 + R_NADDPe  = 0
 M_no2_c:         - R_NO3t7 + R_NO2t2r + R_NO3R2 + R_NO3R1 - R_NTRIR2x  = 0
 M_no2_e:         R_NO3t7 - R_NO2t2r - R_EX_no2_e  = 0
 M_no3_c:         R_NO3t7 - R_NO3R2 - R_NO3R1  = 0
 M_no3_e:         - R_NO3t7 - R_EX_no3_e  = 0
 M_o2_c:          - 0.5 R_OMMBLHX + R_CAT - 0.5 R_OMPHHX + R_SPODM + R_O2t
                  - R_DKMPPD - R_3HPPPNH - R_PDX5PO - R_CINNDO - R_ASPO6
                  - R_DHCINDO - R_3HCINNMH - R_TAUDO - R_CPPPGO - R_PEAMNO
                  - R_PYAM5PO - 0.5 R_CYTBD - 0.5 R_CYTBO3 - R_PPPNDO
                  - R_HPPPNDO - 1.5 R_PPPGO - 0.5 R_OPHHX  = 0
 M_o2_e:          - R_O2t - R_EX_o2_e  = 0
 M_o2s_c:         - 2 R_SPODM  = 0
 M_oaa_c:         R_ASPTA + R_MDH - R_MCITS + R_TARTD + R_PPC + R_MDH2 + R_MDH3
                  - R_CS + R_CITL - R_PPCK  = 0
 M_ocdca_c:       R_OCDCAt2 - R_FAO3  = 0
 M_ocdca_e:       - R_OCDCAt2 - R_EX_ocdca_e  = 0
 M_ocdcea_c:      0.5 R_PLIPA3 - R_AACPS5 + 0.5 R_PLIPA1 + 0.5 R_PLIPA2
                  + 0.5 R_LPLIPA2 + 0.5 R_LPLIPA3 + 0.5 R_LPLIPA1  = 0
 M_octdp_c:       - R_HBZOPT - R_DHNAOT + R_OCTDPS  = 0
 M_octeACP_c:     R_AACPS5 + R_C181SN - R_PASYN_EC  = 0
 M_ohpb_c:        - R_OHPBAT + R_PERD  = 0
 M_op4en_c:       - R_OP4ENH + R_HKNDDH + R_HKNTDH  = 0
 M_orn_c:         R_ORNabc + R_ACODA - R_ORNDC - R_OCBT - R_PTRCORNt7
                  - R_ARGORNt7 - R_ORNTA  = 0
 M_orn_e:         - R_ORNabc - R_EX_orn_e + R_PTRCORNt7 + R_ARGORNt7  = 0
 M_orot5p_c:      - R_OMPDC - R_ORPT  = 0
 M_orot_c:        R_DHORD5 + R_ORPT + R_DHORD2  = 0
 M_pa_EC_c:       - 0.02 R_DASYN_EC - 0.02 R_PAPA_EC + 0.02 R_DAGK_EC
                  + 0.02 R_PASYN_EC + 0.02 R_CDAPPA_EC  = 0
 M_pac_c:         R_ALDD19x - R_PACCOAL  = 0
 M_pacald_c:      R_PEAMNO - R_ALDD19x  = 0
 M_palmACP_c:     R_C160SN + R_AACPS3 - 0.72 R_PASYN_EC  = 0
 M_pan4p_c:       R_PPCDC - R_PTPATi  = 0
 M_pant_R_c:      - R_PANTS + R_DPR  = 0
 M_pap_c:         R_PAPSR + R_ACPS1 - R_BPNT  = 0
 M_paps_c:        - R_PAPSR + R_ADSK  = 0
 M_pc_EC_c:       - 0.02 R_PLIPA3  = 0
 M_pdx5p_c:       R_PYDXNK - R_PDX5PO + R_PDX5PS - R_PDXPP  = 0
 M_pe_EC_c:       0.02 R_PSD_EC - 0.02 R_PLIPA2 - 0.001935 R_BiomassEcoli
                  - 0.02 R_PEPT_EC  = 0
 M_peamn_c:       - R_PEAMNO  = 0
 M_pep_c:         - R_GLCpts - R_FRUpts2 - R_MNLpts - R_SUCpts - R_DHAPT
                  - R_MALTpts - R_UAGCVT - R_MANpts + R_ENO - R_SBTpts
                  - R_KDOPS - R_FRUpts - R_PYK + R_PPS - R_TREpts - R_PPC
                  - R_PSCVT - R_ACGApts - R_GALTpts - R_DDPA - R_GAMpts
                  - R_ACMANApts + R_PPCK  = 0
 M_peptido_EC_c:  R_PPTGS - 0.0276 R_BiomassEcoli  = 0
 M_pg_EC_c:       - 0.04 R_CLPNS_EC - 0.02 R_LPLIPA6 - 0.02 R_PLIPA1
                  + 0.02 R_PGPP_EC - 0.000464 R_BiomassEcoli - 0.02 R_LPLIPA5
                  - 0.02 R_LPLIPA4  = 0
 M_pgp_EC_c:      - 0.02 R_PGPP_EC + 0.02 R_PGSA_EC  = 0
 M_phaccoa_c:     R_PACCOAL  = 0
 M_phe_L_c:       R_PHEt2r - R_PHETA1 - 0.176 R_BiomassEcoli  = 0
 M_phe_L_e:       - R_PHEt2r - R_EX_phe_L_e  = 0
 M_pheme_c:       R_FCLT - R_HEMEOS  = 0
 M_phom_c:        R_HSK - R_THRS  = 0
 M_phpyr_c:       R_PHETA1 + R_PPNDH  = 0
 M_phthr_c:       - R_4HTHRS - R_PDX5PS + R_OHPBAT  = 0
 M_pi_c:          R_CTPS2 - R_GLCP - 2 R_MAN6Pt6_2 + R_PYDXPP + R_ORNabc
                  + R_PRAGSr + R_ARBabc + R_ATPM + R_NTD9 + R_NTD8 + R_NTD1
                  + R_NTD3 + R_NTD2 + R_NTD5 + R_VALabc + R_GALabc + R_CYSabc
                  + R_CBL1abc + R_G5SD + R_MALTHXabc + R_DNMPPA + R_DKMPPD
                  + R_UAGCVT + R_PGPP_EC + R_ACCOACr + R_UDCPDP + R_ACCOAL
                  + R_DAPabc + R_GLNabc + R_TAURabc - R_MLTP3 + R_DKMPPD2
                  + R_NAMNPP - R_PTAr + R_OCBT + R_TSULabc - R_PTA2 + R_BPNT
                  + R_GLUCYS + R_PTRCabc + R_Kabc - R_DURIPP
                  + 45.5628 R_BiomassEcoli - R_AGPR + R_FE2abc - R_MLTP1
                  - R_MLTP2 + R_GLUabc + R_KDOPS + R_GMHEPPA + R_THRS
                  + R_NADPPPS + R_CBLAT_DELETE + R_MALTabc + R_METDabc
                  + R_PROabc + R_SUCCabc - R_ASAD + R_UAMAS - R_PUNP5 - R_PUNP4
                  - R_PUNP6 - R_PUNP1 - R_PUNP2 + R_PAPA_EC + R_SADT2
                  + R_UAMAGS - R_TMDPP - R_PUNP7 + R_ASNabc - R_GAPD + R_METAT
                  + R_MALTTRabc + R_RBFSa + R_PPS + R_4HTHRS + R_PPC + 2 R_PPA
                  + R_PSCVT + R_ADSS + R_TRE6PP + R_CBIAT_DELETE + R_DHQS
                  + R_G1PP + R_MALTPTabc - R_PUNP3 + R_SPMDabc + R_MALTTTRabc
                  + R_AIRC2 + R_MGSA - R_GLYC3Pt6 + R_METabc - 2 R_G6Pt6_2
                  + R_PMDPHT + R_XYLabc + R_FBP + 2 R_PIabc + R_NTD10 + R_DHFS
                  + R_NTD4 + R_NTD7 + R_MAN1PT2 + R_PGLYCP + R_UAAGDS + R_PSP_L
                  + R_PDX5PS + R_RZ5PP + R_ARGabc + R_DDPA + R_QULNS + R_ALAabc
                  - R_ATPS4r + R_SULabc + R_GTHS + R_GSPMDS + R_GLNS + R_UGMDDS
                  + R_SUCOAS + R_NTD11 + R_NTD6 + R_GLYC3Pabc + R_LYSabc
                  + R_MI1PP + R_CBPS + R_HYPOE + R_PRFGS + R_THMabc + R_CHORS
                  + R_RIBabc + R_HISabc + R_ILEabc + R_GLYBabc + R_HISTP
                  + R_PDXPP + R_PIt2r + R_CHLabc + R_ASPCT + R_THRabc + R_GART
                  + R_PRASCS - 2 R_FUCPt6_2 + R_ALAALAr + R_DBTSr + R_LEUabc
                  + R_ASPabc + R_KDOPP + R_SELNPS - R_PYNP2r + R_PRAIS  = 0
 M_pi_e:          2 R_MAN6Pt6_2 + R_GLYC3Pt6 + 2 R_G6Pt6_2 - R_PIabc
                  - R_EXoutput_pi_e - R_PIt2r + R_EXinput_pi_e + 2 R_FUCPt6_2
                   = 0
 M_pmcoa_c:       - R_AOXSr  = 0
 M_pnto_R_c:      R_PANTS - R_PNTK + R_PNTOt4  = 0
 M_pnto_R_e:      - R_EX_pnto_R_e - R_PNTOt4  = 0
 M_ppa_c:         - R_ACCOAL + R_PPCSCT + R_PPAKr  = 0
 M_ppap_c:        R_PTA2 - R_PPAKr  = 0
 M_ppbng_c:       R_PPBNGS - 4 R_HMBS  = 0
 M_ppcoa_c:       R_ACCOAL - R_PTA2 - R_MCITS - R_PPCSCT + R_OBTFL + R_MMCD
                   = 0
 M_pphn_c:        - R_PPND - R_PPNDH + R_CHORM  = 0
 M_ppi_c:         R_DASYN_EC + R_PANTS + R_DHPS2 + R_HBZOPT + R_DHNAOT
                  + R_FMNAT + R_AACPS5 + R_AACPS4 + R_GUAPRT + R_AACPS2
                  + R_AACPS1 + R_GLUTRS + R_GRTT + R_DUTPDP + R_GMHEPAT
                  + R_NMNAT + R_DNTPPA + R_NAMNPP + R_ACBIPGT + R_PPNCL2
                  + R_GLGC + 0.7302 R_BiomassEcoli + R_NTPP1 + R_SUCBZL
                  + R_CBLAT_DELETE + R_NADS1 + R_SERASr + R_UPPRT + R_DHBSr
                  + R_SADT2 + R_ANPRT + R_METAT + R_NTPP8 + R_NTPP6 + R_NTPP7
                  + R_NTPP4 + R_NTPP5 + R_NTPP2 + R_NTPP3 + R_ATPPRT - R_PPA
                  + R_CBIAT_DELETE + R_THZPSN + R_GALU + R_AACPS3 + 5 R_OCTDPS
                  + R_MEPCT + R_TMPPP + R_PACCOAL + R_XPPT + R_KDOCT2 + R_G1PTT
                  + R_HXPRT - R_ORPT + R_ACS + R_ASNS1 + 8 R_UDCPDPS + R_NNDPR
                  + R_ADPT + R_DMATT + R_GMPS2 + R_NNAT + R_GTPCII2 + R_HEMEOS
                  + R_GLUPRT + R_FAO1 + R_FAO2 + R_UAGDP + R_ARGSS + R_PTPATi
                  + R_ASNS2 + R_ADNCYC + R_FAO3 + R_PRATPP  = 0
 M_ppp9_c:        - R_FCLT + R_PPPGO  = 0
 M_pppg9_c:       R_CPPPGO - R_PPPGO  = 0
 M_pppi_c:        R_NTPTP2 + R_NTPTP1  = 0
 M_pppn_c:        - R_PPPNDO + R_PPPNt2r  = 0
 M_pppn_e:        - R_EX_pppn_e - R_PPPNt2r  = 0
 M_pram_c:        - R_PRAGSr + R_GLUPRT  = 0
 M_pran_c:        - R_PRAIi + R_ANPRT  = 0
 M_prbamp_c:      - R_PRAMPC + R_PRATPP  = 0
 M_prbatp_c:      R_ATPPRT - R_PRATPP  = 0
 M_prfp_c:        R_PRAMPC - R_PRMICIi  = 0
 M_prlp_c:        R_PRMICIi - R_IG3PS  = 0
 M_pro_L_c:       R_PROt2r + R_PROt4 - 0.21 R_BiomassEcoli + R_PROabc + R_P5CR
                  - R_PROD2  = 0
 M_pro_L_e:       - R_PROt2r - R_PROt4 - R_PROabc - R_EX_pro_L_e  = 0
 M_prpp_c:        R_PRPPS - R_GUAPRT - R_NAMNPP - R_UPPRT - R_ANPRT - R_ATPPRT
                  - R_XPPT - R_HXPRT + R_ORPT - R_NNDPR - R_ADPT - R_GLUPRT
                   = 0
 M_ps_EC_c:       0.02 R_PSSA_EC - 0.02 R_PSD_EC - 0.000052 R_BiomassEcoli  = 0
 M_pser_L_c:      R_PSERT - R_PSP_L  = 0
 M_ptrc_c:        - R_SPMS + R_ORNDC + R_PTRCabc - 0.035 R_BiomassEcoli
                  + R_AGMT + R_PTRCORNt7 - R_PTRCTA + R_PTRCt2r  = 0
 M_ptrc_e:        - R_PTRCabc - R_PTRCORNt7 - R_EX_ptrc_e - R_PTRCt2r  = 0
 M_pyam5p_c:      R_ALATA_D2 + R_ALATA_L2 - R_PYAM5PO + R_PYDAMK - R_HYPOE  = 0
 M_pydam_c:       - R_PYDAMK + R_HYPOE  = 0
 M_pydx5p_c:      - R_PYDXPP + R_PDX5PO - R_ALATA_D2 - R_ALATA_L2 + R_PYDXK
                  + R_PYAM5PO  = 0
 M_pydx_c:        R_PYDXPP - R_PYDXK  = 0
 M_pydxn_c:       - R_PYDXNK + R_PDXPP  = 0
 M_pyr_c:         R_ACNML + R_GLCpts + R_FRUpts2 + R_ADCL + R_MNLpts + R_SUCpts
                  + R_EDA + R_DHAPT + R_SERD_D + R_CHRPL + R_SERD_L + R_MALTpts
                  - R_PDH - R_ACHBS + R_SHCHCS2 + R_ICHORT + R_MCITL2
                  + R_MANpts + R_TRPAS2 + R_DAAD + R_SBTpts + R_GLCRAL
                  + R_ALATA_D2 + R_LDH_D + R_FRUpts + R_ALATA_L2 + R_PYK
                  - R_PPS + R_TREpts - 2 R_ACLS + R_ACGApts + R_ALATA_L
                  + R_HOPNTAL - R_DXPS - R_POX + R_LDH_D2 + R_GALTpts - R_PFL
                  + R_CYSDS + R_VPAMT + R_ME1 + R_ME2 + R_GAMpts + R_PYRt2r
                  - R_DHDPS + R_DDPGALA + R_ANS + R_ACMANApts + R_CYSTL
                  + R_L_LACD3 + R_L_LACD2  = 0
 M_pyr_e:         - R_EX_pyr_e - R_PYRt2r  = 0
 M_q8_c:          - R_GLYCTO2 - R_HYD1 - R_ASPO3 - R_G3PD5 - R_SUCD4 - R_POX
                  + R_CYTBD - R_LDH_D2 + R_CYTBO3 + R_NO3R1 - R_FDH2 - R_MDH2
                  - R_DHORD2 - R_NADH5 - R_NADH6 - R_GLCDe - R_L_LACD2  = 0
 M_q8h2_c:        R_GLYCTO2 + R_HYD1 + R_DMQMT + R_ASPO3 + R_G3PD5 + R_SUCD4
                  + R_POX - R_CYTBD + R_LDH_D2 - R_CYTBO3 - R_NO3R1 + R_FDH2
                  + R_MDH2 + R_DHORD2 + R_NADH5 + R_NADH6 + R_GLCDe + R_L_LACD2
                   = 0
 M_quln_c:        R_QULNS - R_NNDPR  = 0
 M_r1p_c:         R_PUNP5 + R_PUNP1 + R_PUNP7 - R_PPM + R_PUNP3 + R_PYNP2r  = 0
 M_r5p_c:         - R_PRPPS + R_NMNN + R_AMPN + R_PPM + R_CMPN + R_NMNt7
                  - R_TKT1 + R_RBK - R_RPI  = 0
 M_rbl_L_c:       - R_RBK_L1 + R_ARAI  = 0
 M_rdmbzi_c:      - R_ADOCBLS + R_RZ5PP  = 0
 M_rhcys_c:       - R_RHCCE + R_AHCYSNS  = 0
 M_rib_D_c:       - R_RBK + R_RIBabc  = 0
 M_rib_D_e:       - R_RIBabc - R_EX_rib_D_e  = 0
 M_ribflv_c:      R_RBFSb - R_RBFK  = 0
 M_rml1p_c:       - R_RMPA + R_RMK  = 0
 M_rml_c:         - R_RMK + R_RMI  = 0
 M_rmn_c:         R_RMNt - R_RMI  = 0
 M_rmn_e:         - R_RMNt - R_EX_rmn_e  = 0
 M_ru5p_D_c:      R_GND - R_DB4PS - R_A5PISO - R_RPE + R_RPI  = 0
 M_ru5p_L_c:      - R_RBP4E + R_RBK_L1 + R_X5PL3E  = 0
 M_s7p_c:         - R_TALA - R_S7PI + R_TKT1  = 0
 M_sbt6p_c:       - R_SBTPD + R_SBTpts  = 0
 M_sbt_D_e:       - R_SBTpts - R_EX_sbt_D_e  = 0
 M_sbzcoa_c:      R_SUCBZL - R_NPHS  = 0
 M_scl_c:         R_SHCHD2 - R_SHCHF  = 0
 M_seln_c:        - R_SELNPS  = 0
 M_selnp_c:       R_SELNPS  = 0
 M_ser_D_c:       - R_SERD_D + R_DSERt2r  = 0
 M_ser_D_e:       - R_EX_ser_D_e - R_DSERt2r  = 0
 M_ser_L_c:       R_SERt2r - R_SERD_L - R_PSSA_EC - R_GHMT2 - R_TRPS2 - R_TRPS1
                  - 0.205 R_BiomassEcoli - R_SERASr + R_GPDDA3 + R_SERt4
                  + R_PSP_L - R_SERAT  = 0
 M_ser_L_e:       - R_SERt2r - R_EX_ser_L_e - R_SERt4  = 0
 M_seramp_c:      - 3 R_ENTCS + R_SERASr  = 0
 M_sheme_c:       R_SHCHF  = 0
 M_skm5p_c:       R_SHKK - R_PSCVT  = 0
 M_skm_c:         R_SHK3Dr - R_SHKK  = 0
 M_sl26da_c:      - R_SDPDS - R_SDPTA  = 0
 M_sl2a6o_c:      R_THDPS + R_SDPTA  = 0
 M_so3_c:         R_PAPSR + R_TAUDO + R_CYANST + R_SULR  = 0
 M_so4_c:         - R_SADT2 + R_SULabc  = 0
 M_so4_e:         R_EXinput_so4_e - R_EXoutput_so4_e - R_SULabc  = 0
 M_spmd_c:        R_SPMS - 0.007 R_BiomassEcoli + R_SPMDabc + R_GSPMDA
                  - R_GSPMDS - R_SPMDAT1 - R_SPMDAT2  = 0
 M_spmd_e:        - R_EX_spmd_e - R_SPMDabc  = 0
 M_ssaltpp_c:     - R_SHCHCS2 + R_OXGDC2  = 0
 M_suc6p_c:       R_SUCpts - R_FFSD  = 0
 M_sucarg_c:      R_AST - R_SADH  = 0
 M_sucbz_c:       R_SUCBZS - R_SUCBZL  = 0
 M_succ_c:        R_SUCCt2_2 + R_SDPDS + R_MCITL2 + R_ICL - R_SUCFUMt
                  + R_SUCCt2_3 + R_SSALy + R_SSALx + R_TAUDO + R_SUCCabc
                  - R_SUCCt2b - R_CITt7 - R_TARTRt7 + R_HKNDDH + R_SGDS
                  - R_PPCSCT - R_SUCD1i - R_SUCOAS + R_ASPO5 + R_FRD3 + R_SHSL1
                  + R_FRD2  = 0
 M_succ_e:        R_EXinput_succ_e - R_SUCCt2_2 - R_EXoutput_succ_e + R_SUCFUMt
                  - R_SUCCt2_3 - R_SUCCabc + R_SUCCt2b + R_CITt7 + R_TARTRt7
                   = 0
 M_succoa_c:      - R_MMM2 - R_HSST + R_AKGDH - R_AST - 0.000003 R_BiomassEcoli
                  - R_THDPS + R_PPCSCT + R_SUCOAS  = 0
 M_sucglu_c:      R_SGSAD - R_SGDS  = 0
 M_sucgsa_c:      - R_SGSAD + R_SOTA  = 0
 M_suchms_c:      R_HSST - R_SHSL1  = 0
 M_sucorn_c:      R_SADH - R_SOTA  = 0
 M_sucr_e:        - R_SUCpts - R_EX_sucr_e  = 0
 M_sucsal_c:      R_ABTA - R_SSALy - R_SSALx  = 0
 M_tag6p_D_c:     - R_PFK_2 + R_GLTPD  = 0
 M_tagdp_D_c:     - R_TGBPA + R_PFK_2  = 0
 M_tagur_c:       R_GUI2 + R_TAGURr  = 0
 M_tartr_L_c:     - R_TARTD + R_TARTRt7  = 0
 M_tartr_L_e:     - R_EX_tartr_L_e - R_TARTRt7  = 0
 M_taur_c:        R_TAURabc - R_TAUDO  = 0
 M_taur_e:        - R_TAURabc - R_EX_taur_e  = 0
 M_tcynt_c:       R_CYANST  = 0
 M_tdeACP_c:      R_AACPS2 + R_C141SN - 0.1 R_PASYN_EC  = 0
 M_thdp_c:        R_DHDPRy - R_THDPS  = 0
 M_thf_c:         - R_GLYCL - R_GHMT2 + R_MOHMT + R_AICART + R_METS + R_DHFR
                  + R_FTHFD + R_GARFT  = 0
 M_thm_c:         - R_TMKr + R_THMabc  = 0
 M_thm_e:         - R_EXoutput_thm_e + R_EXinput_thm_e - R_THMabc  = 0
 M_thmmp_c:       R_TMKr - R_TMPKr + R_TMPPP  = 0
 M_thmpp_c:       R_SHCHCS2 + R_TMPKr - R_OXGDC2  = 0
 M_thr_L_c:       - 0.241 R_BiomassEcoli + R_THRS - R_THRD + R_THRt4 - R_THRD_L
                  + R_THRt2r + R_THRabc - R_THRAr  = 0
 M_thr_L_e:       - R_EX_thr_L_e - R_THRt4 - R_THRt2r - R_THRabc  = 0
 M_thym_c:        R_TMDPP  = 0
 M_thymd_c:       - R_TMDK1 + R_NTD5 + R_THMDt2r - R_TMDPP  = 0
 M_thymd_e:       - R_THMDt2r - R_EX_thymd_e  = 0
 M_tma_c:         R_TMAOR2 + R_TMAOR1  = 0
 M_tma_e:         R_TMAOR1e + R_TMAOR2e - R_EX_tma_e  = 0
 M_tmao_c:        - R_TMAOR2 - R_TMAOR1  = 0
 M_tmao_e:        - R_TMAOR1e - R_TMAOR2e - R_EX_tmao_e  = 0
 M_trdox_c:       R_PAPSR + R_RNDR4 + R_RNTR2 + R_RNTR3 - R_TRDR + R_RNTR1
                  + R_RNDR1 + R_RNTR4 + R_RNDR2 + R_RNDR3  = 0
 M_trdrd_c:       - R_PAPSR - R_RNDR4 - R_RNTR2 - R_RNTR3 + R_TRDR - R_RNTR1
                  - R_RNDR1 - R_RNTR4 - R_RNDR2 - R_RNDR3  = 0
 M_tre6p_c:       R_TREpts + R_TRE6PS - R_TRE6PP - R_TRE6PH  = 0
 M_tre_c:         - R_TREH + R_TRE6PP  = 0
 M_tre_e:         - R_TREHe - R_TREpts + R_EXinput_tre_e - R_EXoutput_tre_e
                   = 0
 M_trnaglu_c:     - R_GLUTRS + R_GLUTRR  = 0
 M_trp_L_c:       - R_TRPAS2 + R_TRPt2r + R_TRPS2 + R_TRPS1
                  - 0.054 R_BiomassEcoli  = 0
 M_trp_L_e:       - R_TRPt2r - R_EX_trp_L_e  = 0
 M_tsul_c:        R_TSULabc - R_CYANST  = 0
 M_tsul_e:        - R_TSULabc - R_EX_tsul_e  = 0
 M_ttdca_c:       0.02 R_PLIPA3 - R_AACPS1 + 0.02 R_PLIPA1 + 0.02 R_PLIPA2
                  + 0.02 R_LPLIPA2 + 0.02 R_LPLIPA3 + 0.02 R_LPLIPA1
                  + R_TTDCAt2 - R_FAO1  = 0
 M_ttdca_e:       - R_EX_ttdca_e - R_TTDCAt2  = 0
 M_ttdcea_c:      0.05 R_PLIPA3 - R_AACPS2 + 0.05 R_PLIPA1 + 0.05 R_PLIPA2
                  + 0.05 R_LPLIPA2 + 0.05 R_LPLIPA3 + 0.05 R_LPLIPA1  = 0
 M_tyr_L_c:       - R_TYRTA - 0.131 R_BiomassEcoli - R_THZPSN + R_TYRt2r  = 0
 M_tyr_L_e:       - R_EX_tyr_L_e - R_TYRt2r  = 0
 M_u23ga_c:       R_U23GAAT - R_USHD - R_LPADSS  = 0
 M_u3aga_c:       - R_UHGADA + R_UAGAAT  = 0
 M_u3hga_c:       - R_U23GAAT + R_UHGADA  = 0
 M_uaagmda_c:     R_UAGPT3 - R_PPTGS  = 0
 M_uaccg_c:       R_UAGCVT - R_UAPGR  = 0
 M_uacgam_c:      - R_UAGPT3 - R_ACGAMT - R_UAGCVT - R_UAG2Ei - R_UAGAAT
                  + R_UAGDP  = 0
 M_uacmam_c:      - R_UACMAMO + R_UAG2Ei  = 0
 M_uacmamu_c:     - R_ACMAMUT + R_UACMAMO  = 0
 M_uagmda_c:      - R_UAGPT3 + R_PAPPT3  = 0
 M_uama_c:        R_UAMAS - R_UAMAGS  = 0
 M_uamag_c:       R_UAMAGS - R_UAAGDS  = 0
 M_uamr_c:        - R_UAMAS + R_UAPGR  = 0
 M_udcpdp_c:      - R_UDCPDP + R_PPTGS + R_UDCPDPS + R_ECAP_EC  = 0
 M_udcpp_c:       - R_ACGAMT + R_UDCPDP - R_PAPPT3  = 0
 M_udp_c:         R_UAGPT3 + R_ACMAMUT - R_RNDR4 + R_UMPK - R_NDPK2
                  + 2 R_LPSSYN_EC + R_TRE6PS + R_LPADSS  = 0
 M_udpg_c:        - R_UGLT - R_UDPG4E - 0.003 R_BiomassEcoli - R_UDPGD
                  - 2 R_LPSSYN_EC - R_TRE6PS + R_GALU  = 0
 M_udpgal_c:      R_UGLT + R_UDPG4E - R_UDPGALM  = 0
 M_udpgalfur_c:   R_UDPGALM  = 0
 M_udpglcur_c:    R_UDPGD  = 0
 M_ugmd_c:        R_UAAGDS - R_UGMDDS  = 0
 M_ugmda_c:       - R_PAPPT3 + R_UGMDDS  = 0
 M_ump_c:         R_ACGAMT - R_NTD2 - R_UMPK + R_USHD + R_UPPRT + R_NTPP8
                  + R_PAPPT3 + R_OMPDC + R_URIK2  = 0
 M_unaga_c:       - R_ACMAMUT + R_ACGAMT  = 0
 M_unagamu_c:     R_ACMAMUT - R_AADDGT  = 0
 M_unagamuf_c:    R_AADDGT - R_ECAP_EC  = 0
 M_uppg3_c:       - R_UPP3MT + R_UPP3S - R_UPPDC1  = 0
 M_ura_c:         R_DURIPP + R_URAt2r - R_UPPRT + R_CSND + R_PYNP2r  = 0
 M_ura_e:         R_EXinput_ura_e - R_URAt2r - R_EXoutput_ura_e  = 0
 M_urdglyc_c:     R_ALLTAH - R_UGLYCH  = 0
 M_urea_c:        R_ALLTAH + R_UREAt + R_AGMT  = 0
 M_urea_e:        - R_UREAt - R_EX_urea_e  = 0
 M_uri_c:         R_NTD2 + R_URIt2r + R_CYTD - R_URIK2 - R_PYNP2r  = 0
 M_uri_e:         R_EXinput_uri_e - R_URIt2r - R_EXoutput_uri_e  = 0
 M_utp_c:         - R_CTPS2 + R_NDPK2 - 0.136 R_BiomassEcoli - R_NTPP8 - R_GALU
                  - R_RNTR4 - R_UAGDP  = 0
 M_val_L_c:       R_VALabc - R_VALTA - 0.402 R_BiomassEcoli + R_VPAMT
                  + R_VALt2r  = 0
 M_val_L_e:       - R_VALabc - R_EXoutput_val_L_e - R_VALt2r
                  + R_EXinput_val_L_e  = 0
 M_xan_c:         R_PUNP7 + R_XANt2 - R_XPPT + R_XANt + R_GUAD  = 0
 M_xan_e:         - R_EXoutput_xan_e + R_EXinput_xan_e - R_XANt2 - R_XANt  = 0
 M_xmp_c:         R_IMPD - R_NTD10 + R_XPPT - R_GMPS2  = 0
 M_xtsn_c:        - R_PUNP7 + R_XTSNt2r + R_NTD10  = 0
 M_xtsn_e:        - R_EX_xtsn_e - R_XTSNt2r  = 0
 M_xu5p_D_c:      R_RBP4E + R_XYLK - R_TKT1 - R_TKT2 + R_RPE  = 0
 M_xu5p_L_c:      R_KG6PDC - R_X5PL3E  = 0
 M_xyl_D_c:       R_XYLabc + R_XYLt2 - R_XYLI1  = 0
 M_xyl_D_e:       - R_EX_xyl_D_e - R_XYLabc - R_XYLt2  = 0
 M_xylu_D_c:      - R_XYLK + R_XYLI1  = 0
 ind_R_EXinput_pi_e:      bin_R_EXinput_pi_e = 0 -> R_EXinput_pi_e  = 0
 ind_R_EXinput_nh4_e:     bin_R_EXinput_nh4_e = 0 -> R_EXinput_nh4_e  = 0
 ind_R_EXinput_uri_e:     bin_R_EXinput_uri_e = 0 -> R_EXinput_uri_e  = 0
 ind_R_EXinput_asp_L_e:   bin_R_EXinput_asp_L_e = 0 -> R_EXinput_asp_L_e  = 0
 ind_R_EXinput_glyclt_e:  bin_R_EXinput_glyclt_e = 0 -> R_EXinput_glyclt_e  = 0
 ind_R_EXinput_tre_e:     bin_R_EXinput_tre_e = 0 -> R_EXinput_tre_e  = 0
 ind_R_EXinput_glcn_e:    bin_R_EXinput_glcn_e = 0 -> R_EXinput_glcn_e  = 0
 ind_R_EXinput_k_e:       bin_R_EXinput_k_e = 0 -> R_EXinput_k_e  = 0
 ind_R_EXinput_glc_e:     bin_R_EXinput_glc_e = 0 -> R_EXinput_glc_e  = 0
 ind_R_EXinput_dgsn_e:    bin_R_EXinput_dgsn_e = 0 -> R_EXinput_dgsn_e  = 0
 ind_R_EXinput_gly_e:     bin_R_EXinput_gly_e = 0 -> R_EXinput_gly_e  = 0
 ind_R_EXinput_12ppd_S_e: bin_R_EXinput_12ppd_S_e = 0 -> R_EXinput_12ppd_S_e
                           = 0
 ind_R_EXinput_thm_e:     bin_R_EXinput_thm_e = 0 -> R_EXinput_thm_e  = 0
 ind_R_EXinput_so4_e:     bin_R_EXinput_so4_e = 0 -> R_EXinput_so4_e  = 0
 ind_R_EXinput_gam_e:     bin_R_EXinput_gam_e = 0 -> R_EXinput_gam_e  = 0
 ind_R_EXinput_galct_D_e: bin_R_EXinput_galct_D_e = 0 -> R_EXinput_galct_D_e
                           = 0
 ind_R_EXinput_3hcinnm_e: bin_R_EXinput_3hcinnm_e = 0 -> R_EXinput_3hcinnm_e
                           = 0
 ind_R_EXinput_malttr_e:  bin_R_EXinput_malttr_e = 0 -> R_EXinput_malttr_e  = 0
 ind_R_EXinput_man6p_e:   bin_R_EXinput_man6p_e = 0 -> R_EXinput_man6p_e  = 0
 ind_R_EXinput_ile_L_e:   bin_R_EXinput_ile_L_e = 0 -> R_EXinput_ile_L_e  = 0
 ind_R_EXinput_lys_L_e:   bin_R_EXinput_lys_L_e = 0 -> R_EXinput_lys_L_e  = 0
 ind_R_EXinput_3hpppn_e:  bin_R_EXinput_3hpppn_e = 0 -> R_EXinput_3hpppn_e  = 0
 ind_R_EXinput_ura_e:     bin_R_EXinput_ura_e = 0 -> R_EXinput_ura_e  = 0
 ind_R_EXinput_xan_e:     bin_R_EXinput_xan_e = 0 -> R_EXinput_xan_e  = 0
 ind_R_EXinput_fe2_e:     bin_R_EXinput_fe2_e = 0 -> R_EXinput_fe2_e  = 0
 ind_R_EXinput_succ_e:    bin_R_EXinput_succ_e = 0 -> R_EXinput_succ_e  = 0
 ind_R_EXinput_acnam_e:   bin_R_EXinput_acnam_e = 0 -> R_EXinput_acnam_e  = 0
 ind_R_EXinput_val_L_e:   bin_R_EXinput_val_L_e = 0 -> R_EXinput_val_L_e  = 0
 ind_R_EXinput_ala_L_e:   bin_R_EXinput_ala_L_e = 0 -> R_EXinput_ala_L_e  = 0
Bounds
-999999 <= R_GLYCt <= 999999
 0 <= R_PLIPA3 <= 999999
 0 <= R_TREH <= 999999
 0 <= R_CTPS2 <= 999999
 0 <= R_ACNML <= 999999
 0 <= R_LYSDC <= 999999
 0 <= R_ALDD2x <= 999999
-999999 <= R_MME <= 999999
 0 <= R_GLCpts <= 999999
 0 <= R_NACODA <= 999999
-999999 <= R_UGLT <= 999999
 0 <= R_FRUpts2 <= 999999
-999999 <= R_PRPPS <= 999999
 0 <= R_TMDS <= 999999
 0 <= R_OMMBLHX <= 999999
-999999 <= R_DASYN_EC <= 999999
 0 <= R_IPMD <= 999999
 0 <= R_ADCS <= 999999
 0 <= R_ASNN <= 999999
-999999 <= R_PROt2r <= 999999
-999999 <= R_IDOND <= 999999
 0 <= R_PYDXNK <= 999999
 0 <= R_PANTS <= 999999
 0 <= R_MMM2 <= 999999
 0 <= R_EXoutput_xan_e <= 1000
-999999 <= R_AOXSr <= 999999
 0 <= R_GLYCL <= 999999
 0 <= R_GLYCK <= 999999
 0 <= R_ADCL <= 999999
 0 <= R_GLCP <= 999999
 0 <= R_MOAT <= 999999
 0 <= R_EXinput_k_e <= 100
 0 <= R_MAN6Pt6_2 <= 999999
 0 <= R_RMNt <= 999999
 0 <= R_EXinput_succ_e <= 100
 0 <= R_SPMS <= 999999
-999999 <= R_GLYBt2r <= 999999
 0 <= R_DHPS2 <= 999999
 0 <= R_HBZOPT <= 999999
-999999 <= R_SERt2r <= 999999
-999999 <= R_CYTK2 <= 999999
 0 <= R_PYDXPP <= 999999
 0 <= R_EX_acald_e <= 999999
 0 <= R_DHNAOT <= 999999
-999999 <= R_LCAD <= 999999
 0 <= R_INSK <= 999999
 0 <= R_MNLpts <= 999999
 0 <= R_FCLT <= 999999
 0 <= R_SUCpts <= 999999
 0 <= R_ORNabc <= 999999
-999999 <= R_RMPA <= 999999
 0 <= R_CAT <= 999999
 0 <= R_TMDK1 <= 999999
 0 <= R_EXinput_dgsn_e <= 100
-999999 <= R_Kt2r <= 999999
 0 <= R_FMNAT <= 999999
 0 <= R_EXinput_xan_e <= 100
 0 <= R_EXoutput_lys_L_e <= 1000
 0 <= R_EDD <= 999999
 0 <= R_FCLK <= 999999
 0 <= R_EDA <= 999999
 0 <= R_TDSK <= 999999
 0 <= R_EX_pnto_R_e <= 999999
 0 <= R_DHAPT <= 999999
 0 <= R_SUCCt2_2 <= 999999
 0 <= R_OMPHHX <= 999999
 0 <= R_UAGPT3 <= 999999
-999999 <= R_TMKr <= 999999
 0 <= R_SERD_D <= 999999
-999999 <= R_DGK1 <= 999999
 0 <= R_CHRPL <= 999999
-999999 <= R_SHK3Dr <= 999999
 0 <= R_SERD_L <= 999999
 0 <= R_U23GAAT <= 999999
-999999 <= R_PRAGSr <= 999999
-999999 <= R_ALLTNt2r <= 999999
-999999 <= R_PPM2 <= 999999
 0 <= R_AMMQT8_2 <= 999999
 0 <= R_SPODM <= 999999
 0 <= R_ARBabc <= 999999
 0 <= R_DHCIND <= 999999
 0 <= R_SDPDS <= 999999
      R_ATPM = 7.6
 0 <= R_ACMAMUT <= 999999
 0 <= R_ACGAMT <= 999999
-999999 <= R_SBTPD <= 999999
 0 <= R_GSNt2 <= 999999
 0 <= R_MALTpts <= 999999
 0 <= R_PAPSR <= 999999
 0 <= R_ALLTAH <= 999999
-999999 <= R_CLPNS_EC <= 999999
 0 <= R_PPND <= 999999
      R_EXoutput_glc_e = 0
 0 <= R_PDH <= 999999
 0 <= R_BUTCT <= 999999
 0 <= R_AACPS5 <= 999999
 0 <= R_AACPS4 <= 999999
 0 <= R_GUAPRT <= 999999
 0 <= R_AACPS2 <= 999999
 0 <= R_AACPS1 <= 999999
 0 <= R_NADDP <= 999999
 0 <= R_NTD9 <= 999999
 0 <= R_NTD8 <= 999999
 0 <= R_ACODA <= 999999
-999999 <= R_GUAt <= 999999
 0 <= R_NTD1 <= 999999
 0 <= R_NTD3 <= 999999
 0 <= R_NTD2 <= 999999
 0 <= R_NTD5 <= 999999
 0 <= R_GLUDC <= 999999
 0 <= R_GLYCLTDx <= 999999
 0 <= R_GLYCLTDy <= 999999
-999999 <= R_ALAR <= 999999
 0 <= R_EXoutput_succ_e <= 1000
-999999 <= R_PSSA_EC <= 999999
-999999 <= R_NAt3_1 <= 999999
 0 <= R_NAt3_2 <= 999999
-999999 <= R_GLUDy <= 999999
-999999 <= R_DHAt <= 999999
 0 <= R_ORNDC <= 999999
 0 <= R_AADDGT <= 999999
 0 <= R_OHPHM <= 999999
 0 <= R_RNDR4 <= 999999
 0 <= R_EXoutput_gam_e <= 1000
 0 <= R_LPLIPA6 <= 999999
 0 <= R_PPCDC <= 999999
 0 <= R_PRAMPC <= 999999
-999999 <= R_LYSt2r <= 999999
-999999 <= R_ACOATA <= 999999
 0 <= R_G5SADs <= 999999
 0 <= R_FFSD <= 999999
 0 <= R_DURIt2 <= 999999
-999999 <= R_GLUt2r <= 999999
-999999 <= R_UMPK <= 999999
 0 <= R_VALabc <= 999999
 0 <= R_ASP1DC <= 999999
 0 <= R_CBMK <= 999999
 0 <= R_GALabc <= 999999
 0 <= R_ACHBS <= 999999
-999999 <= R_UDPG4E <= 999999
 0 <= R_DHORD5 <= 999999
 0 <= R_DMPPS <= 999999
 0 <= R_RHCCE <= 999999
 0 <= R_CYSabc <= 999999
 0 <= R_SHKK <= 999999
 0 <= R_GLUTRS <= 999999
 0 <= R_GLUTRR <= 999999
 0 <= R_CBL1abc <= 999999
 0 <= R_EXoutput_dgsn_e <= 1000
 0 <= R_G5SD <= 999999
 0 <= R_DHPPDA2 <= 999999
 0 <= R_ACPS1 <= 999999
 0 <= R_MALTHXabc <= 999999
 0 <= R_EXinput_so4_e <= 100
 0 <= R_HSTPT <= 999999
-999999 <= R_O2t <= 999999
-999999 <= R_12PPDt <= 999999
 0 <= R_DNMPPA <= 999999
-999999 <= R_ADSL2r <= 999999
 0 <= R_NNAM <= 999999
 0 <= R_GLU5K <= 999999
 0 <= R_USHD <= 999999
 0 <= R_SHCHCS2 <= 999999
-999999 <= R_RBP4E <= 999999
 0 <= R_PSD_EC <= 999999
 0 <= R_C141SN <= 999999
 0 <= R_PLIPA1 <= 999999
-999999 <= R_CYTDt2r <= 999999
 0 <= R_EXoutput_val_L_e <= 1000
 0 <= R_PLIPA2 <= 999999
 0 <= R_ICHORT <= 999999
 0 <= R_DKMPPD <= 999999
 0 <= R_DPCOAK <= 999999
 0 <= R_UACMAMO <= 999999
 0 <= R_FRUK <= 999999
 0 <= R_PNTK <= 999999
-999999 <= R_FCLPA <= 999999
 0 <= R_ADNK1 <= 999999
-999999 <= R_MCITL2 <= 999999
 0 <= R_PROt4 <= 999999
 0 <= R_UAGCVT <= 999999
 0 <= R_EXinput_fe2_e <= 100
 0 <= R_DMSOR1 <= 999999
 0 <= R_DMSOR2 <= 999999
 0 <= R_IDOND2 <= 999999
 0 <= R_PGPP_EC <= 999999
 0 <= R_RBK_L1 <= 999999
 0 <= R_2DGULRx <= 999999
 0 <= R_2DGULRy <= 999999
-999999 <= R_ACCOACr <= 999999
 0 <= R_DHDPRy <= 999999
 0 <= R_3HPPPNH <= 999999
 0 <= R_GF6PTA <= 999999
 0 <= R_NO3t7 <= 999999
 0 <= R_GHMT2 <= 999999
 0 <= R_MANpts <= 999999
 0 <= R_UDCPDP <= 999999
 0 <= R_EXoutput_nh4_e <= 1000
 0 <= R_EXoutput_ala_L_e <= 1000
 0 <= R_HSK <= 999999
-999999 <= R_PHEt2r <= 999999
 0 <= R_FAO4 <= 999999
 0 <= R_EXoutput_gly_e <= 1000
 0 <= R_GTPCI <= 999999
-999999 <= R_TPI <= 999999
 0 <= R_ACCOAL <= 999999
 0 <= R_GRTT <= 999999
 0 <= R_DAPabc <= 999999
-999999 <= R_ASPTA <= 999999
-999999 <= R_NDPK8 <= 999999
-999999 <= R_NDPK7 <= 999999
-999999 <= R_NDPK6 <= 999999
-999999 <= R_NDPK5 <= 999999
-999999 <= R_NDPK4 <= 999999
-999999 <= R_NDPK3 <= 999999
-999999 <= R_NDPK2 <= 999999
-999999 <= R_NDPK1 <= 999999
 0 <= R_DUTPDP <= 999999
 0 <= R_DURIK1 <= 999999
 0 <= R_HSST <= 999999
 0 <= R_ASPt2 <= 999999
 0 <= R_EX_tyr_L_e <= 999999
-999999 <= R_M1PD <= 999999
 0 <= R_GLNabc <= 999999
 0 <= R_AKGDH <= 999999
 0 <= R_GALS3 <= 999999
 0 <= R_LEUTAi <= 999999
 0 <= R_TAURabc <= 999999
-999999 <= R_PDX5PO <= 999999
-999999 <= R_MLTP3 <= 999999
 0 <= R_LACZ <= 999999
 0 <= R_G6PDA <= 999999
-999999 <= R_ENO <= 999999
-999999 <= R_THMDt2r <= 999999
 0 <= R_PFK <= 999999
 0 <= R_UHGADA <= 999999
-999999 <= R_TYRTA <= 999999
 0 <= R_GNK <= 999999
-999999 <= R_DDGLCNt2r <= 999999
-999999 <= R_GLUR <= 999999
 0 <= R_KARA1i <= 999999
 0 <= R_SHCHD2 <= 999999
-999999 <= R_AMAOTr <= 999999
 0 <= R_GMHEPAT <= 999999
 0 <= R_RNTR2 <= 999999
 0 <= R_RNTR3 <= 999999
 0 <= R_GND <= 999999
 0 <= R_NMNAT <= 999999
 0 <= R_EXoutput_asp_L_e <= 1000
-999999 <= R_LCARS <= 999999
-999999 <= R_AKGt2r <= 999999
-999999 <= R_MANAO <= 999999
 0 <= R_EXinput_ura_e <= 100
 0 <= R_GSNK <= 999999
 0 <= R_NMNP <= 999999
 0 <= R_GLUN <= 999999
-999999 <= R_ADHEr <= 999999
 0 <= R_C120SN <= 999999
 0 <= R_DNTPPA <= 999999
 0 <= R_DXYLK <= 999999
 0 <= R_GLYATi <= 999999
 0 <= R_DKMPPD2 <= 999999
-999999 <= R_ADSL1r <= 999999
 0 <= R_NAMNPP <= 999999
-999999 <= R_PTAr <= 999999
-999999 <= R_OCBT <= 999999
 0 <= R_EX_thr_L_e <= 999999
 0 <= R_GLYCTO2 <= 999999
 0 <= R_GLYCTO3 <= 999999
 0 <= R_ADOCBLS <= 999999
 0 <= R_EXinput_gly_e <= 100
-999999 <= R_TRPAS2 <= 999999
 0 <= R_CINNDO <= 999999
 0 <= R_EX_glyald_e <= 999999
 0 <= R_DAAD <= 999999
 0 <= R_ACBIPGT <= 999999
-999999 <= R_TRPt2r <= 999999
 0 <= R_GMPR <= 999999
 0 <= R_DGSNt2 <= 999999
 0 <= R_IPDPS <= 999999
-999999 <= R_ILEt2r <= 999999
 0 <= R_SBTpts <= 999999
 0 <= R_EX_gsn_e <= 999999
 0 <= R_ICL <= 999999
 0 <= R_ABTA <= 999999
 0 <= R_HYD2 <= 999999
 0 <= R_HYD3 <= 999999
 0 <= R_HYD1 <= 999999
 0 <= R_PPNCL2 <= 999999
 0 <= R_NTPTP2 <= 999999
 0 <= R_OMBZLM <= 999999
 0 <= R_PPTGS <= 999999
-999999 <= R_TMPKr <= 999999
 0 <= R_ENTCS <= 999999
 0 <= R_TSULabc <= 999999
 0 <= R_UAG2Ei <= 999999
 0 <= R_GLGC <= 999999
 0 <= R_ACACCT <= 999999
 0 <= R_TREHe <= 999999
 0 <= R_PTA2 <= 999999
 0 <= R_DPR <= 999999
 0 <= R_OMCDC <= 999999
 0 <= R_EDTXS4 <= 999999
 0 <= R_EDTXS3 <= 999999
 0 <= R_EDTXS2 <= 999999
 0 <= R_EDTXS1 <= 999999
-999999 <= R_SUCFUMt <= 999999
-999999 <= R_H2Ot <= 999999
 0 <= R_SUCCt2_3 <= 999999
 0 <= R_SSALy <= 999999
 0 <= R_SSALx <= 999999
-999999 <= R_FORt <= 999999
 0 <= R_ASPO4 <= 999999
 0 <= R_ASPO6 <= 999999
 0 <= R_EX_pppn_e <= 999999
-999999 <= R_EX_na1_e <= 999999
-999999 <= R_GLYALDt <= 999999
 0 <= R_DHCINDO <= 999999
-999999 <= R_VALTA <= 999999
 0 <= R_MLTG1 <= 999999
 0 <= R_BPNT <= 999999
 0 <= R_3HCINNMH <= 999999
 0 <= R_GCALDD <= 999999
 0 <= R_EX_spmd_e <= 999999
 0 <= R_EXinput_glcn_e <= 100
-999999 <= R_PHETA1 <= 999999
 0 <= R_DMQMT <= 999999
 0 <= R_EX_nmn_e <= 999999
 0 <= R_IGPS <= 999999
-999999 <= R_UREAt <= 999999
 0 <= R_HETZK <= 999999
 0 <= R_EX_ttdca_e <= 999999
-999999 <= R_MDH <= 999999
 0 <= R_EX_dcyt_e <= 999999
 0 <= R_DOGULNR <= 999999
-999999 <= R_GUI1 <= 999999
-999999 <= R_GUI2 <= 999999
 0 <= R_BSORy <= 999999
 0 <= R_XYLI2i <= 999999
-999999 <= R_LCTSt <= 999999
 0 <= R_ARGDC <= 999999
 0 <= R_GLUCYS <= 999999
 0 <= R_AGDC <= 999999
-999999 <= R_NO2t2r <= 999999
-999999 <= R_ILETA <= 999999
-999999 <= R_PGM <= 999999
 0 <= R_PGL <= 999999
-999999 <= R_PGK <= 999999
-999999 <= R_PGI <= 999999
 0 <= R_PTRCabc <= 999999
 0 <= R_UGLYCH <= 999999
-999999 <= R_BUTt2r <= 999999
 0 <= R_EX_phe_L_e <= 999999
 0 <= R_AST <= 999999
 0 <= R_TRPS3 <= 999999
 0 <= R_TRPS2 <= 999999
 0 <= R_TRPS1 <= 999999
 0 <= R_GLCRAL <= 999999
 0 <= R_EXoutput_12ppd_S_e <= 1000
 0 <= R_Kabc <= 999999
-999999 <= R_F6PA <= 999999
 0 <= R_ALATA_D2 <= 999999
 0 <= R_GUAt2 <= 999999
-999999 <= R_DURIPP <= 999999
 0 <= R_EX_chol_e <= 999999
 0 <= R_TAUDO <= 999999
 0 <= R_EXoutput_thm_e <= 1000
 0 <= R_GALt2 <= 999999
 1 <= R_BiomassEcoli <= 2
 0 <= R_EX_no2_e <= 999999
-999999 <= R_AGPR <= 999999
 0 <= R_ABUTt2 <= 999999
 0 <= R_PRMICIi <= 999999
-999999 <= R_ARGSL <= 999999
-999999 <= R_CHLt2r <= 999999
 0 <= R_KARA2i <= 999999
 0 <= R_TRDR <= 999999
-999999 <= R_HYXNt <= 999999
 0 <= R_EXoutput_so4_e <= 1000
 0 <= R_FE2abc <= 999999
-999999 <= R_MLTP1 <= 999999
 0 <= R_GDMANE <= 999999
 0 <= R_MCITS <= 999999
-999999 <= R_MLTP2 <= 999999
 0 <= R_NADTRHD <= 999999
 0 <= R_CSNt2 <= 999999
 0 <= R_EXinput_12ppd_S_e <= 100
 0 <= R_EX_glyc_e <= 999999
 0 <= R_UDPGD <= 999999
 0 <= R_EX_alltn_e <= 999999
 0 <= R_MCITD <= 999999
 0 <= R_DMSOR1e <= 999999
 0 <= R_2DGLCNRy <= 999999
 0 <= R_2DGLCNRx <= 999999
 0 <= R_GLUabc <= 999999
 0 <= R_EX_nac_e <= 999999
 0 <= R_NTPP1 <= 999999
 0 <= R_KDOPS <= 999999
 0 <= R_ASPO3 <= 999999
-999999 <= R_MAN6PI <= 999999
 0 <= R_SUCBZS <= 999999
 0 <= R_GMHEPPA <= 999999
 0 <= R_THRS <= 999999
 0 <= R_EX_ser_D_e <= 999999
 0 <= R_LGTHL <= 999999
 0 <= R_SUCBZL <= 999999
 0 <= R_EXinput_glyclt_e <= 100
 0 <= R_NADPPPS <= 999999
-999999 <= R_LDH_D <= 999999
-999999 <= R_CBLAT_DELETE <= 999999
 0 <= R_MALTabc <= 999999
 0 <= R_EXinput_thm_e <= 100
 0 <= R_APRAUR <= 999999
 0 <= R_METDabc <= 999999
 0 <= R_EX_crn_e <= 999999
 0 <= R_NMNN <= 999999
-999999 <= R_DALAt2r <= 999999
 0 <= R_FRUpts <= 999999
 0 <= R_PROabc <= 999999
 0 <= R_EXinput_malttr_e <= 100
 0 <= R_EX_pyr_e <= 999999
 0 <= R_NADS1 <= 999999
 0 <= R_EX_asn_L_e <= 999999
 0 <= R_AGMT <= 999999
-999999 <= R_SERASr <= 999999
 0 <= R_PSERT <= 999999
 0 <= R_DB4PS <= 999999
 0 <= R_SUCCabc <= 999999
 0 <= R_EX_met_D_e <= 999999
 0 <= R_MALS <= 999999
-999999 <= R_TALA <= 999999
-999999 <= R_ASAD <= 999999
-999999 <= R_GLUABUTt7 <= 999999
 0 <= R_TMAOR2 <= 999999
 0 <= R_EX_tartr_L_e <= 999999
-999999 <= R_IMPC <= 999999
-999999 <= R_URAt2r <= 999999
 0 <= R_EX_leu_L_e <= 999999
 0 <= R_EX_for_e <= 999999
 0 <= R_UPPRT <= 999999
 0 <= R_CRNt7 <= 999999
-999999 <= R_ADEt2r <= 999999
-999999 <= R_ACONT <= 999999
 0 <= R_AMPN <= 999999
 0 <= R_DDGLK <= 999999
 0 <= R_PPNDH <= 999999
 0 <= R_UAMAS <= 999999
-999999 <= R_PUNP5 <= 999999
-999999 <= R_PUNP4 <= 999999
 0 <= R_GP4GH <= 999999
-999999 <= R_PUNP6 <= 999999
-999999 <= R_PUNP1 <= 999999
 0 <= R_AMPMS <= 999999
-999999 <= R_PUNP2 <= 999999
 0 <= R_NNDMBRT <= 999999
 0 <= R_EX_thymd_e <= 999999
 0 <= R_IPDDIi <= 999999
 0 <= R_OP4ENH <= 999999
 0 <= R_PAPA_EC <= 999999
-999999 <= R_DHBSr <= 999999
 0 <= R_EXoutput_fe2_e <= 1000
 0 <= R_ACALDi <= 999999
 0 <= R_EX_cys_L_e <= 999999
 0 <= R_SADT2 <= 999999
-999999 <= R_DSERt2r <= 999999
 0 <= R_CDPMEK <= 999999
 0 <= R_UAMAGS <= 999999
-999999 <= R_GLCRt2r <= 999999
 0 <= R_PRAIi <= 999999
 0 <= R_IMPD <= 999999
 0 <= R_LPSSYN_EC <= 999999
 0 <= R_ALATA_L2 <= 999999
 0 <= R_MNNH <= 999999
 0 <= R_EX_cit_e <= 999999
-999999 <= R_TMDPP <= 999999
 0 <= R_MTHFR2 <= 999999
 0 <= R_ANPRT <= 999999
 0 <= R_EX_dmso_e <= 999999
 0 <= R_EX_ser_L_e <= 999999
 0 <= R_MECDPS <= 999999
 0 <= R_RNTR1 <= 999999
-999999 <= R_PUNP7 <= 999999
 0 <= R_EXinput_3hcinnm_e <= 100
 0 <= R_ASNabc <= 999999
 0 <= R_EXoutput_man6p_e <= 1000
 0 <= R_TARTD <= 999999
 0 <= R_GPDDA5 <= 999999
 0 <= R_GPDDA4 <= 999999
 0 <= R_GPDDA3 <= 999999
 0 <= R_GPDDA2 <= 999999
 0 <= R_GPDDA1 <= 999999
 0 <= R_GLCt2 <= 999999
 0 <= R_AMANK <= 999999
-999999 <= R_GAPD <= 999999
 0 <= R_EXoutput_ura_e <= 1000
 0 <= R_METAT <= 999999
 0 <= R_NTPP8 <= 999999
 0 <= R_ALLTN <= 999999
 0 <= R_NTPP6 <= 999999
 0 <= R_NTPP7 <= 999999
 0 <= R_NTPP4 <= 999999
 0 <= R_NTPP5 <= 999999
 0 <= R_NTPP2 <= 999999
 0 <= R_NTPP3 <= 999999
 0 <= R_RNDR1 <= 999999
 0 <= R_DHNPA2 <= 999999
 0 <= R_C161SN <= 999999
 0 <= R_PYK <= 999999
 0 <= R_THDPS <= 999999
 0 <= R_MALTTRabc <= 999999
 0 <= R_EX_xyl_D_e <= 999999
 0 <= R_CPPPGO <= 999999
 0 <= R_GLYOX <= 999999
 0 <= R_ACGK <= 999999
 0 <= R_EX_rmn_e <= 999999
-999999 <= R_DHORTS <= 999999
 0 <= R_PYDXK <= 999999
-999999 <= R_HCO3E <= 999999
 0 <= R_ACGS <= 999999
 0 <= R_TDPADGAT <= 999999
 0 <= R_CHORM <= 999999
 0 <= R_ATPPRT <= 999999
 0 <= R_RBFSb <= 999999
 0 <= R_RBFSa <= 999999
 0 <= R_EX_dad_2_e <= 999999
-999999 <= R_PGSA_EC <= 999999
 0 <= R_EXinput_uri_e <= 100
 0 <= R_GOFUCR <= 999999
 0 <= R_PPS <= 999999
 0 <= R_NTPTP1 <= 999999
-999999 <= R_FUM <= 999999
 0 <= R_TREpts <= 999999
 0 <= R_4HTHRS <= 999999
 0 <= R_ACONMT <= 999999
 0 <= R_HEX1 <= 999999
 0 <= R_PPC <= 999999
 0 <= R_PPA <= 999999
 0 <= R_EX_fuc_L_e <= 999999
 0 <= R_EX_ac_e <= 999999
 0 <= R_EXinput_nh4_e <= 100
-999999 <= R_PPM <= 999999
 0 <= R_EX_2ddglcn_e <= 999999
 0 <= R_SUCCt2b <= 999999
 0 <= R_CITt7 <= 999999
 0 <= R_EX_orn_e <= 999999
-999999 <= R_PSCVT <= 999999
 0 <= R_MACPD <= 999999
 0 <= R_EX_lcts_e <= 999999
 0 <= R_PEAMNO <= 999999
 0 <= R_TDPDRR <= 999999
-999999 <= R_FCI <= 999999
 0 <= R_ADSS <= 999999
-999999 <= R_GLYCLTt2r <= 999999
 0 <= R_THD2 <= 999999
 0 <= R_EX_csn_e <= 999999
-999999 <= R_DTMPK <= 999999
 0 <= R_SERt4 <= 999999
 0 <= R_ASPt2_3 <= 999999
 0 <= R_ASPt2_2 <= 999999
-999999 <= R_ACOTA <= 999999
 0 <= R_TRE6PS <= 999999
 0 <= R_TRE6PP <= 999999
 0 <= R_EXinput_ile_L_e <= 100
-999999 <= R_DHQD <= 999999
-999999 <= R_ALAt2r <= 999999
 0 <= R_EX_malthx_e <= 999999
-999999 <= R_PTRCORNt7 <= 999999
-999999 <= R_CBIAT_DELETE <= 999999
 0 <= R_EXoutput_k_e <= 1000
 0 <= R_MOHMT <= 999999
 0 <= R_ACLS <= 999999
-999999 <= R_INDOLEt2r <= 999999
-999999 <= R_ASNt2r <= 999999
 0 <= R_GLUt4 <= 999999
 0 <= R_TRE6PH <= 999999
 0 <= R_DHQS <= 999999
 0 <= R_GLCS1 <= 999999
-999999 <= R_A5PISO <= 999999
 0 <= R_G1PP <= 999999
 0 <= R_EX_15dap_e <= 999999
 0 <= R_MALTPTabc <= 999999
 0 <= R_S7PI <= 999999
 0 <= R_MDRPD <= 999999
 0 <= R_C160SN <= 999999
 0 <= R_TMAOR1e <= 999999
 0 <= R_PTRCTA <= 999999
 0 <= R_THZPSN <= 999999
 0 <= R_EX_no3_e <= 999999
 0 <= R_EX_26dap_M_e <= 999999
 0 <= R_EX_melib_e <= 999999
 0 <= R_EX_gbbtn_e <= 999999
-999999 <= R_PUNP3 <= 999999
 0 <= R_EX_galctn_D_e <= 999999
 0 <= R_PAPPT3 <= 999999
 0 <= R_SPMDabc <= 999999
 0 <= R_PPBNGS <= 999999
-999999 <= R_CYTK1 <= 999999
 0 <= R_CYTDK2 <= 999999
 0 <= R_RBFK <= 999999
 0 <= R_HMBS <= 999999
 0 <= R_DAGK_EC <= 999999
 0 <= R_CMPN <= 999999
 0 <= R_GLYCTO4 <= 999999
 0 <= R_EX_glyb_e <= 999999
 0 <= R_LPLIPA5 <= 999999
 0 <= R_ACGApts <= 999999
 0 <= R_UPP3MT <= 999999
-999999 <= R_DHBD <= 999999
 0 <= R_MALTTTRabc <= 999999
-999999 <= R_PEPT_EC <= 999999
 0 <= R_EX_glyc3p_e <= 999999
 0 <= R_DDGALK <= 999999
 0 <= R_C181SN <= 999999
 0 <= R_P5CR <= 999999
-999999 <= R_GALCTt2r <= 999999
 0 <= R_CYNTt2 <= 999999
-999999 <= R_ALATA_L <= 999999
-999999 <= R_NH4t <= 999999
-999999 <= R_GTHOr <= 999999
 0 <= R_G3PD7 <= 999999
 0 <= R_G3PD5 <= 999999
-999999 <= R_G3PD2 <= 999999
 0 <= R_THRD <= 999999
-999999 <= R_AIRC3 <= 999999
 0 <= R_AIRC2 <= 999999
 0 <= R_CYNTAH <= 999999
 0 <= R_PYAM5PO <= 999999
 0 <= R_THRt4 <= 999999
-999999 <= R_TYRt2r <= 999999
-999999 <= R_ACKr <= 999999
 0 <= R_HOPNTAL <= 999999
 0 <= R_EX_xtsn_e <= 999999
 0 <= R_P5CD <= 999999
-999999 <= R_ADNt2r <= 999999
 0 <= R_SGSAD <= 999999
 0 <= R_EX_gln_L_e <= 999999
 0 <= R_NMNt7 <= 999999
 0 <= R_OCDCAt2 <= 999999
 0 <= R_DXPS <= 999999
-999999 <= R_GALU <= 999999
-999999 <= R_SUCD4 <= 999999
 0 <= R_CYANST <= 999999
 0 <= R_IGPDH <= 999999
 0 <= R_MGSA <= 999999
 0 <= R_POX <= 999999
 0 <= R_LPLIPA2 <= 999999
 0 <= R_LPLIPA3 <= 999999
 0 <= R_DMSOR2e <= 999999
 0 <= R_LPLIPA1 <= 999999
 0 <= R_DAPDC <= 999999
 0 <= R_EX_dms_e <= 999999
 0 <= R_LPLIPA4 <= 999999
-999999 <= R_URIt2r <= 999999
 0 <= R_GLYC3Pt6 <= 999999
 0 <= R_EX_ocdca_e <= 999999
 0 <= R_CYTBD <= 999999
 0 <= R_METabc <= 999999
-999999 <= R_ETOHt2r <= 999999
 0 <= R_AACPS3 <= 999999
 0 <= R_EX_taur_e <= 999999
 0 <= R_G6Pt6_2 <= 999999
 0 <= R_PMDPHT <= 999999
 0 <= R_XYLK <= 999999
 0 <= R_SADH <= 999999
-999999 <= R_CRNCBCT <= 999999
 0 <= R_XYLabc <= 999999
-999999 <= R_TKT1 <= 999999
 0 <= R_LDH_D2 <= 999999
-999999 <= R_EX_co2_e <= 999999
 0 <= R_NO3R2 <= 999999
 0 <= R_CYTBO3 <= 999999
 0 <= R_NO3R1 <= 999999
 0 <= R_FBP <= 999999
 0 <= R_EX_mnl_e <= 999999
 0 <= R_EXinput_ala_L_e <= 100
-999999 <= R_AICART <= 999999
-999999 <= R_GLYt2r <= 999999
 0 <= R_FUMt2_3 <= 999999
 0 <= R_FUMt2_2 <= 999999
 0 <= R_XANt2 <= 999999
-999999 <= R_XTSNt2r <= 999999
 0 <= R_ETHAAL <= 999999
 0 <= R_EXinput_gam_e <= 100
 0 <= R_EXoutput_galct_D_e <= 1000
 0 <= R_RNTR4 <= 999999
 0 <= R_ALDD19x <= 999999
 0 <= R_FDH3 <= 999999
 0 <= R_FDH2 <= 999999
 0 <= R_GALTpts <= 999999
-999999 <= R_TARTRt7 <= 999999
 0 <= R_PIabc <= 999999
 0 <= R_NTD10 <= 999999
 0 <= R_ADSK <= 999999
-999999 <= R_CRNCDH <= 999999
 0 <= R_DHAD2 <= 999999
 0 <= R_TMAOR2e <= 999999
 0 <= R_MTRK <= 999999
 0 <= R_DHFS <= 999999
-999999 <= R_MTRI <= 999999
 0 <= R_BETALDHy <= 999999
 0 <= R_BETALDHx <= 999999
 0 <= R_TDPGDH <= 999999
-999999 <= R_LEUt2r <= 999999
 0 <= R_EX_acgam_e <= 999999
 0 <= R_MECDPDH <= 999999
 0 <= R_PPPNDO <= 999999
 0 <= R_RMK <= 999999
 0 <= R_DKGLCNR2y <= 999999
 0 <= R_NTD4 <= 999999
 0 <= R_EX_duri_e <= 999999
 0 <= R_HPPPNDO <= 999999
 0 <= R_NTD7 <= 999999
 0 <= R_METS <= 999999
-999999 <= R_ACACT1r <= 999999
 0 <= R_MAN1PT2 <= 999999
 0 <= R_OCTDPS <= 999999
 0 <= R_EX_ala_D_e <= 999999
 0 <= R_MEPCT <= 999999
 0 <= R_HKNDDH <= 999999
-999999 <= R_URIDK2r <= 999999
 0 <= R_PGLYCP <= 999999
 0 <= R_PROD2 <= 999999
 0 <= R_UAAGDS <= 999999
 0 <= R_PSP_L <= 999999
 0 <= R_AP5AH <= 999999
-999999 <= R_ICDHyr <= 999999
 0 <= R_EX_fuc1p_L_e <= 999999
 0 <= R_EXoutput_h2o_e <= 1000
 0 <= R_SOTA <= 999999
 0 <= R_EX_arg_L_e <= 999999
 0 <= R_EX_idon_L_e <= 999999
 0 <= R_C140SN <= 999999
 0 <= R_PDX5PS <= 999999
 0 <= R_PFL <= 999999
 0 <= R_CSND <= 999999
-999999 <= R_OHPBAT <= 999999
 0 <= R_RZ5PP <= 999999
 0 <= R_EXoutput_malttr_e <= 1000
 0 <= R_EXoutput_acnam_e <= 1000
 0 <= R_UAPGR <= 999999
 0 <= R_EX_tsul_e <= 999999
-999999 <= R_SDPTA <= 999999
 0 <= R_DADNt2 <= 999999
      R_EXinput_glc_e = 0
 0 <= R_FHL <= 999999
 0 <= R_ARGabc <= 999999
 0 <= R_DDPA <= 999999
 0 <= R_TMPPP <= 999999
 0 <= R_QULNS <= 999999
 0 <= R_DKGLCNR1 <= 999999
 0 <= R_SGDS <= 999999
 0 <= R_GMHEPK <= 999999
 0 <= R_EXinput_asp_L_e <= 100
 0 <= R_TTDCAt2 <= 999999
-999999 <= R_GALURt2r <= 999999
 0 <= R_ALAabc <= 999999
 0 <= R_NADK <= 999999
-999999 <= R_ATPS4r <= 999999
 0 <= R_SULabc <= 999999
 0 <= R_DRPA <= 999999
 0 <= R_CYSDS <= 999999
 0 <= R_PACCOAL <= 999999
 0 <= R_EX_malt_e <= 999999
 0 <= R_EX_maltpt_e <= 999999
 0 <= R_CYTD <= 999999
-999999 <= R_PERD <= 999999
 0 <= R_MDH2 <= 999999
 0 <= R_MDH3 <= 999999
-999999 <= R_TAGURr <= 999999
-999999 <= R_GLCURt2r <= 999999
 0 <= R_UDPGALM <= 999999
 0 <= R_XPPT <= 999999
 0 <= R_ADD <= 999999
-999999 <= R_PPPNt2r <= 999999
-999999 <= R_PMANM <= 999999
-999999 <= R_ARBt2r <= 999999
-999999 <= R_MTHFD <= 999999
-999999 <= R_MTHFC <= 999999
 0 <= R_GTHS <= 999999
 0 <= R_GSPMDA <= 999999
 0 <= R_PGCD <= 999999
 0 <= R_HDCAt2 <= 999999
 0 <= R_KDOCT2 <= 999999
 0 <= R_THRD_L <= 999999
 0 <= R_EX_sbt_D_e <= 999999
-999999 <= R_TGBPA <= 999999
-999999 <= R_GK1 <= 999999
 0 <= R_GSPMDS <= 999999
 0 <= R_OXGDC2 <= 999999
 0 <= R_UNK3 <= 999999
 0 <= R_PPCSCT <= 999999
 0 <= R_VPAMT <= 999999
 0 <= R_EX_din_e <= 999999
 0 <= R_EX_fum_e <= 999999
 0 <= R_MTAN <= 999999
 0 <= R_OPHBDC <= 999999
 0 <= R_EXoutput_uri_e <= 1000
 0 <= R_EX_amp_e <= 999999
 0 <= R_AP4AH <= 999999
-999999 <= R_PGMT <= 999999
 0 <= R_HKNTDH <= 999999
 0 <= R_PMPK <= 999999
-999999 <= R_GALCTNt2r <= 999999
 0 <= R_KG6PDC <= 999999
-999999 <= R_ARGORNt7 <= 999999
 0 <= R_MALt2_3 <= 999999
 0 <= R_MALt2_2 <= 999999
 0 <= R_EX_cytd_e <= 999999
 0 <= R_TMAOR1 <= 999999
 0 <= R_ADOCBIK <= 999999
 0 <= R_G1SATi <= 999999
-999999 <= R_ARAI <= 999999
 0 <= R_MLTG2 <= 999999
 0 <= R_MLTG3 <= 999999
 0 <= R_MLTG4 <= 999999
 0 <= R_MLTG5 <= 999999
 0 <= R_GLNS <= 999999
 0 <= R_GLYCDx <= 999999
-999999 <= R_VALt2r <= 999999
 0 <= R_UGMDDS <= 999999
 0 <= R_NAt3_15 <= 999999
 0 <= R_EX_pro_L_e <= 999999
-999999 <= R_ACt2r <= 999999
 0 <= R_EXoutput_glyclt_e <= 1000
 0 <= R_TDPDRE <= 999999
 0 <= R_EX_etoh_e <= 999999
 0 <= R_ME1 <= 999999
 0 <= R_ME2 <= 999999
 0 <= R_ALTRH <= 999999
 0 <= R_EX_man_e <= 999999
 0 <= R_G1PTT <= 999999
 0 <= R_GAMpts <= 999999
 0 <= R_EX_lac_L_e <= 999999
 0 <= R_EX_hxan_e <= 999999
 0 <= R_EX_sucr_e <= 999999
 0 <= R_DHPTDCs <= 999999
-999999 <= R_IPPMIa <= 999999
-999999 <= R_IPPMIb <= 999999
-999999 <= R_UAGAAT <= 999999
-999999 <= R_HISt2r <= 999999
 0 <= R_EX_urea_e <= 999999
-999999 <= R_DAPE <= 999999
 0 <= R_SUCD1i <= 999999
 0 <= R_EX_trp_L_e <= 999999
 0 <= R_RBK <= 999999
-999999 <= R_TKT2 <= 999999
-999999 <= R_SUCOAS <= 999999
 0 <= R_ASPO5 <= 999999
 0 <= R_DINSt2 <= 999999
 0 <= R_NTD11 <= 999999
-999999 <= R_INSt2r <= 999999
 0 <= R_UPP3S <= 999999
 0 <= R_UPPDC1 <= 999999
 0 <= R_DCYTt2 <= 999999
 0 <= R_NTD6 <= 999999
 0 <= R_ABUTD <= 999999
 0 <= R_GLYC3Pabc <= 999999
 0 <= R_EX_cbl1_e <= 999999
 0 <= R_EX_ptrc_e <= 999999
 0 <= R_PPPGO <= 999999
 0 <= R_OMPDC <= 999999
 0 <= R_EXoutput_pi_e <= 1000
 0 <= R_LYSabc <= 999999
 0 <= R_EXinput_tre_e <= 100
 0 <= R_ACNAMt2 <= 999999
 0 <= R_NTRIR2x <= 999999
-999999 <= R_GLCNt2r <= 999999
 0 <= R_HXPRT <= 999999
-999999 <= R_IDONt2r <= 999999
 0 <= R_CYSS <= 999999
 0 <= R_MI1PP <= 999999
-999999 <= R_SERAT <= 999999
 0 <= R_NPHS <= 999999
 0 <= R_CBPS <= 999999
 0 <= R_PYDAMK <= 999999
 0 <= R_AMANAPE <= 999999
-999999 <= R_ADK1 <= 999999
 0 <= R_RNDR2 <= 999999
-999999 <= R_ADK3 <= 999999
 0 <= R_MOAT2 <= 999999
-999999 <= R_ADK4 <= 999999
 0 <= R_RNDR3 <= 999999
 0 <= R_EXoutput_3hpppn_e <= 1000
 0 <= R_DHPPD <= 999999
-999999 <= R_ACALDt <= 999999
 0 <= R_NADH10 <= 999999
 0 <= R_AB6PGH <= 999999
 0 <= R_HYPOE <= 999999
 0 <= R_KAS14 <= 999999
-999999 <= R_RPE <= 999999
 0 <= R_KAS16 <= 999999
 0 <= R_TDPAGTA <= 999999
-999999 <= R_ORPT <= 999999
-999999 <= R_PYRt2r <= 999999
 0 <= R_EX_g6p_e <= 999999
-999999 <= R_5DGLCNR <= 999999
-999999 <= R_RPI <= 999999
 0 <= R_OBTFL <= 999999
 0 <= R_ACS <= 999999
 0 <= R_SHCHF <= 999999
 0 <= R_IPPS <= 999999
 0 <= R_EXinput_man6p_e <= 100
 0 <= R_XYLt2 <= 999999
 0 <= R_DHDPS <= 999999
 0 <= R_PRFGS <= 999999
 0 <= R_EX_tmao_e <= 999999
-999999 <= R_FUCt <= 999999
 0 <= R_ASNS1 <= 999999
 0 <= R_URIK2 <= 999999
-999999 <= R_ALCD19 <= 999999
-999999 <= R_CRNBTCT <= 999999
 0 <= R_HMPK1 <= 999999
 0 <= R_THMabc <= 999999
-999999 <= R_PGAMT <= 999999
-999999 <= R_DDPGALA <= 999999
 0 <= R_EX_glcr_e <= 999999
 0 <= R_EX_4abut_e <= 999999
 0 <= R_EX_galur_e <= 999999
-999999 <= R_ACACt2 <= 999999
 0 <= R_GLXCL <= 999999
 0 <= R_NADDPe <= 999999
-999999 <= R_AHC <= 999999
 0 <= R_CHORS <= 999999
 0 <= R_RIBabc <= 999999
 0 <= R_PFK_2 <= 999999
-999999 <= R_DHFR <= 999999
 0 <= R_CS <= 999999
 0 <= R_EX_lac_D_e <= 999999
-999999 <= R_HCINNMt2r <= 999999
 0 <= R_HISabc <= 999999
 0 <= R_EX_glcur_e <= 999999
 0 <= R_UDCPDPS <= 999999
 0 <= R_NNDPR <= 999999
 0 <= R_ILEabc <= 999999
 0 <= R_DCTPD <= 999999
 0 <= R_EX_glu_L_e <= 999999
 0 <= R_GLUSy <= 999999
 0 <= R_TRSAR <= 999999
 0 <= R_EX_ade_e <= 999999
 0 <= R_GLYBabc <= 999999
 0 <= R_ADPT <= 999999
-999999 <= R_FBA <= 999999
 0 <= R_DHORD2 <= 999999
 0 <= R_GLYK <= 999999
 0 <= R_HISTP <= 999999
-999999 <= R_GALKr <= 999999
 0 <= R_IG3PS <= 999999
 0 <= R_DMATT <= 999999
 0 <= R_MICITD <= 999999
-999999 <= R_HPPPNt2r <= 999999
 0 <= R_PASYN_EC <= 999999
 0 <= R_EX_fru_e <= 999999
-999999 <= R_THRt2r <= 999999
 0 <= R_PDXPP <= 999999
 0 <= R_HISTD <= 999999
 0 <= R_GLCRD <= 999999
 0 <= R_OPHHX <= 999999
 0 <= R_EX_acmana_e <= 999999
-999999 <= R_BTS2 <= 999999
-20 <= R_EX_o2_e <= 0
 0 <= R_NMNDA <= 999999
 0 <= R_KAS15 <= 999999
 0 <= R_FTHFD <= 999999
 0 <= R_EX_akg_e <= 999999
-999999 <= R_PIt2r <= 999999
 0 <= R_CHLabc <= 999999
 0 <= R_NACUP <= 999999
 0 <= R_EX_mal_L_e <= 999999
 0 <= R_EXinput_pi_e <= 100
 0 <= R_EX_dha_e <= 999999
 0 <= R_ORNTA <= 999999
 0 <= R_GMPS2 <= 999999
 0 <= R_NNAT <= 999999
 0 <= R_DKGLCNR2x <= 999999
 0 <= R_AMALT1 <= 999999
 0 <= R_EXoutput_ile_L_e <= 1000
 0 <= R_AMALT3 <= 999999
 0 <= R_AMALT4 <= 999999
 0 <= R_EX_his_L_e <= 999999
 0 <= R_EX_cynt_e <= 999999
 0 <= R_ASPCT <= 999999
 0 <= R_THRabc <= 999999
 0 <= R_GALCTND <= 999999
-999999 <= R_EX_h_e <= 999999
 0 <= R_GTPCII2 <= 999999
 0 <= R_EX_arab_L_e <= 999999
 0 <= R_MELIBt2 <= 999999
 0 <= R_HEMEOS <= 999999
 0 <= R_GART <= 999999
-999999 <= R_THRAr <= 999999
-999999 <= R_GLTPD <= 999999
 0 <= R_HPYRRx <= 999999
 0 <= R_HPYRRy <= 999999
 0 <= R_BSORx <= 999999
 0 <= R_EX_but_e <= 999999
 0 <= R_NADH5 <= 999999
 0 <= R_ANS <= 999999
 0 <= R_NADH7 <= 999999
 0 <= R_NADH6 <= 999999
 0 <= R_NADH9 <= 999999
 0 <= R_NADH8 <= 999999
 0 <= R_EXinput_3hpppn_e <= 100
-999999 <= R_PRASCS <= 999999
 0 <= R_CDAPPA_EC <= 999999
-999999 <= R_XANt <= 999999
 0 <= R_GMAND <= 999999
 0 <= R_DHAD1 <= 999999
 0 <= R_DADA <= 999999
 0 <= R_G1PACT <= 999999
 0 <= R_EX_acac_e <= 999999
-999999 <= R_DADK <= 999999
 0 <= R_FUCPt6_2 <= 999999
 0 <= R_EX_rib_D_e <= 999999
-999999 <= R_ALAALAr <= 999999
 0 <= R_GLUPRT <= 999999
-999999 <= R_E4PD <= 999999
 0 <= R_GLCDe <= 999999
 0 <= R_CITL <= 999999
 0 <= R_FAO1 <= 999999
 0 <= R_CADVt <= 999999
 0 <= R_FAO2 <= 999999
 0 <= R_FRD3 <= 999999
 0 <= R_UAGDP <= 999999
 0 <= R_SHSL1 <= 999999
 0 <= R_AMALT2 <= 999999
 0 <= R_EX_met_L_e <= 999999
 0 <= R_LPADSS <= 999999
-999999 <= R_PPAKr <= 999999
-999999 <= R_ADMDCr <= 999999
 0 <= R_EXinput_val_L_e <= 100
 0 <= R_ARGSS <= 999999
 0 <= R_HPPK2 <= 999999
 0 <= R_ACMANApts <= 999999
-999999 <= R_DBTSr <= 999999
-999999 <= R_D_LACt2 <= 999999
 0 <= R_EXoutput_tre_e <= 1000
 0 <= R_SPMDAT1 <= 999999
 0 <= R_EX_ins_e <= 999999
 0 <= R_SPMDAT2 <= 999999
 0 <= R_X5PL3E <= 999999
 0 <= R_PPCK <= 999999
-999999 <= R_CO2t <= 999999
 0 <= R_AHCYSNS <= 999999
 0 <= R_LEUabc <= 999999
 0 <= R_ECAP_EC <= 999999
 0 <= R_ICHORSi <= 999999
 0 <= R_PTPATi <= 999999
-999999 <= R_L_LACt2r <= 999999
 0 <= R_GALCTD <= 999999
 0 <= R_EXoutput_3hcinnm_e <= 1000
 0 <= R_ASNS2 <= 999999
 0 <= R_ASPabc <= 999999
 0 <= R_EX_adn_e <= 999999
 0 <= R_EX_indole_e <= 999999
 0 <= R_PNTOt4 <= 999999
-999999 <= R_PTRCt2r <= 999999
 0 <= R_EX_galt_e <= 999999
 0 <= R_CYSTL <= 999999
 0 <= R_EX_maltttr_e <= 999999
 0 <= R_ADA <= 999999
 0 <= R_AGMHE <= 999999
-999999 <= R_HSDy <= 999999
-999999 <= R_SULR <= 999999
 0 <= R_EX_nad_e <= 999999
 0 <= R_GUAD <= 999999
 0 <= R_EX_hdca_e <= 999999
 0 <= R_ADNCYC <= 999999
 0 <= R_L_LACD3 <= 999999
 0 <= R_L_LACD2 <= 999999
 0 <= R_EXinput_acnam_e <= 100
 0 <= R_EX_tma_e <= 999999
 0 <= R_KDOPP <= 999999
 0 <= R_EX_gal_e <= 999999
-999999 <= R_G6PDH2r <= 999999
 0 <= R_FAO3 <= 999999
-999999 <= R_MCOATA <= 999999
 0 <= R_EXinput_galct_D_e <= 100
 0 <= R_SELNPS <= 999999
-999999 <= R_XYLI1 <= 999999
 0 <= R_EXoutput_glcn_e <= 1000
-999999 <= R_PYNP2r <= 999999
 0 <= R_ASPT <= 999999
 0 <= R_MMCD <= 999999
 0 <= R_PRATPP <= 999999
 0 <= R_EX_gua_e <= 999999
 0 <= R_DXPRIi <= 999999
-999999 <= R_RMI <= 999999
-999999 <= R_GARFT <= 999999
-999999 <= R_HPYRI <= 999999
 0 <= R_DCYTD <= 999999
 0 <= R_G3PD6 <= 999999
 0 <= R_PRAIS <= 999999
-999999 <= R_ASPK <= 999999
 0 <= R_EXinput_lys_L_e <= 100
 0 <= R_FRD2 <= 999999
 0 <= bin_R_EXinput_pi_e <= 1
 0 <= bin_R_EXinput_nh4_e <= 1
 0 <= bin_R_EXinput_uri_e <= 1
 0 <= bin_R_EXinput_asp_L_e <= 1
 0 <= bin_R_EXinput_glyclt_e <= 1
 0 <= bin_R_EXinput_tre_e <= 1
 0 <= bin_R_EXinput_glcn_e <= 1
 0 <= bin_R_EXinput_k_e <= 1
 0 <= bin_R_EXinput_glc_e <= 1
 0 <= bin_R_EXinput_dgsn_e <= 1
 0 <= bin_R_EXinput_gly_e <= 1
 0 <= bin_R_EXinput_12ppd_S_e <= 1
 0 <= bin_R_EXinput_thm_e <= 1
 0 <= bin_R_EXinput_so4_e <= 1
 0 <= bin_R_EXinput_gam_e <= 1
 0 <= bin_R_EXinput_galct_D_e <= 1
 0 <= bin_R_EXinput_3hcinnm_e <= 1
 0 <= bin_R_EXinput_malttr_e <= 1
 0 <= bin_R_EXinput_man6p_e <= 1
 0 <= bin_R_EXinput_ile_L_e <= 1
 0 <= bin_R_EXinput_lys_L_e <= 1
 0 <= bin_R_EXinput_3hpppn_e <= 1
 0 <= bin_R_EXinput_ura_e <= 1
 0 <= bin_R_EXinput_xan_e <= 1
 0 <= bin_R_EXinput_fe2_e <= 1
 0 <= bin_R_EXinput_succ_e <= 1
 0 <= bin_R_EXinput_acnam_e <= 1
 0 <= bin_R_EXinput_val_L_e <= 1
 0 <= bin_R_EXinput_ala_L_e <= 1
Binaries
 bin_R_EXinput_pi_e  bin_R_EXinput_nh4_e  bin_R_EXinput_uri_e 
 bin_R_EXinput_asp_L_e  bin_R_EXinput_glyclt_e  bin_R_EXinput_tre_e 
 bin_R_EXinput_glcn_e  bin_R_EXinput_k_e  bin_R_EXinput_glc_e 
 bin_R_EXinput_dgsn_e  bin_R_EXinput_gly_e  bin_R_EXinput_12ppd_S_e 
 bin_R_EXinput_thm_e  bin_R_EXinput_so4_e  bin_R_EXinput_gam_e 
 bin_R_EXinput_galct_D_e  bin_R_EXinput_3hcinnm_e  bin_R_EXinput_malttr_e 
 bin_R_EXinput_man6p_e  bin_R_EXinput_ile_L_e  bin_R_EXinput_lys_L_e 
 bin_R_EXinput_3hpppn_e  bin_R_EXinput_ura_e  bin_R_EXinput_xan_e 
 bin_R_EXinput_fe2_e  bin_R_EXinput_succ_e  bin_R_EXinput_acnam_e 
 bin_R_EXinput_val_L_e  bin_R_EXinput_ala_L_e 
End
