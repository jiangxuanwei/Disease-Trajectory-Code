# ============================================================
# Multiomics Signatures Identification
# Author: YGR
# Last Modified: Feb 2025
# ============================================================

# Load required libraries
library("mstate")
library(dplyr)
library(openxlsx)
library(survival)
library(dplyr)
library(openxlsx)
data<-read.csv("CMD_CA.csv")
df6 <- read.csv("covariates_impute.csv")
protein <- read.csv("protein_UKB_filled.csv")
data1 <-merge(data,df6,by="n_eid")
data2 <-merge(data1,protein,by.x="n_eid",by.y = "eid")
data1 <- data2
# Standardize protein measurements (columns 55-2974)
data1[, 55:2974] <- scale(data1[, 55:2974])

# Create training/test splits
train_data_birth <- data1[data1$birthplace == 1, ]
test_data_birth <- data1[is.na(data1$birthplace) | data1$birthplace != 1, ]

# Covariates for adjustment
covariate <- c("Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", 
               "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", 
               "cancer_treat", "CMD_treat", "cancer_screen")
covs_P <-   c("Age", "Sex", "Ethnicity", "Edu", "Employed", "Tdi", "Smoke", "Drink", 
              "Mets", "Dietscore", "Sleeptime", "BMI", "T2dhis", "Cvdhis", "hiscancer", 
              "cancer_treat", "CMD_treat", "cancer_screen")
# List of proteins to analyze (2920 proteins)
traits <- c("aarsd1", "abhd14b","abl1", "acaa1", "acan", "ace2", "acox1", "acp5", "acp6", "acta2", 
                "actn4", "acvrl1", "acy1", "ada", "ada2", "adam15", "adam22", 
                "adam23", "adam8", "adamts13", "adamts15", "adamts16", "adamts8", 
                "adcyap1r1", "adgrb3", "adgre2", "adgre5", "adgrg1", "adgrg2", 
                "adh4", "adm", "afp", "ager", "agr2", "agr3", "agrn", "agrp", 
                "agxt", "ahcy", "ahsp", "aif1", "aifm1", "ak1", "akr1b1", "akr1c4", 
                "akt1s1", "akt3", "alcam", "aldh1a1", "aldh3a1", "alpp", "ambn", 
                "ambp", "amfr", "amigo2", "amn", "amy2a", "amy2b", "ang", "angpt1", 
                "angpt2", "angptl1", "angptl2", "angptl3", "angptl4", "angptl7", 
                "ankrd54", "anpep", "anxa10", "anxa11", "anxa3", "anxa4", "anxa5", 
                "aoc1", "aoc3", "apbb1ip", "apex1", "aplp1", "apoh", "apom", 
                "app", "aprt", "areg", "arg1", "arhgap1", "arhgap25", "arhgef12", 
                "arid4b", "arnt", "arsa", "arsb", "art3", "artn", "asah2", "asgr1", 
                "atf2", "atg4a", "atox1", "atp5if1", "atp5po", "atp6ap2", "atp6v1d", 
                "atp6v1f", "atxn10", "axin1", "axl", "azu1", "b4galt1", "b4gat1", 
                "bach1", "bag3", "bag6", "baiap2", "bambi", "bank1", "bax", "bcam", 
                "bcan", "bcl2l11", "bcr", "bgn", "bid", "bin2", "birc2", "blmh", 
                "blvrb", "bmp4", "bmp6", "boc", "bpifb1", "brk1", "bsg", "bst1", 
                "bst2", "btc", "btn2a1", "btn3a2", "c19orf12", "c1qa", "c1qtnf1", 
                "c2", "c2cd2l", "c4bpb", "ca1", "ca11", "ca12", "ca13", "ca14", 
                "ca2", "ca3", "ca4", "ca5a", "ca6", "ca9", "calb1", "calb2", 
                "calca", "calcoco1", "camkk1", "cant1", "capg", "carhsp1", "casp1", 
                "casp10", "casp2", "casp3", "casp8", "cblif", "cbln4", "cc2d1a", 
                "ccdc80", "ccl11", "ccl13", "ccl14", "ccl15", "ccl16", "ccl17", 
                "ccl18", "ccl19", "ccl2", "ccl20", "ccl21", "ccl22", "ccl23", 
                "ccl24", "ccl25", "ccl26", "ccl27", "ccl28", "ccl3", "ccl4", 
                "ccl5", "ccl7", "ccl8", "ccn1", "ccn2", "ccn3", "ccn4", "ccn5", 
                "ccs", "cct5", "cd109", "cd14", "cd160", "cd163", "cd164", "cd177", 
                "cd1c", "cd200", "cd200r1", "cd207", "cd209", "cd22", "cd244", 
                "cd27", "cd274", "cd276", "cd28", "cd2ap", "cd300c", "cd300e", 
                "cd300lf", "cd300lg", "cd302", "cd33", "cd34", "cd38", "cd4", 
                "cd40", "cd40lg", "cd46", "cd48", "cd5", "cd55", "cd58", "cd59", 
                "cd6", "cd63", "cd69", "cd70", "cd74", "cd79b", "cd83", "cd84", 
                "cd8a", "cd93", "cd99", "cd99l2", "cdc27", "cdc37", "cdcp1", 
                "cdh1", "cdh15", "cdh17", "cdh2", "cdh3", "cdh5", "cdh6", "cdhr1", 
                "cdhr2", "cdhr5", "cdkn1a", "cdkn2d", "cdnf", "cdon", "cdsn", 
                "ceacam1", "ceacam21", "ceacam3", "ceacam5", "ceacam8", "cebpb", 
                "cela3a", "cep164", "cep20", "cep43", "cep85", "cert", "ces1", 
                "ces2", "ces3", "cetn2", "cfc1", "cga", "cgref1", "chac2", "chek2", 
                "chgb", "chi3l1", "chit1", "chl1", "chmp1a", "chrdl1", "chrdl2", 
                "ciapin1", "ckap4", "ckmt1a_ckmt1b", "clc", "clec10a", "clec11a", 
                "clec14a", "clec1a", "clec1b", "clec4a", "clec4c", "clec4d", 
                "clec4g", "clec5a", "clec6a", "clec7a", "clip2", "clmp", "clpp", 
                "clps", "clspn", "clstn1", "clstn2", "clta", "clul1", "cndp1", 
                "cnpy2", "cnpy4", "cnst", "cntn1", "cntn2", "cntn3", "cntn4", 
                "cntn5", "cntnap2", "col18a1", "col1a1", "col4a1", "col6a3", 
                "col9a1", "colec12", "comp", "comt", "cope", "coro1a", "cox5b", 
                "cpa1", "cpa2", "cpb1", "cpe", "cpm", "cpped1", "cpvl", "cpxm1", 
                "cr2", "cracr2a", "cradd", "creg1", "creld2", "crh", "crhbp", 
                "crhr1", "crim1", "crip2", "crisp2", "crkl", "crlf1", "crnn", 
                "crtac1", "crtam", "crx", "csf1", "csf2ra", "csf3", "cst3", "cst5", 
                "cst6", "cst7", "cstb", "ctf1", "ctrb1", "ctrc", "ctsb", "ctsc", 
                "ctsd", "ctsf", "ctsh", "ctsl", "ctso", "ctss", "ctsv", "ctsz", 
                "cx3cl1", "cxadr", "cxcl1", "cxcl10", "cxcl11", "cxcl12", "cxcl13", 
                "cxcl14", "cxcl16", "cxcl17", "cxcl3", "cxcl5", "cxcl6", "cxcl8", 
                "cxcl9", "dab2", "dag1", "dapp1", "dars1", "dbi", "dbnl", "dcbld2", 
                "dcn", "dctn1", "dctn2", "dctn6", "dctpp1", "dcxr", "ddah1", 
                "ddc", "ddr1", "ddx58", "decr1", "defa1_defa1b", "defb4a_defb4b", 
                "dffa", "dgkz", "diablo", "dkk1", "dkk3", "dkk4", "dkkl1", "dlk1", 
                "dll1", "dnaja2", "dnajb1", "dnajb8", "dner", "dnmbp", "dnph1", 
                "dok2", "dpep1", "dpep2", "dpp10", "dpp4", "dpp6", "dpp7", "dpt", 
                "dpy30", "draxin", "drg2", "dsc2", "dsg2", "dsg3", "dsg4", "dtx3", 
                "duox2", "dusp3", "ebag9", "ebi3_il27", "ece1", "eda2r", "edar", 
                "edil3", "efemp1", "efna1", "efna4", "egf", "egfl7", "egfr", 
                "egln1", "eif4b", "eif4ebp1", "eif4g1", "eif5a", "eloa", "enah", 
                "eng", "eno1", "eno2", "enpp2", "enpp5", "enpp7", "entpd2", "entpd5", 
                "entpd6", "epcam", "epha1", "epha10", "epha2", "ephb4", "ephb6", 
                "ephx2", "epo", "eps8l2", "erbb2", "erbb3", "erbb4", "erbin", 
                "ereg", "erp44", "esam", "esm1", "ezr", "f11r", "f2r", "f3", 
                "f7", "f9", "fabp1", "fabp2", "fabp4", "fabp5", "fabp6", "fabp9", 
                "fadd", "fam3b", "fam3c", "fap", "fas", "faslg", "fbp1", "fcar", 
                "fcer2", "fcgr2a", "fcgr2b", "fcgr3b", "fcn2", "fcrl1", "fcrl2", 
                "fcrl3", "fcrl5", "fcrl6", "fcrlb", "fen1", "fes", "fetub", "fgf19", 
                "fgf2", "fgf21", "fgf23", "fgf5", "fgfbp1", "fgfr2", "fgr", "fhit", 
                "fis1", "fkbp1b", "fkbp4", "fkbp5", "fkbp7", "fli1", "flrt2", 
                "flt1", "flt3", "flt3lg", "flt4", "fmnl1", "fmr1", "folr1", "folr2", 
                "folr3", "fosb", "foxo1", "foxo3", "frzb", "fst", "fstl3", "fuca1", 
                "furin", "fus", "fut3_fut5", "fut8", "fxn", "fxyd5", "fyb1", 
                "gal", "galnt10", "galnt2", "galnt3", "galnt7", "gas6", "gbp2", 
                "gbp4", "gcg", "gcnt1", "gdf15", "gdf2", "gdnf", "gfap", "gfer", 
                "gfod2", "gfra1", "gfra2", "gfra3", "gga1", "ggh", "ggt1", "ggt5", 
                "gh1", "gh2", "ghrhr", "ghrl", "gkn1", "glb1", "glo1", "glod4", 
                "glrx", "glt8d2", "gmpr", "gne", "gnly", "golm2", "gopc", "gp1ba", 
                "gp2", "gp6", "gpa33", "gpc1", "gpc5", "gpkow", "gpnmb", "gpr37", 
                "grap2", "grk5", "grn", "grpel1", "gsap", "gsta1", "gsta3", "gstp1", 
                "guca2a", "gusb", "gys1", "gzma", "gzmb", "gzmh", "hagh", "hao1", 
                "hars1", "havcr1", "havcr2", "hbegf", "hbq1", "hcls1", "hdgf", 
                "hebp1", "hexim1", "hgf", "hgs", "hk2", "hla_dra", "hla_e", "hmbs", 
                "hmox1", "hmox2", "hnmt", "hnrnpk", "hpcal1", "hpgds", "hs3st3b1", 
                "hs6st1", "hsd11b1", "hsp90b1", "hspa1a", "hspb1", "hspb6", "hspg2", 
                "htra2", "hyal1", "hyou1", "ica1", "icam1", "icam2", "icam3", 
                "icam4", "icam5", "icoslg", "idi2", "ids", "idua", "ifng", "ifngr1", 
                "ifngr2", "ifnl1", "ifnlr1", "igf1r", "igf2r", "igfbp1", "igfbp2", 
                "igfbp3", "igfbp4", "igfbp6", "igfbp7", "igfbpl1", "igsf3", "igsf8", 
                "ikbkg", "ikzf2", "il10", "il10ra", "il10rb", "il11", "il12a_il12b", 
                "il12b", "il12rb1", "il13", "il13ra1", "il15", "il15ra", "il16", 
                "il17a", "il17c", "il17d", "il17f", "il17ra", "il17rb", "il18", 
                "il18bp", "il18r1", "il18rap", "il19", "il1a", "il1b", "il1r1", 
                "il1r2", "il1rap", "il1rl1", "il1rl2", "il1rn", "il2", "il20", 
                "il20ra", "il22ra1", "il24", "il2ra", "il2rb", "il32", "il33", 
                "il34", "il3ra", "il4", "il4r", "il5", "il5ra", "il6", "il6r", 
                "il6st", "il7", "il7r", "ilkap", "impa1", "ing1", "inhbc", "inpp1", 
                "inppl1", "ipcef1", "iqgap2", "irag2", "irak1", "irak4", "islr2", 
                "ism1", "itga11", "itga5", "itga6", "itgam", "itgav", "itgb1", 
                "itgb1bp1", "itgb1bp2", "itgb2", "itgb5", "itgb6", "itgb7", "itih3", 
                "itm2a", "ivd", "jam2", "jchain", "jun", "kazald1", "kcnip4", 
                "kdr", "kel", "kifbp", "kir2dl3", "kir3dl1", "kirrel2", "kit", 
                "kitlg", "klb", "klk1", "klk10", "klk11", "klk12", "klk13", "klk14", 
                "klk4", "klk6", "klk8", "klrb1", "klrd1", "krt14", "krt18", "krt19", 
                "krt5", "kyat1", "kynu", "l1cam", "lactb2", "lag3", "lair1", 
                "lair2", "lama4", "lamp2", "lamp3", "lap3", "lat", "lat2", "layn", 
                "lbp", "lbr", "lcn2", "ldlr", "lefty2", "lep", "lepr", "lgals1", 
                "lgals3", "lgals4", "lgals7_lgals7b", "lgals8", "lgals9", "lgmn", 
                "lhb", "lhpp", "lif", "lifr", "lilra2", "lilra5", "lilrb1", "lilrb2", 
                "lilrb4", "lilrb5", "lpcat2", "lpl", "lpo", "lrig1", "lrp1", 
                "lrp11", "lrpap1", "lrrc25", "lrrn1", "lsm1", "lsp1", "lta", 
                "lta4h", "ltbp2", "ltbp3", "ltbr", "lto1", "lxn", "ly6d", "ly75", 
                "ly9", "ly96", "lyar", "lyn", "lypd1", "lypd3", "lypd8", "mad1l1", 
                "maea", "maged1", "manf", "mansc1", "map2k6", "map3k5", "map4k5", 
                "mapk9", "mapt", "marco", "masp1", "matn2", "matn3", "mavs", 
                "max", "mb", "mcam", "mcfd2", "mdga1", "mdk", "med18", "megf10", 
                "megf9", "mep1b", "mepe", "mertk", "mesd", "met", "metap1", "metap1d", 
                "metap2", "mfap3", "mfap5", "mfge8", "mgll", "mgmt", "mia", "micb_mica", 
                "mif", "milr1", "mitd1", "mln", "mme", "mmp1", "mmp10", "mmp12", 
                "mmp13", "mmp3", "mmp7", "mmp8", "mmp9", "mnda", "mog", "mphosph8", 
                "mpi", "mpig6b", "mpo", "mrpl46", "msln", "msmb", "msr1", "msra", 
                "mstn", "mtpn", "muc13", "muc16", "mvk", "myo9b", "myoc", "mzb1", 
                "mzt1", "naaa", "nadk", "nampt", "nbl1", "nbn", "ncam1", "ncam2", 
                "ncan", "ncf2", "nck2", "ncln", "ncr1", "ncs1", "ndrg1", "ndufs6", 
                "nectin2", "nectin4", "nefl", "nell1", "nell2", "nfasc", "nfatc1", 
                "nfatc3", "nfkbie", "ngf", "nid1", "nid2", "ninj1", "nme3", "nmnat1", 
                "nomo1", "nos1", "nos3", "notch1", "notch3", "npdc1", "nppb", 
                "nppc", "nptn", "nptx1", "nptxr", "npy", "nrcam", "nrp1", "nrp2", 
                "nrtn", "nsfl1c", "nt5c3a", "nt5e", "ntf3", "ntf4", "ntprobnp", 
                "ntrk2", "ntrk3", "nub1", "nucb2", "nudc", "nudt2", "nudt5", 
                "nxph1", "obp2b", "odam", "ogfr", "ogn", "olr1", "omd", "omg", 
                "optc", "oscar", "osm", "osmr", "oxt", "p4hb", "padi2", "padi4", 
                "paep", "pag1", "pak4", "pam", "pamr1", "pappa", "park7", "parp1", 
                "pbld", "pcdh1", "pcdh17", "pcsk9", "pdcd1", "pdcd1lg2", "pdcd5", 
                "pdcd6", "pdgfa", "pdgfb", "pdgfc", "pdgfra", "pdgfrb", "pdlim7", 
                "pdp1", "pear1", "pebp1", "pecam1", "pfdn2", "pfkfb2", "pgf", 
                "pglyrp1", "phospho1", "pi3", "pigr", "pik3ap1", "pik3ip1", "pilra", 
                "pilrb", "pklr", "pla2g10", "pla2g15", "pla2g1b", "pla2g2a", 
                "pla2g4a", "pla2g7", "plat", "plau", "plaur", "plin1", "plin3", 
                "plpbp", "pltp", "plxdc1", "plxna4", "plxnb2", "plxnb3", "pm20d1", 
                "pmvk", "pnliprp2", "pnpt1", "podxl", "podxl2", "polr2f", "pon2", 
                "pon3", "ppcdc", "ppib", "ppm1a", "ppme1", "ppp1r12a", "ppp1r2", 
                "ppp1r9b", "ppp3r1", "ppy", "pqbp1", "prcp", "prdx1", "prdx3", 
                "prdx5", "prdx6", "preb", "prelp", "prkab1", "prkar1a", "prkcq", 
                "prkra", "prl", "proc", "prok1", "prss2", "prss27", "prss8", 
                "prtfdc1", "prtg", "prtn3", "psg1", "psip1", "psma1", "psmd9", 
                "psme1", "psme2", "psmg3", "pspn", "psrc1", "pten", "ptgds", 
                "pth1r", "ptk7", "ptn", "ptpn1", "ptpn6", "ptprf", "ptprm", "ptprn2", 
                "ptprs", "pts", "ptx3", "pvalb", "pvr", "pxn", "qdpr", "qpct", 
                "rab37", "rab6a", "rab6b", "rabepk", "rabgap1l", "rad23b", "rangap1", 
                "rarres1", "rarres2", "rasa1", "rassf2", "rbks", "rbp2", "rbp5", 
                "rcor1", "reg1a", "reg1b", "reg3a", "reg4", "relt", "ren", "ret", 
                "retn", "rgma", "rgmb", "rgs8", "rhoc", "rilp", "rnase3", "rnaset2", 
                "rnf41", "robo1", "robo2", "ror1", "rp2", "rrm2", "rrm2b", "rspo1", 
                "rspo3", "rtbdn", "rtn4r", "ruvbl1", "rwdd1", "s100a11", "s100a12", 
                "s100a16", "s100a4", "s100p", "samd9l", "scamp3", "scara5", "scarb1", 
                "scarb2", "scarf1", "scarf2", "scg2", "scg3", "scgb1a1", "scgb3a2", 
                "scgn", "scly", "scp2", "scrn1", "sdc1", "sdc4", "sele", "selp", 
                "selplg", "sema3f", "sema4c", "sema4d", "sema7a", "septin9", 
                "serpina11", "serpina12", "serpina9", "serpinb1", "serpinb5", 
                "serpinb6", "serpinb8", "serpinb9", "serpine1", "sestd1", "setmar", 
                "sez6l", "sez6l2", "sf3b4", "sfrp1", "sftpa1", "sftpa2", "sftpd", 
                "sh2b3", "sh2d1a", "shmt1", "siae", "siglec1", "siglec10", "siglec15", 
                "siglec5", "siglec6", "siglec7", "siglec9", "sirpa", "sirpb1", 
                "sirt2", "sirt5", "sit1", "skap1", "skap2", "slamf1", "slamf6", 
                "a1bg", "aamdc", "abca2", "abo", "abraxas2", "acadm", "acadsb", 
                "ace", "ache", "acot13", "acp1", "acrbp", "acrv1", "acsl1", "actn2", 
                "acy3", "acyp1", "adam12", "adam9", "adamts1", "adamts4", "adamtsl2", 
                "adamtsl4", "adamtsl5", "add1", "adgrd1", "adgre1", "adgrf5", 
                "adgrv1", "adh1b", "adipoq", "adra2a", "afap1", "afm", "agbl2", 
                "agt", "slamf7", "slamf8", "slc16a1", "slc27a4", "slc39a14", 
                "slc39a5", "slit2", "slitrk2", "slitrk6", "smad1", "smad5", "smarca2", 
                "smoc1", "smoc2", "smpd1", "smpdl3a", "snap23", "snap29", "sncg", 
                "snx9", "sod1", "sod2", "sorcs2", "sord", "sort1", "sost", "sparc", 
                "sparcl1", "spink1", "spink4", "spink5", "spink6", "spint1", 
                "spint2", "spock1", "spon1", "spon2", "spp1", "spry2", "src", 
                "srp14", "srpk2", "ssb", "ssc4d", "ssc5d", "st3gal1", "st6gal1", 
                "stambp", "stat5b", "stc1", "stc2", "stip1", "stk11", "stk24", 
                "stk4", "stx16", "stx4", "stx6", "stx8", "stxbp3", "sugt1", "sult1a1", 
                "sult2a1", "sumf2", "susd1", "susd2", "tacc3", "tacstd2", "tafa5", 
                "tank", "tarbp2", "tbc1d17", "tbc1d23", "tbc1d5", "tbcb", "tbcc", 
                "tbl1x", "tcl1a", "tcl1b", "tcn2", "tdgf1", "tdrkh", "tek", "tff1", 
                "tff2", "tff3", "tfpi", "tfpi2", "tfrc", "tgfa", "tgfb1", "tgfbi", 
                "tgfbr2", "tgfbr3", "tgm2", "thbd", "thbs2", "thbs4", "thop1", 
                "thpo", "thy1", "tia1", "tie1", "tigar", "timd4", "timp1", "timp3", 
                "timp4", "tinagl1", "tjap1", "tlr3", "tmprss15", "tmprss5", "tmsb10", 
                "tnc", "tnf", "tnfaip8", "tnfrsf10a", "tnfrsf10b", "tnfrsf10c", 
                "tnfrsf11a", "tnfrsf11b", "tnfrsf12a", "tnfrsf13b", "tnfrsf13c", 
                "tnfrsf14", "tnfrsf19", "tnfrsf1a", "tnfrsf1b", "tnfrsf21", "tnfrsf4", 
                "tnfrsf6b", "tnfrsf8", "tnfrsf9", "tnfsf10", "tnfsf11", "tnfsf12", 
                "tnfsf13", "tnfsf13b", "tnfsf14", "tnni3", "tnr", "tnxb", "tp53", 
                "tp53inp1", "tpmt", "tpp1", "tppp3", "tpsab1", "tpt1", "traf2", 
                "trem2", "treml2", "triap1", "trim21", "trim5", "tshb", "tslp", 
                "tspan1", "tst", "txlna", "txndc15", "txndc5", "txnrd1", "tymp", 
                "tyro3", "ubac1", "ulbp2", "umod", "uso1", "usp8", "uxs1", "vamp5", 
                "vash1", "vasn", "vat1", "vcam1", "vcan", "vegfa", "vegfc", "vegfd", 
                "vim", "vmo1", "vnn2", "vps37a", "vps53", "vsig4", "vsir", "vstm1", 
                "vstm2l", "vta1", "vtcn1", "vwa1", "vwc2", "vwf", "wars", "was", 
                "wasf1", "wasf3", "wfdc12", "wfdc2", "wfikkn1", "wfikkn2", "wif1", 
                "wnt9a", "wwp2", "xcl1", "xg", "xpnpep2", "xrcc4", "yes1", "ythdf3", 
                "zbtb16", "zbtb17", "ahnak", "ahnak2", "ahsa1", "ahsg", "aida", 
                "aif1l", "ak2", "akap12", "akr1b10", "akr7l", "akt2", "aldh2", 
                "aldh5a1", "alms1", "alpi", "amdhd2", "amigo1", "amot", "amotl2", 
                "ampd3", "amy1a_amy1b_amy1c", "ank2", "ankmy2", "ankra2", "anp32c", 
                "anxa1", "anxa2", "ap1g2", "ap2b1", "ap3b1", "ap3s2", "apcs", 
                "apoa1", "apoa2", "apoa4", "apob", "apobr", "apoc1", "apod", 
                "apoe", "apof", "apol1", "appl2", "araf", "arf6", "arfip1", "arg2", 
                "arhgap30", "arhgap45", "arhgap5", "arhgef1", "arhgef10", "arhgef5", 
                "arid3a", "arl13b", "arl2bp", "armcx2", "arntl", "art5", "asah1", 
                "asgr2", "aspn", "aspscr1", "asrgl1", "ass1", "atf4", "atg16l1", 
                "atp1b1", "atp1b2", "atp1b3", "atp1b4", "atp2b4", "atp5f1d", 
                "atp6v1g1", "atp6v1g2", "atraid", "atrn", "atxn2", "atxn2l", 
                "atxn3", "azi2", "b2m", "b3gat3", "b3gnt7", "babam1", "bag4", 
                "bap18", "batf", "bcat1", "bcat2", "bche", "bcl2", "bcl2l1", 
                "bcl2l15", "bcl7a", "bcl7b", "bdnf", "becn1", "bex3", "bglap", 
                "bhlhe40", "bhmt2", "blnk", "bloc1s2", "bloc1s3", "bmp10", "bmper", 
                "bnip2", "bnip3l", "bola1", "bola2_bola2b", "bpifa2", "bpifb2", 
                "brap", "brd1", "brd2", "brd3", "brdt", "brme1", "brsk2", "bsnd", 
                "btd", "btla", "btn1a1", "btnl10", "c1galt1c1", "c1qbp", "c1ql2", 
                "c1qtnf5", "c1qtnf6", "c1qtnf9", "c1r", "c1rl", "c1s", "c2orf69", 
                "c3", "c5", "c7", "c7orf50", "c8b", "c9", "c9orf40", "ca7", "ca8", 
                "cabp2", "cacna1c", "cacna1h", "cacnb1", "cacnb3", "cacybp", 
                "cadps", "calcb", "calcoco2", "caly", "camlg", "camsap1", "capn3", 
                "caps", "casc3", "casp4", "casp7", "casp9", "casq2", "cat", "cbln1", 
                "cbs", "cbx2", "ccar2", "ccdc134", "ccdc28a", "ccdc50", "ccer2", 
                "ccnd2", "ccne1", "cd101", "cd164l2", "cd2", "cd226", "cd248", 
                "cd300a", "cd36", "cd3d", "cd3e", "cd3g", "cd5l", "cd7", "cd72", 
                "cd80", "cd82", "cd86", "cda", "cdan1", "cdc123", "cdc25a", "cdc26", 
                "cdc42bpb", "cdh22", "cdh23", "cdh4", "cdk1", "cdk5rap3", "cdkl5", 
                "ceacam16", "ceacam18", "ceacam19", "ceacam20", "ceacam6", "cebpa", 
                "cela2a", "celsr2", "cemip2", "cend1", "cenpf", "cenpj", "cep112", 
                "cep152", "cep170", "cep290", "cep350", "cetn3", "cfb", "cfd", 
                "cfh", "cfhr2", "cfhr4", "cfhr5", "cfi", "cfp", "cgb3_cgb5_cgb8", 
                "cgn", "chad", "chchd10", "chchd6", "chga", "chm", "chmp6", "chp1", 
                "chrm1", "cilp", "cinp", "cirbp", "cit", "ckb", "clasp1", "clec12a", 
                "clec2l", "clec3b", "clec4m", "clgn", "clic5", "clint1", "clns1a", 
                "clstn3", "clu", "cmc1", "cmip", "cngb3", "cnp", "cntf", "cntnap4", 
                "coch", "col15a1", "col24a1", "col28a1", "col2a1", "col3a1", 
                "col4a4", "col5a1", "col9a2", "commd1", "commd9", "copb2", "coq7", 
                "coro6", "cox6b1", "cpa4", "cpb2", "cplx2", "cpox", "cpq", "cptp", 
                "cpxm2", "cr1", "creb3", "crebzf", "creld1", "crisp3", "crtap", 
                "crybb1", "crybb2", "crygd", "crym", "cryzl1", "csde1", "csf1r", 
                "csf2", "csf2rb", "csf3r", "csh1", "csnk1d", "csnk2a1", "cspg4", 
                "cspg5", "csrp3", "cst1", "ctag1a_ctag1b", "ctbs", "cthrc1", 
                "ctla4", "ctnna1", "ctrl", "ctse", "cuzd1", "cwc15", "cyb5a", 
                "cyb5r2", "cyp24a1", "cyth3", "cytl1", "daam1", "dand5", "dapk2", 
                "dbh", "dbn1", "dcc", "dcdc2c", "dclre1c", "dctd", "dcun1d1", 
                "dcun1d2", "dda1", "ddhd2", "ddi2", "ddt", "ddx1", "ddx25", "ddx39a", 
                "ddx4", "ddx53", "defb103a_defb103b", "defb104a_defb104b", "defb116", 
                "defb118", "dennd2b", "denr", "dgcr6", "dgka", "dhodh", "dhps", 
                "dhrs4l2", "dipk1c", "dipk2b", "dlg4", "dlgap5", "dll4", "dmd", 
                "dmp1", "dnaja1", "dnaja4", "dnajb14", "dnajb2", "dnajb6", "dnajc21", 
                "dnajc6", "dnajc9", "dnlz", "dnm1", "dnm3", "dnpep", "doc2b", 
                "dock9", "dok1", "dscam", "dtd1", "dtnb", "dtx2", "dtymk", "dusp13", 
                "dusp29", "dut", "dxo", "dync1h1", "dynlt1", "dynlt3", "echdc3", 
                "echs1", "eci2", "ecm1", "ecscr", "eddm3b", "edem2", "edf1", 
                "edn1", "ednrb", "eef1d", "efcab14", "efcab2", "efhd1", "efnb2", 
                "egflam", "ehbp1", "ehd3", "eif1ax", "eif2ak2", "eif2ak3", "eif2s2", 
                "eif4e", "eif4g3", "eif5", "elac1", "elavl4", "eln", "elob", 
                "endou", "eno3", "enoph1", "enox2", "enpep", "enpp6", "ensa", 
                "entr1", "ep300", "epb41l5", "epgn", "epha4", "epn1", "eppk1", 
                "erc2", "ercc1", "eri1", "ermap", "ern1", "erp29", "ervv_1", 
                "espl1", "esr1", "esyt2", "evi2b", "evi5", "evpl", "exosc10", 
                "extl1", "f10", "f11", "f12", "f13b", "f2", "fabp3", "fam13a", 
                "fam171a2", "fam171b", "fam172a", "fam20a", "fam3d", "farsa", 
                "fbln2", "fbn2", "fcamr", "fcer1a", "fcn1", "fdx1", "fdx2", "fga", 
                "fgd3", "fgf12", "fgf16", "fgf20", "fgf3", "fgf6", "fgf7", "fgf9", 
                "fgfbp2", "fgfbp3", "fgfr4", "fgl1", "fh", "fhip2a", "fkbp14", 
                "fkbpl", "fn1", "fndc1", "fnta", "folh1", "fos", "foxj3", "frmd4b", 
                "frmd7", "fshb", "fstl1", "ftcd", "fuom", "fut1", "fzd10", "fzd8", 
                "gabarap", "gabarapl1", "gabra4", "gad1", "gad2", "gadd45b", 
                "gadd45gip1", "gage2a", "galnt5", "gamt", "gapdh", "gart", "gas2", 
                "gask1a", "gast", "gata3", "gatd3", "gba", "gbp1", "gbp6", "gc", 
                "gcc1", "gchfr", "gclm", "get3", "gfral", "ggact", "ggct", "ghr", 
                "gid8", "gigyf2", "gimap7", "gimap8", "gip", "gipc2", "gipc3", 
                "gipr", "git1", "gja8", "gla", "gli2", "glp1r", "glrx5", "glyr1", 
                "gm2a", "gmfg", "gmpr2", "gnas", "gngt1", "gnpda1", "gnpda2", 
                "golga3", "gorasp2", "got1", "gp1bb", "gp5", "gpd1", "gpha2", 
                "gpi", "gpihbp1", "gpr101", "gpr158", "gpr15l", "gprc5c", "grhpr", 
                "grik2", "grin2b", "grp", "grsf1", "gsn", "gsr", "gstm4", "gstt2b", 
                "gtf2ird1", "gtpbp2", "gucy2c", "guk1", "h2ap", "hadh", "hbz", 
                "hcg22", "hdac8", "hdac9", "hddc2", "hdgfl2", "heg1", "hepacam2", 
                "heph", "hgfac", "hhex", "hif1a", "hip1", "hip1r", "hjv", "hla_a", 
                "hmcn2", "hmgcl", "hmgcs1", "hmmr", "hnf1a", "hnrnpul1", "hpse", 
                "hras", "hrc", "hrg", "hs1bp3", "hs6st2", "hsbp1", "hsd17b14", 
                "hsd17b3", "hsdl2", "hspa2", "htr1a", "htr1b", "id4", "ido1", 
                "ifi30", "ifit1", "ifit3", "ifnar1", "ifnl2", "ifnw1", "ift20", 
                "igbp1", "igdcc3", "igdcc4", "igf2bp3", "igfl4", "ighmbp2", "iglc2", 
                "iglon5", "igsf21", "igsf9", "il12rb2", "il13ra2", "il20rb", 
                "il21r", "il22", "il25", "il2rg", "il3", "il31", "il31ra", "il36a", 
                "il36g", "il9", "immt", "impact", "impg1", "inhbb", "inpp5d", 
                "inpp5j", "insl3", "insl4", "insl5", "insr", "ism2", "ist1", 
                "itga2", "itgal", "itgax", "itgbl1", "itih1", "itih4", "itih5", 
                "itpa", "itpr1", "itprip", "izumo1", "jam3", "jmjd1c", "jpt2", 
                "kazn", "kcnc4", "kcnh2", "kctd5", "kdm3a", "khdc3l", "khk", 
                "kiaa0319", "kiaa1549", "kiaa1549l", "kiaa2013", "kif1c", "kif20b", 
                "kif22", "kir2dl2", "kir2ds4", "kir3dl2", "kirrel1", "klf4", 
                "klhl41", "klk15", "klk3", "klk7", "klkb1", "klrc1", "klrf1", 
                "klrk1", "krt17", "krt6c", "krt8", "l3hypdh", "lacrt", "lama1", 
                "lamb1", "lamp1", "lamtor5", "larp1", "lats1", "lcat", "lcn15", 
                "lcp1", "ldlrap1", "lect2", "leg1", "lelp1", "leo1", "letm1", 
                "lgals3bp", "lilra3", "lilra4", "lilra6", "lipf", "lmnb1", "lmnb2", 
                "lmod1", "lmod2", "lonp1", "lpa", "lpp", "lrch4", "lrfn2", "lrg1", 
                "lrig3", "lrp2", "lrp2bp", "lrrc37a2", "lrrc38", "lrrc59", "lrrfip1", 
                "lrtm1", "lrtm2", "lsm8", "ltb", "luzp2", "lypla2", "lysmd3", 
                "lyve1", "lyzl2", "lztfl1", "m6pr", "mag", "magea3", "mamdc2", 
                "mamdc4", "man1a2", "man2b2", "maneal", "mansc4", "map1lc3a", 
                "map1lc3b2", "map2", "map2k1", "mapk13", "mapkapk2", "mapre3", 
                "mars1", "mbl2", "mcee", "mcemp1", "mcts1", "mdh1", "mdm1", "mecr", 
                "med21", "megf11", "meltf", "ment", "mep1a", "mfap3l", "mfap4", 
                "micall2", "mindy1", "mink1", "mki67", "mllt1", "mmp15", "mmut", 
                "mn1", "mnat1", "mocs2", "morc3", "morf4l1", "morf4l2", "morn4", 
                "mprip", "mrc1", "mri1", "mrpl24", "mrpl28", "mrpl52", "mrpl58", 
                "mrps16", "mslnl", "mst1", "mtdh", "mthfd2", "mthfsd", "mtif3", 
                "mtr", "mtss1", "mtss2", "mtus1", "muc2", "mucl3", "mxra8", "mybpc1", 
                "mybpc2", "mycbp2", "mydgf", "myh4", "myh7b", "myh9", "myl1", 
                "myl3", "myl4", "myl6b", "mylpf", "myo6", "myom1", "myom2", "myom3", 
                "naa10", "naa80", "nacc1", "naga", "nagk", "nagpa", "nap1l4", 
                "naprt", "nars1", "ncr3lg1", "ndst1", "ndufa5", "ndufb7", "neb", 
                "necap2", "nectin1", "nedd4l", "nedd9", "nek7", "nenf", "neo1", 
                "nexn", "nfat5", "nfe2", "nfic", "nfkb1", "nfkb2", "nfu1", "nfx1", 
                "nfya", "ngfr", "ngrn", "nhlrc3", "nit1", "nit2", "nlgn1", "nlgn2", 
                "nme1", "nmi", "nmrk2", "nmt1", "nop56", "nos2", "notch2", "npc2", 
                "nphs1", "nphs2", "npl", "npr1", "nptx2", "nrgn", "nrn1", "nrxn3", 
                "nt5c", "nt5c1a", "nubp1", "nudt10", "nudt15", "nudt16", "numb", 
                "nup50", "nxpe4", "nxph3", "ocln", "ofd1", "oga", "ogt", "olfm4", 
                "omp", "ophn1", "oplah", "orm1", "osbpl2", "ostn", "otoa", "otud6b", 
                "otud7b", "oxct1", "pacs2", "pafah1b3", "pafah2", "pagr1", "paip2b", 
                "palld", "palm", "palm2", "palm3", "pard3", "paxx", "pbk", "pbxip1", 
                "pcare", "pcbd1", "pcbp2", "pcdh12", "pcdh7", "pcdh9", "pcdhb15", 
                "pcna", "pcsk7", "pcyt2", "pdap1", "pdcl2", "pde1c", "pde4d", 
                "pde5a", "pdia2", "pdia3", "pdia4", "pdia5", "pdlim5", "pdrg1", 
                "pdxdc1", "pdzd2", "pdzk1", "pecr", "penk", "pepd", "per3", "pf4", 
                "pfdn4", "pfdn6", "pga4", "pgd", "pglyrp2", "pglyrp4", "pgm2", 
                "pgr", "phactr2", "phldb1", "phldb2", "phykpl", "pi16", "pibf1", 
                "pikfyve", "pinlyp", "pithd1", "pkd1", "pkd2", "pkn3", "plb1", 
                "plcb1", "plcb2", "plekho1", "plg", "plscr3", "plxdc2", "pmch", 
                "pmm2", "pms1", "pnlip", "pnliprp1", "pnma1", "pnma2", "pof1b", 
                "polr2a", "pomc", "pon1", "postn", "ppbp", "ppie", "ppif", "ppl", 
                "ppm1b", "ppm1f", "ppp1cc", "ppp1r12b", "ppp1r14a", "ppp1r14d", 
                "ppp2r5a", "ppt1", "prame", "prap1", "prc1", "prdx2", "prg2", 
                "prg3", "prkag3", "prkar2a", "prkd2", "prkg1", "prnd", "procr", 
                "pros1", "prr4", "prr5", "prrt3", "prss22", "prss53", "prune2", 
                "psap", "psapl1", "psca", "psmc3", "psmd1", "psmd5", "psmg4", 
                "pstpip2", "ptges2", "ptgr1", "pth", "ptp4a3", "ptpn9", "ptprb", 
                "ptprc", "ptprh", "ptprk", "ptprr", "ptprz1", "ptrhd1", "pttg1", 
                "pxdnl", "pydc1", "pyy", "pzp", "qsox1", "rab10", "rab11fip3", 
                "rab27b", "rab2b", "rab33a", "rab39b", "rab3gap1", "rab44", "rabep1", 
                "rac3", "rad51", "ralb", "raly", "ranbp1", "ranbp2", "rap1a", 
                "rapgef2", "rasgrf1", "rbfox3", "rbm17", "rbm19", "rbm25", "rbp1", 
                "rbp7", "rbpms", "rbpms2", "rcc1", "reck", "reep4", "reg3g", 
                "relb", "reps1", "rest", "rexo2", "rfc4", "rgcc", "rgl2", "rgs10", 
                "rictor", "rida", "rilpl2", "ripk4", "rln1", "rln2", "rnase1", 
                "rnase10", "rnase4", "rnase6", "rnaseh2a", "rnf149", "rnf168", 
                "rnf31", "rnf4", "rnf43", "rnf5", "robo4", "rpa2", "rpe", "rpgr", 
                "rpl14", "rps10", "rras", "rrp15", "rtkn2", "rtn4ip1", "ryr1", 
                "s100a13", "s100a14", "s100a3", "s100g", "saa4", "safb2", "sag", 
                "sap18", "sarg", "sart1", "sat1", "sat2", "satb1", "sbsn", "scgb2a2", 
                "scgb3a1", "scin", "scn2a", "scn2b", "scn3a", "scn3b", "scn4b", 
                "scpep1", "scrg1", "scrib", "sct", "sdccag8", "sdhb", "sdk2", 
                "sec31a", "sel1l", "selenop", "sell", "sema3g", "sema6c", "septin3", 
                "septin7", "septin8", "serpina1", "serpina3", "serpina4", "serpina5", 
                "serpina6", "serpina7", "serpinc1", "serpind1", "serpine2", "serpinf1", 
                "serpinf2", "serping1", "serpinh1", "serpini1", "serpini2", "sez6", 
                "sfrp4", "sgsh", "sh3bgrl2", "sh3bp1", "sh3gl3", "sh3glb2", "shbg", 
                "shc1", "shd", "shh", "shisa5", "shpk", "siglec8", "sil1", "sirt1", 
                "skiv2l", "sla2", "slc12a2", "slc13a1", "slc1a4", "slc28a1", 
                "slc34a3", "slc44a4", "slc4a1", "slc51b", "slc9a3r1", "slc9a3r2", 
                "slirp", "slitrk1", "slk", "slmap", "slurp1", "smad2", "smad3", 
                "smc3", "smndc1", "smpd3", "smpdl3b", "sms", "smtn", "snap25", 
                "snapin", "snca", "sned1", "snrpb2", "snu13", "snx15", "snx18", 
                "snx2", "snx5", "sod3", "sorbs1", "sowaha", "sox2", "sox9", "spaca5_spaca5b", 
                "spag1", "spart", "spesp1", "spink2", "spink8", "spint3", "spred2", 
                "spring1", "sprr1b", "sprr3", "sptbn2", "sptlc1", "srpx", "ssbp1", 
                "ssh3", "ssna1", "st13", "st8sia1", "stab2", "stam", "stat2", 
                "stau1", "steap4", "stoml2", "stx1b", "stx3", "stx5", "stx7", 
                "stxbp1", "sugp1", "sumf1", "suox", "susd4", "susd5", "sv2a", 
                "swap70", "syap1", "syngap1", "syt1", "sytl4", "tab2", "tada3", 
                "tagln3", "taldo1", "tap1", "tarm1", "tars1", "tax1bp1", "tbca", 
                "tbr1", "tcn1", "tcof1", "tcp11", "tctn3", "tdo2", "tdp1", "tef", 
                "terf1", "tet2", "tex101", "tex33", "tf", "tfap2a", "tg", "tgfb2", 
                "tgfbr1", "tgoln2", "thap12", "thrap3", "thsd1", "thtpa", "tigit", 
                "timm10", "timm8a", "timp2", "tjp3", "tk1", "tlr1", "tlr2", "tlr4", 
                "tmco5a", "tmed1", "tmed10", "tmed4", "tmed8", "tmem106a", "tmem132a", 
                "tmem25", "tmod4", "tmprss11b", "tmprss11d", "tnfaip2", "tnfaip8l2", 
                "tnfrsf17", "tnfsf8", "tnip1", "tnn", "tnpo1", "tomm20", "top1", 
                "top1mt", "top2b", "tor1aip1", "tp53bp1", "tp53i3", "tp73", "tpbgl", 
                "tpd52l2", "tpk1", "tpm3", "tppp2", "tpr", "tprkb", "tpsd1", 
                "tpsg1", "traf3", "traf3ip2", "trdmt1", "treh", "treml1", "trim24", 
                "trim25", "trim26", "trim40", "trim58", "trpv3", "tsc1", "tsc22d1", 
                "tsnax", "tspan15", "tspan7", "tspan8", "tspyl1", "ttf2", "ttn", 
                "ttr", "tubb3", "twf2", "txk", "txn", "txndc9", "txnl1", "tyrp1", 
                "ube2b", "ube2l6", "ube2z", "ubqln3", "ubxn1", "ufd1", "ugdh", 
                "uhrf2", "unc5d", "unc79", "ung", "upb1", "upk3a", "upk3bl1", 
                "urod", "uros", "usp25", "usp28", "usp47", "vamp8", "vasp", "vav3", 
                "vcpkmt", "vegfb", "vgf", "vipr1", "vit", "vnn1", "vps28", "vps4b", 
                "vsig10", "vsig10l", "vsig2", "vsnl1", "vstm2b", "vti1a", "vwa5a", 
                "vwc2l", "washc3", "wasl", "wdr46", "wfdc1", "xiap", "yap1", 
                "yars1", "yju2", "yod1", "ywhaq", "yy1", "zbp1", "zcchc8", "zfyve19", 
                "zhx2", "znf174", "znf75d", "znf830", "znrd2", "znrf4", "zp3", 
                "zp4", "zpr1") 
# New column names for restructured data
new_names <- c("n_eid","status","timef","timeto",covs_P,traits,"failcode")

# Function to prepare data for multistate analysis
run_multistate_analysis <- function(data1, new_names, covs_P, traits) {
  prepare_data <- function(data1, new_names) {
    datatrans_list <- list()
    # Create 9 transition patterns
    for(i in 1:9) {
      time_cols <- c(1, (2+(i-1)*3):(4+(i-1)*3))
      common_cols <- c(29:46, 55:2974)
      
      datatrans_list[[i]] <- data1[, c(time_cols, common_cols)]
      datatrans_list[[i]]$transpattern <- i
      
      colnames(datatrans_list[[i]]) <- new_names
    }
    # Combine all transitions
    dataall1 <- do.call(rbind, datatrans_list)
    dataall1$trans <- dataall1$failcode
    # Convert categorical variables to factors
    dataall1 <- dataall1 %>%
      mutate_at(vars(Sex, Ethnicity, Edu, Employed, Smoke, Drink, T2dhis, Cvdhis, 
                     hiscancer, cancer_treat, CMD_treat, cancer_screen), as.factor)
    return(dataall1)
  }
  # Function to run multi-state model for each protein
  run_cox_analysis <- function(msebmt_allname, traits) {
    results <- data.frame(Proteins = character(),
                          HR = character(),
                          P_value = numeric(),
                          stringsAsFactors = FALSE)
    # Create covariate adjustment string for model formula
    covariate_str <- paste(
      "Age.1 + Age.2 + Age.3 + Age.4 + Age.5 + Age.6 + Age.7 + Age.8 + Age.9 + ",
      "Sex.1 + Sex.2 + Sex.3 + Sex.4 + Sex.5 + Sex.6 + Sex.7 + Sex.8 + Sex.9 +",
      "Age.1 + Age.2 + Age.3 + Age.4 + Age.5 + Age.6 + Age.7 + Age.8 + Age.9 + ",
      "Sex.1 + Sex.2 + Sex.3 + Sex.4 + Sex.5 + Sex.6 + Sex.7 + Sex.8 + Sex.9+" ,
      "Ethnicity.1 + Ethnicity.2 + Ethnicity.3 + Ethnicity.4 + Ethnicity.5 + Ethnicity.6 + Ethnicity.7 + Ethnicity.8 + Ethnicity.9 +",
      "Edu1.1 + Edu1.2 + Edu1.3 + Edu1.4 + Edu1.5 + Edu1.6 + Edu1.7 + Edu1.8 + Edu1.9 + ",
      "Edu2.1 + Edu2.2 + Edu2.3 + Edu2.4 + Edu2.5 + Edu2.6 + Edu2.7 + Edu2.8 + Edu2.9 +",
      "Edu3.1 + Edu3.2 + Edu3.3 + Edu3.4 + Edu3.5 + Edu3.6 + Edu3.7 + Edu3.8 + Edu3.9 +", 
      "Employed.1 + Employed.2 + Employed.3 + Employed.4 + Employed.5 + Employed.6 + Employed.7 + Employed.8 + Employed.9 +",
      "Tdi.1 + Tdi.2 + Tdi.3 + Tdi.4 + Tdi.5 + Tdi.6 + Tdi.7 + Tdi.8 + Tdi.9+",
      "Smoke1.1 + Smoke1.2 + Smoke1.3 + Smoke1.4 + Smoke1.5 + Smoke1.6 + Smoke1.7 + Smoke1.8 + Smoke1.9 +",
      "Smoke2.1 + Smoke2.2 + Smoke2.3 + Smoke2.4 + Smoke2.5 + Smoke2.6 + Smoke2.7 + Smoke2.8 + Smoke2.9 +",
      "Drink1.1 + Drink1.2 + Drink1.3 + Drink1.4 + Drink1.5 + Drink1.6 + Drink1.7 + Drink1.8 + Drink1.9 +", 
      "Drink2.1 + Drink2.2 + Drink2.3 + Drink2.4 + Drink2.5 + Drink2.6 + Drink2.7 + Drink2.8 + Drink2.9 + ",
      "Mets.1 + Mets.2 + Mets.3 + Mets.4 + Mets.5 + Mets.6 + Mets.7 + Mets.8 + Mets.9 +",
      "Dietscore.1 + Dietscore.2 + Dietscore.3 + Dietscore.4 + Dietscore.5 + Dietscore.6 + Dietscore.7 + Dietscore.8 + Dietscore.9 +",
      "Sleeptime.1 + Sleeptime.2 + Sleeptime.3 + Sleeptime.4 + Sleeptime.5 + Sleeptime.6 + Sleeptime.7 + Sleeptime.8 + Sleeptime.9 +",
      "BMI.1 + BMI.2 + BMI.3 + BMI.4 + BMI.5 + BMI.6 + BMI.7 + BMI.8 + BMI.9 +",
      "T2dhis.1 + T2dhis.2 + T2dhis.3 + T2dhis.4 + T2dhis.5  + T2dhis.6 + T2dhis.7 + T2dhis.8 + T2dhis.9 +", 
      "Cvdhis.1 + Cvdhis.2 + Cvdhis.3 + Cvdhis.4 + Cvdhis.5 + Cvdhis.6 + Cvdhis.7 + Cvdhis.8 + Cvdhis.9 +",
      "hiscancer.1 + hiscancer.2+ hiscancer.3+ hiscancer.4+ hiscancer.5+ hiscancer.6+ hiscancer.7+ hiscancer.8+ hiscancer.9+", 
      "cancer_treat.1 + cancer_treat.2+ cancer_treat.3+ cancer_treat.4+ cancer_treat.5+ cancer_treat.6+ cancer_treat.7+ cancer_treat.8+ cancer_treat.9+",
      "CMD_treat.1 + CMD_treat.2+ CMD_treat.3+ CMD_treat.4+ CMD_treat.5+ CMD_treat.6+ CMD_treat.7+ CMD_treat.8+ CMD_treat.9+",
      "cancer_screen.1 + cancer_screen.2+ cancer_screen.3+ cancer_screen.4+ cancer_screen.5+ cancer_screen.6+ cancer_screen.7+ cancer_screen.8+ cancer_screen.9+strata(trans)", 
      sep = ""
    )
    # Analyze each protein
    for (met in traits) {
      tryCatch({
        # Create protein terms for all transitions
        met_str <- paste(met, ".1 +", met, ".2 +", met, ".3 +", met, ".4 +", 
                         met, ".5 +", met, ".6 +", met, ".7 +", met, ".8 +", met, ".9", sep = "")
        formula_str <- paste("Surv(timef, timeto, status) ~", met_str, "+", covariate_str)
        model_formula <- as.formula(formula_str)
        # Extract results for each transition
        model <- coxph(model_formula, data = msebmt_allname, method = 'breslow')
        for (i in 1:9) {
          met_var <- paste(met, ".", i, sep = "")
          if (met_var %in% rownames(summary(model)$coefficients)) {
            esti <- summary(model)$conf.int[met_var, "exp(coef)"]
            lower <- summary(model)$conf.int[met_var, "lower .95"]
            upper <- summary(model)$conf.int[met_var, "upper .95"]
            HR <- sprintf("%.2f (%.2f to %.2f)", esti, lower, upper)
            P_value <- format(summary(model)$coefficients[met_var, "Pr(>|z|)"], 
                              scientific = TRUE, digits = 2)
            
            results <- rbind(results, data.frame(Proteins = met_var, 
                                                 HR = HR, 
                                                 P_value = P_value, 
                                                 stringsAsFactors = FALSE))
          }
        }
      }, error = function(e) {
        message(paste("Error in processing", met, ":", e$message))
      })
    }
    # Adjust for multiple testing
    results$Fdr_pvalue <- format(p.adjust(as.numeric(results$P_value), 
                                          method = "fdr"), 
                                 scientific = TRUE, digits = 2)
    results$Bonferroni_pvalue <- format(p.adjust(as.numeric(results$P_value), 
                                                 method = "bonferroni"), 
                                        scientific = TRUE, digits = 2)
    
    return(results)
  }
  dataall1 <- prepare_data(data1, new_names)
  covs_allname <- c(covs_P, traits)
  msebmt_allname <- expand.covs(dataall1, covs_allname, append = TRUE, longnames = FALSE)
  results <- run_cox_analysis(msebmt_allname, traits)
  return(results)
}

# Create formatted Excel output
create_excel_output <- function(train_results, test_results, output_file) {
  wb <- createWorkbook()
  process_and_write_results <- function(results, sheet_name) {
    addWorksheet(wb, sheet_name)
    transition_cols <- list(
      Transition1 = 1,    
      Transition2 = 5,   
      Transition3 = 9,   
      Transition4 = 13,   
      Transition5 = 17,   
      Transition6 = 21,   
      Transition7 = 25,  
      Transition8 = 29,   
      Transition9 = 33    
    )
    
    for(i in 1:9) {
      transition_data <- results[grep(paste0("\\.", i, "$"), results$Proteins), ]
      transition_df <- data.frame(
        HR = transition_data$HR,
        P_value = transition_data$P_value,
        Fdr_pvalue = transition_data$Fdr_pvalue,
        Bonferroni_pvalue = transition_data$Bonferroni_pvalue
      )
      rownames(transition_df) <- gsub(paste0("\\.", i, "$"), "", transition_data$Proteins)
      
      start_col <- transition_cols[[paste0("Transition", i)]]
      writeData(wb, sheet_name, paste0("Transition", i), 
                startCol = start_col, startRow = 1)
      writeData(wb, sheet_name, 
                c("HR", "P_value", "Fdr_pvalue", "Bonferroni_pvalue"), 
                startCol = start_col, startRow = 2, colNames = FALSE)
      writeData(wb, sheet_name, transition_df, 
                startCol = start_col, startRow = 3, rowNames = TRUE, colNames = FALSE)
    }
    setColWidths(wb, sheet_name, cols = 1:40, widths = 15)
    style <- createStyle(
      halign = "center",
      borderStyle = "thin",
      borderColour = "black"
    )
    headerStyle <- createStyle(
      halign = "center",
      textDecoration = "bold",
      borderStyle = "thin",
      borderColour = "black"
    )
    addStyle(wb, sheet_name, style, rows = 1:100, cols = 1:40, gridExpand = TRUE)
    addStyle(wb, sheet_name, headerStyle, rows = 1:2, cols = 1:40, gridExpand = TRUE)
  }
  
  process_and_write_results(train_results, "Training Set")
  process_and_write_results(test_results, "Test Set")
  saveWorkbook(wb, output_file, overwrite = TRUE)
}

# Complete analysis pipeline
run_complete_analysis(train_data_birth, test_data_birth, new_names, covs_P, traits, 
                      output_file = "Protein_birth_results.xlsx")

### Genomics and metabolomics data follow the same analysis pipeline as described above