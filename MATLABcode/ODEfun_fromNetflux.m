function dydt=ODEfun(t,y,params) 
% Assign names for parameters 
[rpar,tau,ymax,speciesNames]=params{:}; 
aAR = 1; 
AC = 2; 
Akt = 3; 
aMHC = 4; 
AngII = 5; 
ANP = 6; 
ANPi = 7; 
AT1R = 8; 
ATF2 = 9; 
BAR = 10; 
bMHC = 11; 
BNP = 12; 
BNPi = 13; 
Calcium = 14; 
CaM = 15; 
CaMK = 16; 
cAMP = 17; 
CaN = 18; 
CellArea = 19; 
cFos = 20; 
cGMP = 21; 
cJun = 22; 
CREB = 23; 
CT1 = 24; 
DAG = 25; 
EGF = 26; 
EGFR = 27; 
eIF2B = 28; 
eIF4E = 29; 
ELK1 = 30; 
ERBB = 31; 
ERK12 = 32; 
ERK5 = 33; 
ET1 = 34; 
ET1R = 35; 
FAK = 36; 
FGF = 37; 
FGFR = 38; 
foxo = 39; 
Gaq11 = 40; 
GATA4 = 41; 
GBG = 42; 
GCA = 43; 
gp130LIFR = 44; 
Gsa = 45; 
GSK3B = 46; 
HDAC = 47; 
IGF1 = 48; 
IGF1R = 49; 
IkB = 50; 
IKK = 51; 
IL6 = 52; 
IL6R = 53; 
Integrins = 54; 
IP3 = 55; 
ISO = 56; 
JAK = 57; 
JNK = 58; 
LIF = 59; 
MAP3K11 = 60; 
MAP3K23 = 61; 
MAP3K4 = 62; 
MAPKAPK = 63; 
MEF2 = 64; 
MEK12 = 65; 
MEK36 = 66; 
MEK4 = 67; 
MEK5 = 68; 
MEK7 = 69; 
MEKK1 = 70; 
MSK1 = 71; 
mTor = 72; 
NE = 73; 
NFAT = 74; 
NFkB = 75; 
NIK = 76; 
NOS = 77; 
NRG1 = 78; 
p38 = 79; 
p70s6k = 80; 
PDK1 = 81; 
PE = 82; 
PI3K = 83; 
PKA = 84; 
PKC = 85; 
PKD = 86; 
PKG1 = 87; 
PLCB = 88; 
PLCG = 89; 
Rac1 = 90; 
Raf1 = 91; 
Raf1A = 92; 
Ras = 93; 
RhoA = 94; 
sACT = 95; 
SERCA = 96; 
sGC = 97; 
SHP2 = 98; 
SRF = 99; 
STAT = 100; 
Stretch = 101; 
TAK1 = 102; 
TGFB = 103; 
TGFR = 104; 
TNFa = 105; 
TNFR = 106; 
PGR = 107; 
NR3C1 = 108; 
CDK2 = 109; 
CEBPB = 110; 
DUSP1 = 111; 
HDAC1 = 112; 
LCK = 113; 
lig = 114; 
dydt = zeros(114,1); 
dydt(aAR) = (OR(act(y(NE),rpar(:,142)),act(y(PE),rpar(:,159)))*ymax(aAR) - y(aAR))/tau(aAR); 
dydt(AC) = (act(y(Gsa),rpar(:,98))*ymax(AC) - y(AC))/tau(AC); 
dydt(Akt) = (act(y(PDK1),rpar(:,158))*ymax(Akt) - y(Akt))/tau(Akt); 
dydt(aMHC) = (AND(rpar(:,21),inhib(y(cFos),rpar(:,21)),inhib(y(cJun),rpar(:,21)),inhib(y(MEF2),rpar(:,21)),inhib(y(NFAT),rpar(:,21)))*ymax(aMHC) - y(aMHC))/tau(aMHC); 
dydt(AngII) = (rpar(1,1)*ymax(AngII) - y(AngII))/tau(AngII); 
dydt(ANP) = (OR(act(y(ATF2),rpar(:,41)),OR(act(y(cFos),rpar(:,53)),OR(act(y(cJun),rpar(:,58)),OR(act(y(CREB),rpar(:,62)),OR(act(y(MEF2),rpar(:,125)),OR(AND(rpar(:,144),act(y(GATA4),rpar(:,144)),act(y(NFAT),rpar(:,144))),act(y(STAT),rpar(:,183))))))))*ymax(ANP) - y(ANP))/tau(ANP); 
dydt(ANPi) = (rpar(1,2)*ymax(ANPi) - y(ANPi))/tau(ANPi); 
dydt(AT1R) = (act(y(AngII),rpar(:,37))*ymax(AT1R) - y(AT1R))/tau(AT1R); 
dydt(ATF2) = (OR(act(y(JNK),rpar(:,112)),act(y(p38),rpar(:,151)))*ymax(ATF2) - y(ATF2))/tau(ATF2); 
dydt(BAR) = (OR(act(y(ISO),rpar(:,108)),act(y(NE),rpar(:,143)))*ymax(BAR) - y(BAR))/tau(BAR); 
dydt(bMHC) = (OR(act(y(ATF2),rpar(:,42)),OR(act(y(cFos),rpar(:,54)),OR(act(y(cJun),rpar(:,59)),OR(act(y(GATA4),rpar(:,90)),OR(act(y(MEF2),rpar(:,126)),OR(act(y(NFAT),rpar(:,146)),act(y(STAT),rpar(:,184))))))))*ymax(bMHC) - y(bMHC))/tau(bMHC); 
dydt(BNP) = (OR(act(y(ATF2),rpar(:,43)),OR(act(y(cFos),rpar(:,55)),OR(act(y(cJun),rpar(:,60)),OR(act(y(ELK1),rpar(:,70)),OR(act(y(MEF2),rpar(:,127)),AND(rpar(:,145),act(y(GATA4),rpar(:,145)),act(y(NFAT),rpar(:,145))))))))*ymax(BNP) - y(BNP))/tau(BNP); 
dydt(BNPi) = (rpar(1,3)*ymax(BNPi) - y(BNPi))/tau(BNPi); 
dydt(Calcium) = (OR(act(y(IP3),rpar(:,107)),OR(act(y(PKA),rpar(:,161)),act(y(CaMK),rpar(:,192))))*ymax(Calcium) - y(Calcium))/tau(Calcium); 
dydt(CaM) = (act(y(Calcium),rpar(:,47))*ymax(CaM) - y(CaM))/tau(CaM); 
dydt(CaMK) = (act(y(CaM),rpar(:,48))*ymax(CaMK) - y(CaMK))/tau(CaMK); 
dydt(cAMP) = (act(y(AC),rpar(:,33))*ymax(cAMP) - y(cAMP))/tau(cAMP); 
dydt(CaN) = (act(y(CaM),rpar(:,49))*ymax(CaN) - y(CaN))/tau(CaN); 
dydt(CellArea) = (OR(inhib(y(foxo),rpar(:,23)),OR(act(y(ATF2),rpar(:,44)),OR(act(y(cJun),rpar(:,61)),OR(act(y(CREB),rpar(:,63)),OR(act(y(GATA4),rpar(:,91)),act(y(MEF2),rpar(:,128)))))))*ymax(CellArea) - y(CellArea))/tau(CellArea); 
dydt(cFos) = (act(y(ERK12),rpar(:,75))*ymax(cFos) - y(cFos))/tau(cFos); 
dydt(cGMP) = (OR(act(y(GCA),rpar(:,95)),act(y(sGC),rpar(:,181)))*ymax(cGMP) - y(cGMP))/tau(cGMP); 
dydt(cJun) = (OR(act(y(ERK12),rpar(:,76)),act(y(JNK),rpar(:,113)))*ymax(cJun) - y(cJun))/tau(cJun); 
dydt(CREB) = (OR(inhib(y(GSK3B),rpar(:,25)),OR(act(y(MAPKAPK),rpar(:,124)),OR(act(y(MSK1),rpar(:,139)),OR(act(y(PKA),rpar(:,162)),act(y(CEBPB),rpar(:,203))))))*ymax(CREB) - y(CREB))/tau(CREB); 
dydt(CT1) = (rpar(1,4)*ymax(CT1) - y(CT1))/tau(CT1); 
dydt(DAG) = (OR(act(y(PLCB),rpar(:,167)),act(y(PLCG),rpar(:,169)))*ymax(DAG) - y(DAG))/tau(DAG); 
dydt(EGF) = (rpar(1,5)*ymax(EGF) - y(EGF))/tau(EGF); 
dydt(EGFR) = (act(y(EGF),rpar(:,66))*ymax(EGFR) - y(EGFR))/tau(EGFR); 
dydt(eIF2B) = (inhib(y(GSK3B),rpar(:,26))*ymax(eIF2B) - y(eIF2B))/tau(eIF2B); 
dydt(eIF4E) = (act(y(mTor),rpar(:,140))*ymax(eIF4E) - y(eIF4E))/tau(eIF4E); 
dydt(ELK1) = (OR(act(y(ERK12),rpar(:,77)),OR(act(y(JNK),rpar(:,114)),act(y(p38),rpar(:,152))))*ymax(ELK1) - y(ELK1))/tau(ELK1); 
dydt(ERBB) = (act(y(NRG1),rpar(:,150))*ymax(ERBB) - y(ERBB))/tau(ERBB); 
dydt(ERK12) = (OR(act(y(MEK12),rpar(:,130)),OR(act(y(LCK),rpar(:,200)),act(y(DUSP1),rpar(:,205))))*ymax(ERK12) - y(ERK12))/tau(ERK12); 
dydt(ERK5) = (act(y(MEK5),rpar(:,134))*ymax(ERK5) - y(ERK5))/tau(ERK5); 
dydt(ET1) = (rpar(1,6)*ymax(ET1) - y(ET1))/tau(ET1); 
dydt(ET1R) = (act(y(ET1),rpar(:,82))*ymax(ET1R) - y(ET1R))/tau(ET1R); 
dydt(FAK) = (act(y(Integrins),rpar(:,106))*ymax(FAK) - y(FAK))/tau(FAK); 
dydt(FGF) = (rpar(1,7)*ymax(FGF) - y(FGF))/tau(FGF); 
dydt(FGFR) = (act(y(FGF),rpar(:,85))*ymax(FGFR) - y(FGFR))/tau(FGFR); 
dydt(foxo) = (OR(AND(rpar(:,18),inhib(y(Akt),rpar(:,18)),inhib(y(CDK2),rpar(:,18))),OR(act(y(HDAC1),rpar(:,198)),act(y(CDK2),rpar(:,204))))*ymax(foxo) - y(foxo))/tau(foxo); 
dydt(Gaq11) = (OR(act(y(aAR),rpar(:,32)),OR(act(y(AT1R),rpar(:,39)),act(y(ET1R),rpar(:,83))))*ymax(Gaq11) - y(Gaq11))/tau(Gaq11); 
dydt(GATA4) = (OR(inhib(y(GSK3B),rpar(:,27)),OR(act(y(ERK12),rpar(:,78)),act(y(p38),rpar(:,153))))*ymax(GATA4) - y(GATA4))/tau(GATA4); 
dydt(GBG) = (OR(act(y(Gaq11),rpar(:,87)),act(y(Gsa),rpar(:,99)))*ymax(GBG) - y(GBG))/tau(GBG); 
dydt(GCA) = (OR(act(y(ANPi),rpar(:,38)),act(y(BNPi),rpar(:,46)))*ymax(GCA) - y(GCA))/tau(GCA); 
dydt(gp130LIFR) = (OR(act(y(CT1),rpar(:,64)),act(y(LIF),rpar(:,115)))*ymax(gp130LIFR) - y(gp130LIFR))/tau(gp130LIFR); 
dydt(Gsa) = (act(y(BAR),rpar(:,45))*ymax(Gsa) - y(Gsa))/tau(Gsa); 
dydt(GSK3B) = (inhib(y(Akt),rpar(:,19))*ymax(GSK3B) - y(GSK3B))/tau(GSK3B); 
dydt(HDAC) = (AND(rpar(:,20),inhib(y(CaMK),rpar(:,20)),inhib(y(PKC),rpar(:,20)),inhib(y(PKD),rpar(:,20)))*ymax(HDAC) - y(HDAC))/tau(HDAC); 
dydt(IGF1) = (rpar(1,8)*ymax(IGF1) - y(IGF1))/tau(IGF1); 
dydt(IGF1R) = (act(y(IGF1),rpar(:,100))*ymax(IGF1R) - y(IGF1R))/tau(IGF1R); 
dydt(IkB) = (inhib(y(IKK),rpar(:,30))*ymax(IkB) - y(IkB))/tau(IkB); 
dydt(IKK) = (OR(act(y(Akt),rpar(:,34)),OR(act(y(NIK),rpar(:,148)),act(y(p38),rpar(:,154))))*ymax(IKK) - y(IKK))/tau(IKK); 
dydt(IL6) = (rpar(1,9)*ymax(IL6) - y(IL6))/tau(IL6); 
dydt(IL6R) = (act(y(IL6),rpar(:,104))*ymax(IL6R) - y(IL6R))/tau(IL6R); 
dydt(Integrins) = (act(y(Stretch),rpar(:,185))*ymax(Integrins) - y(Integrins))/tau(Integrins); 
dydt(IP3) = (OR(act(y(PLCB),rpar(:,168)),act(y(PLCG),rpar(:,170)))*ymax(IP3) - y(IP3))/tau(IP3); 
dydt(ISO) = (rpar(1,10)*ymax(ISO) - y(ISO))/tau(ISO); 
dydt(JAK) = (OR(act(y(AT1R),rpar(:,40)),OR(act(y(gp130LIFR),rpar(:,96)),act(y(IL6R),rpar(:,105))))*ymax(JAK) - y(JAK))/tau(JAK); 
dydt(JNK) = (OR(act(y(MEK4),rpar(:,132)),OR(act(y(MEK7),rpar(:,135)),act(y(DUSP1),rpar(:,207))))*ymax(JNK) - y(JNK))/tau(JNK); 
dydt(LIF) = (rpar(1,11)*ymax(LIF) - y(LIF))/tau(LIF); 
dydt(MAP3K11) = (act(y(Rac1),rpar(:,171))*ymax(MAP3K11) - y(MAP3K11))/tau(MAP3K11); 
dydt(MAP3K23) = (act(y(Ras),rpar(:,175))*ymax(MAP3K23) - y(MAP3K23))/tau(MAP3K23); 
dydt(MAP3K4) = (act(y(Rac1),rpar(:,172))*ymax(MAP3K4) - y(MAP3K4))/tau(MAP3K4); 
dydt(MAPKAPK) = (act(y(p38),rpar(:,155))*ymax(MAPKAPK) - y(MAPKAPK))/tau(MAPKAPK); 
dydt(MEF2) = (OR(inhib(y(HDAC),rpar(:,28)),OR(act(y(ERK5),rpar(:,81)),act(y(p38),rpar(:,156))))*ymax(MEF2) - y(MEF2))/tau(MEF2); 
dydt(MEK12) = (OR(act(y(MAP3K23),rpar(:,118)),OR(act(y(MEKK1),rpar(:,136)),act(y(Raf1),rpar(:,173))))*ymax(MEK12) - y(MEK12))/tau(MEK12); 
dydt(MEK36) = (OR(act(y(MAP3K11),rpar(:,116)),act(y(TAK1),rpar(:,186)))*ymax(MEK36) - y(MEK36))/tau(MEK36); 
dydt(MEK4) = (OR(act(y(MAP3K23),rpar(:,119)),OR(act(y(MAP3K4),rpar(:,122)),act(y(MEKK1),rpar(:,137))))*ymax(MEK4) - y(MEK4))/tau(MEK4); 
dydt(MEK5) = (OR(act(y(MAP3K23),rpar(:,120)),act(y(SHP2),rpar(:,182)))*ymax(MEK5) - y(MEK5))/tau(MEK5); 
dydt(MEK7) = (OR(act(y(MAP3K11),rpar(:,117)),OR(act(y(MAP3K23),rpar(:,121)),OR(act(y(MAP3K4),rpar(:,123)),act(y(MEKK1),rpar(:,138)))))*ymax(MEK7) - y(MEK7))/tau(MEK7); 
dydt(MEKK1) = (act(y(Ras),rpar(:,176))*ymax(MEKK1) - y(MEKK1))/tau(MEKK1); 
dydt(MSK1) = (OR(act(y(ERK12),rpar(:,79)),act(y(p38),rpar(:,157)))*ymax(MSK1) - y(MSK1))/tau(MSK1); 
dydt(mTor) = (act(y(Akt),rpar(:,35))*ymax(mTor) - y(mTor))/tau(mTor); 
dydt(NE) = (rpar(1,12)*ymax(NE) - y(NE))/tau(NE); 
dydt(NFAT) = (OR(AND(rpar(:,24),inhib(y(GSK3B),rpar(:,24)),inhib(y(JNK),rpar(:,24)),inhib(y(p38),rpar(:,24)),inhib(y(PKA),rpar(:,24)),inhib(y(PKG1),rpar(:,24))),OR(act(y(CaN),rpar(:,51)),AND(rpar(:,74),act(y(CaN),rpar(:,74)),act(y(ERK12),rpar(:,74)))))*ymax(NFAT) - y(NFAT))/tau(NFAT); 
dydt(NFkB) = (OR(inhib(y(IkB),rpar(:,29)),act(y(ERK12),rpar(:,80)))*ymax(NFkB) - y(NFkB))/tau(NFkB); 
dydt(NIK) = (act(y(TNFR),rpar(:,190))*ymax(NIK) - y(NIK))/tau(NIK); 
dydt(NOS) = (act(y(Akt),rpar(:,36))*ymax(NOS) - y(NOS))/tau(NOS); 
dydt(NRG1) = (rpar(1,13)*ymax(NRG1) - y(NRG1))/tau(NRG1); 
dydt(p38) = (OR(act(y(MEK36),rpar(:,131)),OR(act(y(MEK4),rpar(:,133)),act(y(DUSP1),rpar(:,206))))*ymax(p38) - y(p38))/tau(p38); 
dydt(p70s6k) = (act(y(mTor),rpar(:,141))*ymax(p70s6k) - y(p70s6k))/tau(p70s6k); 
dydt(PDK1) = (act(y(PI3K),rpar(:,160))*ymax(PDK1) - y(PDK1))/tau(PDK1); 
dydt(PE) = (rpar(1,14)*ymax(PE) - y(PE))/tau(PE); 
dydt(PI3K) = (OR(act(y(EGFR),rpar(:,67)),OR(act(y(ERBB),rpar(:,71)),OR(act(y(GBG),rpar(:,92)),OR(act(y(IGF1R),rpar(:,101)),OR(act(y(JAK),rpar(:,109)),OR(act(y(Ras),rpar(:,177)),act(y(TNFR),rpar(:,191))))))))*ymax(PI3K) - y(PI3K))/tau(PI3K); 
dydt(PKA) = (act(y(cAMP),rpar(:,50))*ymax(PKA) - y(PKA))/tau(PKA); 
dydt(PKC) = (OR(AND(rpar(:,65),act(y(Calcium),rpar(:,65)),act(y(DAG),rpar(:,65))),act(y(TGFR),rpar(:,188)))*ymax(PKC) - y(PKC))/tau(PKC); 
dydt(PKD) = (act(y(PKC),rpar(:,163))*ymax(PKD) - y(PKD))/tau(PKD); 
dydt(PKG1) = (act(y(cGMP),rpar(:,56))*ymax(PKG1) - y(PKG1))/tau(PKG1); 
dydt(PLCB) = (OR(act(y(Gaq11),rpar(:,88)),act(y(IGF1R),rpar(:,102)))*ymax(PLCB) - y(PLCB))/tau(PLCB); 
dydt(PLCG) = (OR(act(y(EGFR),rpar(:,68)),act(y(ERBB),rpar(:,72)))*ymax(PLCG) - y(PLCG))/tau(PLCG); 
dydt(Rac1) = (act(y(Ras),rpar(:,178))*ymax(Rac1) - y(Rac1))/tau(Rac1); 
dydt(Raf1) = (AND(rpar(:,31),inhib(y(PKA),rpar(:,31)),act(y(Raf1A),rpar(:,31)))*ymax(Raf1) - y(Raf1))/tau(Raf1); 
dydt(Raf1A) = (OR(act(y(GBG),rpar(:,93)),OR(act(y(PKC),rpar(:,164)),act(y(Ras),rpar(:,179))))*ymax(Raf1A) - y(Raf1A))/tau(Raf1A); 
dydt(Ras) = (OR(act(y(EGFR),rpar(:,69)),OR(act(y(ERBB),rpar(:,73)),OR(act(y(FAK),rpar(:,84)),OR(act(y(FGFR),rpar(:,86)),OR(AND(rpar(:,94),act(y(CaMK),rpar(:,94)),act(y(CaN),rpar(:,94)),act(y(GBG),rpar(:,94))),OR(act(y(IGF1R),rpar(:,103)),OR(act(y(JAK),rpar(:,110)),act(y(PKC),rpar(:,165)))))))))*ymax(Ras) - y(Ras))/tau(Ras); 
dydt(RhoA) = (AND(rpar(:,174),act(y(Ras),rpar(:,174)),inhib(y(SHP2),rpar(:,174)))*ymax(RhoA) - y(RhoA))/tau(RhoA); 
dydt(sACT) = (OR(AND(rpar(:,52),act(y(cFos),rpar(:,52)),act(y(cJun),rpar(:,52)),act(y(SRF),rpar(:,52))),OR(AND(rpar(:,57),act(y(cJun),rpar(:,57)),act(y(SRF),rpar(:,57))),OR(AND(rpar(:,89),act(y(GATA4),rpar(:,89)),act(y(SRF),rpar(:,89))),OR(act(y(MEF2),rpar(:,129)),act(y(NFAT),rpar(:,147))))))*ymax(sACT) - y(sACT))/tau(sACT); 
dydt(SERCA) = (AND(rpar(:,22),inhib(y(cFos),rpar(:,22)),inhib(y(cJun),rpar(:,22)),inhib(y(NFAT),rpar(:,22)))*ymax(SERCA) - y(SERCA))/tau(SERCA); 
dydt(sGC) = (act(y(NOS),rpar(:,149))*ymax(sGC) - y(sGC))/tau(sGC); 
dydt(SHP2) = (act(y(gp130LIFR),rpar(:,97))*ymax(SHP2) - y(SHP2))/tau(SHP2); 
dydt(SRF) = (act(y(RhoA),rpar(:,180))*ymax(SRF) - y(SRF))/tau(SRF); 
dydt(STAT) = (act(y(JAK),rpar(:,111))*ymax(STAT) - y(STAT))/tau(STAT); 
dydt(Stretch) = (rpar(1,15)*ymax(Stretch) - y(Stretch))/tau(Stretch); 
dydt(TAK1) = (act(y(PKC),rpar(:,166))*ymax(TAK1) - y(TAK1))/tau(TAK1); 
dydt(TGFB) = (rpar(1,16)*ymax(TGFB) - y(TGFB))/tau(TGFB); 
dydt(TGFR) = (act(y(TGFB),rpar(:,187))*ymax(TGFR) - y(TGFR))/tau(TGFR); 
dydt(TNFa) = (rpar(1,17)*ymax(TNFa) - y(TNFa))/tau(TNFa); 
dydt(TNFR) = (act(y(TNFa),rpar(:,189))*ymax(TNFR) - y(TNFR))/tau(TNFR); 
dydt(PGR) = (act(y(lig),rpar(:,193))*ymax(PGR) - y(PGR))/tau(PGR); 
dydt(NR3C1) = (act(y(lig),rpar(:,194))*ymax(NR3C1) - y(NR3C1))/tau(NR3C1); 
dydt(CDK2) = (act(y(PGR),rpar(:,196))*ymax(CDK2) - y(CDK2))/tau(CDK2); 
dydt(CEBPB) = (act(y(NR3C1),rpar(:,202))*ymax(CEBPB) - y(CEBPB))/tau(CEBPB); 
dydt(DUSP1) = (act(y(NR3C1),rpar(:,197))*ymax(DUSP1) - y(DUSP1))/tau(DUSP1); 
dydt(HDAC1) = (act(y(NR3C1),rpar(:,199))*ymax(HDAC1) - y(HDAC1))/tau(HDAC1); 
dydt(LCK) = (act(y(NR3C1),rpar(:,201))*ymax(LCK) - y(LCK))/tau(LCK); 
dydt(lig) = (rpar(1,195)*ymax(lig) - y(lig))/tau(lig); 

% utility functions 
function fact = act(x,rpar) 
% hill activation function with parameters w (weight), n (Hill coeff), EC50 
    w = rpar(1); 
    n = rpar(2); 
    EC50 = rpar(3); 
    beta = (EC50.^n - 1)./(2*EC50.^n - 1); 
    K = (beta - 1).^(1./n); 
    fact = w.*(beta.*x.^n)./(K.^n + x.^n); 
    if fact>w,                 % cap fact(x)<= 1 
        fact = w; 
    end
 
function finhib = inhib(x,rpar) 
% inverse hill function with parameters w (weight), n (Hill coeff), EC50 
    finhib = 1 - act(x,rpar);
 
function z = OR(x,y) 
% OR logic gate 
    z = 1-(1-x)*(1-y);
 
function z = AND(rpar,varargin) 
% AND logic gate, multiplying all of the reactants together 
    w = rpar(1); 
        v = real(cell2mat(varargin)); 
        h = (sum(v)/(nargin-1))^(nargin-2); 
    if w == 0 
        z = 0; 
    elseif h == 0 
        z = 0; 
    else 
        z = prod(v)/h; 
    end 
